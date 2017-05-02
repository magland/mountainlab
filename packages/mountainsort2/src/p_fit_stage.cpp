#include "p_fit_stage.h"

#include <QTime>
#include <diskreadmda.h>
#include <diskreadmda32.h>
#include <mda.h>
#include <mda32.h>
#include "get_sort_indices.h"
#include "mlcommon.h"

typedef QList<bigint> IntList;

namespace P_fit_stage {
Mda sort_firings_by_time(const Mda& firings);
void compute_templates(Mda32& templates_out, Mda32& templates_stdevs_out, const DiskReadMda32& X, const QVector<double>& times, const QVector<bigint>& labels, bigint clip_size);
QVector<bigint> fit_stage_kernel(Mda32& X, Mda32& templates, QVector<double>& times, QVector<bigint>& labels, const Fit_stage_opts& opts, const QList<IntList>& time_channel_mask);
QList<bigint> get_time_channel_mask(const Mda32& template0, const Mda32& template0_stdev, double thresh);
}

bool p_fit_stage(QString timeseries_path, QString firings_path, QString firings_out_path, Fit_stage_opts opts)
{
    //The timeseries data and the dimensions
    DiskReadMda32 X(timeseries_path);
    bigint M = X.N1();
    bigint N = X.N2();
    bigint T = opts.clip_size;

    //Read in the firings array
    qDebug().noquote() << "Reading firings...";
    Mda firingsA;
    firingsA.read(firings_path);
    qDebug().noquote() << "Sorting events...";
    Mda firings = P_fit_stage::sort_firings_by_time(firingsA);

    qDebug().noquote() << "Extracting times and labels...";
    //L is the number of events. Accumulate vectors of times and labels for convenience
    bigint L = firings.N2();
    QVector<double> times;
    QVector<bigint> labels;
    for (bigint j = 0; j < L; j++) {
        times << firings.value(1, j);
        labels << (bigint)firings.value(2, j); //these are labels of clusters
    }

    qDebug().noquote() << "Computing templates...";
    //These are the templates corresponding to the clusters
    Mda32 templates;
    Mda32 templates_stdev;
    P_fit_stage::compute_templates(templates, templates_stdev, X, times, labels, T); //MxTxK

    bigint processing_chunk_size = 1e5;
    bigint processing_chunk_overlap_size = 1e3;

    //Now we do the processing in chunks
    bigint chunk_size = processing_chunk_size;
    bigint overlap_size = processing_chunk_overlap_size;
    if (N < processing_chunk_size) {
        chunk_size = N;
        overlap_size = 0;
    }

    double time_channel_mask_thresh = 0.5;

    bigint K = MLCompute::max(labels);
    QList<IntList> time_channel_mask;
    for (bigint i = 0; i < K; i++) {
        Mda32 template0, template0_stdev;
        templates.getChunk(template0, 0, 0, i, M, T, 1);
        templates_stdev.getChunk(template0_stdev, 0, 0, i, M, T, 1);
        time_channel_mask << P_fit_stage::get_time_channel_mask(template0, template0_stdev, time_channel_mask_thresh); //use only the channels with highest maxval
        qDebug().noquote() << QString("k=%1, mask.size=%2").arg(i + 1).arg(time_channel_mask[i].count());
    }

    QList<bigint> inds_to_use;
    printf("Starting fit stage...\n");
    QTime timer;
    timer.start();
    {
        bigint num_timepoints_handled = 0;
#pragma omp parallel for
        for (bigint timepoint = 0; timepoint < N; timepoint += chunk_size) {
            //for (bigint timepoint = 0; timepoint < N; timepoint += N) { //for debugging
            //QMap<QString, bigint> elapsed_times_local;
            Mda32 chunk; //this will be the chunk we are working on
            Mda32 local_templates; //just a local copy of the templates
            QVector<double> local_times; //the times that fall in this time range
            QVector<bigint> local_labels; //the corresponding labels
            QList<bigint> local_inds; //the corresponding event indices
            Fit_stage_opts local_opts; //a local copy of the opts
            QList<IntList> local_time_channel_mask;
#pragma omp critical(lock1)
            {
                //build the variables above
                local_templates = templates;
                local_opts = opts;
                local_time_channel_mask = time_channel_mask;
                if (!X.readChunk(chunk, 0, timepoint - overlap_size, M, chunk_size + 2 * overlap_size)) {
                    qWarning() << "Problem reading chunk in fit_stage";
                }
                for (bigint jj = 0; jj < L; jj++) {
                    if ((timepoint - overlap_size <= times[jj]) && (times[jj] < timepoint - overlap_size + chunk_size + 2 * overlap_size)) {
                        local_times << times[jj] - (timepoint - overlap_size);
                        local_labels << labels[jj];
                        local_inds << jj;
                    }
                }
            }
            //Our real task is to decide which of these events to keep. Those will be stored in local_inds_to_use
            //"Local" means this chunk in this thread
            QVector<bigint> local_inds_to_use;
            {
                //This is the main kernel operation!!
                local_inds_to_use = P_fit_stage::fit_stage_kernel(chunk, local_templates, local_times, local_labels, local_opts, time_channel_mask);
            }
#pragma omp critical(lock1)
            {
                {
                    for (bigint ii = 0; ii < local_inds_to_use.count(); ii++) {
                        bigint ind0 = local_inds[local_inds_to_use[ii]];
                        double t0 = times[ind0];
                        if ((timepoint <= t0) && (t0 < timepoint + chunk_size)) {
                            inds_to_use << ind0;
                        }
                    }
                }

                num_timepoints_handled += qMin(chunk_size, N - timepoint);
                if (timer.elapsed() > 5000) {
                    qDebug().noquote() << QString("--- Handled %1% of timepoints").arg((int)(num_timepoints_handled * 100.0 / N)) << timer.elapsed();
                    timer.restart();
                }
            }
        }
    }

    qDebug().noquote() << "Setting firings_out...";
    qSort(inds_to_use);
    bigint num_to_use = inds_to_use.count();
    if (times.count()) {
        qDebug().noquote() << QString("using %1/%2 events (%3%)").arg(num_to_use).arg((bigint)times.count()).arg(num_to_use * 100.0 / times.count());
    }
    Mda firings_out(firings.N1(), num_to_use);
    for (bigint i = 0; i < num_to_use; i++) {
        for (bigint j = 0; j < firings.N1(); j++) {
            firings_out.set(firings.get(j, inds_to_use[i]), j, i);
        }
    }

    return firings_out.write64(firings_out_path);
}

namespace P_fit_stage {

Mda sort_firings_by_time(const Mda& firings)
{
    QVector<double> times;
    for (bigint i = 0; i < firings.N2(); i++) {
        times << firings.value(1, i);
    }
    QList<bigint> sort_inds = get_sort_indices_bigint(times);

    Mda F(firings.N1(), firings.N2());
    for (bigint i = 0; i < firings.N2(); i++) {
        for (bigint j = 0; j < firings.N1(); j++) {
            F.setValue(firings.value(j, sort_inds[i]), j, i);
        }
    }

    return F;
}

void compute_templates(Mda32& templates_out, Mda32& templates_stdevs_out, const DiskReadMda32& X, const QVector<double>& times, const QVector<bigint>& labels, bigint clip_size)
{
    bigint M = X.N1();
    bigint T = clip_size;
    bigint L = times.count();

    bigint K = MLCompute::max<bigint>(labels);

    bigint Tmid = (bigint)((T + 1) / 2) - 1;

    Mda sums(M, T, K);
    Mda sumsqrs(M, T, K);
    QList<bigint> counts;
    for (bigint k = 0; k < K; k++)
        counts << 0;
    for (bigint i = 0; i < L; i++) {
        bigint k = labels[i];
        bigint t0 = (bigint)(times[i] + 0.5);
        if (k >= 1) {
            Mda32 X0;
            if (!X.readChunk(X0, 0, t0 - Tmid, M, T)) {
                qWarning() << "Problem reading chunk in compute_templates of fit_stage";
            }
            dtype32* Xptr = X0.dataPtr();
            double* sums_ptr = sums.dataPtr(0, 0, k - 1);
            double* sumsqrs_ptr = sumsqrs.dataPtr(0, 0, k - 1);
            for (bigint i = 0; i < M * T; i++) {
                sums_ptr[i] += Xptr[i];
                sumsqrs_ptr[i] += Xptr[i] * Xptr[i];
            }
            counts[k - 1]++;
        }
    }
    templates_out.allocate(M, T, K);
    templates_stdevs_out.allocate(M, T, K);
    for (bigint k = 0; k < K; k++) {
        for (bigint t = 0; t < T; t++) {
            for (bigint m = 0; m < M; m++) {
                double n = counts[k];
                double s1 = sums.get(m, t, k);
                double s2 = sumsqrs.get(m, t, k);
                double mu = 0;
                if (n) {
                    mu = s1 / n;
                    templates_out.set(mu, m, t, k);
                }
                if (n > 1) {
                    double val = sqrt((s2 - 2 * mu * s1 + n * mu * mu) / (n - 1));
                    templates_stdevs_out.set(val, m, t, k);
                }
                //finish
            }
        }
    }
}

QList<bigint> get_time_channel_mask(const Mda32& template0, const Mda32& template0_stdev, double thresh)
{
    bigint M = template0.N1();
    bigint T = template0.N2();
    QList<bigint> ret;
    int ii = 0;
    for (bigint t = 0; t < T; t++) {
        for (bigint m = 0; m < M; m++) {
            double val = qAbs(template0.value(m, t));
            double val1 = template0_stdev.value(m, t);
            if (val1)
                val = val / val1;
            if (val > thresh)
                ret << ii;
            ii++;
        }
    }
    return ret;
}

bool is_dirty(float* dirty_ptr, const QList<bigint>& tchmask)
{
    for (bigint i = 0; i < tchmask.count(); i++) {
        if (dirty_ptr[i])
            return true;
    }
    return false;
}

/*
double compute_score(bigint N, float* X, float* template0)
{
    Mda resid(1, N);
    double* resid_ptr = resid.dataPtr();
    for (bigint i = 0; i < N; i++)
        resid_ptr[i] = X[i] - template0[i];
    double norm1 = MLCompute::norm(N, X);
    double norm2 = MLCompute::norm(N, resid_ptr);
    return norm1 * norm1 - norm2 * norm2;
}
*/

double compute_score(bigint M, bigint T, float* X_ptr, float* template0, const QList<bigint>& tchmask)
{
    (void)M;
    (void)T;
    double before_sumsqr = 0;
    double after_sumsqr = 0;
    for (bigint j = 0; j < tchmask.count(); j++) {
        bigint i = tchmask[j];
        double val = X_ptr[i];
        before_sumsqr += val * val;
        val -= template0[i];
        after_sumsqr += val * val;
    }
    return before_sumsqr - after_sumsqr;
}

void subtract_scaled_template(bigint N, double* X, double* template0, double scale_min, double scale_max)
{
    double S12 = 0, S22 = 0;
    for (bigint i = 0; i < N; i++) {
        S22 += template0[i] * template0[i];
        S12 += X[i] * template0[i];
    }
    double alpha = 1;
    if (S22)
        alpha = S12 / S22;
    alpha = qMin(scale_max, qMax(scale_min, alpha));
    for (bigint i = 0; i < N; i++) {
        X[i] -= alpha * template0[i];
    }
}

void subtract_scaled_template(bigint M, bigint T, float* X_ptr, float* dirty_ptr, float* template0, const QList<bigint>& tchmask, double scale_min, double scale_max)
{
    (void)M;
    (void)T;
    double S12 = 0, S22 = 0;
    for (bigint j = 0; j < tchmask.count(); j++) {
        bigint i = tchmask[j];
        S22 += template0[i] * template0[i];
        S12 += X_ptr[i] * template0[i];
    }
    double alpha = 1;
    if (S22)
        alpha = S12 / S22;
    alpha = qMin(scale_max, qMax(scale_min, alpha));
    for (bigint j = 0; j < tchmask.count(); j++) {
        bigint i = tchmask[j];
        X_ptr[i] -= alpha * template0[i];
        dirty_ptr[i] = 1;
    }
}

QVector<bigint> find_events_to_use(const QVector<double>& times, const QVector<double>& scores, const Fit_stage_opts& opts)
{
    QVector<bigint> to_use;
    bigint L = times.count();
    for (bigint i = 0; i < L; i++)
        to_use << 0; //start out not using any
    for (bigint i = 0; i < L; i++) {
        if (scores[i] > 0) { //score has to at least be positive
            to_use[i] = 1; //for now we say we are using it
            {
                // but let's check nearby events that may have a larger score

                bigint j = i;
                while ((j >= 0) && (times[j] >= times[i] - opts.clip_size)) {
                    if ((i != j) && (scores[j] >= scores[i])) {
                        to_use[i] = 0; //actually not using it because there is something bigger to the left
                    }
                    j--;
                }
            }
            {
                bigint j = i;
                while ((j < times.count()) && (times[j] <= times[i] + opts.clip_size)) {
                    if (scores[j] > scores[i]) {
                        to_use[i] = 0; //actually not using it because there is something bigger to the right
                    }
                    j++;
                }
            }
        }
    }
    return to_use;
}

QVector<bigint> fit_stage_kernel(Mda32& X, Mda32& templates, QVector<double>& times, QVector<bigint>& labels, const Fit_stage_opts& opts, const QList<IntList>& time_channel_mask)
{
    bigint M = X.N1(); //the number of dimensions
    bigint T = opts.clip_size; //the clip size
    bigint Tmid = (bigint)((T + 1) / 2) - 1; //the center timepoint in a clip (zero-indexed)
    bigint L = times.count(); //number of events we are looking at
    bigint K = MLCompute::max<bigint>(labels); //the maximum label number

    //compute the L2-norms of the templates ahead of time
    QVector<double> template_norms;
    template_norms << 0;
    for (bigint k = 1; k <= K; k++) {
        template_norms << MLCompute::norm(M * T, templates.dataPtr(0, 0, k - 1));
    }

    //keep passing through the data until nothing changes anymore
    bool something_changed = true;
    bigint num_passes = 0;
    //while ((something_changed)&&(num_passes<2)) {

    QVector<double> scores(L);
    Mda32 dirty(X.N1(), X.N2()); //timepoints/channels for which we need to recompute the score if an event overlaps
    for (bigint ii = 0; ii < dirty.totalSize(); ii++) {
        dirty.setValue(1, ii);
    }

    QVector<bigint> event_inds_to_consider;
    for (bigint i = 0; i < L; i++)
        event_inds_to_consider << i;
    QVector<bigint> event_inds_to_use;

    bigint num_score_computes = 0;
    bigint num_to_use = 0;
    while (something_changed) {
        num_passes++;

        QList<bigint> inds_to_try; //indices of the events to try on this pass
        QVector<double> scores_to_try;
        QVector<double> times_to_try;
        QVector<bigint> labels_to_try;

        //QVector<double> template_norms_to_try;
        QVector<bigint> new_event_inds_to_consider;
        for (bigint kk = 0; kk < event_inds_to_consider.count(); kk++) { //loop through the events
            bigint i = event_inds_to_consider[kk];

            double t0 = times[i];
            bigint k0 = labels[i];
            if (k0 > 0) { //make sure we have a positive label (don't know why we wouldn't)
                bigint tt = (bigint)(t0 - Tmid + 0.5); //start time of clip
                double score0 = 0;
                if ((tt >= 0) && (tt + T <= X.N2())) { //make sure we are in range
                    IntList tchmask = time_channel_mask[k0 - 1];
                    if (!is_dirty(dirty.dataPtr(0, tt), tchmask)) {
                        // we don't need to recompute the score
                        score0 = scores[i];
                    }
                    else {
                        //we do need to recompute it.

                        //The score will be how much something like the L2-norm is decreased
                        score0 = compute_score(M, T, X.dataPtr(0, tt), templates.dataPtr(0, 0, k0 - 1), tchmask);
                        num_score_computes++;
                        /*
                        if (score0 < template_norms[k0] * template_norms[k0] * 0.1)
                            score0 = 0; //the norm of the improvement needs to be at least 0.5 times the norm of the template
                            */
                        scores[i] = score0;
                    }
                }
                //the score needs to be at least as large as neglogprior in order to accept the spike
                //double neglogprior = 30;
                double neglogprior = 0;
                if (score0 > neglogprior) {
                    //we are not committing to using this event yet... see below for next step
                    scores_to_try << score0;
                    times_to_try << t0;
                    labels_to_try << k0;
                    inds_to_try << i;
                }
                else {
                    //not gonna use it
                }
            }
        }
        //Look at those events to try and see if we should use them
        QVector<bigint> to_use = find_events_to_use(times_to_try, scores_to_try, opts);

        //at this point, nothing is dirty
        for (bigint i = 0; i < dirty.totalSize(); i++) {
            dirty.set(0, i);
        }

        //amplitude scaling for template subtraction
        double scale_min = 1;
        double scale_max = 1;

        //for all those we are going to "use", we want to subtract out the corresponding templates from the timeseries data
        something_changed = false;
        bigint num_added = 0;
        for (bigint aa = 0; aa < to_use.count(); aa++) {
            if (to_use[aa] == 1) {
                IntList tchmask = time_channel_mask[labels_to_try[aa] - 1];
                something_changed = true;
                num_added++;
                bigint tt = (bigint)(times_to_try[aa] - Tmid + 0.5);
                subtract_scaled_template(M, T, X.dataPtr(0, tt), dirty.dataPtr(0, tt), templates.dataPtr(0, 0, labels_to_try[aa] - 1), tchmask, scale_min, scale_max);
                event_inds_to_use << inds_to_try[aa];
                num_to_use++;
            }
            else {
                //consider it for next time
                new_event_inds_to_consider << inds_to_try[aa];
            }
        }
        event_inds_to_consider = new_event_inds_to_consider;
    }

    //qDebug().noquote() << QString("Processed %1 events in time chunk (using %2%), %3 score computes per event, %4 passes").arg(L).arg(num_to_use * 100.0 / L).arg(num_score_computes * 1.0 / L).arg(num_passes);

    return event_inds_to_use;
}
}
