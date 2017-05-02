#include "p_consolidate_clusters.h"
#include <QFile>
#include <diskreadmda32.h>
#include <mda.h>
#include <mda32.h>
#include "mlcommon.h"

namespace P_consolidate_clusters {
QVector<int> read_labels(QString path);
QVector<double> read_times(QString path);
bool write_labels(QString path, const QVector<int>& labels);
//Mda32 compute_template(QString clips_path, const QVector<int>& labels, int k);
bool should_use_template(const Mda32& template0, Consolidate_clusters_opts opts);
Mda32 compute_templates(const DiskReadMda32& X, const QVector<double>& times, const QVector<int>& labels, int clip_size);
}

bool p_consolidate_clusters(QString timeseries, QString event_times_path, QString labels_path, QString labels_out, Consolidate_clusters_opts opts)
{
    QVector<int> labels = P_consolidate_clusters::read_labels(labels_path);
    QVector<double> times = P_consolidate_clusters::read_times(event_times_path);
    bigint L = labels.count();

    int K = MLCompute::max(labels);

    Mda32 templates = P_consolidate_clusters::compute_templates(timeseries, times, labels, opts.clip_size);

    QVector<int> to_use(K + 1);
    to_use.fill(0);

    for (int k = 1; k <= K; k++) {
        Mda32 template0;
        templates.getChunk(template0, 0, 0, k - 1, templates.N1(), templates.N2(), 1);
        if (P_consolidate_clusters::should_use_template(template0, opts)) {
            to_use[k] = 1;
        }
    }

    QMap<int, int> label_map;
    label_map[0] = 0;
    int knext = 1;
    for (int k = 1; k <= K; k++) {
        if (to_use[k]) {
            label_map[k] = knext;
            knext++;
        }
        else {
            label_map[k] = 0;
        }
    };

    QVector<int> new_labels(L);
    for (int i = 0; i < L; i++) {
        new_labels[i] = label_map.value(labels[i]);
    }

    return P_consolidate_clusters::write_labels(labels_out, new_labels);
}

namespace P_consolidate_clusters {
QVector<int> read_labels(QString path)
{
    Mda X(path);
    bigint L = X.totalSize();
    QVector<int> ret(L);
    for (int i = 0; i < L; i++)
        ret[i] = X.value(i);
    return ret;
}
QVector<double> read_times(QString path)
{
    Mda X(path);
    bigint L = X.totalSize();
    QVector<double> ret(L);
    for (int i = 0; i < L; i++)
        ret[i] = X.value(i);
    return ret;
}
bool write_labels(QString path, const QVector<int>& labels)
{
    Mda X(1, labels.count());
    for (int i = 0; i < labels.count(); i++) {
        X.setValue(labels[i], i);
    }
    return X.write64(path);
}

/*
Mda32 compute_template(QString clips_path, const QVector<int>& labels, int k)
{
    DiskReadMda32 clips(clips_path);
    bigint M = clips.N1();
    bigint T = clips.N2();
    Mda sum(M, T);
    double ct = 0;
    for (bigint i = 0; i < labels.count(); i++) {
        if (labels[i] == k) {
            Mda32 tmp;
            if (!clips.readChunk(tmp, 0, 0, i, M, T, 1)) {
                qWarning() << "Problem reading chunk in compute_template of consolidate_clusters";
                return Mda32();
            }
            for (bigint j = 0; j < M * T; j++) {
                sum.set(sum.get(j) + tmp.get(j), j);
            }
            ct++;
        }
    }
    Mda32 ret(M, T);
    if (ct) {
        for (bigint j = 0; j < M * T; j++) {
            ret.set(sum.get(j) / ct, j);
        }
    }
    return ret;
}
*/

Mda32 compute_templates(const DiskReadMda32& X, const QVector<double>& times, const QVector<int>& labels, int clip_size)
{
    bigint M = X.N1();
    //bigint N = X.N2();
    bigint T = clip_size;

    int Kmax = 0;
    for (bigint i = 0; i < labels.count(); i++) {
        int label0 = labels[i];
        if (label0 > Kmax)
            Kmax = label0;
    }

    Mda sums(M, T, Kmax);
    QVector<double> counts(Kmax);
    counts.fill(0);

    double* sums_ptr = sums.dataPtr();

    bigint Tmid = (bigint)((T + 1) / 2) - 1;
    printf("computing templates (M=%ld,T=%ld,K=%d,L=%d)...\n", M, T, Kmax, times.count());
    for (bigint i = 0; i < times.count(); i++) {
        bigint t1 = times[i] - Tmid;
        //bigint t2 = t1 + T - 1;
        int k = labels[i];
        if (k > 0) {
            Mda32 tmp;
            tmp.allocate(M, T);
            if (!X.readChunk(tmp, 0, t1, M, T)) {
                qWarning() << "Problem reading chunk of timeseries list:" << t1;
                return false;
            }
            float* tmp_ptr = tmp.dataPtr();
            bigint offset = M * T * (k - 1);
            for (bigint aa = 0; aa < M * T; aa++) {
                sums_ptr[offset + aa] += tmp_ptr[aa];
            }
            counts[k - 1]++;
        }
    }

    Mda32 templates(M, T, Kmax);
    float* templates_ptr = templates.dataPtr();

    bigint bb = 0;
    for (bigint k = 0; k < Kmax; k++) {
        for (bigint aa = 0; aa < M * T; aa++) {
            if (counts[k])
                templates_ptr[bb] = sums_ptr[bb] / counts[k];
            bb++;
        }
    }

    return templates;
}

bool should_use_template(const Mda32& template0, Consolidate_clusters_opts opts)
{
    int peak_location_tolerance = 10;

    bigint M = template0.N1();
    bigint T = template0.N2();
    bigint Tmid = (int)((T + 1) / 2) - 1;

    double abs_peak_on_central_channel = 0;
    bigint abs_peak_t_on_central_channel = Tmid;
    if (opts.central_channel > 0) {
        for (bigint t = 0; t < T; t++) {
            double val = template0.value(opts.central_channel - 1, t);
            if (fabs(val) > abs_peak_on_central_channel) {
                abs_peak_on_central_channel = fabs(val);
                abs_peak_t_on_central_channel = t;
            }
        }
    }

    double abs_peak = 0;
    bigint abs_peak_t = Tmid;
    for (bigint t = 0; t < T; t++) {
        for (bigint m = 0; m < M; m++) {
            double val = template0.value(m, t);
            if (fabs(val) > abs_peak) {
                abs_peak = fabs(val);
                abs_peak_t = t;
            }
        }
    }
    if (opts.central_channel > 0) {
        if (abs_peak_on_central_channel < abs_peak * opts.consolidation_factor)
            return false;
        if (fabs(abs_peak_t_on_central_channel - Tmid) > peak_location_tolerance)
            return false;
    }
    else {
        if (fabs(abs_peak_t - Tmid) > peak_location_tolerance)
            return false;
    }

    return true;

    /*
    int peak_location_tolerance = 10;

    bigint M = template0.N1();
    bigint T = template0.N2();
    bigint Tmid = (int)((T + 1) / 2) - 1;


    if (opts.central_channel > 0) {
        QVector<double> energies(M);
        energies.fill(0);
        double max_energy = 0;
        for (bigint t = 0; t < T; t++) {
            for (bigint m = 0; m < M; m++) {
                double val = template0.value(m, t);
                energies[m] += val * val;
                if ((m != opts.central_channel - 1) && (energies[m] > max_energy))
                    max_energy = energies[m];
            }
        }
        if (energies[opts.central_channel - 1] < max_energy * opts.consolidation_factor)
            return false;
    }

    double abs_peak_val = 0;
    bigint abs_peak_ind = 0;
    bigint abs_peak_chan = 0;
    for (bigint t = 0; t < T; t++) {
        for (bigint m = 0; m < M; m++) {
            double value = template0.value(m, t);
            if (fabs(value) > abs_peak_val) {
                abs_peak_val = fabs(value);
                abs_peak_ind = t;
                abs_peak_chan = m + 1;
            }
        }
    }

    if (fabs(abs_peak_ind - Tmid) > peak_location_tolerance) {
        return false;
    }
    if ((opts.central_channel > 0) && (abs_peak_chan != opts.central_channel))
        return false;
    return true;
    */
}
}
