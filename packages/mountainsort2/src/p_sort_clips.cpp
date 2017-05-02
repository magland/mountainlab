#include "p_sort_clips.h"

#include <QTime>
#include <mda32.h>
#include "pca.h"
#include "isosplit5.h"
#include "mlcommon.h"
#include <QCoreApplication>

namespace P_sort_clips {
QVector<int> sort_clips_subset(const Mda32& clips, const QVector<bigint>& indices, Sort_clips_opts opts);
Mda32 dimension_reduce_clips(Mda32& clips, bigint num_features_per_channel, bigint max_samples);
Mda32 compute_templates(Mda32& clips, const QVector<int>& labels);
}

bool p_sort_clips(QString clips_path, QString labels_out, Sort_clips_opts opts)
{
    qDebug().noquote() << "p_sort_clips";
    Mda32 clips(clips_path);
    {
        QTime timer;
        timer.start();
        clips = P_sort_clips::dimension_reduce_clips(clips, opts.num_features, opts.max_samples);
        qDebug().noquote() << QString("Time elapsed for dimension reduction per channel (%1x%2x%3): %4 sec").arg(clips.N1()).arg(clips.N2()).arg(clips.N3()).arg(timer.elapsed() * 1.0 / 1000);
    }
    bigint M = clips.N1();
    bigint T = clips.N2();
    bigint L = clips.N3();

    QVector<bigint> indices(L);
    for (bigint i = 0; i < L; i++) {
        indices[i] = i;
    }

    qDebug().noquote() << "Sorting clips...";
    QVector<int> labels = P_sort_clips::sort_clips_subset(clips, indices, opts);

    if (opts.remove_outliers) {
        qDebug().noquote() << "Computing templates...";
        Mda32 templates0 = P_sort_clips::compute_templates(clips, labels);
        float* tptr = templates0.dataPtr();
        float* cptr = clips.dataPtr();
        qDebug().noquote() << "Removing outliers...";
        QVector<float> diff(M * T);
        bigint aa = 0;
        bigint num_outliers_removed = 0;
        QVector<double> clip_norms(L);
        QVector<double> template_norms(L);
        QVector<double> resid_norms(L);
        for (bigint i = 0; i < L; i++) {
            int k = labels[i];
            if (k > 0) {
                float* tptr0 = &tptr[(k - 1) * M * T];
                float* cptr0 = &cptr[aa];
                for (int bb = 0; bb < M * T; bb++) {
                    diff[bb] = cptr0[bb] - tptr0[bb];
                }
                clip_norms[i] = sqrt(MLCompute::dotProduct(M * T, &cptr[aa], &cptr[aa]));
                template_norms[i] = sqrt(MLCompute::dotProduct(M * T, tptr0, tptr0));
                resid_norms[i] = sqrt(MLCompute::dotProduct(M * T, diff.data(), diff.data()));
            }
            aa += M * T;
        }
        double mean_resid_norms = MLCompute::mean(resid_norms);
        double stdev_resid_norms = MLCompute::stdev(resid_norms);
        for (bigint i = 0; i < L; i++) {
            if (resid_norms[i] > clip_norms[i]) {
                //subtracting off the template made the L2-norm increase
                labels[i] = 0;
            }
            if (resid_norms[i] > mean_resid_norms + 4 * stdev_resid_norms) {
                //the resid is larger than most.
                labels[i] = 0;
            }
            if (!labels[i])
                num_outliers_removed++;
        }
        qDebug().noquote() << QString("Removed %1 outliers out of %2 (%3%)").arg(num_outliers_removed).arg(L).arg(num_outliers_removed * 100.0 / L);
    }

    Mda ret(1, L);
    for (bigint i = 0; i < L; i++) {
        ret.setValue(labels[i], i);
    }

    return ret.write64(labels_out);
}

namespace P_sort_clips {

Mda32 dimension_reduce_clips(Mda32& clips, bigint num_features_per_channel, bigint max_samples)
{
    bigint M = clips.N1();
    bigint T = clips.N2();
    bigint L = clips.N3();
    float* clips_ptr = clips.dataPtr();

    qDebug().noquote() << QString("Dimension reduce clips %1x%2x%3").arg(M).arg(T).arg(L);

    Mda32 ret(M, num_features_per_channel, L);
    float* retptr = ret.dataPtr();
    for (bigint m = 0; m < M; m++) {
        Mda32 reshaped(T, L);
        float* reshaped_ptr = reshaped.dataPtr();
        bigint aa = 0;
        bigint bb = m;
        for (bigint i = 0; i < L; i++) {
            for (bigint t = 0; t < T; t++) {
                //reshaped.set(clips.get(bb),aa);
                reshaped_ptr[aa] = clips_ptr[bb];
                aa++;
                bb += M;
                //reshaped.setValue(clips.value(m, t, i), t, i);
            }
        }
        Mda32 CC, FF, sigma;
        pca_subsampled(CC, FF, sigma, reshaped, num_features_per_channel, false, max_samples);
        float* FF_ptr = FF.dataPtr();
        aa = 0;
        bb = m;
        for (bigint i = 0; i < L; i++) {
            for (bigint a = 0; a < num_features_per_channel; a++) {
                //ret.setValue(FF.value(a, i), m, a, i);
                //ret.set(FF.get(aa),bb);
                retptr[bb] = FF_ptr[aa];
                aa++;
                bb += M;
            }
        }
    }
    return ret;
}

QVector<int> sort_clips_subset(const Mda32& clips, const QVector<bigint>& indices, Sort_clips_opts opts)
{
    bigint M = clips.N1();
    bigint T = clips.N2();
    bigint L0 = indices.count();

    qDebug().noquote() << QString("Sorting clips %1x%2x%3").arg(M).arg(T).arg(L0);

    Mda32 FF;
    {
        // do this inside a code block so memory gets released
        Mda32 clips_reshaped(M * T, L0);
        bigint iii = 0;
        for (bigint ii = 0; ii < L0; ii++) {
            for (bigint t = 0; t < T; t++) {
                for (bigint m = 0; m < M; m++) {
                    clips_reshaped.set(clips.value(m, t, indices[ii]), iii);
                    iii++;
                }
            }
        }

        Mda32 CC, sigma;
        QTime timer;
        timer.start();
        pca_subsampled(CC, FF, sigma, clips_reshaped, opts.num_features, false, opts.max_samples); //should we subtract the mean?
        qDebug().noquote() << QString("Time elapsed for pca (%1x%2x%3): %4 sec").arg(M).arg(T).arg(L0).arg(timer.elapsed() * 1.0 / 1000);
    }

    isosplit5_opts i5_opts;
    i5_opts.isocut_threshold = opts.isocut_threshold;
    i5_opts.K_init = opts.K_init;

    QVector<int> labels0(L0);
    {
        QTime timer;
        timer.start();
        isosplit5(labels0.data(), opts.num_features, L0, FF.dataPtr(), i5_opts);
        qDebug().noquote() << QString("Time elapsed for isosplit (%1x%2) - K=%3: %4 sec").arg(FF.N1()).arg(FF.N2()).arg(MLCompute::max(labels0)).arg(timer.elapsed() * 1.0 / 1000);
    }

    bigint K0 = MLCompute::max(labels0);
    if (K0 <= 1) {
        //only 1 cluster
        return labels0;
    }
    else {
        //branch method
        QVector<int> labels_new(L0);
        bigint k_offset = 0;
        for (bigint k = 1; k <= K0; k++) {
            QVector<bigint> indices2, indices2b;
            for (bigint a = 0; a < L0; a++) {
                if (labels0[a] == k) {
                    indices2 << indices[a];
                    indices2b << a;
                }
            }
            if (indices2.count() > 0) {
                QVector<int> labels2 = sort_clips_subset(clips, indices2, opts);
                bigint K2 = MLCompute::max(labels2);
                for (bigint bb = 0; bb < indices2.count(); bb++) {
                    labels_new[indices2b[bb]] = k_offset + labels2[bb];
                }
                k_offset += K2;
            }
        }

        return labels_new;
    }
}
Mda32 compute_templates(Mda32& clips, const QVector<int>& labels)
{
    int M = clips.N1();
    int T = clips.N2();
    bigint L = clips.N3();
    int Kmax = MLCompute::max(labels);
    Mda sum(M, T, Kmax);
    QVector<double> counts(Kmax);
    double* sumptr = sum.dataPtr();
    float* cptr = clips.dataPtr();
    bigint aa = 0;
    for (bigint i = 0; i < L; i++) {
        int k = labels[i];
        if (k > 0) {
            bigint bb = M * T * (k - 1);
            for (bigint jj = 0; jj < M * T; jj++) {
                sumptr[bb] += cptr[aa];
                aa++;
                bb++;
            }
            counts[k - 1]++;
        }
    }
    Mda32 ret(M, T, Kmax);
    for (int k = 0; k < Kmax; k++) {
        double factor = 0;
        if (counts[k])
            factor = 1.0 / counts[k];
        for (int t = 0; t < T; t++) {
            for (int m = 0; m < M; m++) {
                ret.setValue(sum.value(m, t, k) * factor, m, t, k);
            }
        }
    }
    return ret;
}
}

namespace P_reorder_labels {
struct template_comparer_struct {
    int channel;
    double template_peak;
    int index;
};
struct template_comparer {
    bool operator()(const template_comparer_struct& a, const template_comparer_struct& b) const
    {
        //put the zeros at the end
        if ((a.template_peak == 0) && (b.template_peak > 0))
            return false;
        if ((a.template_peak > 0) && (b.template_peak == 0))
            return true;

        if (a.channel < b.channel)
            return true;
        else if (a.channel == b.channel) {
            if (a.template_peak > b.template_peak)
                return true;
            else if (a.template_peak == b.template_peak)
                return (a.index < b.index);
            else
                return false;
        }
        else
            return false;
    }
};
template_comparer_struct compute_comparer(const Mda32& template0, int index)
{
    //QList<double> abs_peak_values;
    template_comparer_struct ret;
    int peak_channel = 1;
    double abs_peak_value = 0;
    for (int t = 0; t < template0.N2(); t++) {
        for (int m = 0; m < template0.N1(); m++) {
            double val = template0.get(m, t);
            if (fabs(val) > abs_peak_value) {
                abs_peak_value = fabs(val);
                peak_channel = m + 1;
            }
        }
    }
    ret.index = index;
    ret.channel = peak_channel;
    ret.template_peak = abs_peak_value;
    return ret;
}
}

bool p_reorder_labels(QString templates_path, QString firings_path, QString firings_out_path)
{
    Mda32 templates(templates_path);
    Mda firings(firings_path);
    qDebug().noquote() << "p_reorder_labels" << templates.N1() << templates.N2() << templates.N3() << firings.N1() << firings.N2();
    QList<P_reorder_labels::template_comparer_struct> list;
    for (int i = 0; i < templates.N3(); i++) {
        Mda32 template0;
        templates.getChunk(template0, 0, 0, i, templates.N1(), templates.N2(), 1);
        list << P_reorder_labels::compute_comparer(template0, i);
    }
    qSort(list.begin(), list.end(), P_reorder_labels::template_comparer());
    //QList<int> sort_indices;
    QMap<int, int> label_map;
    for (int i = 0; i < list.count(); i++) {
        label_map[list[i].index + 1] = i + 1;
    }
    qDebug().noquote() << "label_map" << label_map;
    for (bigint i = 0; i < firings.N2(); i++) {
        int label0 = firings.value(2, i);
        if (label_map.contains(label0)) {
            label0 = label_map[label0];
        }
        firings.setValue(label0, 2, i);
    }
    return firings.write64(firings_out_path);
}
