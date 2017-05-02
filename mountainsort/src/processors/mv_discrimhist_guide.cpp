#include "mv_discrimhist_guide.h"

#include "diskreadmda.h"
#include "extract_clips.h"
#include "mlcommon.h"
#include "compute_templates_0.h"
#include "msmisc.h"
#include "mv_discrimhist.h"

struct discrimhist_guide_data {
    int k1 = 0, k2 = 0;
    QVector<double> data1;
    QVector<double> data2;
    double separation_score = 0;
};

struct discrimhist_guide_data_comparer {
    bool operator()(const discrimhist_guide_data& a, const discrimhist_guide_data& b) const
    {
        return (a.separation_score < b.separation_score);
    }
};

/// TODO parallelize mv_distrimhist

//void get_discrimhist_guide_data(QVector<double>& ret1, QVector<double>& ret2, const DiskReadMda& timeseries, const DiskReadMda& firings, int k1, int k2, int clip_size);
Mda compute_distance_matrix(DiskReadMda32 timeseries, DiskReadMda firings, mv_discrimhist_guide_opts opts);
void get_pairs_to_compare(QVector<int>& ret_k1, QVector<int>& ret_k2, const Mda& distance_matrix, int num_histograms, QSet<int> clusters_to_exclude);
double compute_separation_score(const QVector<double>& data1, const QVector<double>& data2);

bool mv_discrimhist_guide(QString timeseries_path, QString firings_path, QString output_path, mv_discrimhist_guide_opts opts)
{
    DiskReadMda32 timeseries(timeseries_path);
    DiskReadMda firings(firings_path);

    Mda distance_matrix = compute_distance_matrix(DiskReadMda32(timeseries_path), DiskReadMda(firings_path), opts);
    QVector<int> k1s, k2s;
    get_pairs_to_compare(k1s, k2s, distance_matrix, opts.num_histograms * 2, opts.clusters_to_exclude); //get more than we need to be reduced later after sorting

    QList<discrimhist_guide_data> datas;
    for (int ii = 0; ii < k1s.count(); ii++) {
        int k1 = k1s[ii];
        int k2 = k2s[ii];
        discrimhist_guide_data DD;
        DD.k1 = k1;
        DD.k2 = k2;
        get_discrimhist_data(DD.data1, DD.data2, timeseries, firings, k1, k2, opts.clip_size, opts.method);
        DD.separation_score = compute_separation_score(DD.data1, DD.data2);
        datas << DD;
    }

    qSort(datas.begin(), datas.end(), discrimhist_guide_data_comparer());

    datas = datas.mid(0, opts.num_histograms);

    int total_count = 0;
    for (int i = 0; i < datas.count(); i++) {
        total_count += datas[i].data1.count();
        total_count += datas[i].data2.count();
    }

    Mda output(4, total_count);
    //first two rows are k1/k2, third row is k1 or k2, fourth row is the value
    int jj = 0;
    for (int i = 0; i < datas.count(); i++) {
        int k1 = datas[i].k1;
        int k2 = datas[i].k2;
        for (int k = 0; k < datas[i].data1.count(); k++) {
            output.setValue(k1, 0, jj);
            output.setValue(k2, 1, jj);
            output.setValue(k1, 2, jj);
            output.setValue(datas[i].data1[k], 3, jj);
            jj++;
        }
        for (int k = 0; k < datas[i].data2.count(); k++) {
            output.setValue(k1, 0, jj);
            output.setValue(k2, 1, jj);
            output.setValue(k2, 2, jj);
            output.setValue(datas[i].data2[k], 3, jj);
            jj++;
        }
    }

    output.write32(output_path);

    return true;
}

/*
void get_discrimhist_guide_data(QVector<double>& ret1, QVector<double>& ret2, const DiskReadMda& timeseries, const DiskReadMda& firings, int k1, int k2, int clip_size)
{
    QVector<double> times1, times2;
    for (int i = 0; i < firings.N2(); i++) {
        int label = (int)firings.value(2, i);
        if (label == k1) {
            times1 << firings.value(1, i);
        }
        if (label == k2) {
            times2 << firings.value(1, i);
        }
    }
    Mda clips1 = extract_clips(timeseries, times1, clip_size);
    Mda clips2 = extract_clips(timeseries, times2, clip_size);

    Mda centroid1 = compute_mean_clip(clips1);
    Mda centroid2 = compute_mean_clip(clips2);

    Mda diff(centroid1.N1(), centroid1.N2());
    for (int i2 = 0; i2 < centroid1.N2(); i2++) {
        for (int i1 = 0; i1 < centroid1.N1(); i1++) {
            diff.setValue(centroid2.value(i1, i2) - centroid1.value(i1, i2), i1, i2);
        }
    }

    double* ptr_clips1 = clips1.dataPtr();
    double* ptr_clips2 = clips2.dataPtr();

    int N = centroid1.N1() * centroid1.N2();
    double* ptr_diff = diff.dataPtr();
    double norm0 = MLCompute::norm(N, ptr_diff);
    if (!norm0)
        norm0 = 1;

    ret1.clear();
    for (int i = 0; i < clips1.N3(); i++) {
        ret1 << MLCompute::dotProduct(N, ptr_diff, &ptr_clips1[N * i]) / (norm0 * norm0);
    }
    ret2.clear();
    for (int i = 0; i < clips2.N3(); i++) {
        ret2 << MLCompute::dotProduct(N, ptr_diff, &ptr_clips2[N * i]) / (norm0 * norm0);
    }
}
*/

double compute_distance(int N, float* ptr1, float* ptr2)
{
    double sumsqr = 0;
    for (int i = 0; i < N; i++) {
        double tmp = ptr1[i] - ptr2[i];
        sumsqr += tmp * tmp;
    }
    return sqrt(sumsqr);
}

bool is_zero(const Mda32& X)
{
    return ((X.minimum() == 0) && (X.maximum() == 0));
}

Mda compute_distance_matrix(DiskReadMda32 timeseries, DiskReadMda firings, mv_discrimhist_guide_opts opts)
{
    QVector<double> times;
    QVector<int> labels;
    for (int i = 0; i < firings.N2(); i++) {
        times << firings.value(1, i);
        labels << (int)firings.value(2, i);
    }
    Mda32 X = compute_templates_0(timeseries, times, labels, opts.clip_size);
    int M = X.N1();
    int T = X.N2();
    int K = X.N3();

    Mda ret(K, K);
    for (int k1 = 0; k1 < K; k1++) {
        Mda32 tmp1;
        X.getChunk(tmp1, 0, 0, k1, M, T, 1);
        for (int k2 = 0; k2 < K; k2++) {
            Mda32 tmp2;
            X.getChunk(tmp2, 0, 0, k2, M, T, 1);
            double dist = compute_distance(M * T, tmp1.dataPtr(), tmp2.dataPtr());
            if ((is_zero(tmp1)) || (is_zero(tmp2))) {
                ret.setValue(0, k1, k2);
            }
            else {
                ret.setValue(dist, k1, k2);
            }
        }
    }
    return ret;
}

struct pair_struct {
    double dist;
    int k1, k2;
};

struct pair_comparer {
    bool operator()(const pair_struct& a, const pair_struct& b) const
    {
        if (a.dist < b.dist)
            return true;
        else
            return false;
    }
};

void get_pairs_to_compare(QVector<int>& ret_k1, QVector<int>& ret_k2, const Mda& distance_matrix, int num_histograms, QSet<int> clusters_to_exclude)
{
    QList<pair_struct> pairs;
    for (int i1 = 0; i1 < distance_matrix.N1(); i1++) {
        for (int i2 = i1 + 1; i2 < distance_matrix.N2(); i2++) {
            pair_struct pp;
            pp.dist = distance_matrix.value(i1, i2);
            pp.k1 = i1 + 1;
            pp.k2 = i2 + 1;
            if ((!clusters_to_exclude.contains(pp.k1)) && (!clusters_to_exclude.contains(pp.k2))) {
                if (pp.dist) {
                    pairs << pp;
                }
            }
        }
    }

    qSort(pairs.begin(), pairs.end(), pair_comparer());

    ret_k1.clear();
    ret_k2.clear();
    for (int ii = 0; ii < qMin(num_histograms, pairs.count()); ii++) {
        ret_k1 << pairs[ii].k1;
        ret_k2 << pairs[ii].k2;
    }
}

double compute_separation_score(const QVector<double>& data1, const QVector<double>& data2)
{
    //return data1.count();
    double mean1 = MLCompute::mean(data1);
    double mean2 = MLCompute::mean(data2);
    QVector<double> mean_subtracted_vals;
    for (int i = 0; i < data1.count(); i++) {
        mean_subtracted_vals << data1[i] - mean1;
    }
    for (int i = 0; i < data2.count(); i++) {
        mean_subtracted_vals << data2[i] - mean2;
    }

    double stdev = MLCompute::stdev(mean_subtracted_vals);
    if (stdev) {
        return (qAbs(mean2 - mean1) / stdev);
    }
    else {
        return 0;
    }
}
