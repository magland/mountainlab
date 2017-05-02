/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
** Created: 4/7/2016
*******************************************************/

#include "compute_detectability_scores.h"

#include <diskreadmda.h>
#include "extract_clips.h"
#include "mlcommon.h"
#include "mlcommon.h"
#include <QTime>
#include <math.h>
#include "get_pca_features.h"
#include "msmisc.h"

bool compute_detectability_scores(QString timeseries_path, QString firings_path, QString firings_out_path, const compute_detectability_scores_opts& opts)
{

    /// TODO figure out what is the bottleneck and try to speed it up

    DiskReadMda32 timeseries;
    timeseries.setPath(timeseries_path);
    printf("reading firings array.\n");
    Mda firings;
    firings.read(firings_path);

    int T = opts.clip_size;
    int L = firings.N2();

    //prepare the channels, labels, times, and allocate the scores
    printf("Preparing data arrays\n");
    QVector<int> channels, labels;
    QVector<double> times;
    QVector<double> scores;
    for (int i = 0; i < L; i++) {
        channels << (int)firings.get(0, i) - 1; //convert to zero-based indexing
        times << firings.get(1, i) - 1; //convert to zero-based indexing
        labels << (int)firings.get(2, i);
        scores << 0;
    }

    //we used to extract the clips here... that was taking too much RAM
    //Mda clips=extract_clips(timeseries,times,T);

    int K = MLCompute::max<int>(labels);

    QList<Subcluster> subclusters;

#pragma omp parallel for
    for (int k = 1; k <= K; k++) { //iterate through all clusters
        int count1 = 0;
        int count2 = 0;

        QVector<int> inds_k;
#pragma omp critical
        inds_k = find_label_inds(labels, k); //find the indices corresonding to this cluster

        if (inds_k.count() > 0) { //if the cluster is non-empty
            Mda32 clips_k;
            Mda32 noise_shape;
            int channel;
#pragma omp critical
            {
                printf("k=%d/%d\n", k, K);
                QVector<double> times_k;
                for (int i = 0; i < inds_k.count(); i++) {
                    times_k << times[inds_k[i]];
                }
                channel = channels[inds_k[0]]; //the channel should be the same for all events in the cluster
                clips_k = extract_clips(timeseries, times_k, T);
                noise_shape = estimate_noise_shape(timeseries, T, channel); //estimate the noise shape corresponding to detection on this particular channel
            }
            Define_Shells_Opts opts2;
            opts2.min_shell_size = opts.min_shell_size;
            opts2.shell_increment = opts.shell_increment;
            QList<Subcluster> subclusters0 = compute_subcluster_detectability_scores(noise_shape, clips_k, channel, opts2); //actually compute the scores
#pragma omp critical
            {
                for (int ii = 0; ii < subclusters0.count(); ii++) {
                    Subcluster SC = subclusters0[ii];
                    QList<int> new_inds; //we want the inds to correspond to the original inds
                    for (int j = 0; j < SC.inds.count(); j++) {
                        new_inds << inds_k[SC.inds[j]];
                    }
                    SC.inds = new_inds;
                    SC.label = k;
                    subclusters << SC;
                    count1 += SC.inds.count();
                    if (SC.inds.count() > 0)
                        count2++;
                }
            }
        }
    }

    int num_subclusters = subclusters.count();

    //normalize by setting slope of peak vs score equal to 1
    QVector<double> subcluster_abs_peaks, subcluster_scores;
    for (int ii = 0; ii < num_subclusters; ii++) {
        subcluster_abs_peaks << fabs(subclusters[ii].peak);
        subcluster_scores << subclusters[ii].detectability_score;
    }
    double slope = compute_slope(subcluster_abs_peaks, subcluster_scores);
    if (slope) {
        for (int ii = 0; ii < num_subclusters; ii++) {
            subcluster_scores[ii] /= slope;
            subclusters[ii].detectability_score /= slope;
        }
    }

    for (int ii = 0; ii < num_subclusters; ii++) {
        QList<int> inds0 = subclusters[ii].inds;
        double score = subclusters[ii].detectability_score;
        for (int jj = 0; jj < inds0.count(); jj++) {
            scores[inds0[jj]] = score;
        }
    }

    Mda firings_out;
    int N1_new = qMax((long long)firings.N1(), 6LL);
    firings_out.allocate(N1_new, firings.N2());
    for (int j = 0; j < firings.N2(); j++) {
        for (int i = 0; i < firings.N1(); i++) {
            firings_out.setValue(firings.value(i, j), i, j);
        }
        firings_out.setValue(scores[j], 5, j);
    }

    firings_out.write64(firings_out_path);

    return true;
}

QVector<int> find_label_inds(const QVector<int>& labels, int k)
{
    QVector<int> ret;
    for (int i = 0; i < labels.count(); i++) {
        if (labels[i] == k)
            ret << i;
    }
    return ret;
}

Mda32 get_subclips(Mda32& clips, const QList<int>& inds)
{
    int M = clips.N1();
    int T = clips.N2();
    int L2 = inds.count();
    Mda32 ret;
    ret.allocate(M, T, L2);
    for (int i = 0; i < L2; i++) {
        int aaa = M * T * i;
        int bbb = M * T * inds[i];
        for (int t = 0; t < T; t++) {
            for (int m = 0; m < M; m++) {
                ret.set(clips.get(bbb), aaa);
                aaa++;
                bbb++;
            }
        }
    }
    return ret;
}

QList<Shell> define_shells(const QVector<double>& peaks, const Define_Shells_Opts& opts)
{
    QList<Shell> shells;
    if (peaks.count() == 0)
        return shells;

    //negatives and positives
    if (MLCompute::min(peaks) < 0) {
        QList<int> inds_neg, inds_pos;
        for (int i = 0; i < peaks.count(); i++) {
            if (peaks[i] < 0)
                inds_neg << i;
            else
                inds_pos << i;
        }
        QVector<double> neg_peaks_neg;
        for (int i = 0; i < inds_neg.count(); i++)
            neg_peaks_neg << -peaks[inds_neg[i]];
        QVector<double> peaks_pos;
        for (int i = 0; i < inds_pos.count(); i++)
            peaks_pos << peaks[inds_pos[i]];
        QList<Shell> shells_neg = define_shells(neg_peaks_neg, opts);
        QList<Shell> shells_pos = define_shells(peaks_pos, opts);
        for (int i = shells_neg.count() - 1; i >= 0; i--) {
            Shell SS = shells_neg[i];
            QList<int> new_inds;
            for (int i = 0; i < SS.inds.count(); i++)
                new_inds << inds_neg[SS.inds[i]];
            SS.inds = new_inds;
            shells << SS;
        }
        for (int i = 0; i < shells_pos.count(); i++) {
            Shell SS = shells_pos[i];
            QList<int> new_inds;
            for (int i = 0; i < SS.inds.count(); i++)
                new_inds << inds_pos[SS.inds[i]];
            SS.inds = new_inds;
            shells << SS;
        }
        return shells;
    }

    //only positives
    double max_peak = MLCompute::max(peaks);
    double peak_lower = 0;
    double peak_upper = peak_lower + opts.shell_increment;
    while (peak_lower <= max_peak) {
        QList<int> inds1;
        QList<int> inds2;
        for (int i = 0; i < peaks.count(); i++) {
            if ((peak_lower <= peaks[i]) && (peaks[i] < peak_upper)) {
                inds1 << i;
            }
            if (peaks[i] >= peak_upper) {
                inds2 << i;
            }
        }
        int ct1 = inds1.count();
        int ct2 = inds2.count();
        if (peak_upper > max_peak) {
            Shell SSS;
            SSS.inds = inds1;
            peak_lower = peak_upper;
            shells << SSS;
        }
        else if ((ct1 >= opts.min_shell_size) && (ct2 >= opts.min_shell_size)) {
            Shell SSS;
            SSS.inds = inds1;
            peak_lower = peak_upper;
            peak_upper = peak_lower + opts.shell_increment;
            shells << SSS;
        }
        else {
            peak_upper = peak_upper + opts.shell_increment;
        }
    }
    return shells;
}

QVector<double> randsample_with_replacement(int N, int K)
{
    QVector<double> ret;
    for (int i = 0; i < K; i++) {
        ret << qrand() % N;
    }
    return ret;
}

Mda32 estimate_noise_shape(DiskReadMda32& X, int T, int ch)
{
    // Computes the expected shape of the template in a noise cluster
    // which may be considered as a row in the noise covariance matrix
    // X is the MxN array of raw or pre-processed data
    // T is the clip size
    // ch is the channel where detection takes place
    QTime timer;
    timer.start();
    int M = X.N1();
    int N = X.N2();
    int Tmid = (int)((T + 1) / 2) - 1;
    int num_rand_times = 10000;
    QVector<double> rand_times = randsample_with_replacement(N - 2 * T, num_rand_times);
    for (int i = 0; i < rand_times.count(); i++)
        rand_times[i] += T;
    Mda32 rand_clips = extract_clips(X, rand_times, T);
    QVector<double> peaks;
    for (int i = 0; i < rand_times.count(); i++) {
        peaks << rand_clips.get(ch, Tmid, i);
    }

    Mda noise_shape_0(M, T);
    double* noise_shape_0_ptr = noise_shape_0.dataPtr();
    float* rand_clips_ptr = rand_clips.dataPtr();
    int bbb = 0;
    for (int i = 0; i < rand_times.count(); i++) {
        if (fabs(peaks[i]) <= 2) { //focus on noise clips (where amplitude is low)
            int aaa = 0;
            for (int t = 0; t < T; t++) {
                for (int m = 0; m < M; m++) {
                    noise_shape_0_ptr[aaa] += peaks[i] * rand_clips_ptr[bbb];
                    aaa++;
                    bbb++;
                }
            }
        }
        else
            bbb += M * T; //important and tricky!!!
    }
    double noise_shape_0_norm = 0;
    for (int t = 0; t < T; t++) {
        for (int m = 0; m < M; m++) {
            double val = noise_shape_0.get(m, t);
            noise_shape_0_norm += val * val;
        }
    }
    noise_shape_0_norm = sqrt(noise_shape_0_norm);
    if (!noise_shape_0_norm)
        noise_shape_0_norm = 1;

    Mda32 noise_shape;
    noise_shape.allocate(M, T);
    {
        for (int t = 0; t < T; t++) {
            for (int m = 0; m < M; m++) {
                noise_shape.set(noise_shape_0.get(m, t) / noise_shape_0_norm, m, t);
            }
        }
    }
    return noise_shape;
}

Mda compute_features(Mda32& clips, int num_features)
{
    int M = clips.N1();
    int T = clips.N2();
    int L = clips.N3();

    Mda clips_reshaped(M * T, L);
    int NNN = M * T * L;
    for (int iii = 0; iii < NNN; iii++) {
        clips_reshaped.set(clips.get(iii), iii);
    }

    Mda FF, CC, sigma;
    pca(CC, FF, sigma, clips_reshaped, num_features, false);

    return FF;

    //Mda ret;
    //ret.allocate(num_features, L);
    //get_pca_features(M * T, L, num_features, ret.dataPtr(), clips.dataPtr());
    //return ret;
}

void compute_geometric_median(int M, int N, double* output, double* input, int num_iterations)
{
    double* weights = (double*)malloc(sizeof(double) * N);
    double* dists = (double*)malloc(sizeof(double) * N);
    for (int j = 0; j < N; j++)
        weights[j] = 1;
    for (int it = 1; it <= num_iterations; it++) {
        float sumweights = 0;
        for (int j = 0; j < N; j++)
            sumweights += weights[j];
        if (sumweights)
            for (int j = 0; j < N; j++)
                weights[j] /= sumweights;
        for (int i = 0; i < M; i++)
            output[i] = 0;
        {
            //compute output
            int kkk = 0;
            for (int j = 0; j < N; j++) {
                int i = 0;
                for (int m = 0; m < M; m++) {
                    output[i] += weights[j] * input[kkk];
                    i++;
                    kkk++;
                }
            }
        }
        {
            //compute dists
            int kkk = 0;
            for (int j = 0; j < N; j++) {
                int i = 0;
                double sumsqr = 0;
                for (int m = 0; m < M; m++) {
                    double diff = output[i] - input[kkk];
                    i++;
                    kkk++;
                    sumsqr += diff * diff;
                }
                dists[j] = sqrt(sumsqr);
            }
        }
        {
            //compute weights
            for (int j = 0; j < N; j++) {
                if (dists[j])
                    weights[j] = 1 / dists[j];
                else
                    weights[j] = 0;
            }
        }
    }
    free(dists);
    free(weights);
}

Mda32 compute_geometric_median_template(Mda32& clips)
{
    int M = clips.N1();
    int T = clips.N2();
    int L = clips.N3();
    Mda32 ret;
    ret.allocate(M, T);
    if (L == 0)
        return ret;
    int num_features = 18;
    Mda FF = compute_features(clips, num_features);
    Mda FFmm;
    FFmm.allocate(num_features, 1);
    compute_geometric_median(num_features, L, FFmm.dataPtr(), FF.dataPtr());
    QVector<double> dists;
    for (int i = 0; i < L; i++) {
        double tmp = 0;
        for (int a = 0; a < num_features; a++) {
            double val = FF.get(a, i) - FFmm.get(a, 0L);
            tmp += val * val;
        }
        tmp = sqrt(tmp);
        dists << tmp;
    }
    QVector<double> sorted_dists = dists;
    qSort(sorted_dists);
    double dist_cutoff = sorted_dists[(int)(sorted_dists.count() * 0.3)];
    QList<int> inds;
    for (int i = 0; i < dists.count(); i++) {
        if (dists[i] <= dist_cutoff)
            inds << i;
    }
    Mda32 subclips = get_subclips(clips, inds);
    return compute_mean_clip(subclips);
}

double compute_template_ip(Mda32& T1, Mda32& T2)
{
    double ip = 0;
    for (int t = 0; t < T1.N2(); t++) {
        for (int m = 0; m < T1.N1(); m++) {
            ip += T1.get(m, t) * T2.get(m, t);
        }
    }
    return ip;
}

double compute_template_norm(Mda32& T)
{
    return sqrt(compute_template_ip(T, T));
}

QList<Subcluster> compute_subcluster_detectability_scores(Mda32& noise_shape, Mda32& clips, int channel, const Define_Shells_Opts& opts)
{
    QTime timer;
    timer.start();
    int M = clips.N1();
    int T = clips.N2();
    int L = clips.N3();
    int Tmid = (int)((T + 1) / 2) - 1;

    //get the list of peaks
    QVector<double> peaks;
    for (int i = 0; i < L; i++) {
        peaks << clips.get(channel, Tmid, i);
    }
    //define the shells
    QList<Shell> shells = define_shells(peaks, opts);

    QList<Subcluster> subclusters;
    for (int s = 0; s < shells.count(); s++) {
        QList<int> inds_s = shells[s].inds;
        Mda32 clips_s = get_subclips(clips, inds_s);
        //Mda subtemplate=compute_geometric_median_template(clips_s);
        Mda32 subtemplate = compute_mean_clip(clips_s);
        double ip = compute_template_ip(noise_shape, subtemplate);
        Mda32 subtemplate_resid;
        subtemplate_resid.allocate(M, T);
        for (int t = 0; t < T; t++) {
            for (int m = 0; m < M; m++) {
                subtemplate_resid.set(subtemplate.get(m, t) - ip * noise_shape.get(m, t), m, t);
            }
        }
        double subtemplate_resid_norm = compute_template_norm(subtemplate_resid);
        Subcluster SC;
        SC.inds = inds_s;
        SC.detectability_score = subtemplate_resid_norm; // here's the score!
        SC.peak = subtemplate.get(channel, Tmid);
        subclusters << SC;
    }
    return subclusters;
}

double compute_slope(const QVector<double>& X, const QVector<double>& Y)
{
    if (X.count() == 0)
        return 0;
    double mu_x = 0, mu_y = 0;
    for (int i = 0; i < X.count(); i++) {
        mu_x += X.value(i);
        mu_y += Y.value(i);
    }
    mu_x /= X.count();
    mu_y /= Y.count();
    double numer = 0, denom = 0;
    for (int i = 0; i < X.count(); i++) {
        numer += (X.value(i) - mu_x) * (Y.value(i) - mu_y);
        denom += (X.value(i) - mu_x) * (X.value(i) - mu_x);
    }
    if (denom)
        return numer / denom;
    else
        return 0;
}
