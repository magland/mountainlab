#include "branch_cluster_v2.h"
#include "diskreadmda.h"
#include <stdio.h>
#include "get_pca_features.h"
#include <math.h>
#include "isosplit2.h"
#include <QDebug>
#include <QTime>
#include "extract_clips.h"
#include "msmisc.h"
#include "compute_templates_0.h"
#include "get_sort_indices.h"

QList<int> do_branch_cluster_v2(Mda& clips, const Branch_Cluster_V2_Opts& opts, int channel_for_display);
QList<double> compute_peaks_v2(Mda& clips, int ch);
QList<int> consolidate_labels_v2(DiskReadMda& X, const QList<double>& times, const QList<int>& labels, int ch, int clip_size, int detect_interval, double consolidation_factor);
QList<long> get_sort_indices(const QList<int>& channels, const QList<double>& template_peaks);

bool branch_cluster_v2(const QString& timeseries_path, const QString& detect_path, const QString& adjacency_matrix_path, const QString& output_firings_path, const Branch_Cluster_V2_Opts& opts)
{
    /// TODO display total number of clusters at the end

    printf("Starting branch_cluster_v2\n");
    DiskReadMda X;
    X.setPath(timeseries_path);
    int M = X.N1();

    DiskReadMda detect;
    detect.setPath(detect_path);
    int L = detect.N2();

    Mda AM;
    if (!adjacency_matrix_path.isEmpty()) {
        AM.read(adjacency_matrix_path);
    } else {
        AM.allocate(M, M);
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < M; j++) {
                AM.set(1, i, j);
            }
        }
    }

    if ((AM.N1() != M) || (AM.N2() != M)) {
        printf("Error: incompatible dimensions between AM and X.\n");
        return false;
    }

    Mda firings0;
    firings0.allocate(5, L); //L is the max it could be

    int jjjj = 0;
    int k_offset = 0;
#pragma omp parallel for
    for (int m = 0; m < M; m++) {
        Mda clips;
        QList<double> times;
#pragma omp critical
        {
            QList<int> neighborhood;
            neighborhood << m;
            for (int a = 0; a < M; a++)
                if ((AM.value(m, a)) && (a != m))
                    neighborhood << a;
            for (int i = 0; i < L; i++) {
                if (detect.value(0, i) == (m + 1)) {
                    times << detect.value(1, i) - 1; //convert to 0-based indexing
                }
            }
            clips = extract_clips(X, times, neighborhood, opts.clip_size);
        }
        QList<int> labels = do_branch_cluster_v2(clips, opts, m);
#pragma omp critical
        {
            labels = consolidate_labels_v2(X, times, labels, m, opts.clip_size, opts.detect_interval, opts.consolidation_factor);
            QList<double> peaks = compute_peaks_v2(clips, 0);

            for (int i = 0; i < times.count(); i++) {
                if (labels[i]) {
                    firings0.setValue(m + 1, 0, jjjj); //channel
                    firings0.setValue(times[i] + 1, 1, jjjj); //times //convert back to 1-based indexing
                    firings0.setValue(labels[i] + k_offset, 2, jjjj); //labels
                    firings0.setValue(peaks[i], 3, jjjj); //peaks
                    jjjj++;
                }
            }
            k_offset += compute_max(labels);
        }
    }

    int L_true = jjjj;
    Mda firings;
    firings.allocate(firings0.N1(), L_true);
    for (int i = 0; i < L_true; i++) {
        for (int j = 0; j < firings0.N1(); j++) {
            firings.setValue(firings0.value(j, i), j, i);
        }
    }

    //Now reorder the labels
    {
        printf("Reordering labels...\n");
        QList<int> labels;
        for (long i = 0; i < L; i++) {
            int k = (int)firings.value(2, i);
            labels << k;
        }
        int K = compute_max(labels);
        QList<int> channels;
        for (int k = 0; k < K; k++)
            channels << 0;
        for (long i = 0; i < L; i++) {
            int k = (int)firings.value(2, i);
            if (k >= 1) {
                channels[k - 1] = (int)firings.value(0, i);
            }
        }
        int T_for_peaks = 3;
        int Tmid_for_peaks = (int)((T_for_peaks + 1) / 2) - 1;
        Mda templates = compute_templates_0(X, firings, T_for_peaks); //MxTxK
        QList<double> template_peaks;
        for (int k = 0; k < K; k++) {
            if (channels[k] >= 1) {
                template_peaks << templates.value(channels[k] - 1, Tmid_for_peaks, k);
            } else {
                template_peaks << 0;
            }
        }
        QList<long> sort_inds = get_sort_indices(channels, template_peaks);
        QList<long> label_map;
        for (int k = 0; k <= K; k++)
            label_map << 0;
        for (int j = 0; j < sort_inds.count(); j++)
            label_map[sort_inds[j] + 1] = j + 1;
        for (long i = 0; i < L; i++) {
            int k = (int)firings.value(2, i);
            if (k >= 1) {
                k = label_map[k];
                firings.setValue(k, 2, i);
            }
        }
    }

    firings.write64(output_firings_path);

    return true;
}

struct template_comparer_struct {
    int channel;
    double template_peak;
    long index;
};
struct template_comparer {
    bool operator()(const template_comparer_struct& a, const template_comparer_struct& b) const
    {
        if (a.channel < b.channel)
            return true;
        else if (a.channel == b.channel) {
            if (a.template_peak < b.template_peak)
                return true;
            else if (a.template_peak == b.template_peak)
                return (a.index < b.index);
            else
                return false;
        } else
            return false;
    }
};

QList<long> get_sort_indices(const QList<int>& channels, const QList<double>& template_peaks)
{
    QList<template_comparer_struct> list;
    for (long i = 0; i < channels.count(); i++) {
        template_comparer_struct tmp;
        tmp.channel = channels[i];
        tmp.template_peak = template_peaks[i];
        tmp.index = i;
        list << tmp;
    }
    qSort(list.begin(), list.end(), template_comparer());
    QList<long> ret;
    for (int i = 0; i < list.count(); i++) {
        ret << list[i].index;
    }
    return ret;
}

QList<int> consolidate_labels_v2(DiskReadMda& X, const QList<double>& times, const QList<int>& labels, int ch, int clip_size, int detect_interval, double consolidation_factor)
{
    printf("Consolidation factor = %g\n",consolidation_factor);
    int M = X.N1();
    int T = clip_size;
    int K = compute_max(labels);
    int Tmid = (int)((T + 1) / 2) - 1;
    QList<int> all_channels;
    for (int m = 0; m < M; m++)
        all_channels << m;
    int label_mapping[K + 1];
    label_mapping[0] = 0;
    int kk = 1;
    for (int k = 1; k <= K; k++) {
        QList<double> times_k;
        for (int i = 0; i < times.count(); i++) {
            if (labels[i] == k)
                times_k << times[i];
        }
        Mda clips_k = extract_clips(X, times_k, all_channels, clip_size);
        Mda template_k = compute_mean_clip(clips_k);
        QList<double> energies;
        for (int m = 0; m < M; m++)
            energies << 0;
        double max_energy=0;
        for (int t = 0; t < T; t++) {
            for (int m = 0; m < M; m++) {
                double val = template_k.value(m, t);
                energies[m] += val * val;
                if ((m!=ch)&&(energies[m]>max_energy)) max_energy=energies[m];
            }
        }
        //double max_energy = compute_max(energies);
        bool okay = true;
        if (energies[ch] < max_energy * consolidation_factor)
            okay = false;
        double abs_peak_val = 0;
        int abs_peak_ind = 0;
        for (int t = 0; t < T; t++) {
            double value = template_k.value(ch, t);
            if (fabs(value) > abs_peak_val) {
                abs_peak_val = fabs(value);
                abs_peak_ind = t;
            }
        }
        if (fabs(abs_peak_ind - Tmid) > detect_interval) {
            okay = false;
        }
        if (okay) {
            label_mapping[k] = kk;
            kk++;
        } else
            label_mapping[k] = 0;
    }
    QList<int> ret;
    for (int i = 0; i < labels.count(); i++) {
        ret << label_mapping[labels[i]];
    }
    printf("Channel %d: Using %d of %d clusters.\n", ch + 1, compute_max(ret), K);
    return ret;
}

QList<double> compute_peaks_v2(Mda& clips, int ch)
{
    int T = clips.N2();
    int L = clips.N3();
    int t0 = (T + 1) / 2 - 1;
    QList<double> ret;
    for (int i = 0; i < L; i++) {
        ret << clips.value(ch, t0, i);
    }
    return ret;
}

QList<double> compute_abs_peaks_v2(Mda& clips, int ch)
{
    int T = clips.N2();
    int L = clips.N3();
    int t0 = (T + 1) / 2 - 1;
    QList<double> ret;
    for (int i = 0; i < L; i++) {
        ret << fabs(clips.value(ch, t0, i));
    }
    return ret;
}

QList<int> find_peaks_below_threshold_v2(QList<double>& peaks, double threshold)
{
    QList<int> ret;
    for (int i = 0; i < peaks.count(); i++) {
        if (peaks[i] < threshold)
            ret << i;
    }
    return ret;
}

QList<int> find_peaks_above_threshold_v2(QList<double>& peaks, double threshold)
{
    QList<int> ret;
    for (int i = 0; i < peaks.count(); i++) {
        if (peaks[i] >= threshold)
            ret << i;
    }
    return ret;
}

void normalize_features_v2(Mda& F)
{
    int M = F.N1();
    int N = F.N2();
    QList<double> norms;
    int aa = 0;
    for (int i = 0; i < N; i++) {
        double sumsqr = 0;
        for (int j = 0; j < M; j++) {
            double val = F.get(aa);
            sumsqr += val * val;
            aa++;
        }
        norms << sqrt(sumsqr);
    }
    aa = 0;
    for (int i = 0; i < N; i++) {
        double factor = 1;
        if (norms[i])
            factor = 1 / norms[i];
        for (int j = 0; j < M; j++) {
            F.set(F.get(aa) * factor, aa);
            aa++;
        }
    }
}

QList<int> do_cluster_with_normalized_features(Mda& clips, const Branch_Cluster_V2_Opts& opts)
{
    int M = clips.N1();
    int T = clips.N2();
    int L = clips.N3();
    int nF = opts.num_features;
    Mda FF;
    FF.allocate(nF, L);
    get_pca_features(M * T, L, nF, FF.dataPtr(), clips.dataPtr());
    normalize_features_v2(FF);
    return isosplit2(FF);
}

QList<int> do_cluster_without_normalized_features(Mda& clips, const Branch_Cluster_V2_Opts& opts)
{
    QTime timer;
    timer.start();
    int M = clips.N1();
    int T = clips.N2();
    int L = clips.N3();
    int nF = opts.num_features;
    Mda FF;
    FF.allocate(nF, L);
    get_pca_features(M * T, L, nF, FF.dataPtr(), clips.dataPtr(), opts.num_pca_representatives);
    //normalize_features(FF);
    QList<int> ret = isosplit2(FF);
    return ret;
}

QList<double> compute_dists_from_template(Mda& clips, Mda& template0)
{
    int M = clips.N1();
    int T = clips.N2();
    int L = clips.N3();
    double* ptr1 = clips.dataPtr();
    double* ptr2 = template0.dataPtr();
    QList<double> ret;
    int aaa = 0;
    for (int i = 0; i < L; i++) {
        int bbb = 0;
        double sumsqr = 0;
        for (int t = 0; t < T; t++) {
            for (int m = 0; m < M; m++) {
                double diff0 = ptr1[aaa] - ptr2[bbb];
                sumsqr += diff0 * diff0;
                aaa++;
                bbb++;
            }
        }
        ret << sqrt(sumsqr);
    }
    return ret;
}

QList<int> do_branch_cluster_v2(Mda& clips, const Branch_Cluster_V2_Opts& opts, int channel_for_display)
{
    printf("do_branch_cluster_v2 %ldx%ldx%ld (channel %d)\n",clips.N1(),clips.N2(),clips.N3(),channel_for_display+1);
    int M = clips.N1();
    int T = clips.N2();
    int L = clips.N3();
    QList<double> peaks = compute_peaks_v2(clips, 0);
    QList<double> abs_peaks = compute_abs_peaks_v2(clips, 0);

    //In the case we have both positive and negative peaks, just split into two tasks!
    double min0 = compute_min(peaks);
    double max0 = compute_max(peaks);
    if ((min0 < 0) && (max0 >= 0)) {
        //find the event inds corresponding to negative and positive peaks
        QList<int> inds_neg, inds_pos;
        for (int i = 0; i < L; i++) {
            if (peaks[i] < 0)
                inds_neg << i;
            else
                inds_pos << i;
        }

        //grab the negative and positive clips
        Mda clips_neg = grab_clips_subset(clips, inds_neg);
        Mda clips_pos = grab_clips_subset(clips, inds_pos);

        //cluster the negatives and positives separately
        printf("Channel %d: NEGATIVES (%d)\n", channel_for_display + 1, inds_neg.count());
        QList<int> labels_neg = do_branch_cluster_v2(clips_neg, opts, channel_for_display);
        printf("Channel %d: POSITIVES (%d)\n", channel_for_display + 1, inds_pos.count());
        QList<int> labels_pos = do_branch_cluster_v2(clips_pos, opts, channel_for_display);

        //Combine them together
        int K_neg = compute_max(labels_neg);
        QList<int> labels;
        for (int i = 0; i < L; i++)
            labels << 0;
        for (int i = 0; i < inds_neg.count(); i++) {
            labels[inds_neg[i]] = labels_neg[i];
        }
        for (int i = 0; i < inds_pos.count(); i++) {
            if (labels_pos[i])
                labels[inds_pos[i]] = labels_pos[i] + K_neg;
            else
                labels[inds_pos[i]] = 0;
        }
        return labels;
    }

    //First we simply cluster all of the events
    //QList<int> labels0=do_cluster_with_normalized_features(clips,opts);
    QTime timer;
    timer.start();
    QList<int> labels0 = do_cluster_without_normalized_features(clips, opts);
    int K0 = compute_max(labels0);

    if (K0 > 1) {
        //if we found more than one cluster, then we should divide and conquer
        //we apply the same procedure to each cluster and then combine all of the clusters together.
        printf("Channel %d: K=%d\n", channel_for_display + 1, K0);
        QList<int> labels;
        for (int i = 0; i < L; i++)
            labels << 0;
        int kk_offset = 0;
        for (int k = 1; k <= K0; k++) {
            QList<int> inds_k;
            for (int a = 0; a < L; a++) {
                if (labels0[a] == k)
                    inds_k << a;
            }
            Mda clips_k = grab_clips_subset(clips, inds_k);
            QList<int> labels_k = do_branch_cluster_v2(clips_k, opts, channel_for_display);
            for (int a = 0; a < inds_k.count(); a++) {
                labels[inds_k[a]] = labels_k[a] + kk_offset;
            }
            kk_offset += compute_max(labels_k);
        }
        return labels;
    } else {
        //otherwise, we have only one cluster
        //so we need to increase the threshold to see if we can get things to split at higher amplitude
        double abs_peak_threshold = 0;
        double max_abs_peak = compute_max(abs_peaks);

        //increase abs_peak_threshold by opts.shell_increment until we have at least opts.min_shell_size below and above the threshold
        while (true) {
            QList<int> inds_below = find_peaks_below_threshold_v2(abs_peaks, abs_peak_threshold);
            if ((inds_below.count() >= opts.min_shell_size) && (L - inds_below.count() >= opts.min_shell_size)) {
                break;
            }
            if (abs_peak_threshold > max_abs_peak) {
                break;
            }
            abs_peak_threshold += opts.shell_increment;
        }
        if (abs_peak_threshold > max_abs_peak) {
            //we couldn't split it. So fine, we'll just say there is only one cluster
            QList<int> labels;
            for (int i = 0; i < L; i++)
                labels << 1;
            return labels;
        } else {
            //we now split things into two categories based on abs_peak_threshold
            QList<int> inds_below = find_peaks_below_threshold_v2(abs_peaks, abs_peak_threshold);
            QList<int> inds_above = find_peaks_above_threshold_v2(abs_peaks, abs_peak_threshold);
            Mda clips_above = grab_clips_subset(clips, inds_above);

            //Apply the procedure to the events above the threshold
            QList<int> labels_above = do_branch_cluster_v2(clips_above, opts, channel_for_display);
            int K_above = compute_max(labels_above);

            if (K_above <= 1) {
                //there is really only one cluster
                QList<int> labels;
                for (int i = 0; i < L; i++)
                    labels << 1;
                return labels;
            } else {
                //there is more than one cluster. Let's divide up the based on the nearest
                //let's consider only the next shell above
                QList<double> abs_peaks_above;
                for (int i = 0; i < inds_above.count(); i++)
                    abs_peaks_above << abs_peaks[inds_above[i]];
                QList<int> inds_next_shell = find_peaks_below_threshold_v2(abs_peaks_above, abs_peak_threshold + opts.shell_increment);
                Mda clips_next_shell = grab_clips_subset(clips_above, inds_next_shell);
                QList<int> labels_next_shell;
                for (int i = 0; i < inds_next_shell.count(); i++)
                    labels_next_shell << labels_above[inds_next_shell[i]];

                //compute the centroids for the next shell above
                Mda centroids;
                centroids.allocate(M, T, K_above);
                for (int kk = 1; kk <= K_above; kk++) {
                    QList<int> inds_kk;
                    for (int i = 0; i < labels_next_shell.count(); i++) {
                        if (labels_next_shell[i] == kk)
                            inds_kk << i;
                    }
                    Mda clips_kk = grab_clips_subset(clips_next_shell, inds_kk);
                    Mda centroid0 = compute_mean_clip(clips_kk);
                    for (int t = 0; t < T; t++) {
                        for (int m = 0; m < M; m++) {
                            centroids.setValue(centroid0.value(m, t), m, t, kk - 1);
                        }
                    }
                }

                //set the labels for all of the inds above
                QList<int> labels;
                for (int i = 0; i < L; i++)
                    labels << 0;
                for (int i = 0; i < inds_above.count(); i++) {
                    labels[inds_above[i]] = labels_above[i];
                }

                //grab the clips for the events below and compute the distances to all the centroids of the next shell above
                Mda clips_below = grab_clips_subset(clips, inds_below);
                Mda distances;
                distances.allocate(inds_below.count(), K_above);
                for (int k = 1; k <= K_above; k++) {
                    QList<int> tmp;
                    tmp << k - 1;
                    Mda centroid0 = grab_clips_subset(centroids, tmp);
                    QList<double> dists = compute_dists_from_template(clips_below, centroid0);
                    for (int i = 0; i < inds_below.count(); i++) {
                        distances.setValue(dists[i], i, k - 1);
                    }
                }

                //label the events below based on distance to threshold
                for (int i = 0; i < inds_below.count(); i++) {
                    int best_k = 0;
                    double best_dist = distances.value(i, 0L);
                    for (int k = 0; k < K_above; k++) {
                        double dist0 = distances.value(i, k);
                        if (dist0 < best_dist) {
                            best_dist = dist0;
                            best_k = k;
                        }
                    }
                    labels[inds_below[i]] = best_k + 1; //convert back to 1-based indexing
                }
                return labels;
            }
        }
    }
}
