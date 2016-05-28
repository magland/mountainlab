/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
** Created: 3/29/2016
*******************************************************/

#include "compute_templates_0.h"
#include "msmisc.h"
#include <math.h>

Mda compute_templates_0(DiskReadMda& X, Mda& firings, int clip_size)
{
    QList<double> times;
    QList<int> labels;
    long L = firings.N2();
    for (long i = 0; i < L; i++) {
        times << firings.value(1, i);
        labels << (int)firings.value(2, i);
    }
    return compute_templates_0(X, times, labels, clip_size);
}

Mda compute_templates_0(DiskReadMda& X, const QList<double>& times, const QList<int>& labels, int clip_size)
{
    int M = X.N1();
    int T = clip_size;
    long L = times.count();

    int K = compute_max(labels);

    int Tmid = (int)((T + 1) / 2) - 1;

    Mda templates(M, T, K);
    QList<long> counts;
    for (int k = 0; k < K; k++)
        counts << k;
    for (long i = 0; i < L; i++) {
        int k = labels[i];
        long t0 = (long)(times[i] + 0.5);
        if (k >= 1) {
            Mda X0;
            X.readChunk(X0, 0, t0 - Tmid, M, T);
            double* Xptr = X0.dataPtr();
            double* Tptr = templates.dataPtr(0, 0, k - 1);
            for (int i = 0; i < M * T; i++) {
                Tptr[i] += Xptr[i];
            }
            counts[k - 1]++;
        }
    }
    for (int k = 0; k < K; k++) {
        for (int t = 0; t < T; t++) {
            for (int m = 0; m < M; m++) {
                if (counts[k]) {
                    templates.set(templates.get(m, t, k) / counts[k], m, t, k);
                }
            }
        }
    }

    return templates;
}

void compute_templates_stdevs(Mda& templates, Mda& stdevs, DiskReadMda& X, const QList<double>& times, const QList<int>& labels, int clip_size)
{
    int M = X.N1();
    int T = clip_size;
    long L = times.count();

    int K = compute_max(labels);

    int Tmid = (int)((T + 1) / 2) - 1;

    Mda sums(M, T, K);
    Mda sumsqrs(M, T, K);
    QList<long> counts;
    for (int k = 0; k < K; k++)
        counts << k;
    for (long i = 0; i < L; i++) {
        int k = labels[i];
        long t0 = (long)(times[i] + 0.5);
        if (k >= 1) {
            Mda X0;
            X.readChunk(X0, 0, t0 - Tmid, M, T);
            double* Xptr = X0.dataPtr();
            double* sum_ptr = sums.dataPtr(0, 0, k - 1);
            double* sumsqr_ptr = sumsqrs.dataPtr(0, 0, k - 1);
            for (int i = 0; i < M * T; i++) {
                sum_ptr[i] += Xptr[i];
                sumsqr_ptr[i] += Xptr[i] * Xptr[i];
            }
            counts[k - 1]++;
        }
    }

    templates.allocate(M, T, K);
    stdevs.allocate(M, T, K);
    for (int k = 0; k < K; k++) {
        for (int t = 0; t < T; t++) {
            for (int m = 0; m < M; m++) {
                if (counts[k] >= 2) {
                    double sum0 = sums.get(m, t, k);
                    double sumsqr0 = sumsqrs.get(m, t, k);
                    templates.set(sum0 / counts[k], m, t, k);
                    stdevs.set(sqrt(sumsqr0 / counts[k] - (sum0 * sum0) / (counts[k] * counts[k])), m, t, k);
                }
            }
        }
    }
}
