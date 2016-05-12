/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
** Created: 4/7/2016
*******************************************************/

#ifndef COMPUTE_DETECTABILITY_SCORES_H
#define COMPUTE_DETECTABILITY_SCORES_H

#include <QList>
#include <QString>
#include <diskreadmda.h>
#include <mda.h>

struct compute_detectability_scores_opts {
    int clip_size;
    double shell_increment;
    int min_shell_size;
};

bool compute_detectability_scores(QString timeseries_path, QString firings_path, QString firings_out_path, const compute_detectability_scores_opts& opts);

struct Subcluster {
    QList<long> inds;
    double detectability_score;
    double peak;
    int label;
};

struct Shell {
    QList<long> inds;
};
struct Define_Shells_Opts {
    double shell_increment;
    int min_shell_size;
};

QList<long> find_label_inds(const QList<long>& labels, int k);
Mda get_subclips(Mda& clips, const QList<long>& inds);
QList<Shell> define_shells(const QList<double>& peaks, const Define_Shells_Opts& opts);
QList<double> randsample_with_replacement(long N, long K);
Mda estimate_noise_shape(DiskReadMda& X, int T, int ch);
Mda compute_features(Mda& clips, int num_features);
void compute_geometric_median(int M, int N, double* output, double* input, int num_iterations = 10);
Mda compute_geometric_median_template(Mda& clips);
double compute_template_ip(Mda& T1, Mda& T2);
double compute_template_norm(Mda& T);
QList<Subcluster> compute_subcluster_detectability_scores(Mda& noise_shape, Mda& clips, int channel, const Define_Shells_Opts& opts);
double compute_slope(const QList<double>& X, const QList<double>& Y);

#endif // COMPUTE_DETECTABILITY_SCORES_H
