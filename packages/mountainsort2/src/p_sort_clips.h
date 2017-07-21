#ifndef P_SORT_CLIPS_H
#define P_SORT_CLIPS_H

#include <QString>
#include "mlcommon.h"

struct Sort_clips_opts {
    int num_features = 8;
    double isocut_threshold = 1;
    double K_init = 200;
    bigint max_samples = 10000; //for subsampled pca
    bool remove_outliers = false;
    int weighted_pca = 0;
};

bool p_sort_clips(QString clips, QString firings_out, Sort_clips_opts opts);
bool p_reorder_labels(QString templates, QString firings, QString firings_out);

#endif // P_SORT_CLIPS_H
