/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
*******************************************************/

#ifndef MVUTILS_H
#define MVUTILS_H

#include "diskarraymodel_new.h"
#include "diskreadmda.h"
#include <QList>
#include <QColor>

///A firing event defined by a time and a label
struct MVEvent {
    double time;
    int label;
};

///An epoch or time interval within a timeseries. t_begin and t_end are in timepoint units
struct Epoch {
    QString name;
    double t_begin;
    double t_end;
};
///Read a set of epochs from a text file (special format)
QList<Epoch> read_epochs(const QString& path);

///Utility
Mda compute_mean_waveform(DiskArrayModel_New* C);
///Utility
Mda compute_mean_stdev_waveform(DiskArrayModel_New* C);
///Utility
Mda compute_features(DiskArrayModel_New* C);
///Utility
Mda compute_features(const QList<DiskArrayModel_New*>& C);
///TODO: Low-priority Handle this properly
QColor get_heat_map_color(double val);

void user_save_image(const QImage& img);

#endif // MVUTILS_H
