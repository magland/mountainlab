/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
** Created: 7/4/2016
*******************************************************/

#ifndef MDACONVERT_H
#define MDACONVERT_H

#include <QDebug>
#include "mlcommon.h"

struct mdaconvert_opts {
    QString input_path;
    QString input_dtype; // uint16, float32, ...
    QString input_format; // mda, raw, ...
    int input_num_header_rows = 0; //for csv
    int input_num_header_cols = 0; //for csv
    int num_channels = 0; //for nrd

    QString output_path;
    QString output_dtype; // uint16, float32, ...
    QString output_format; // mda, raw, ...

    QList<bigint> dims;

    bool check_input_file_size = true;
};
bool mdaconvert(const mdaconvert_opts& opts);

#endif // MDACONVERT_H
