/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
** Created: 2/22/2017
*******************************************************/
#ifndef P_EXTRACT_NEIGHBORHOOD_TIMESERIES_H
#define P_EXTRACT_NEIGHBORHOOD_TIMESERIES_H

#include <QString>
#include "mlcommon.h"

bool p_extract_neighborhood_timeseries(QString timeseries, QString timeseries_out, QList<int> channels);
bool p_extract_geom_channels(QString geom, QString geom_out, QList<int> channels);

#endif // P_EXTRACT_NEIGHBORHOOD_TIMESERIES_H
