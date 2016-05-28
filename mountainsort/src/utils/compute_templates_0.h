/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
** Created: 3/29/2016
*******************************************************/

#ifndef COMPUTE_TEMPLATES_0_H
#define COMPUTE_TEMPLATES_0_H

#include "diskreadmda.h"

Mda compute_templates_0(DiskReadMda& X, Mda& firings, int clip_size);
Mda compute_templates_0(DiskReadMda& X, const QList<double>& times, const QList<int>& labels, int clip_size);
void compute_templates_stdevs(Mda& ret_templates, Mda& ret_stdevs, DiskReadMda& X, const QList<double>& times, const QList<int>& labels, int clip_size);

#endif // COMPUTE_TEMPLATES_0_H
