/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
** Created: 5/11/2016
*******************************************************/

#ifndef MLUTILS_H
#define MLUTILS_H

#include <QString>
#include <QDebug>

QString cfp(const QString &path); //cannonical file path
QString compute_checksum_of_file(const QString &path);

#endif // MLUTILS_H

