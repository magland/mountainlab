/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
** Created: 5/11/2016
*******************************************************/

#include "mlutils.h"

#include <QCryptographicHash>
#include <QFileInfo>

QString cfp(const QString &path)
{
    /// Witold. How can I return a can a canonical file path when the file does not yet exist? This can be important for consistency!
    //for now I just return the path
    return path;
    /*
    if (QFile::exists(path)) {
        return QFileInfo(path).canonicalFilePath();
    }
    else {
        /// Witold. How can I return a can a canonical file path when the file does not yet exist?
        return path;
    }
    */
}


QString compute_checksum_of_file(const QString &path)
{
    QFile file(path);
    if (!file.open(QIODevice::ReadOnly))
        return "";
    QCryptographicHash hash(QCryptographicHash::Sha1);
    hash.addData(&file);
    file.close();
    QString ret = QString(hash.result().toHex());
    return ret;
}
