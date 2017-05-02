/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
** Created: 4/5/2016
*******************************************************/

#ifndef CACHEMANAGER_H
#define CACHEMANAGER_H

#include <QString>
#include <QVariant>

class CacheManagerPrivate;
class CacheManager {
public:
    enum Duration {
        ShortTerm,
        LongTerm
    };

    friend class CacheManagerPrivate;
    CacheManager();
    virtual ~CacheManager();

    void setLocalBasePath(const QString& path);
    void setIntermediateFileFolder(const QString& folder);
    QString makeRemoteFile(const QString& mlproxy_url, const QString& file_name = "", Duration duration = ShortTerm);
    QString makeLocalFile(const QString& file_name = "", Duration duration = ShortTerm);
    QString makeIntermediateFile(const QString& file_name = "");
    QString localTempPath();

    void setTemporaryFileDuration(QString path, qint64 duration_sec);
    void setTemporaryFileExpirePid(QString path, qint64 pid);
    QString makeExpiringFile(QString file_name, qint64 duration_sec);
    void removeExpiredFiles();

    void cleanUp();

    static CacheManager* globalInstance();

    //private slots:
    //    void slot_remove_on_delete();

private:
    CacheManagerPrivate* d;
};

#endif // CACHEMANAGER_H
