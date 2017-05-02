/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
** Created: 4/5/2016
*******************************************************/

#include "cachemanager.h"

#include <QDir>
#include <QDebug>
#include <QDateTime>
#include <QThread>
#include <QJsonDocument>
#include "mlcommon.h"
#include <signal.h>

#define DEFAULT_LOCAL_BASE_PATH MLUtil::tempPath()

class CacheManagerPrivate {
public:
    CacheManager* q;
    QString m_local_base_path;
    QString m_intermediate_file_folder;

    QString create_random_file_name();
};

CacheManager::CacheManager()
{
    d = new CacheManagerPrivate;
    d->q = this;
}

CacheManager::~CacheManager()
{
    delete d;
}

void CacheManager::setLocalBasePath(const QString& path)
{
    d->m_local_base_path = path;
    if (!QDir(path).exists()) {
        QString parent = QFileInfo(path).path();
        QString name = QFileInfo(path).fileName();
        if (!QDir(parent).mkdir(name)) {
            qWarning() << "Unable to create local base path" << path;
        }
    }
    if (!QDir(path + "/tmp_short_term").exists()) {
        QDir(path).mkdir("tmp_short_term");
    }
    if (!QDir(path + "/tmp_long_term").exists()) {
        QDir(path).mkdir("tmp_long_term");
    }
    QFile::Permissions perm = QFileDevice::ReadUser | QFileDevice::WriteUser | QFileDevice::ExeUser | QFileDevice::ReadGroup | QFileDevice::WriteGroup | QFileDevice::ExeGroup | QFileDevice::ReadOther | QFileDevice::WriteOther | QFileDevice::ExeOther;
    QFile::setPermissions(path, perm);
    QFile::setPermissions(path + "/tmp_short_term", perm);
    QFile::setPermissions(path + "/tmp_long_term", perm);
}

void CacheManager::setIntermediateFileFolder(const QString& folder)
{
    d->m_intermediate_file_folder = folder;
}

QString CacheManager::makeRemoteFile(const QString& mlproxy_url, const QString& file_name_in, CacheManager::Duration duration)
{
    if (mlproxy_url.isEmpty()) {
        return makeLocalFile(file_name_in, duration);
    }

    QString file_name = file_name_in;
    if (file_name.isEmpty())
        file_name = d->create_random_file_name();
    QString str;
    if (duration == ShortTerm)
        str = "tmp_short_term";
    else if (duration == LongTerm)
        str = "tmp_long_term";
    else {
        qWarning() << "Unexpected problem" << __FUNCTION__ << __FILE__ << __LINE__;
        return "";
    }
    return QString("%1/mdaserver/%2/%3").arg(mlproxy_url).arg(str).arg(file_name);
}

QString CacheManager::makeLocalFile(const QString& file_name_in, CacheManager::Duration duration)
{
    QString file_name = file_name_in;
    if (file_name.isEmpty())
        file_name = d->create_random_file_name();

    QString str;
    if (duration == ShortTerm)
        str = "tmp_short_term";
    else if (duration == LongTerm)
        str = "tmp_long_term";
    else {
        qWarning() << "Unexpected problem" << __FUNCTION__ << __FILE__ << __LINE__;
        return "";
    }
    QString ret = QString("%1/%2/%3").arg(localTempPath()).arg(str).arg(file_name);

    return ret;
}

QString CacheManager::makeIntermediateFile(const QString& file_name_in)
{
    QString file_name = file_name_in;
    if (file_name.isEmpty())
        file_name = d->create_random_file_name();

    QString dirname;
    if (d->m_intermediate_file_folder.isEmpty()) {
        dirname = localTempPath() + "/tmp_long_term";
    }
    else if (QFileInfo(d->m_intermediate_file_folder).isAbsolute())
        dirname = d->m_intermediate_file_folder;
    else {
        QDir(localTempPath()).mkdir("/tmp_long_term");
        dirname = localTempPath() + "/tmp_long_term/" + d->m_intermediate_file_folder;
    }
    QDir(QFileInfo(dirname).path()).mkdir(QFileInfo(dirname).fileName());

    QString ret = QString("%1/%2").arg(dirname).arg(file_name);

    return ret;
}

QString CacheManager::localTempPath()
{
    if (d->m_local_base_path.isEmpty()) {
        //qWarning() << "Local base path has not been set. Using default: " + QString(DEFAULT_LOCAL_BASE_PATH);
        this->setLocalBasePath(DEFAULT_LOCAL_BASE_PATH);
    }
    return d->m_local_base_path;
}

void CacheManager::setTemporaryFileDuration(QString path, qint64 duration_sec)
{
    //duration_sec = 10;

    QString str = MLUtil::computeSha1SumOfString(path);
    if (!QFile::exists(QString("%1/expiration_records").arg(localTempPath())))
        QDir(localTempPath()).mkdir("expiration_records");
    QString fname = QString("%1/expiration_records/%2.json").arg(CacheManager::localTempPath()).arg(str);
    QJsonObject obj;
    if (QFile::exists(fname)) {
        obj = QJsonDocument::fromJson(TextFile::read(fname).toUtf8()).object();
    }
    obj["path"] = path;
    if (QFile::exists(path)) {
        obj["creation_timestamp"] = QFileInfo(path).created().toString("yyyy-MM-dd hh:mm:ss");
    }
    else {
        obj["creation_timestamp"] = "";
    }
    QDateTime timestamp = QFileInfo(path).created().addSecs(duration_sec);
    obj["expiration_timestamp"] = timestamp.toString("yyyy-MM-dd hh:mm:ss");
    ;

    QString json = QJsonDocument(obj).toJson();
    TextFile::write(fname, json);
}

void CacheManager::setTemporaryFileExpirePid(QString path, qint64 pid)
{
    QString str = MLUtil::computeSha1SumOfString(path);
    if (!QFile::exists(QString("%1/expiration_records").arg(localTempPath())))
        QDir(localTempPath()).mkdir("expiration_records");
    QString fname = QString("%1/expiration_records/%2.json").arg(CacheManager::localTempPath()).arg(str);
    QJsonObject obj;
    if (QFile::exists(fname)) {
        obj = QJsonDocument::fromJson(TextFile::read(fname).toUtf8()).object();
    }
    obj["path"] = path;
    if (QFile::exists(path)) {
        obj["creation_timestamp"] = QFileInfo(path).created().toString("yyyy-MM-dd hh:mm:ss");
    }
    else {
        obj["creation_timestamp"] = "";
    }
    obj["expiration_pid"] = pid;

    QString json = QJsonDocument(obj).toJson();
    TextFile::write(fname, json);
}

QString CacheManager::makeExpiringFile(QString file_name, qint64 duration_sec)
{
    QString ret = this->makeLocalFile(file_name, ShortTerm);
    this->setTemporaryFileDuration(ret, duration_sec);
    return ret;
}

struct CMFileRec {
    QString path;
    int elapsed_sec;
    double size_gb;
};

QList<CMFileRec> get_file_records(const QString& path)
{
    QStringList fnames = QDir(path).entryList(QStringList("*"), QDir::Files, QDir::Name);
    QList<CMFileRec> records;
    foreach (QString fname, fnames) {
        CMFileRec rec;
        rec.path = path + "/" + fname;
        rec.elapsed_sec = QFileInfo(rec.path).lastModified().secsTo(QDateTime::currentDateTime());
        rec.size_gb = QFileInfo(rec.path).size() * 1.0 / 1e9;
        records << rec;
    }

    QStringList dirnames = QDir(path).entryList(QStringList("*"), QDir::Dirs | QDir::NoDotAndDotDot | QDir::NoSymLinks, QDir::Name);
    foreach (QString dirname, dirnames) {
        QList<CMFileRec> RR = get_file_records(path + "/" + dirname);
        records.append(RR);
    }
    return records;
}

struct CMFileRec_comparer {
    bool operator()(const CMFileRec& a, const CMFileRec& b) const
    {
        if (a.elapsed_sec > b.elapsed_sec)
            return true;
        else
            return false;
    }
};

void sort_by_elapsed(QList<CMFileRec>& records)
{
    qSort(records.begin(), records.end(), CMFileRec_comparer());
}

bool pid_exists(int pid)
{
    return (kill(pid, 0) == 0);
}

bool is_expired(CMFileRec rec)
{
    QString fname = QFileInfo(rec.path).fileName();
    QStringList list = fname.split("--");
    if (list.count() < 2)
        return false;
    QStringList list2 = list.value(0).split(".");
    foreach (QString str, list2) {
        QStringList vals = str.split("=");
        if (vals.count() == 2) {
            if (vals.value(0) == "dursec") {
                bigint numsec = vals.value(1).toLongLong();
                if (numsec < rec.elapsed_sec)
                    return true;
            }
            else if (vals.value(0) == "pid") {
                int pid = vals.value(1).toInt();
                if (!pid_exists(pid))
                    return true;
            }
        }
    }
    return false;
}

void CacheManager::removeExpiredFiles()
{
    //get a list of the expiration records
    QString dirname = QString("%1/expiration_records").arg(CacheManager::localTempPath());
    QStringList list = QDir(dirname).entryList(QStringList("*.json"), QDir::Files, QDir::Name);
    foreach (QString str, list) {
        QString fname = QString("%1/%2").arg(dirname).arg(str);
        //make sure it was created at least a few seconds ago
        if (QFileInfo(fname).lastModified().secsTo(QDateTime::currentDateTime()) > 3) {
            QString json = TextFile::read(fname);
            QJsonObject obj = QJsonDocument::fromJson(json.toUtf8()).object();
            QString path0 = obj["path"].toString();
            QString creation_timestamp_str = obj["creation_timestamp"].toString();
            //QDateTime creation_timestamp=QDateTime::fromString(creation_timestamp_str,"yyyy-MM-dd hh:mm:ss");
            QString expiration_timestamp_str = obj["expiration_timestamp"].toString();
            QDateTime expiration_timestamp = QDateTime::fromString(expiration_timestamp_str, "yyyy-MM-dd hh:mm:ss");
            int expiration_pid = obj["expiration_pid"].toInt();
            if (QFile(path0).exists()) {
                //the temporary file exists
                bool ok = true;
                if (!creation_timestamp_str.isEmpty()) {
                    //there was a creation timestamp. let's see if it matches.
                    if (QFileInfo(path0).created().toString("yyyy-MM-dd hh:mm:ss") != creation_timestamp_str) {
                        qWarning() << "Temporary file creation timestamp does not match that in the expiration record. Removing expiration record.";
                        QFile::remove(fname);
                        ok = false;
                    }
                }
                if ((ok) && (!expiration_timestamp_str.isEmpty())) {
                    //an expiration time exists
                    if (expiration_timestamp.secsTo(QDateTime::currentDateTime()) > 0) {
                        //the file is expired
                        if (QFile::remove(path0)) {
                            QFile::remove(fname);
                        }
                        else {
                            qWarning() << "Unable to remove expired temporary file: " + path0;
                        }
                    }
                }
                if ((ok) && (expiration_pid)) {
                    if (!pid_exists(expiration_pid)) {
                        //the pid no longer exists
                        if (QFile::remove(path0)) {
                            QFile::remove(fname);
                        }
                        else {
                            qWarning() << "Unable to remove expired temporary file: " + path0;
                        }
                    }
                }
            }
            else {
                //the temporary file does not exist
                if (QFileInfo(fname).lastModified().secsTo(QDateTime::currentDateTime()) > 10) {
                    //the expiration record has been around for at least 10 seconds
                    QFile::remove(fname);
                }
            }
        }
    }

    /*
    QList<CMFileRec> recs=get_file_records(this->localTempPath()+"/expiring");
    for (int i=0; i<recs.count(); i++) {
        CMFileRec rec=recs[i];
        if (is_expired(rec)) {
            QFile::remove(rec.path);
        }
    }
    */
}

void CacheManager::cleanUp()
{
    double max_gb = MLUtil::configValue("general", "max_cache_size_gb").toDouble();
    if (!max_gb) {
        qWarning() << "max_gb is zero. You probably need to adjust the mountainlab configuration files.";
        return;
    }
    if (!d->m_local_base_path.endsWith("/mountainlab")) {
        qWarning() << "For safety, the temporary path must end with /mountainlab";
        return;
    }
    QList<CMFileRec> records = get_file_records(d->m_local_base_path);
    double total_size_gb = 0;
    for (int i = 0; i < records.count(); i++) {
        total_size_gb += records[i].size_gb;
    }
    int num_files_removed = 0;
    double amount_removed = 0;
    if (total_size_gb > max_gb) {
        double amount_to_remove = total_size_gb - 0.75 * max_gb; //let's get it down to 75% of the max allowed
        sort_by_elapsed(records);
        for (int i = 0; i < records.count(); i++) {
            if (amount_removed >= amount_to_remove) {
                break;
            }
            if (!QFile::remove(records[i].path)) {
                qWarning() << "Unable to remove file while cleaning up cache: " + records[i].path;
                return;
            }
            amount_removed += records[i].size_gb;
            num_files_removed++;
        }
    }
    if (num_files_removed) {
        qWarning() << QString("CacheManager removed %1 GB and %2 files").arg(amount_removed).arg(num_files_removed);
    }
}

Q_GLOBAL_STATIC(CacheManager, theInstance)
CacheManager* CacheManager::globalInstance()
{
    return theInstance;
}

/*
void CacheManager::slot_remove_on_delete()
{
    QString fname=sender()->property("CacheManager_file_to_remove").toString();
    if (!fname.isEmpty()) {
        if (!QFile::remove(fname)) {
            qWarning() << "Unable to remove local cached file:" << fname;
        }
    }
    else {
        qWarning() << "Unexpected problem" << __FUNCTION__ << __FILE__ << __LINE__;
    }
}
*/

QString CacheManagerPrivate::create_random_file_name()
{

    int num1 = QDateTime::currentMSecsSinceEpoch();
    long num2 = (long)QThread::currentThreadId();
    int num3 = qrand();
    return QString("ms.%1.%2.%3.tmp").arg(num1).arg(num2).arg(num3);
}
