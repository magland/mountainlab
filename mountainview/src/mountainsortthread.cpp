/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
** Created: 4/4/2016
*******************************************************/

#include "mountainsortthread.h"
#include <QCoreApplication>
#include <QMap>
#include <QProcess>
#include <QVariant>
#include <QDebug>
#include "msmisc.h"
#include "cachemanager.h"
#include "mlutils.h"

class MountainsortThreadPrivate {
public:
    MountainsortThread* q;
    QString m_processor_name;
    QMap<QString, QVariant> m_parameters;
    QString m_mscmdserver_url;

    QString get_remote_url_from_parameters();
    QString remote_url_of_path(const QString& path);
    QString create_temporary_output_file_name(const QString& remote_url, const QString& processor_name, const QMap<QString, QVariant>& params, const QString& parameter_name);
};

MountainsortThread::MountainsortThread()
{
    d = new MountainsortThreadPrivate;
    d->q = this;
}

MountainsortThread::~MountainsortThread()
{
    delete d;
}

void MountainsortThread::setProcessorName(const QString& pname)
{
    d->m_processor_name = pname;
    this->setStatus(pname);
}

QString MountainsortThread::makeOutputFilePath(const QString& pname)
{
    QString ret = d->create_temporary_output_file_name(d->get_remote_url_from_parameters(), d->m_processor_name, d->m_parameters, pname);
    d->m_parameters[pname] = ret;
    return ret;
}

void MountainsortThread::setInputParameters(const QMap<QString, QVariant>& parameters)
{
    d->m_parameters = parameters;
}

void MountainsortThread::setMscmdServerUrl(const QString& url)
{
    d->m_mscmdserver_url = url;
}

void MountainsortThread::compute()
{

    if (d->m_mscmdserver_url.isEmpty()) {
	QString mountainsort_exe = cfp(qApp->applicationDirPath() + "/../../mountainsort/bin/mountainsort");
        QStringList args;
        args << d->m_processor_name;
        QStringList keys = d->m_parameters.keys();
        foreach (QString key, keys) {
            args << QString("--%1=%2").arg(key).arg(d->m_parameters.value(key).toString());
        }
        this->setStatus("Local "+d->m_processor_name,"Executing locally: "+mountainsort_exe,0.5);
        if (QProcess::execute(mountainsort_exe, args) != 0) {
            qWarning() << "Problem running mountainsort" << mountainsort_exe << args;
        }
        this->setStatus("","Done.",1);
    }
    else {
        QString url = d->m_mscmdserver_url + "/?";
        url += "processor=" + d->m_processor_name + "&";
        QStringList keys = d->m_parameters.keys();
        foreach (QString key, keys) {
            url += QString("%1=%2&").arg(key).arg(d->m_parameters.value(key).toString());
        }
        this->setStatus("Remote "+d->m_processor_name,"http_get_text: "+url,0.5);
        http_get_text(url);
        this->setStatus("","",1);
    }
}

QString MountainsortThreadPrivate::create_temporary_output_file_name(const QString& remote_url, const QString& processor_name, const QMap<QString, QVariant>& params, const QString& parameter_name)
{
    QString str = processor_name + ":";
    QStringList keys = params.keys();
    qSort(keys);
    foreach (QString key, keys) {
        str += key + "=" + params.value(key).toString() + "&";
    }

    QString file_name = QString("%1_%2.tmp").arg(compute_hash(str)).arg(parameter_name);
    QString ret = CacheManager::globalInstance()->makeRemoteFile(remote_url, file_name, CacheManager::LongTerm);
    return ret;
}

QString MountainsortThreadPrivate::get_remote_url_from_parameters()
{
    QString ret;
    QStringList keys = m_parameters.keys();
    foreach (QString key, keys) {
        QString val = m_parameters.value(key).toString();
        if (val.startsWith("http://")) {
            QString remote_url = remote_url_of_path(val);
            if (!remote_url.isEmpty()) {
                if (ret.isEmpty()) {
                    ret = remote_url;
                }
            }
            else {
                qWarning() << QString("More than one remote name in parameters for %1. Using %2.").arg(m_processor_name).arg(ret);
                qWarning() << m_parameters;
            }
        }
    }
    return ret;
}

QString MountainsortThreadPrivate::remote_url_of_path(const QString& path)
{
    if (path.startsWith("http://")) {
        int ind = path.indexOf("/", QString("http://").count());
        if (ind < 0)
            return "";
        return path.mid(0, ind);
    }
    else {
        return "";
    }
}
