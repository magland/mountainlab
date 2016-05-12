/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
** Created: 4/27/2016
*******************************************************/

#include "processmanager.h"
#include <QDir>
#include <QDebug>
#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonArray>
#include <QProcess>
#include <QTime>
#include <QEventLoop>
#include <QCryptographicHash>
#include "mpdaemon.h"
#include "textfile.h"
#include <QCoreApplication>
#include "mpdaemon.h"

struct PMProcess {
    MLProcessInfo info;
    QProcess* qprocess;
};

class ProcessManagerPrivate {
public:
    ProcessManager* q;

    QMap<QString, MLProcessor> m_processors;
    QMap<QString, PMProcess> m_processes;

    void clear_all_processes();
    void update_process_info(QString id);

    static MLProcessor create_processor_from_json_object(QJsonObject obj);
    static MLParameter create_parameter_from_json_object(QJsonObject obj);
    static QJsonObject compute_unique_process_object(MLProcessor P, const QVariantMap& parameters);
    static bool all_input_and_output_files_exist(MLProcessor P, const QVariantMap& parameters);
    static QString compute_unique_object_code(QJsonObject obj);
    static QJsonObject create_file_object(const QString& fname);
};

ProcessManager::ProcessManager()
{
    d = new ProcessManagerPrivate;
    d->q = this;
}

ProcessManager::~ProcessManager()
{
    d->clear_all_processes();
    delete d;
}

bool ProcessManager::loadProcessors(const QString& path, bool recursive)
{
    QStringList fnames = QDir(path).entryList(QStringList("*.mp"), QDir::Files, QDir::Name);
    foreach(QString fname, fnames)
    {
        if (!this->loadProcessorFile(path + "/" + fname))
            return false;
    }
    if (recursive) {
        QStringList subdirs = QDir(path).entryList(QStringList("*"), QDir::Dirs | QDir::NoDotAndDotDot, QDir::Name);
        foreach(QString subdir, subdirs)
        {
            if (!this->loadProcessors(path + "/" + subdir))
                return false;
        }
    }
    return true;
}

bool ProcessManager::loadProcessorFile(const QString& path)
{
    QString json;
    if (QFileInfo(path).isExecutable()) {
        QProcess pp;
        pp.start(path, QStringList("spec"));
        if (!pp.waitForFinished()) {
            qWarning() << "Problem with executable processor file, waiting for finish: " + path;
            return false;
        }
        pp.waitForReadyRead();
        QString output = pp.readAll();
        json = output;
        if (json.isEmpty()) {
            qWarning() << "Executable processor file did not return output for spec: " + path;
            return false;
        }
    } else {
        json = read_text_file(path);
        if (json.isEmpty()) {
            qWarning() << "Processor file is empty: " + path;
            return false;
        }
    }
    QJsonParseError error;
    QJsonObject obj = QJsonDocument::fromJson(json.toLatin1(), &error).object();
    if (error.error != QJsonParseError::NoError) {
        qWarning() << "Json parse error: " << error.errorString();
        return false;
    }
    if (!obj["processors"].isArray()) {
        qWarning() << "Problem with processor file: processors field is missing or not any array: " + path;
        return false;
    }
    QJsonArray processors = obj["processors"].toArray();
    for (int i = 0; i < processors.count(); i++) {
        if (!processors[i].isObject()) {
            qWarning() << "Problem with processor file: processor is not an object: " + path;
            return false;
        }
        MLProcessor P = d->create_processor_from_json_object(processors[i].toObject());
        P.basepath = QFileInfo(path).path();
        if (P.name.isEmpty()) {
            qWarning() << "Problem with processor file: processor error: " + path;
            return false;
        }
        d->m_processors[P.name] = P;
    }
    return true;
}

QStringList ProcessManager::processorNames() const
{
    return d->m_processors.keys();
}

MLProcessor ProcessManager::processor(const QString& name)
{
    return d->m_processors.value(name);
}

QString ProcessManager::startProcess(const QString& processor_name, const QVariantMap& parameters)
{
    if (!this->checkParameters(processor_name, parameters))
        return "";

    if (!d->m_processors.contains(processor_name)) {
        qWarning() << "Unable to find processor: " + processor_name;
        return "";
    }
    MLProcessor P = d->m_processors[processor_name];
    QString exe_command = P.exe_command;
    exe_command.replace(QRegExp("\\$\\(basepath\\)"), P.basepath);
    {
        QString ppp;
        {
            QStringList keys = P.inputs.keys();
            foreach(QString key, keys)
            {
                exe_command.replace(QRegExp(QString("\\$%1").arg(key)), parameters[key].toString());
                ppp += QString("--%1=%2 ").arg(key).arg(parameters[key].toString());
            }
        }
        {
            QStringList keys = P.outputs.keys();
            foreach(QString key, keys)
            {
                exe_command.replace(QRegExp(QString("\\$%1").arg(key)), parameters[key].toString());
                ppp += QString("--%1=%2 ").arg(key).arg(parameters[key].toString());
            }
        }
        {
            QStringList keys = P.parameters.keys();
            foreach(QString key, keys)
            {
                exe_command.replace(QRegExp(QString("$%1").arg(key)), parameters[key].toString());
                ppp += QString("--%1=%2 ").arg(key).arg(parameters[key].toString());
            }
        }

        exe_command.replace(QRegExp("\\$\\(arguments\\)"), ppp);
    }

    QString id = make_random_id();
    PMProcess PP;
    PP.info.exe_command = exe_command;
    PP.info.parameters = parameters;
    PP.info.processor_name = processor_name;
    PP.info.finished = false;
    PP.info.exit_code = 0;
    PP.info.exit_status = QProcess::NormalExit;
    PP.qprocess = new QProcess;
    PP.qprocess->setProcessChannelMode(QProcess::MergedChannels);
    //connect(PP.qprocess,SIGNAL(readyRead()),this,SLOT(slot_qprocess_output()));
    QObject::connect(PP.qprocess, SIGNAL(finished(int, QProcess::ExitStatus)), this, SLOT(slot_process_finished()));
    printf("STARTING: %s.\n", PP.info.exe_command.toLatin1().data());
    PP.qprocess->start(PP.info.exe_command);
    PP.qprocess->setProperty("pp_id", id);
    if (!PP.qprocess->waitForStarted(2000)) {
        qWarning() << "Problem starting process: " + exe_command;
        delete PP.qprocess;
        return "";
    }
    d->m_processes[id] = PP;
    return id;
}

bool ProcessManager::waitForFinished(const QString& process_id, int msecs)
{
    if (!d->m_processes.contains(process_id))
        return false;
    QProcess* qprocess = d->m_processes[process_id].qprocess;
    return MPDaemon::waitForFinishedAndWriteOutput(qprocess);
}

bool ProcessManager::checkParameters(const QString& processor_name, const QVariantMap& parameters)
{
    if (!d->m_processors.contains(processor_name)) {
        qWarning() << "checkProcess: Unable to find processor: " + processor_name;
        return false;
    }
    MLProcessor P = d->m_processors[processor_name];
    {
        QStringList keys = P.inputs.keys();
        foreach(QString key, keys)
        {
            if (!parameters.contains(key)) {
                qWarning() << QString("checkProcess: Missing input in %1: %2").arg(processor_name).arg(key);
                return false;
            }
        }
    }
    {
        QStringList keys = P.outputs.keys();
        foreach(QString key, keys)
        {
            if (!parameters.contains(key)) {
                qWarning() << QString("checkProcess: Missing output in %1: %2").arg(processor_name).arg(key);
                return false;
            }
        }
    }
    {
        QStringList keys = P.parameters.keys();
        foreach(QString key, keys)
        {
            if (!P.parameters[key].optional) {
                if (!parameters.contains(key)) {
                    qWarning() << QString("checkProcess: Missing required parameter in %1: %2").arg(processor_name).arg(key);
                    return false;
                }
            }
        }
    }
    return true;
}

bool ProcessManager::processAlreadyCompleted(const QString& processor_name, const QVariantMap& parameters)
{
    if (!d->m_processors.contains(processor_name))
        return false;

    MLProcessor P = d->m_processors[processor_name];

    if (!d->all_input_and_output_files_exist(P, parameters))
        return false;

    QJsonObject obj = d->compute_unique_process_object(P, parameters);

    QString code = d->compute_unique_object_code(obj);

    QString path = MPDaemon::daemonPath() + "/completed_processes/" + code + ".json";

    return QFile::exists(path);
}

QStringList ProcessManager::allProcessIds() const
{
    return d->m_processes.keys();
}

void ProcessManager::clearProcess(const QString& id)
{
    if (!d->m_processes.contains(id))
        return;
    QProcess* qprocess = d->m_processes[id].qprocess;
    d->m_processes.remove(id);
    if (qprocess->state() == QProcess::Running) {
        d->m_processes[id].qprocess->kill();
        ;
    }
    delete qprocess;
}

void ProcessManager::clearAllProcesses()
{
    d->clear_all_processes();
}

MLProcessInfo ProcessManager::processInfo(const QString& id)
{
    if (!d->m_processes.contains(id)) {
        MLProcessInfo dummy;
        dummy.finished = true;
        dummy.exit_code = 0;
        dummy.exit_status = QProcess::NormalExit;
        return dummy;
    }
    d->update_process_info(id);
    return d->m_processes[id].info;
}

bool ProcessManager::isFinished(const QString& id)
{
    return processInfo(id).finished;
}

Q_GLOBAL_STATIC(ProcessManager, theInstance)
ProcessManager* ProcessManager::globalInstance()
{
    return theInstance;
}

void ProcessManager::slot_process_finished()
{
    QProcess* qprocess = qobject_cast<QProcess*>(sender());
    if (!qprocess) {
        qWarning() << "Unexpected problem in slot_process_finished: qprocess is null.";
        return;
    }
    {
        qprocess->waitForReadyRead();
        QByteArray str1 = qprocess->readAll();
        if (!str1.isEmpty()) {
            printf("%s", str1.data());
        }
    }
    QString id = qprocess->property("pp_id").toString();
    if (!d->m_processes.contains(id)) {
        qWarning() << "Unexpected problem in slot_process_finished. id not found in m_processes: " + id;
        return;
    }
    d->update_process_info(id);
    if ((qprocess->exitCode() == 0) && (qprocess->exitStatus() == QProcess::NormalExit)) {
        QString processor_name = d->m_processes[id].info.processor_name;
        QVariantMap parameters = d->m_processes[id].info.parameters;
        if (!d->m_processors.contains(processor_name)) {
            qWarning() << "Unexpected problem in slot_process_finished. processor not found!!! " + processor_name;
        } else {
            MLProcessor processor = d->m_processors[processor_name];
            QJsonObject obj = d->compute_unique_process_object(processor, parameters);
            QString code = d->compute_unique_object_code(obj);
            QString fname = MPDaemon::daemonPath() + "/completed_processes/" + code + ".json";
            QString json = QJsonDocument(obj).toJson();
            if (QFile::exists(fname))
                QFile::remove(fname); //shouldn't be needed
            if (write_text_file(fname + ".tmp", json)) {
                QFile::rename(fname + ".tmp", fname);
            }
        }
    }
    emit this->processFinished(id);
}

void ProcessManager::slot_qprocess_output()
{
    QProcess* P = qobject_cast<QProcess*>(sender());
    if (!P)
        return;
    QByteArray str = P->readAll();
    if (!str.isEmpty()) {
        printf("%s", str.data());
    }
}

void ProcessManagerPrivate::clear_all_processes()
{
    foreach(PMProcess P, m_processes)
    {
        delete P.qprocess;
    }
    m_processes.clear();
}

void ProcessManagerPrivate::update_process_info(QString id)
{
    if (!m_processes.contains(id))
        return;
    PMProcess* PP = &m_processes[id];
    QProcess* qprocess = PP->qprocess;
    if (qprocess->state() == QProcess::NotRunning) {
        PP->info.finished = true;
        PP->info.exit_code = qprocess->exitCode();
        PP->info.exit_status = qprocess->exitStatus();
    }
    PP->info.standard_output += qprocess->readAll();
}

MLProcessor ProcessManagerPrivate::create_processor_from_json_object(QJsonObject obj)
{
    MLProcessor P;
    P.name = obj["name"].toString();
    P.version = obj["version"].toString();
    P.description = obj["description"].toString();

    QJsonArray inputs = obj["inputs"].toArray();
    for (int i = 0; i < inputs.count(); i++) {
        MLParameter param = create_parameter_from_json_object(inputs[i].toObject());
        P.inputs[param.name] = param;
    }

    QJsonArray outputs = obj["outputs"].toArray();
    for (int i = 0; i < outputs.count(); i++) {
        MLParameter param = create_parameter_from_json_object(outputs[i].toObject());
        P.outputs[param.name] = param;
    }

    QJsonArray parameters = obj["parameters"].toArray();
    for (int i = 0; i < parameters.count(); i++) {
        MLParameter param = create_parameter_from_json_object(parameters[i].toObject());
        P.parameters[param.name] = param;
    }

    P.exe_command = obj["exe_command"].toString();

    return P;
}

MLParameter ProcessManagerPrivate::create_parameter_from_json_object(QJsonObject obj)
{
    MLParameter param;
    param.name = obj["name"].toString();
    param.ptype = obj["ptype"].toString();
    param.description = obj["description"].toString();
    param.optional = obj["optional"].toBool();
    param.default_value = obj["default_value"].toVariant();
    return param;
}

QJsonObject ProcessManagerPrivate::compute_unique_process_object(MLProcessor P, const QVariantMap& parameters)
{
    /*
     * Returns an object that depends uniquely on the following:
     *   1. Processor name and version
     *   2. The paths, sizes, and modification times of the input files (together with their parameter names)
     *   3. Same for the output files
     *   4. The parameters converted to strings
     */

    QJsonObject obj;

    obj["processor_name"] = P.name;
    obj["processor_version"] = P.version;
    {
        QJsonObject inputs;
        QStringList input_pnames = P.inputs.keys();
        qSort(input_pnames);
        foreach(QString input_pname, input_pnames)
        {
            QString fname = parameters[input_pname].toString();
            inputs[input_pname] = create_file_object(fname);
        }
        obj["inputs"] = inputs;
    }
    {
        QJsonObject outputs;
        QStringList output_pnames = P.outputs.keys();
        qSort(output_pnames);
        foreach(QString output_pname, output_pnames)
        {
            QString fname = parameters[output_pname].toString();
            outputs[output_pname] = create_file_object(fname);
        }
        obj["outputs"] = outputs;
    }
    {
        QJsonObject parameters0;
        QStringList pnames = P.parameters.keys();
        qSort(pnames);
        foreach(QString pname, pnames)
        {
            parameters0[pname] = parameters[pname].toString();
        }
        obj["parameters"] = parameters0;
    }
    return obj;
}

bool ProcessManagerPrivate::all_input_and_output_files_exist(MLProcessor P, const QVariantMap& parameters)
{
    QStringList file_pnames = P.inputs.keys();
    file_pnames.append(P.outputs.keys());
    foreach(QString pname, file_pnames)
    {
        QString fname = parameters.value(pname).toString();
        if (!fname.isEmpty()) {
            if (!QFile::exists(fname))
                return false;
        }
    }
    return true;
}

QString ProcessManagerPrivate::compute_unique_object_code(QJsonObject obj)
{
    /// Witold I need a string that depends on the json object. However I am worried about the order of the fields. Is there a way to make this canonical?
    QByteArray json = QJsonDocument(obj).toJson();
    QCryptographicHash hash(QCryptographicHash::Sha1);
    hash.addData(json);
    return QString(hash.result().toHex());
}

QJsonObject ProcessManagerPrivate::create_file_object(const QString& fname)
{
    QJsonObject obj;
    if (fname.isEmpty())
        return obj;
    obj["path"] = fname;
    if (!QFile::exists(fname)) {
        obj["size"] = 0;
        return obj;
    }
    obj["size"] = QFileInfo(fname).size();
    obj["last_modified"] = QFileInfo(fname).lastModified().toString("yyyy-MM-dd-hh-mm-ss-zzz");
    return obj;
}

QChar make_random_alphanumeric()
{
    static int val = 0;
    val++;
    QTime time = QTime::currentTime();
    int num = qHash(time.toString("hh:mm:ss:zzz") + QString::number(qrand() + val));
    if (num < 0)
        num = -num;
    num = num % 36;
    if (num < 26)
        return QChar('A' + num);
    else
        return QChar('0' + num - 26);
}
QString make_random_id(int numchars)
{
    QString ret;
    for (int i = 0; i < numchars; i++) {
        ret.append(make_random_alphanumeric());
    }
    return ret;
}
