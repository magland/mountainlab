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
#include "mlcommon.h"

#include <QCoreApplication>
#include <QThread>
#include <QTimer>
#include <cachemanager.h>
#include "mpdaemon.h"

struct PMProcess {
    MLProcessInfo info;
    QString tempdir = "";
    bool preserve_tempdir = true; // for debugging set to true, otherwise will clean up the tempdir when process has finished
    bool exec_mode = false;
    QProcess* qprocess;
};

class ProcessManagerPrivate {
public:
    ProcessManager* q;

    QStringList m_processor_paths;
    QMap<QString, MLProcessor> m_processors;
    QMap<QString, PMProcess> m_processes;
    //QStringList m_server_urls;
    //QString m_server_base_path;

    void clear_all_processes();
    void update_process_info(QString id);
    QString resolve_file_name_p(QString fname);
    QVariantMap resolve_file_names_in_parameters(QString processor_name, const QVariantMap& parameters);
    QString compute_unique_object_code(QJsonObject obj);
    QJsonObject compute_unique_process_object(MLProcessor P, const QVariantMap& parameters, bool rprv_inputs);
    bool all_input_and_output_files_exist(MLProcessor P, const QVariantMap& parameters, bool allow_rprv_inputs, bool allow_rprv_outputs);
    QJsonObject create_file_object(const QString& fname, bool allow_rprv_inputs);
    void reload_processors();

    static MLProcessor create_processor_from_json_object(QJsonObject obj);
    static MLParameter create_parameter_from_json_object(QJsonObject obj);
};

ProcessManager::ProcessManager()
{
    d = new ProcessManagerPrivate;
    d->q = this;

    QTimer::singleShot(1000, this, SLOT(slot_monitor()));
}

ProcessManager::~ProcessManager()
{
    d->clear_all_processes();
    delete d;
}

void ProcessManager::setProcessorPaths(const QStringList& paths)
{
    d->m_processor_paths = paths;
    d->reload_processors();
}

void ProcessManager::reloadProcessors()
{
    d->reload_processors();
}

/*
void ProcessManager::setServerUrls(const QStringList& urls)
{
    d->m_server_urls = urls;
}

void ProcessManager::setServerBasePath(const QString& path)
{
    d->m_server_base_path = path;
}
*/

bool ProcessManager::loadProcessors(const QString& path, bool recursive)
{
    QStringList fnames = QDir(path).entryList(QStringList("*.mp"), QDir::Files, QDir::Name);
    foreach (QString fname, fnames) {
        if (!this->loadProcessorFile(path + "/" + fname)) {
            return false;
        }
    }
    if (recursive) {
        QStringList subdirs = QDir(path).entryList(QStringList("*"), QDir::Dirs | QDir::NoDotAndDotDot, QDir::Name);
        foreach (QString subdir, subdirs) {
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
        if ((QFile::exists(path + ".spec")) && (QFileInfo(path).lastModified().secsTo(QFileInfo(path + ".spec").lastModified()) >= 0) && (QFileInfo(path + ".spec").lastModified().secsTo(QDateTime::currentDateTime()) <= 60)) {
            json = TextFile::read(path + ".spec"); // read the saved spec so we don't need to make the system call next time, but let's be careful about it... regenerate every 60 seconds
        }
        else {
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
                qWarning() << "Potential problem with executable processor file: " + path + ". Expected json output but got empty string.";
                if (QFileInfo(path).size() < 1e6) {
                    json = TextFile::read(path);
                    //now test it, since it is executable we are suspicious...
                    QJsonParseError error;
                    QJsonObject obj = QJsonDocument::fromJson(json.toLatin1(), &error).object();
                    if (error.error != QJsonParseError::NoError) {
                        qWarning() << "Executable processor file did not return output for spec: " + path;
                        return false;
                    }
                    else {
                        //we are okay -- apparently the text file .mp got marked as executable by the user, so let's proceed
                    }
                }
                else {
                    qWarning() << "File is too large to be a text file. Executable processor file did not return output for spec: " + path;
                    return false;
                }
            }
            else {
                TextFile::write(path + ".spec", json); // so we don't need to make the system call this time
            }
        }
    }
    else {
        json = TextFile::read(path);
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

void delete_tempdir(QString tempdir)
{
    //return; // Uncomment this line for debugging -- but remember to put it back!!!
    // to be safe!
    if (!tempdir.startsWith(CacheManager::globalInstance()->localTempPath())) {
        qWarning() << "Unexpected tempdir in delete_tempdir: " + tempdir;
        return;
    }

    {
        QStringList list = QDir(tempdir).entryList(QStringList("*"), QDir::Files, QDir::Name);
        foreach (QString fname, list) {
            QFile::remove(tempdir + "/" + fname);
        }
    }
    {
        QStringList list = QDir(tempdir).entryList(QStringList("*"), QDir::Dirs | QDir::NoDotAndDotDot, QDir::Name);
        foreach (QString foldername, list) {
            delete_tempdir(tempdir + "/" + foldername);
        }
    }
    QDir(QFileInfo(tempdir).path()).rmdir(QFileInfo(tempdir).fileName());
}

bool do_mkdir(QString path, int num_tries = 1)
{
    if (!QDir(QFileInfo(path).path()).mkdir(QFileInfo(path).fileName())) {
        if (num_tries > 1) {
            do_mkdir(path, num_tries - 1);
        }
    }
    return true;
}

QString ProcessManager::startProcess(const QString& processor_name, const QVariantMap& parameters_in, const RequestProcessResources& RPR, bool exec_mode, bool preserve_tempdir)
{
    QVariantMap parameters = d->resolve_file_names_in_parameters(processor_name, parameters_in);

    if (!this->checkParameters(processor_name, parameters)) {
        qWarning() << "Problem checking parameters";
        return "";
    }

    if (!d->m_processors.contains(processor_name)) {
        qWarning() << "Unable to find processor (165): " + processor_name;
        return "";
    }
    MLProcessor P = d->m_processors[processor_name];

    {
        // check that we have the required parameters
        QStringList keys = P.parameters.keys();
        foreach (QString key, keys) {
            if (!parameters.contains(key)) {
                if (!P.parameters[key].optional) {
                    qWarning() << QString("Missing required parameter in processor %1: %2").arg(P.name).arg(key);
                    return "";
                }
            }
        }
    }
    this->setDefaultParameters(processor_name, parameters);

    QString id = MLUtil::makeRandomId();

    QString tempdir = CacheManager::globalInstance()->makeLocalFile("tempdir_" + id, CacheManager::ShortTerm);
    if (!do_mkdir(tempdir, 3)) {
        qWarning() << "Error creating temporary directory for process: " + tempdir;
        QCoreApplication::exit(-1);
    }

    QString exe_command = P.exe_command;
    exe_command.replace(QRegExp("\\$\\(basepath\\)"), P.basepath);
    exe_command.replace(QRegExp("\\$\\(tempdir\\)"), tempdir);
    {
        QString ppp;
        {
            QStringList keys = P.inputs.keys();
            foreach (QString key, keys) {
                QStringList list = MLUtil::toStringList(parameters[key]);
                exe_command.replace(QRegExp(QString("\\$%1\\$").arg(key)), list.value(0)); //note that only the first file name is inserted here
                foreach (QString str, list) {
                    ppp += QString("--%1=%2 ").arg(key).arg(str);
                }
            }
        }
        {
            QStringList keys = P.outputs.keys();
            foreach (QString key, keys) {
                QStringList list = MLUtil::toStringList(parameters[key]);
                exe_command.replace(QRegExp(QString("\\$%1\\$").arg(key)), list.value(0)); //note that only the first file name is inserted here
                foreach (QString str, list) {
                    ppp += QString("--%1=%2 ").arg(key).arg(str);
                }
            }
        }
        {
            QStringList keys = P.parameters.keys();
            foreach (QString key, keys) {
                exe_command.replace(QRegExp(QString("\\$%1\\$").arg(key)), parameters[key].toString());
                ppp += QString("--%1=%2 ").arg(key).arg(parameters[key].toString());
            }
        }

        if (RPR.request_num_threads) {
            ppp += QString("--_request_num_threads=%1 ").arg(RPR.request_num_threads);
        }
        ppp += QString("--_tempdir=%1 ").arg(tempdir);

        exe_command.replace(QRegExp("\\$\\(arguments\\)"), ppp);
    }

    PMProcess PP;
    PP.info.exe_command = exe_command;
    PP.info.parameters = parameters;
    PP.info.processor_name = processor_name;
    PP.info.finished = false;
    PP.info.exit_code = 0;
    PP.info.exit_status = QProcess::NormalExit;
    PP.tempdir = tempdir;
    PP.preserve_tempdir = preserve_tempdir;
    PP.exec_mode = exec_mode;
    PP.qprocess = new QProcess;
    PP.qprocess->setProcessChannelMode(QProcess::MergedChannels);
    //connect(PP.qprocess,SIGNAL(readyRead()),this,SLOT(slot_qprocess_output()));
    QObject::connect(PP.qprocess, SIGNAL(finished(int, QProcess::ExitStatus)), this, SLOT(slot_process_finished()));
    printf("STARTING: %s.\n", PP.info.exe_command.toLatin1().data());
    PP.info.start_time = QDateTime::currentDateTime();

    // Do it this way so that special characters are handled exactly like a system call
    MPDaemon::start_bash_command_and_kill_when_pid_is_gone(PP.qprocess, PP.info.exe_command, QCoreApplication::applicationPid());
    //PP.qprocess->start(PP.info.exe_command); //not this way

    PP.qprocess->setProperty("pp_id", id);
    if (!PP.qprocess->waitForStarted(2000)) {
        qWarning() << "Problem starting process: " + exe_command;
        delete_tempdir(tempdir);
        delete PP.qprocess;
        return "";
    }
    d->m_processes[id] = PP;
    return id;
}

bool ProcessManager::waitForFinished(const QString& process_id, int parent_pid)
{
    if (!d->m_processes.contains(process_id))
        return false;
    QProcess* qprocess = d->m_processes[process_id].qprocess;
    return MPDaemon::waitForFinishedAndWriteOutput(qprocess, parent_pid);
}

bool ProcessManager::checkParameters(const QString& processor_name, const QVariantMap& parameters)
{
    if (processor_name.isEmpty())
        return true; //empty processor name is for handling .prv files
    if (!d->m_processors.contains(processor_name)) {
        qWarning() << "checkProcess: Unable to find processor (233): " + processor_name;
        return false;
    }
    if (processor_name.isEmpty())
        return true; //empty processor name is for handling .prv files
    MLProcessor P = d->m_processors[processor_name];
    {
        QStringList keys = P.inputs.keys();
        foreach (QString key, keys) {
            if (!P.inputs[key].optional) {
                if (!parameters.contains(key)) {
                    qWarning() << QString("checkProcess: Missing input in %1: %2").arg(processor_name).arg(key);
                    return false;
                }
            }
        }
    }
    {
        QStringList keys = P.outputs.keys();
        foreach (QString key, keys) {
            if (!parameters.contains(key)) {
                if (!P.outputs[key].optional) {
                    qWarning() << QString("checkProcess: Missing required output in %1: %2").arg(processor_name).arg(key);
                    return false;
                }
            }
        }
    }
    {
        QStringList keys = P.parameters.keys();
        foreach (QString key, keys) {
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

void ProcessManager::setDefaultParameters(const QString& processor_name, QVariantMap& parameters)
{
    if (processor_name.isEmpty())
        return;
    if (!d->m_processors.contains(processor_name)) {
        return;
    }
    MLProcessor P = d->m_processors[processor_name];
    {
        QStringList keys = P.parameters.keys();
        foreach (QString key, keys) {
            if (P.parameters[key].optional) {
                if (!parameters.contains(key)) {
                    parameters[key] = P.parameters[key].default_value;
                }
            }
        }
    }
}

bool ProcessManager::processAlreadyCompleted(const QString& processor_name, const QVariantMap& parameters, bool allow_rprv_inputs, bool allow_rprv_outputs)
{
    if (!d->m_processors.contains(processor_name))
        return false;

    MLProcessor P = d->m_processors[processor_name];

    if (!d->all_input_and_output_files_exist(P, parameters, allow_rprv_inputs, allow_rprv_outputs))
        return false;

    QJsonObject obj = d->compute_unique_process_object(P, parameters, allow_rprv_inputs);

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
        qprocess->kill();
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

bool output_object_still_exists(const QJsonObject& obj)
{
    QString path = obj["path"].toString();
    QDateTime last_modified = QDateTime::fromString(obj["last_modified"].toString(), "yyyy-MM-dd-hh-mm-ss-zzz");
    bigint size = obj["size"].toDouble();
    if (QFile::exists(path)) {
        if (QFileInfo(path).size() == size) {
            if (QFileInfo(path).lastModified().secsTo(last_modified) == 0) { //round to nearest second is important i think
                return true;
            }
        }
    }
    else {
        if (QFile::exists(path + ".rprv")) {
            QJsonObject rprv = QJsonDocument::fromJson(TextFile::read(path + ".rprv").toUtf8()).object();
            QDateTime last_modified0 = QDateTime::fromString(rprv["original_last_modified"].toString(), "yyyy-MM-dd-hh-mm-ss-zzz");
            bigint size0 = rprv["original_size"].toDouble();
            if (size0 == size) {
                if (last_modified0.secsTo(last_modified) == 0) { //round to nearest second is important i think
                    return true;
                }
            }
        }
    }
    return false;
}

void ProcessManager::cleanUpCompletedProcessRecords()
{
    qDebug().noquote() << "Cleaning up completed process records...";
    QString path = MPDaemon::daemonPath() + "/completed_processes";
    QStringList list = QDir(path).entryList(QStringList("*.json"), QDir::Files, QDir::Name);
    foreach (QString fname, list) {
        QString path2 = path + "/" + fname;
        QJsonObject obj = QJsonDocument::fromJson(TextFile::read(path2).toUtf8()).object();
        QJsonObject outputs = obj["outputs"].toObject();
        QStringList output_pnames = outputs.keys();
        foreach (QString pname, output_pnames) {
            QJsonObject A = outputs[pname].toObject();
            if (A["size"].toDouble()) { //maybe there is no such output
                if (!output_object_still_exists(A)) {
                    qDebug().noquote() << "Removing completed process record: " + path2;
                    QFile::remove(path2);
                }
            }
        }
    }
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
        qCritical() << "Unexpected problem in slot_process_finished: qprocess is null.";
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
        qWarning() << "Process not found in slot_process_finished. It was probably cleared beforehand: " + id;
        return;
    }
    d->update_process_info(id);
    if ((qprocess->exitCode() == 0) && (qprocess->exitStatus() == QProcess::NormalExit)) {
        QString processor_name = d->m_processes[id].info.processor_name;
        QVariantMap parameters = d->m_processes[id].info.parameters;
        if (!d->m_processors.contains(processor_name)) {
            qCritical() << "Unexpected problem in slot_process_finished. processor not found!!! " + processor_name;
        }
        else {
            MLProcessor processor = d->m_processors[processor_name];
            if (!d->m_processes[id].exec_mode) { //in exec_mode we don't keep track of which processes have already completed
                QJsonObject obj = d->compute_unique_process_object(processor, parameters, false);
                QString code = d->compute_unique_object_code(obj);
                QString fname = MPDaemon::daemonPath() + "/completed_processes/" + code + ".json";
                QString json = QJsonDocument(obj).toJson();
                if (QFile::exists(fname))
                    QFile::remove(fname); //shouldn't be needed
                if (TextFile::write(fname + ".tmp", json)) {
                    QFile::rename(fname + ".tmp", fname);
                    QFile::Permissions perm = QFileDevice::ReadUser | QFileDevice::WriteUser | QFileDevice::ExeUser | QFileDevice::ReadGroup | QFileDevice::WriteGroup | QFileDevice::ExeGroup | QFileDevice::ReadOther | QFileDevice::WriteOther | QFileDevice::ExeOther;
                    QFile::setPermissions(fname, perm);
                }
            }
        }
    }
    if (!d->m_processes[id].preserve_tempdir)
        delete_tempdir(d->m_processes[id].tempdir);
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

QString execute_and_read_stdout(QString cmd)
{
    QProcess P;
    P.start(cmd);
    P.waitForStarted();
    P.waitForFinished();
    return P.readAllStandardOutput();
}

void ProcessManager::slot_monitor()
{
    QStringList ids = d->m_processes.keys();
    foreach (QString id, ids) {
        PMProcess* PP = &d->m_processes[id];
        if (PP->qprocess) {
            QString cmd = QString("ps -p %1 -o rss,%cpu --noheader").arg(PP->qprocess->pid());
            QString str = execute_and_read_stdout(cmd);
            QStringList list = str.split(" ", QString::SkipEmptyParts);
            MonitorStats MS;
            MS.timestamp = QDateTime::currentDateTime();
            MS.mem_bytes = list.value(0).toLong() * 1000;
            MS.cpu_pct = list.value(1).toDouble();
            PP->info.monitor_stats << MS;
        }
    }

    QTimer::singleShot(1000, this, SLOT(slot_monitor()));
}

void ProcessManagerPrivate::clear_all_processes()
{
    foreach (PMProcess P, m_processes) {
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
        PP->info.finish_time = QDateTime::currentDateTime();
        PP->info.finished = true;
        PP->info.exit_code = qprocess->exitCode();
        PP->info.exit_status = qprocess->exitStatus();
    }
    PP->info.standard_output += qprocess->readAll();
}

QString ProcessManagerPrivate::resolve_file_name_p(QString fname)
{
    //for now I have completely removed the client/server functionality

    if (fname.endsWith(".prv")) {
        QString txt = TextFile::read(fname);
        QJsonParseError err;
        QJsonObject obj = QJsonDocument::fromJson(txt.toUtf8(), &err).object();
        if (err.error != QJsonParseError::NoError) {
            qWarning() << "Error parsing .prv file: " + fname;
            return "";
        }
        QString path0 = locate_prv(obj);
        if (path0.isEmpty()) {
            qWarning() << "Unable to locate prv file originally at: " + obj["original_path"].toString();
            return "";
        }
        return path0;
    }
    else
        return fname;

    /*
    //This is terrible, we need to fix it!
    QString fname = fname_in;
    foreach (QString str, server_urls) {
        if (fname.startsWith(str + "/mdaserver")) {
            fname = server_base_path + "/" + fname.mid((str + "/mdaserver").count());
            if (fname.mid(server_base_path.count()).contains("..")) {
                qWarning() << "Illegal .. in file path: " + fname;
                fname = "";
            }
        }
    }
    */

    return fname;
}

QVariantMap ProcessManagerPrivate::resolve_file_names_in_parameters(QString processor_name, const QVariantMap& parameters_in)
{
    QVariantMap parameters = parameters_in;

    if (!m_processors.contains(processor_name)) {
        qWarning() << "Unable to find processor in resolve_file_names_in_parameters(): " + processor_name;
        return parameters_in;
    }
    MLProcessor MLP = m_processors[processor_name];
    foreach (MLParameter P, MLP.inputs) {
        QStringList list = MLUtil::toStringList(parameters[P.name]);
        if (list.count() == 1) {
            parameters[P.name] = resolve_file_name_p(list[0]);
        }
        else {
            QVariantList list2;
            foreach (QString str, list) {
                list2 << resolve_file_name_p(str);
            }
            parameters[P.name] = list2;
        }
    }
    foreach (MLParameter P, MLP.outputs) {
        QStringList list = MLUtil::toStringList(parameters[P.name]);
        if (list.count() == 1) {
            parameters[P.name] = resolve_file_name_p(list[0]);
        }
        else {
            QVariantList list2;
            foreach (QString str, list) {
                list2 << resolve_file_name_p(str);
            }
            parameters[P.name] = list2;
        }
    }
    return parameters;
}

MLProcessor ProcessManagerPrivate::create_processor_from_json_object(QJsonObject obj)
{
    MLProcessor P;
    P.spec = obj;
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

QJsonObject ProcessManagerPrivate::compute_unique_process_object(MLProcessor P, const QVariantMap& parameters, bool allow_rprv_inputs)
{
    /*
     * Returns an object that depends uniquely on the following:
     *   1. Version of mountainprocess
     *   2. Processor name and version
     *   3. The paths, sizes, and modification times of the input files (together with their parameter names)
     *   4. Same for the output files
     *   5. The parameters converted to strings
     */

    QJsonObject obj;

    obj["mountainprocess_version"] = "0.1";
    obj["processor_name"] = P.name;
    obj["processor_version"] = P.version;
    {
        QJsonObject inputs;
        QStringList input_pnames = P.inputs.keys();
        qSort(input_pnames);
        foreach (QString input_pname, input_pnames) {
            QStringList fnames = MLUtil::toStringList(parameters[input_pname]);
            if (fnames.count() == 1) {
                inputs[input_pname] = create_file_object(resolve_file_name_p(fnames[0]), allow_rprv_inputs);
            }
            else {
                QJsonArray array;
                foreach (QString fname0, fnames) {
                    QString fname = resolve_file_name_p(fname0);
                    array.append(create_file_object(fname, allow_rprv_inputs));
                }
                inputs[input_pname] = array;
            }
        }
        obj["inputs"] = inputs;
    }
    {
        QJsonObject outputs;
        QStringList output_pnames = P.outputs.keys();
        qSort(output_pnames);
        foreach (QString output_pname, output_pnames) {
            QStringList fnames = MLUtil::toStringList(parameters[output_pname]);
            if (fnames.count() == 1) {
                outputs[output_pname] = create_file_object(resolve_file_name_p(fnames[0]), true);
            }
            else {
                QJsonArray array;
                foreach (QString fname0, fnames) {
                    QString fname = resolve_file_name_p(fname0);
                    array.append(create_file_object(fname, true));
                }
                outputs[output_pname] = array;
            }
        }
        obj["outputs"] = outputs;
    }
    {
        QJsonObject parameters0;
        QStringList pnames = P.parameters.keys();
        qSort(pnames);
        foreach (QString pname, pnames) {
            parameters0[pname] = parameters[pname].toString();
        }
        obj["parameters"] = parameters0;
    }
    return obj;
}

bool ProcessManagerPrivate::all_input_and_output_files_exist(MLProcessor P, const QVariantMap& parameters, bool allow_rprv_input_files, bool allow_rprv_output_files)
{
    QStringList input_file_pnames = P.inputs.keys();
    QStringList output_file_pnames = P.outputs.keys();

    foreach (QString pname, input_file_pnames) {
        QStringList fnames = MLUtil::toStringList(parameters.value(pname));
        foreach (QString fname0, fnames) {
            if (!fname0.isEmpty()) {
                QString fname = resolve_file_name_p(fname0);
                if (!fname.isEmpty()) {
                    if (!QFile::exists(fname)) {
                        if (allow_rprv_input_files) {
                            if (!QFile::exists(fname + ".rprv")) {
                                return false;
                            }
                        }
                        else {
                            return false;
                        }
                    }
                }
            }
        }
    }

    foreach (QString pname, output_file_pnames) {
        QStringList fnames = MLUtil::toStringList(parameters.value(pname));
        foreach (QString fname0, fnames) {
            if (!fname0.isEmpty()) {
                QString fname = resolve_file_name_p(fname0);
                if (!fname.isEmpty()) {
                    if (!QFile::exists(fname)) {
                        if (allow_rprv_output_files) {
                            if (!QFile::exists(fname + ".rprv")) {
                                return false;
                            }
                        }
                        else {
                            return false;
                        }
                    }
                }
            }
        }
    }

    return true;
}

QString ProcessManagerPrivate::compute_unique_object_code(QJsonObject obj)
{
    /// Witold I need a string that depends on the json object. However I am worried about the order of the fields. Is there a way to make this canonical?
    /// Jeremy: You can sort all keys in the dictionary, convert that to string and
    ///         calculate hash of that. However this is going to be CPU consuming
    QByteArray json = QJsonDocument(obj).toJson();
    QCryptographicHash hash(QCryptographicHash::Sha1);
    hash.addData(json);
    return QString(hash.result().toHex());
}

QJsonObject ProcessManagerPrivate::create_file_object(const QString& fname_in, bool allow_rprv_inputs)
{
    QString fname = resolve_file_name_p(fname_in);
    QJsonObject obj;
    if (fname.isEmpty())
        return obj;
    obj["path"] = fname;
    if (!QFile::exists(fname)) {
        if ((allow_rprv_inputs) && (QFile::exists(fname_in + ".rprv"))) {
            QString json0 = TextFile::read(fname_in + ".rprv");
            QJsonObject obj0 = QJsonDocument::fromJson(json0.toUtf8()).object();
            obj["size"] = obj0["original_size"];
            obj["last_modified"] = obj0["original_last_modified"];
        }
        else {
            obj["size"] = 0;
            return obj;
        }
    }
    else {
        obj["size"] = QFileInfo(fname).size();
        obj["last_modified"] = QFileInfo(fname).lastModified().toString("yyyy-MM-dd-hh-mm-ss-zzz");
    }
    if (QFileInfo(fname).isDir()) {
        QStringList fnames = QDir(fname).entryList(QDir::Files, QDir::Name);
        QJsonArray files_array;
        foreach (QString fname2, fnames) {
            QJsonObject obj0;
            obj0["name"] = fname2;
            obj0["object"] = create_file_object(fname + "/" + fname2, false);
            files_array.push_back(obj0);
        }
        obj["files"] = files_array;

        QStringList dirnames = QDir(fname).entryList(QDir::Dirs | QDir::NoDotAndDotDot, QDir::Name);
        QJsonArray dirs_array;
        foreach (QString dirname2, dirnames) {
            QJsonObject obj0;
            obj0["name"] = dirname2;
            obj0["object"] = create_file_object(fname + "/" + dirname2, false);
            dirs_array.push_back(obj0);
        }
        obj["directories"] = dirs_array;
    }
    return obj;
}

void ProcessManagerPrivate::reload_processors()
{
    QMap<QString, MLProcessor> saved = m_processors;
    m_processors.clear();
    foreach (QString processor_path, m_processor_paths) {
        if (!q->loadProcessors(processor_path)) {
            //something happened, let's revert to saved version
            m_processors = saved;
            return;
        }
    }
}
