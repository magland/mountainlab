/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland, Witold Wysota
** Created: 5/2/2016
*******************************************************/

#include "mpdaemon.h"
#include <QEventLoop>
#include <QCoreApplication>
#include <QTimer>
#include <QDir>
#include <QJsonDocument>
#include <QJsonObject>
#include <QFileSystemWatcher>
#include <QJsonArray>
#include <QProcess>

#include <QDebug>
#include "cachemanager.h"
#include "unistd.h" //for usleep
#include <sys/stat.h> //for mkfifo
#include "processmanager.h"
#include "mlcommon.h"
#include <QSettings>
#include <QSharedMemory>
#include <objectregistry.h>
#include <qprocessmanager.h>
#include <signal.h>

static bool stopDaemon = false;

void sighandler(int num)
{
    if (num == SIGINT || num == SIGTERM) {
        stopDaemon = true;
        qApp->quit();
    }
}

class MountainProcessServerClient : public LocalServer::Client {
public:
    MountainProcessServerClient(QLocalSocket* sock, LocalServer::Server* parent)
        : LocalServer::Client(sock, parent)
    {
    }
    ~MountainProcessServerClient()
    {
    }

    void close();

protected:
    bool handleMessage(const QByteArray& ba) Q_DECL_OVERRIDE;
    bool getState();
};

void append_line_to_file(QString fname, QString line)
{
    QFile ff(fname);
    if (ff.open(QFile::Append | QFile::WriteOnly)) {
        ff.write(line.toLatin1());
        ff.write(QByteArray("\n"));
        ff.close();
    }
    else {
        static bool reported = false;
        if (!reported) {
            reported = true;
            qCritical() << "Unable to write to log file" << fname;
        }
    }
}

void debug_log(const char* function, const char* file, int line)
{
    QString fname = CacheManager::globalInstance()->localTempPath() + QString("/mpdaemon_debug_%1.log").arg(QString(qgetenv("USER")));
    QString line0 = QString("%1: %2 %3:%4").arg(QDateTime::currentDateTime().toString("yy-MM-dd:hh:mm:ss.zzz")).arg(function).arg(file).arg(line);
    line0 += " ARGS: ";
    foreach (QString arg, qApp->arguments()) {
        line0 += arg + " ";
    }
    append_line_to_file(fname, line0);
}

#if 0
MPDaemon::MPDaemon()
{
    debug_log(__FUNCTION__, __FILE__, __LINE__);

    d = new MPDaemonPrivate;
    d->q = this;
    d->m_is_running = false;

    d->m_total_resources_available.num_threads = 12;
    d->m_total_resources_available.memory_gb = 8;
}
#endif

void kill_process_and_children(int pid)
{
    /// Witold, do we need to worry about making this cross-platform?
    QString cmd = QString("CPIDS=$(pgrep -P %1); (sleep 33 && kill -KILL $CPIDS &); kill -TERM $CPIDS").arg(pid);
    int ret = system(cmd.toUtf8().data());
    Q_UNUSED(ret);
}

void kill_process_and_children(QProcess* P)
{
    kill_process_and_children(P->processId());
}

#include "signal.h"
bool pidExists(qint64 pid)
{
    // check whether process exists (works on Linux)
    //return (kill(pid, 0) == 0);

    // the following is supposed to work even when the process is owned by a different user
    return (getpgid(pid) >= 0);
}

QStringList json_array_to_stringlist(QJsonArray X)
{
    QStringList ret;
    for (int i = 0; i < X.count(); i++) {
        ret << X[i].toString();
    }
    return ret;
}

QJsonObject variantmap_to_json_obj(QVariantMap map)
{
    return QJsonObject::fromVariantMap(map);
}

QStringList paths_to_file_names(const QStringList& paths)
{
    QStringList ret;
    foreach (QString path, paths) {
        ret << QFileInfo(path).fileName();
    }
    return ret;
}

QJsonObject pript_struct_to_obj(MPDaemonPript S, RecordType rt)
{
    QJsonObject ret;
    ret["is_finished"] = S.is_finished;
    ret["is_running"] = S.is_running;
    if (rt != AbbreviatedRecord) {
        ret["parameters"] = variantmap_to_json_obj(S.parameters);
        ret["output_fname"] = S.output_fname;
        ret["preserve_tempdir"] = S.preserve_tempdir;
        ret["stdout_fname"] = S.stdout_fname;
        ret["parent_pid"] = QString("%1").arg(S.parent_pid);
    }
    ret["force_run"] = S.force_run;
    ret["working_path"] = S.working_path;
    ret["id"] = S.id;
    ret["success"] = S.success;
    ret["error"] = S.error;
    ret["timestamp_queued"] = S.timestamp_queued.toString("yyyy-MM-dd|hh:mm:ss.zzz");
    ret["timestamp_started"] = S.timestamp_started.toString("yyyy-MM-dd|hh:mm:ss.zzz");
    ret["timestamp_finished"] = S.timestamp_finished.toString("yyyy-MM-dd|hh:mm:ss.zzz");
    ret["request_num_threads"] = S.RPR.request_num_threads;
    if (S.prtype == ScriptType) {
        ret["prtype"] = "script";
        if (rt != AbbreviatedRecord) {
            ret["script_paths"] = QJsonArray::fromStringList(S.script_paths);
            ret["script_path_checksums"] = QJsonArray::fromStringList(S.script_path_checksums);
        }
        else {
            ret["script_names"] = QJsonArray::fromStringList(paths_to_file_names(S.script_paths));
        }
    }
    else {
        ret["prtype"] = "process";
        ret["processor_name"] = S.processor_name;
        ret["processor_spec"] = S.processor_spec;
    }
    if (rt == RuntimeRecord) {
        ret["runtime_opts"] = runtime_opts_struct_to_obj(S.runtime_opts);
        ret["runtime_results"] = S.runtime_results;
    }
    return ret;
}

MPDaemonPript pript_obj_to_struct(QJsonObject obj)
{
    MPDaemonPript ret;
    ret.is_finished = obj.value("is_finished").toBool();
    ret.is_running = obj.value("is_running").toBool();
    ret.parameters = obj.value("parameters").toObject().toVariantMap();
    ret.id = obj.value("id").toString();
    ret.output_fname = obj.value("output_fname").toString();
    ret.preserve_tempdir = obj.value("preserve_tempdir").toBool();
    ret.stdout_fname = obj.value("stdout_fname").toString();
    ret.success = obj.value("success").toBool();
    ret.error = obj.value("error").toString();
    ret.parent_pid = obj.value("parent_pid").toString().toLongLong();
    ret.force_run = obj.value("force_run").toBool();
    ret.working_path = obj.value("working_path").toString();
    ret.timestamp_queued = QDateTime::fromString(obj.value("timestamp_queued").toString(), "yyyy-MM-dd|hh:mm:ss.zzz");
    ret.timestamp_started = QDateTime::fromString(obj.value("timestamp_started").toString(), "yyyy-MM-dd|hh:mm:ss.zzz");
    ret.timestamp_finished = QDateTime::fromString(obj.value("timestamp_finished").toString(), "yyyy-MM-dd|hh:mm:ss.zzz");
    ret.RPR.request_num_threads = obj.value("request_num_threads").toInt();
    if (obj.value("prtype").toString() == "script") {
        ret.prtype = ScriptType;
        ret.script_paths = json_array_to_stringlist(obj.value("script_paths").toArray());
        ret.script_path_checksums = json_array_to_stringlist(obj.value("script_path_checksums").toArray());
    }
    else {
        ret.prtype = ProcessType;
        ret.processor_name = obj.value("processor_name").toString();
        ret.processor_spec = obj.value("processor_spec").toObject();
    }
    return ret;
}

bool is_at_most(ProcessResources PR1, ProcessResources PR2, ProcessResources total_available)
{
    if (total_available.memory_gb != 0) {
        if (PR1.memory_gb > PR2.memory_gb)
            return false;
    }
    if (total_available.num_threads != 0) {
        if (PR1.num_threads > PR2.num_threads)
            return false;
    }
    if (total_available.num_processes != 0) {
        if (PR1.num_processes > PR2.num_processes)
            return false;
    }
    return true;
}

QJsonObject runtime_opts_struct_to_obj(ProcessRuntimeOpts opts)
{
    QJsonObject ret;
    ret["memory_gb_allotted"] = opts.memory_gb_allotted;
    ret["num_threads_allotted"] = opts.num_threads_allotted;
    return ret;
}

bool MountainProcessServerClient::getState()
{
    MountainProcessServer* s = static_cast<MountainProcessServer*>(server());
    writeMessage(QJsonDocument(s->state()).toJson());
    return true;
}

void MountainProcessServerClient::close()
{
    MountainProcessServer* srvr = static_cast<MountainProcessServer*>(server());
    srvr->unregisterLogListener(this);
    LocalServer::Client::close();
}

bool MountainProcessServerClient::handleMessage(const QByteArray& ba)
{
    QJsonParseError error;
    QJsonObject obj = QJsonDocument::fromJson(ba, &error).object();
    if (error.error != QJsonParseError::NoError) {
        qCritical() << "Error in slot_commands_directory_changed parsing json file";
    }
    printf("Received commmand: %s\n", obj["command"].toString().toUtf8().constData());
    if (obj["command"] == "stop") {
        MountainProcessServer* s = static_cast<MountainProcessServer*>(server());
        s->stop();
    }
    if (obj["command"] == "get-daemon-state") {
        return getState();
    }
    if (obj["command"] == "get-log") {
        MountainProcessServer* s = static_cast<MountainProcessServer*>(server());
        QJsonArray log = s->log();
        writeMessage(QJsonDocument(log).toJson());
        return true;
    }
    if (obj["command"] == "log-listener") {
        MountainProcessServer* srvr = static_cast<MountainProcessServer*>(server());
        srvr->registerLogListener(this);
        writeMessage(QJsonDocument(srvr->log()).toJson());
        return true;
    }
    if (obj["command"] == "queue-script") {
        MountainProcessServer* srvr = static_cast<MountainProcessServer*>(server());
        MPDaemonPript S = pript_obj_to_struct(obj);
        S.prtype = ScriptType;
        S.timestamp_queued = QDateTime::currentDateTime();
        if (!srvr->queueScript(S)) {
            // TODO: Send back info to the client?
        }
        return true;
    }
    if (obj["command"] == "queue-process") {
        MountainProcessServer* srvr = static_cast<MountainProcessServer*>(server());
        MPDaemonPript P = pript_obj_to_struct(obj);
        P.prtype = ProcessType;
        P.timestamp_queued = QDateTime::currentDateTime();
        if (!srvr->queueProcess(P)) {
            // TODO: Send back info to the client?
        }
        return true;
    }
    if (obj["command"] == "clear-processing") {
        MountainProcessServer* srvr = static_cast<MountainProcessServer*>(server());
        srvr->clearProcessing();
        return true;
    }
    writeMessage("BAD COMMAND");
    return true;
}

MountainProcessServer::MountainProcessServer(QObject* parent)
    : LocalServer::Server(parent)
{
}

MountainProcessServer::~MountainProcessServer()
{
}

void MountainProcessServer::distributeLogMessage(const QJsonObject& msg)
{
    foreach (LocalServer::Client* client, m_listeners) {
        client->writeMessage(QJsonDocument(msg).toJson());
    }
}

void MountainProcessServer::registerLogListener(LocalServer::Client* listener)
{
    m_listeners.append(listener);
}

void MountainProcessServer::unregisterLogListener(LocalServer::Client* listener)
{
    m_listeners.removeOne(listener);
}

QJsonObject MountainProcessServer::state()
{
    QJsonObject ret;

    ret["is_running"] = m_is_running;
    QJsonObject scripts;
    QJsonObject processes;
    QStringList keys = m_pripts.keys();
    foreach (QString key, keys) {
        if (m_pripts[key].prtype == ScriptType)
            scripts[key] = pript_struct_to_obj(m_pripts[key], AbbreviatedRecord);
        else
            processes[key] = pript_struct_to_obj(m_pripts[key], AbbreviatedRecord);
    }
    ret["scripts"] = scripts;
    ret["processes"] = processes;
    return ret;
}

QJsonArray MountainProcessServer::log()
{
    return m_log;
}

void MountainProcessServer::contignousLog()
{
}

bool MountainProcessServer::queueScript(const MPDaemonPript& script)
{
    if (m_pripts.contains(script.id)) {
        writeLogRecord("error",
            "message",
            QString("Unable to queue script. "
                    "Process or script with id %1 already exists: ").arg(script.id));
        return false;
    }
    writeLogRecord("queue-script", "pript_id", script.id);
    m_pripts[script.id] = script;
    write_pript_file(script);
    return true;
}

bool MountainProcessServer::queueProcess(const MPDaemonPript& process)
{
    if (m_pripts.contains(process.id)) {
        writeLogRecord("error",
            "message",
            QString("Unable to queue process. "
                    "Process or script with id %1 already exists: ").arg(process.id));
        return false;
    }
    writeLogRecord("queue-process", process.id);
    m_pripts[process.id] = process;
    write_pript_file(process);
    return true;
}

bool MountainProcessServer::clearProcessing()
{
    QStringList keys = m_pripts.keys();
    foreach (QString key, keys) {
        stop_or_remove_pript(key);
    }
    return true;
}

bool MountainProcessServer::start()
{
    if (m_is_running)
        return false;

    if (!startServer())
        return false;

    signal(SIGINT, sighandler);
    signal(SIGTERM, sighandler);

    m_is_running = true;

    writeLogRecord("start-daemon");
    QTimer timer;
    connect(&timer, &QTimer::timeout, [this]() {
        iterate();
    });
    timer.start(100);
    qApp->exec();
    m_is_running = false;
    writeLogRecord("stop-daemon");
    signal(SIGINT, SIG_DFL);
    signal(SIGTERM, SIG_DFL);
    stopServer();
    return true;
}

bool MountainProcessServer::stop()
{
    m_is_running = false;
    qApp->quit();
    return true;
}

void MountainProcessServer::setLogPath(const QString& lp)
{
    if (m_logPath == lp)
        return;
    m_logPath = lp;
}

void MountainProcessServer::setTotalResourcesAvailable(ProcessResources PR)
{
    m_total_resources_available = PR;
}

LocalServer::Client* MountainProcessServer::createClient(QLocalSocket* sock)
{
    MountainProcessServerClient* client = new MountainProcessServerClient(sock, this);
    return client;
}

void MountainProcessServer::clientAboutToBeDestroyed(LocalServer::Client* client)
{
    unregisterLogListener(client);
}

bool MountainProcessServer::startServer()
{
    if (!acquireServer()) {
        printf("Another daemon seems to be running. Closing.\n");
        return false;
    }
    return true;
}

void MountainProcessServer::stopServer()
{
    releaseServer();
}

bool MountainProcessServer::acquireServer()
{
    if (!shm)
        shm = new QSharedMemory(shmName(), this);
    else if (shm->isAttached())
        shm->detach();
    // algorithm:
    // create a shared memory segment
    // if successful, there is no active server, so we become one
    //   - write your own PID into the segment
    // otherwise try to attach to the segment
    //   - check the pid stored in the segment
    //   - see if process with that pid exists
    //   - if yes, bail out
    //   - if not, overwrite the pid with your own
    // repeat the process until either create or attach is successfull

    while (true) {
        if (shm->create(sizeof(MountainProcessDescriptor))) {
            shm->lock(); // there is potentially a race condition here -> someone might have locked the segment
            // before we did and written its own pid there
            // we assume that by default memory is zeroed (it seems to be on Linux)
            // so we can check the version
            MountainProcessDescriptor* desc = reinterpret_cast<MountainProcessDescriptor*>(shm->data());

            // on Linux the memory seems to be zeroed by default
            if (desc->version != 0) {
                // someone has hijacked our segment
                shm->unlock();
                shm->detach();
                continue; // try again
            }
            desc->version = 1;
            desc->pid = (pid_t)QCoreApplication::applicationPid();
            acquireSocket();
            shm->unlock();
            return true;
        }
        if (shm->attach()) {
            shm->lock();
            MountainProcessDescriptor* desc = reinterpret_cast<MountainProcessDescriptor*>(shm->data());
            // on Linux the memory seems to be zeroed by default
            if (desc->version != 0) {
                if (kill(desc->pid, 0) == 0) {
                    // pid exists, server is likely active
                    shm->unlock();
                    shm->detach();
                    return false;
                }
            }
            // server has crashed or we have hijacked the segment
            desc->version = 1;
            desc->pid = (pid_t)QCoreApplication::applicationPid();
            acquireSocket();
            shm->unlock();
            return true;
        }
    }
}

bool MountainProcessServer::releaseServer()
{
    if (shm && shm->isAttached()) {
        shm->lock();
        MountainProcessDescriptor* desc = reinterpret_cast<MountainProcessDescriptor*>(shm->data());
        // just to be sure we're the reigning server:
        if (desc->pid == (pid_t)QCoreApplication::applicationPid()) {
            // to avoid a race condition (released the server but not ended the process just yet), empty the block
            desc->version = 0;
            desc->pid = 0;
            releaseSocket();
        }
        shm->unlock();
        shm->detach();
        return true;
    }
    return false;
}

bool MountainProcessServer::acquireSocket()
{
    return listen(socketName());
}

bool MountainProcessServer::releaseSocket()
{
    shutdown();
    return true;
}

void MountainProcessServer::iterate()
{
    // hack by jfm to temporarily implement mp-list-daemons
    {
        QString daemon_id = qgetenv("MP_DAEMON_ID");
        QSettings settings(QSettings::UserScope, "Magland", "MountainLab");
        QStringList list = settings.value("mp-list-daemons-candidates").toStringList();
        if (!list.contains(daemon_id)) {
            list.append(daemon_id);
            settings.setValue("mp-list-daemons-candidates", list);
        }
    }

    stop_orphan_processes_and_scripts();
    handle_scripts();
    handle_processes();
}

void MountainProcessServer::writeLogRecord(QString record_type, QString key1, QVariant val1, QString key2, QVariant val2, QString key3, QVariant val3)
{
    QVariantMap map;
    if (!key1.isEmpty()) {
        map[key1] = val1;
    }
    if (!key2.isEmpty()) {
        map[key2] = val2;
    }
    if (!key3.isEmpty()) {
        map[key3] = val3;
    }
    QJsonObject obj = variantmap_to_json_obj(map);
    writeLogRecord(record_type, obj);
}

void MountainProcessServer::writeLogRecord(QString record_type, const QJsonObject& obj)
{
    QJsonObject X;
    X["record_type"] = record_type;
    X["timestamp"] = QDateTime::currentDateTime().toString("yyyy-MM-dd|hh:mm:ss.zzz");
    X["data"] = obj;
    QString line = QJsonDocument(X).toJson(QJsonDocument::Compact);
    while (m_log.size() >= 10)
        m_log.pop_front();
    m_log.push_back(X);
    distributeLogMessage(X);
}

void MountainProcessServer::write_pript_file(const MPDaemonPript& P)
{
    if (m_logPath.isEmpty())
        return;
    if (P.id.isEmpty())
        return;
    QString fname;
    if (P.prtype == ScriptType) {
        fname = QString("%1/scripts/%2.json").arg(m_logPath).arg(P.id);
    }
    else if (P.prtype == ProcessType) {
        fname = QString("%1/processes/%2.json").arg(m_logPath).arg(P.id);
    }
    QJsonObject obj = pript_struct_to_obj(P, RuntimeRecord);
    QString json = QJsonDocument(obj).toJson();
    TextFile::write(fname, json);
}

bool MountainProcessServer::stop_or_remove_pript(const QString& key)
{
    if (!m_pripts.contains(key))
        return false;
    MPDaemonPript* PP = &m_pripts[key];
    if ((PP->is_running)) {
        if (PP->qprocess) {
            qWarning() << "Terminating qprocess: " + key;
            //PP->qprocess->terminate(); // I think it's okay to terminate a process. It won't cause this program to crash.
            if (PP->parent_pid) {
                //if we can, kill the parent pid instead so that the qprocess will exit gracefully when it discovers it has no parent
                kill_process_and_children(PP->parent_pid);
            }
            else {
                kill_process_and_children(PP->qprocess);
                delete PP->qprocess;
            }
        }
        finish_and_finalize(*PP);
        m_pripts.remove(key);
        if (PP->prtype == ScriptType)
            writeLogRecord("stop-script", "pript_id", key, "reason", "requested");
        else
            writeLogRecord("stop-process", "pript_id", key, "reason", "requested");
    }
    else {
        m_pripts.remove(key);
        if (PP->prtype == ScriptType)
            writeLogRecord("unqueue-script", "pript_id", key, "reason", "requested");
        else
            writeLogRecord("unqueue-process", "pript_id", key, "reason", "requested");
    }
    return true;
}

void MountainProcessServer::finish_and_finalize(MPDaemonPript& P)
{
    P.is_finished = true;
    P.is_running = false;
    P.timestamp_finished = QDateTime::currentDateTime();
    if (!P.stdout_fname.isEmpty()) {
        P.runtime_results["stdout"] = TextFile::read(P.stdout_fname);
    }
    write_pript_file(P);
}

void MountainProcessServer::stop_orphan_processes_and_scripts()
{
    QStringList keys = m_pripts.keys();
    foreach (QString key, keys) {
        if (!m_pripts[key].is_finished) {
            if ((m_pripts[key].parent_pid) && (!pidExists(m_pripts[key].parent_pid))) {
                debug_log(__FUNCTION__, __FILE__, __LINE__);
                if (m_pripts[key].qprocess) {
                    if (m_pripts[key].prtype == ScriptType) {
                        writeLogRecord("stop-script", "pript_id", key, "reason", "orphan", "parent_pid", m_pripts[key].parent_pid);
                        qWarning() << "Terminating orphan script qprocess: " + key;
                    }
                    else {
                        writeLogRecord("stop-process", "pript_id", key, "reason", "orphan", "parent_pid", m_pripts[key].parent_pid);
                        qWarning() << "Terminating orphan process qprocess: " + key;
                    }

                    m_pripts[key].qprocess->disconnect(); //so we don't go into the finished slot
                    //m_pripts[key].qprocess->terminate();
                    kill_process_and_children(m_pripts[key].qprocess);
                    delete m_pripts[key].qprocess;
                    finish_and_finalize(m_pripts[key]);
                    m_pripts.remove(key);
                }
                else {
                    if (m_pripts[key].prtype == ScriptType) {
                        writeLogRecord("unqueue-script", "pript_id", key, "reason", "orphan", "parent_pid", m_pripts[key].parent_pid);
                        qWarning() << "Removing orphan script: " + key + " " + m_pripts[key].script_paths.value(0);
                    }
                    else {
                        writeLogRecord("unqueue-process", "pript_id", key, "reason", "orphan", "parent_pid", m_pripts[key].parent_pid);
                        qWarning() << "Removing orphan process: " + key + " " + m_pripts[key].processor_name;
                    }
                    finish_and_finalize(m_pripts[key]);
                    m_pripts.remove(key);
                }
            }
        }
    }
}

bool MountainProcessServer::handle_scripts()
{
    int max_simultaneous_scripts = 100;
    if (num_running_scripts() < max_simultaneous_scripts) {
        int old_num_pending_scripts = num_pending_scripts();
        if (num_pending_scripts() > 0) {
            if (launch_next_script()) {
                printf("%d scripts running.\n", num_running_scripts());
            }
            else {
                if (num_pending_scripts() == old_num_pending_scripts) {
                    qCritical() << "Failed to launch_next_script and the number of pending scripts has not decreased (unexpected). This has potential for infinite loop. So we are aborting.";
                    abort();
                }
            }
        }
    }
    return true;
}

bool MountainProcessServer::handle_processes()
{
    ProcessManager* PM = ProcessManager::globalInstance();
    PM->reloadProcessors();

    ProcessResources pr_available = compute_process_resources_available();
    QStringList keys = m_pripts.keys();
    foreach (QString key, keys) {
        if (m_pripts[key].prtype == ProcessType) {
            if ((!m_pripts[key].is_running) && (!m_pripts[key].is_finished)) {
                ProcessResources pr_needed = compute_process_resources_needed(m_pripts[key]);
                if (is_at_most(pr_needed, pr_available, m_total_resources_available)) {
                    if (process_parameters_are_okay(key)) {
                        if (okay_to_run_process(key)) { //check whether there are io file conflicts at the moment
                            if (launch_pript(key)) {
                                pr_available.num_threads -= m_pripts[key].runtime_opts.num_threads_allotted;
                                pr_available.memory_gb -= m_pripts[key].runtime_opts.memory_gb_allotted;
                                pr_available.num_processes -= 1;
                            }
                        }
                    }
                    else {
                        writeLogRecord("unqueue-process", "pript_id", key, "reason", "processor not found or parameters are incorrect.");
                        m_pripts.remove(key);
                    }
                }
            }
        }
    }
    return true;
}

int MountainProcessServer::num_running_pripts(PriptType prtype) const
{
    int ret = 0;
    QStringList keys = m_pripts.keys();
    foreach (QString key, keys) {
        if (m_pripts[key].is_running) {
            if (m_pripts[key].prtype == prtype)
                ret++;
        }
    }
    return ret;
}

int MountainProcessServer::num_pending_pripts(PriptType prtype) const
{
    int ret = 0;
    QStringList keys = m_pripts.keys();
    foreach (QString key, keys) {
        if ((!m_pripts[key].is_running) && (!m_pripts[key].is_finished)) {
            if (m_pripts[key].prtype == prtype)
                ret++;
        }
    }
    return ret;
}

bool MountainProcessServer::launch_next_script()
{
    QStringList keys = m_pripts.keys();
    foreach (QString key, keys) {
        if (m_pripts[key].prtype == ScriptType) {
            if ((!m_pripts[key].is_running) && (!m_pripts[key].is_finished)) {
                if (launch_pript(key)) {
                    return true;
                }
            }
        }
    }
    return false;
}

bool MountainProcessServer::launch_pript(QString pript_id)
{
    MPDaemonPript* S;
    if (!m_pripts.contains(pript_id)) {
        qCritical() << "Unexpected problem. No pript exists with id: " + pript_id;
        abort();
        return false;
    }
    S = &m_pripts[pript_id];
    if (S->is_running) {
        qCritical() << "Unexpected problem. Pript is already running: " + pript_id;
        abort();
        return false;
    }
    if (S->is_finished) {
        qCritical() << "Unexpected problem. Pript is already finished: " + pript_id;
        abort();
        return false;
    }

    QString exe = qApp->applicationFilePath();
    QStringList args;

    if (S->prtype == ScriptType) {
        debug_log(__FUNCTION__, __FILE__, __LINE__);
        args << "run-script";
        if (!S->output_fname.isEmpty()) {
            args << "--_script_output=" + S->output_fname;
        }
        args << "--_working_path=" + S->working_path;
        for (int ii = 0; ii < S->script_paths.count(); ii++) {
            QString fname = S->script_paths[ii];
            if (!QFile::exists(fname)) {
                QString message = "Script file does not exist: " + fname;
                qWarning() << message;
                writeLogRecord("error", "message", message);
                writeLogRecord("unqueue-script", "pript_id", pript_id, "reason", message);
                m_pripts.remove(pript_id);
                return false;
            }
            if (MLUtil::computeSha1SumOfFile(fname) != S->script_path_checksums.value(ii)) {
                QString message = "Script checksums do not match. Script file has changed since queueing: " + fname + " Not launching process: " + pript_id;
                qWarning() << message;
                qWarning() << MLUtil::computeSha1SumOfFile(fname) << "<>" << S->script_path_checksums.value(ii);
                writeLogRecord("error", "message", message);
                writeLogRecord("unqueue-script", "pript_id", pript_id, "reason", "Script file has changed: " + fname);
                m_pripts.remove(pript_id);
                return false;
            }
            args << fname;
        }
        QJsonObject parameters = variantmap_to_json_obj(S->parameters);
        QString parameters_json = QJsonDocument(parameters).toJson();
        QString par_fname = CacheManager::globalInstance()->makeLocalFile(S->id + ".par", CacheManager::ShortTerm);
        TextFile::write(par_fname, parameters_json);
        CacheManager::globalInstance()->setTemporaryFileDuration(par_fname, 600);
        args << par_fname;
    }
    else if (S->prtype == ProcessType) {
        debug_log(__FUNCTION__, __FILE__, __LINE__);
        args << "run-process";
        args << S->processor_name;
        if (!S->output_fname.isEmpty())
            args << "--_process_output=" + S->output_fname;
        if (S->preserve_tempdir)
            args << "--_preserve_tempdir";
        args << QString("--_parent_pid=%1").arg(S->parent_pid);
        QStringList pkeys = S->parameters.keys();
        foreach (QString pkey, pkeys) {
            QStringList list = MLUtil::toStringList(S->parameters[pkey]);
            foreach (QString str, list) {
                args << QString("--%1=%2").arg(pkey).arg(str);
            }
        }
        if (S->RPR.request_num_threads) {
            args << QString("--_request_num_threads=%1").arg(S->RPR.request_num_threads);
        }
        //S->runtime_opts.num_threads_allotted = S->num_threads_requested;
        //S->runtime_opts.memory_gb_allotted = S->memory_gb_requested;
    }
    if (S->force_run) {
        args << "--_force_run";
    }
    args << "--_working_path=" + S->working_path;
    debug_log(__FUNCTION__, __FILE__, __LINE__);
    QProcess* qprocess = new QProcess;
    qprocess->setProperty("pript_id", pript_id);
    qprocess->setProcessChannelMode(QProcess::MergedChannels);
    QObject::connect(qprocess, SIGNAL(readyRead()), this, SLOT(slot_qprocess_output()));
    if (S->prtype == ScriptType) {
        printf("   Launching script %s: ", pript_id.toLatin1().data());
        writeLogRecord("start-script", "pript_id", pript_id);
        foreach (QString fname, S->script_paths) {
            QString str = QFileInfo(fname).fileName();
            printf("%s ", str.toLatin1().data());
        }
        printf("\n");
    }
    else {
        printf("   Launching process %s %s: ", S->processor_name.toLatin1().data(), pript_id.toLatin1().data());
        writeLogRecord("start-process", "pript_id", pript_id);
        QString cmd = args.join(" ");
        printf("%s\n", cmd.toLatin1().data());
    }

    QObject::connect(qprocess, SIGNAL(finished(int)), this, SLOT(slot_pript_qprocess_finished()));

    debug_log(__FUNCTION__, __FILE__, __LINE__);

    if (S->parent_pid) {
        MPDaemon::start_bash_command_and_kill_when_pid_is_gone(qprocess, exe, args, S->parent_pid);
    }
    else {
        qprocess->start(exe, args);
        /*
        qprocess->setProgram(exe);
        qprocess->setArguments(args);
        QProcessManager *manager=ObjectRegistry::getObject<QProcessManager>();
        if (!manager) {
            qCritical() << "No QProcessManager object found in object registry. Aborting.";
            abort();
        }
        manager->start(qprocess);
        */
    }
    if (qprocess->waitForStarted()) {
        if (S->prtype == ScriptType) {
            writeLogRecord("started-script", "pript_id", pript_id, "pid", (int)qprocess->processId());
        }
        else {
            writeLogRecord("started-process", "pript_id", pript_id, "pid", (int)qprocess->processId());
        }
        S->qprocess = qprocess;
        if (!S->stdout_fname.isEmpty()) {
            S->stdout_file = new QFile(S->stdout_fname);
            if (!S->stdout_file->open(QFile::WriteOnly)) {
                qCritical() << "Unable to open stdout file for writing: " + S->stdout_fname;
                writeLogRecord("error", "message", "Unable to open stdout file for writing: " + S->stdout_fname);
                delete S->stdout_file;
                S->stdout_file = 0;
            }
        }
        S->is_running = true;
        S->timestamp_started = QDateTime::currentDateTime();
        write_pript_file(*S);
        return true;
    }
    else {
        debug_log(__FUNCTION__, __FILE__, __LINE__);
        if (S->prtype == ScriptType) {
            writeLogRecord("stop-script", "pript_id", pript_id, "reason", "Unable to start script.");
            qCritical() << "Unable to start script: " + S->id;
        }
        else {
            writeLogRecord("stop-process", "pript_id", pript_id, "reason", "Unable to start process.");
            qCritical() << "Unable to start process: " + S->processor_name + " " + S->id;
        }
        qprocess->disconnect();
        delete qprocess;
        write_pript_file(*S);
        m_pripts.remove(pript_id);
        return false;
    }
}

ProcessResources MountainProcessServer::compute_process_resources_available() const
{
    ProcessResources ret = m_total_resources_available;
    QStringList keys = m_pripts.keys();
    foreach (QString key, keys) {
        if (m_pripts[key].prtype == ProcessType) {
            if (m_pripts[key].is_running) {
                ProcessRuntimeOpts rtopts = m_pripts[key].runtime_opts;
                ret.num_threads -= rtopts.num_threads_allotted;
                ret.memory_gb -= rtopts.memory_gb_allotted;
                ret.num_processes -= 1;
            }
        }
    }
    return ret;
}

ProcessResources MountainProcessServer::compute_process_resources_needed(MPDaemonPript P) const
{
    ProcessResources ret;
    ret.num_threads = P.RPR.request_num_threads;
    if (ret.num_threads < 1)
        ret.num_threads = 1;
    //ret.memory_gb = P.memory_gb_requested;
    ret.memory_gb = 1;
    ret.num_processes = 1;
    return ret;
}

bool MountainProcessServer::process_parameters_are_okay(const QString& key) const
{
    //check that the processor is registered and that the parameters are okay
    ProcessManager* PM = ProcessManager::globalInstance();
    if (!m_pripts.contains(key))
        return false;
    MPDaemonPript P0 = m_pripts[key];
    MLProcessor MLP = PM->processor(P0.processor_name);
    if (QJsonDocument(MLP.spec).toJson() != QJsonDocument(P0.processor_spec).toJson()) {
        qWarning() << "Processor spec does not match. Perhaps the processor has changed and the daemon needs to be restarted.";
        return false;
    }
    if (MLP.name != P0.processor_name) {
        qWarning() << "Unable to find processor **: " + P0.processor_name;
        return false;
    }
    if (!PM->checkParameters(P0.processor_name, P0.parameters)) {
        qWarning() << "Failure in checkParameters: " + P0.processor_name;
        return false;
    }
    return true;
}

bool MountainProcessServer::okay_to_run_process(const QString& key) const
{
    //next we check that there are no running processes where the input or output files conflict with the proposed input or output files
    //but note that it is okay for the same input file to be used in more than one simultaneous process
    QSet<QString> pending_input_paths;
    QSet<QString> pending_output_paths;
    QStringList pripts_keys = m_pripts.keys();
    foreach (QString key0, pripts_keys) {
        if (m_pripts[key0].prtype == ProcessType) {
            if (m_pripts[key0].is_running) {
                QStringList input_paths = get_input_paths(m_pripts[key0]);
                foreach (QString path0, input_paths) {
                    pending_input_paths.insert(path0);
                }
                QStringList output_paths = get_output_paths(m_pripts[key0]);
                foreach (QString path0, output_paths) {
                    pending_output_paths.insert(path0);
                }
            }
        }
    }
    {
        QStringList proposed_input_paths = get_input_paths(m_pripts[key]);
        QStringList proposed_output_paths = get_output_paths(m_pripts[key]);
        foreach (QString path0, proposed_input_paths) {
            if (pending_output_paths.contains(path0))
                return false;
        }
        foreach (QString path0, proposed_output_paths) {
            if (pending_input_paths.contains(path0))
                return false;
            if (pending_output_paths.contains(path0))
                return false;
        }
    }
    return true;
}

QStringList MountainProcessServer::get_input_paths(MPDaemonPript P) const
{
    QStringList ret;
    if (P.prtype != ProcessType)
        return ret;
    ProcessManager* PM = ProcessManager::globalInstance();
    MLProcessor MLP = PM->processor(P.processor_name);
    QStringList pnames = MLP.inputs.keys();
    foreach (QString pname, pnames) {
        QStringList paths0 = MLUtil::toStringList(P.parameters[pname]); //is this right?
        foreach (QString path0, paths0) {
            if (!path0.isEmpty()) {
                ret << path0;
            }
        }
    }
    return ret;
}

QStringList MountainProcessServer::get_output_paths(MPDaemonPript P) const
{
    QStringList ret;
    if (P.prtype != ProcessType)
        return ret;
    ProcessManager* PM = ProcessManager::globalInstance();
    MLProcessor MLP = PM->processor(P.processor_name);
    QStringList pnames = MLP.outputs.keys();
    foreach (QString pname, pnames) {
        QStringList paths0 = MLUtil::toStringList(P.parameters[pname]); //is this right?
        foreach (QString path0, paths0) {
            if (!path0.isEmpty()) {
                ret << path0;
            }
        }
    }
    return ret;
}

void MountainProcessServer::slot_pript_qprocess_finished()
{
    debug_log(__FUNCTION__, __FILE__, __LINE__);

    QProcess* P = qobject_cast<QProcess*>(sender());
    if (!P)
        return;
    QString pript_id = P->property("pript_id").toString();
    MPDaemonPript* S;
    if (m_pripts.contains(pript_id)) {
        S = &m_pripts[pript_id];
    }
    else {
        writeLogRecord("error", "message", "Unexpected problem in slot_pript_qprocess_finished. Unable to find script or process with id: " + pript_id);
        qCritical() << "Unexpected problem in slot_pript_qprocess_finished. Unable to find script or process with id: " + pript_id;
        return;
    }
    if (!S->output_fname.isEmpty()) {
        debug_log(__FUNCTION__, __FILE__, __LINE__);
        QString runtime_results_json = TextFile::read(S->output_fname);
        if (runtime_results_json.isEmpty()) {
            S->success = false;
            S->error = "Could not read results file ****: " + S->output_fname;
        }
        else {
            QJsonParseError error;
            S->runtime_results = QJsonDocument::fromJson(runtime_results_json.toLatin1(), &error).object();
            if (error.error != QJsonParseError::NoError) {
                S->success = false;
                S->error = "Error parsing json of runtime results.";
            }
            else {
                S->success = S->runtime_results["success"].toBool();
                S->error = S->runtime_results["error"].toString();
            }
        }
    }
    else {
        S->success = true;
    }
    finish_and_finalize(*S);

    QJsonObject obj0;
    obj0["pript_id"] = pript_id;
    obj0["reason"] = "finished";
    obj0["success"] = S->success;
    obj0["error"] = S->error;
    if (S->prtype == ScriptType) {
        writeLogRecord("stop-script", obj0);
        printf("  Script %s finished ", pript_id.toLatin1().data());
    }
    else {
        writeLogRecord("stop-process", obj0);
        printf("  Process %s %s finished ", S->processor_name.toLatin1().data(), pript_id.toLatin1().data());
    }
    if (S->success)
        printf("successfully\n");
    else
        printf("with error: %s\n", S->error.toLatin1().data());
    if (S->qprocess) {
        if (S->qprocess->state() == QProcess::Running) {
            if (S->qprocess->waitForFinished(1000)) {
                delete S->qprocess;
            }
            else {
                writeLogRecord("error", "pript_id", pript_id, "message", "Process did not finish after waiting even though we are in the slot for finished!!");
                qCritical() << "Process did not finish after waiting even though we are in the slot for finished!!";
            }
        }
        S->qprocess = 0;
    }
    if (S->stdout_file) {
        if (S->stdout_file->isOpen()) {
            S->stdout_file->close();
        }
        delete S->stdout_file;
        S->stdout_file = 0;
    }
}

void MountainProcessServer::slot_qprocess_output()
{
    QProcess* P = qobject_cast<QProcess*>(sender());
    if (!P)
        return;
    QByteArray str = P->readAll();
    QString pript_id = P->property("pript_id").toString();
    if ((m_pripts.contains(pript_id)) && (m_pripts[pript_id].stdout_file)) {
        if (m_pripts[pript_id].stdout_file->isOpen()) {
            m_pripts[pript_id].stdout_file->write(str);
            m_pripts[pript_id].stdout_file->flush();
        }
    }
    else {
        printf("%s", str.data());
    }
}

void MPDaemon::wait(qint64 msec)
{
    usleep(msec * 1000);
}

bool MPDaemon::waitForFileToAppear(QString fname, qint64 timeout_ms, bool remove_on_appear, qint64 parent_pid, QString stdout_fname)
{
    debug_log(__FUNCTION__, __FILE__, __LINE__);

    QTime timer;
    timer.start();
    QFile stdout_file(stdout_fname);
    bool failed_to_open_stdout_file = false;

    /*
    QString i_am_alive_fname;
    if (parent_pid) {
        i_am_alive_fname=CacheManager::globalInstance()->makeLocalFile(QString("i_am_alive.%1.txt").arg(parent_pid));
    }
    QTime timer_i_am_alive; timer_i_am_alive.start();
    */

    QTime debug_timer;
    debug_timer.start();
    while (1) {
        wait(200);

        bool terminate_file_exists = QFile::exists(fname); //do this before we check other things, like the stdout

        if ((timeout_ms >= 0) && (timer.elapsed() > timeout_ms))
            return false;
        if ((parent_pid) && (!MPDaemon::pidExists(parent_pid))) {
            qWarning() << "Exiting waitForFileToAppear because parent process is gone.";
            break;
        }
        if ((!stdout_fname.isEmpty()) && (!failed_to_open_stdout_file)) {
            if (QFile::exists(stdout_fname)) {
                if (!stdout_file.isOpen()) {
                    if (!stdout_file.open(QFile::ReadOnly)) {
                        qCritical() << "Unable to open stdout file for reading: " + stdout_fname;
                        failed_to_open_stdout_file = true;
                    }
                }
                if (stdout_file.isOpen()) {
                    QByteArray str = stdout_file.readAll();
                    if (!str.isEmpty()) {
                        printf("%s", str.data());
                    }
                }
            }
        }
        if (terminate_file_exists)
            break;

        /*
        if (debug_timer.elapsed() > 20000) {
            qWarning() << QString("Still waiting for file to appear after %1 sec: %2").arg(timer.elapsed() * 1.0 / 1000).arg(fname);
            debug_timer.restart();
        }
        */
    }
    if (stdout_file.isOpen())
        stdout_file.close();
    if (remove_on_appear) {
        QFile::remove(fname);
    }
    return true;
}

QString MPDaemon::daemonPath()
{
    // Witold, it turns out we don't want a separate path for each daemon id. This causes problems.
    QString ret = CacheManager::globalInstance()->localTempPath() + "/mpdaemon";
    MLUtil::mkdirIfNeeded(ret);
    MLUtil::mkdirIfNeeded(ret + "/completed_processes");
    QFile::Permissions perm = QFileDevice::ReadUser | QFileDevice::WriteUser | QFileDevice::ExeUser | QFileDevice::ReadGroup | QFileDevice::WriteGroup | QFileDevice::ExeGroup | QFileDevice::ReadOther | QFileDevice::WriteOther | QFileDevice::ExeOther;
    QFile::setPermissions(ret, perm);
    QFile::setPermissions(ret + "/completed_processes", perm);
    return ret;
}

bool MPDaemon::waitForFinishedAndWriteOutput(QProcess* P, int parent_pid)
{
    debug_log(__FUNCTION__, __FILE__, __LINE__);

    P->waitForStarted();
    while (P->state() == QProcess::Running) {
        P->waitForReadyRead(100);
        QByteArray str = P->readAll();
        if (str.count() > 0) {
            printf("%s", str.data());
        }
        qApp->processEvents();
        if (parent_pid) {
            if (!pidExists(parent_pid)) {
                qWarning() << "Parent pid does not exist. Terminating process...." << P->pid();
                kill_process_and_children(P);
                qWarning() << "Process killed.";
                return false;
            }
        }
    }
    {
        P->waitForReadyRead();
        QByteArray str = P->readAll();
        if (str.count() > 0) {
            printf("%s", str.data());
        }
    }
    return (P->state() != QProcess::Running);
}

bool MPDaemon::pidExists(qint64 pid)
{
    return (kill(pid, 0) == 0);
}

void MPDaemon::start_bash_command_and_kill_when_pid_is_gone(QProcess* qprocess, QString exe_command, int pid)
{
    QString bash_script_fname = CacheManager::globalInstance()->makeLocalFile();

    QString script;
    script += "#!/bin/bash\n\n";
    script += exe_command + " &\n";
    script += "cmdpid=$!\n"; //get the pid of the exe_command
    script += "trap \"kill $cmdpid; exit 255;\" SIGINT SIGTERM\n"; //capture the terminate signal and pass it on
    script += QString("while kill -0 %1 >/dev/null 2>&1; do\n").arg(pid); //while the (parent) pid still exists
    script += "    if kill -0 $cmdpid > /dev/null 2>&1; then\n";
    script += "        sleep 1;\n"; //we are still going
    script += "    else\n"; //else the exe process is done
    script += "        wait $cmdpid\n"; //get the return code for the process that has already completed
    script += "        exit $?\n";
    script += "    fi\n";
    script += "done ;\n";
    script += "kill $cmdpid\n"; //the parent pid is gone
    script += "exit 255\n"; //return error exit code

    TextFile::write(bash_script_fname, script);
    QProcessManager* manager = ObjectRegistry::getObject<QProcessManager>();
    if (!manager) {
        qCritical() << "No QProcessManager object found in object registry. Aborting.";
        abort();
    }
    qprocess->setProgram("/bin/bash");
    qprocess->setArguments(QStringList(bash_script_fname));
    QSharedPointer<QProcess> pp = manager->start(qprocess);

    CacheManager::globalInstance()->setTemporaryFileExpirePid(bash_script_fname, pp->processId());

    //qprocess->start("/bin/bash", QStringList(bash_script_fname));
}

void MPDaemon::start_bash_command_and_kill_when_pid_is_gone(QProcess* qprocess, QString exe, QStringList args, int pid)
{
    QString exe_command = exe + " " + args.join(" ");
    start_bash_command_and_kill_when_pid_is_gone(qprocess, exe_command, pid);
}
