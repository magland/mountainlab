/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
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
#include "textfile.h"
#include <QDebug>
#include "cachemanager.h"
#include "unistd.h" //for usleep
#include <sys/stat.h> //for mkfifo
#include "processmanager.h"
#include "mlutils.h"

class MPDaemonPrivate {
public:
    MPDaemon* q;
    bool m_is_running;
    QFileSystemWatcher m_watcher;
    QMap<QString, MPDaemonPript> m_pripts;
    ProcessResources m_total_resources_available;
    QString m_log_path;

    void process_command(QJsonObject obj);
    void writeLogRecord(QString record_type, QString key1 = "", QVariant val1 = QVariant(), QString key2 = "", QVariant val2 = QVariant(), QString key3 = "", QVariant val3 = QVariant());
    void writeLogRecord(QString record_type, const QJsonObject& obj);

    bool handle_scripts();
    bool handle_processes();
    bool process_parameters_are_okay(const QString& key);
    bool okay_to_run_process(const QString& key);
    QStringList get_input_paths(MPDaemonPript P);
    QStringList get_output_paths(MPDaemonPript P);
    bool write_running_file();
    void write_daemon_state();
    bool stop_or_remove_pript(const QString& key);
    void write_pript_file(const MPDaemonPript& P);
    void finish_and_finalize(MPDaemonPript& P);
    int num_running_scripts()
    {
        return num_running_pripts(ScriptType);
    }
    int num_running_processes()
    {
        return num_running_pripts(ProcessType);
    }
    int num_pending_scripts()
    {
        return num_pending_pripts(ScriptType);
    }
    int num_pending_processes()
    {
        return num_pending_pripts(ProcessType);
    }

    void stop_orphan_processes_and_scripts();

    /////////////////////////////////
    int num_running_pripts(PriptType prtype);
    int num_pending_pripts(PriptType prtype);
    bool launch_next_script();
    bool launch_pript(QString id);
    ProcessResources compute_process_resources_available();
    ProcessResources compute_process_resources_needed(MPDaemonPript P);
};

void append_line_to_file(QString fname, QString line)
{
    QFile ff(fname);
    if (ff.open(QFile::Append | QFile::WriteOnly)) {
        ff.write(line.toLatin1());
        ff.write(QByteArray("\n"));
        ff.close();
    } else {
        static bool reported = false;
        if (!reported) {
            reported = true;
            qCritical() << "Unable to write to log file" << fname;
        }
    }
}

void debug_log(const char* function, const char* file, int line)
{
    QString fname = CacheManager::globalInstance()->localTempPath() + "/mpdaemon_debug.log";
    QString line0 = QString("%1: %2 %3:%4").arg(QDateTime::currentDateTime().toString("yy-MM-dd:hh:mm:ss.zzz")).arg(function).arg(file).arg(line);
    line0 += " ARGS: ";
    foreach(QString arg, qApp->arguments())
    {
        line0 += arg + " ";
    }
    append_line_to_file(fname, line0);
}

MPDaemon::MPDaemon()
{
    debug_log(__FUNCTION__, __FILE__, __LINE__);

    d = new MPDaemonPrivate;
    d->q = this;
    d->m_is_running = false;

    d->m_total_resources_available.num_threads = 12;
    d->m_total_resources_available.memory_gb = 8;

    connect(&d->m_watcher, SIGNAL(directoryChanged(QString)), this, SLOT(slot_commands_directory_changed()));
}

MPDaemon::~MPDaemon()
{
    debug_log(__FUNCTION__, __FILE__, __LINE__);

    foreach(MPDaemonPript P, d->m_pripts)
    {
        if (P.qprocess) {
            if (P.qprocess->state() == QProcess::Running)
                P.qprocess->terminate(); // I think it's okay to terminate a process. It won't cause this program to crash.
        }
    }

    delete d;
}

void MPDaemon::setTotalResourcesAvailable(ProcessResources PR)
{
    d->m_total_resources_available = PR;
}

void MPDaemon::setLogPath(const QString& path)
{
    d->m_log_path = path;
    mkdir_if_doesnt_exist(path);
    mkdir_if_doesnt_exist(path + "/scripts");
    mkdir_if_doesnt_exist(path + "/processes");
}

bool MPDaemon::run()
{
    debug_log(__FUNCTION__, __FILE__, __LINE__);

    if (d->m_is_running) {
        qWarning() << "Unexpected problem in MPDaemon:run(). Daemon is already running.";
        return false;
    }

    if (!d->write_running_file()) { //this also checks whether another daemon is running,
        return false;
    }
    d->write_daemon_state();

    d->m_is_running = true;

    d->writeLogRecord("start-daemon");

    if (!d->m_watcher.directories().isEmpty()) {
        d->m_watcher.removePaths(d->m_watcher.directories());
    }
    d->m_watcher.addPath(MPDaemon::daemonPath() + "/daemon_commands");

    QTime timer0;
    timer0.start();
    QTime timer1;
    timer1.start();
    QTime timer2;
    timer2.start();
    long num_cycles = 0;
    while (d->m_is_running) {
        if (timer0.elapsed() > 2000) {
            if (!d->write_running_file())
                return false;
        }
        if (timer1.elapsed() > 5000) {
            d->writeLogRecord("timer1", "num_cycles", (long long)num_cycles);
            num_cycles = 0;
            timer1.restart();
            d->write_daemon_state();
            printf(".");
        }
        if (timer2.elapsed() > 10 * 60000) {
            d->writeLogRecord("timer2");
            timer2.restart();
            printf("\n");
        }
        d->stop_orphan_processes_and_scripts();
        d->handle_scripts();
        d->handle_processes();
        qApp->processEvents();
        MPDaemon::wait(100);
        num_cycles++;
    }

    d->writeLogRecord("stop-daemon");

    return true;
}

void MPDaemon::clearProcessing()
{
    QStringList keys = d->m_pripts.keys();
    foreach(QString key, keys)
    {
        d->stop_or_remove_pript(key);
    }
}

QString MPDaemon::daemonPath()
{
    QString ret = CacheManager::globalInstance()->localTempPath() + "/mpdaemon";
    mkdir_if_doesnt_exist(ret);
    mkdir_if_doesnt_exist(ret + "/daemon_state");
    mkdir_if_doesnt_exist(ret + "/daemon_commands");
    mkdir_if_doesnt_exist(ret + "/completed_processes");
    return ret;
}

QString MPDaemon::makeTimestamp(const QDateTime& dt)
{
    return dt.toString("yyyy-MM-dd-hh-mm-ss-zzz");
}

QDateTime MPDaemon::parseTimestamp(const QString& timestamp)
{
    return QDateTime::fromString(timestamp, "yyyy-MM-dd-hh-mm-ss-zzz");
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

void MPDaemon::wait(qint64 msec)
{
    usleep(msec * 1000);
}

void MPDaemon::slot_commands_directory_changed()
{
    QString path = MPDaemon::daemonPath() + "/daemon_commands";
    QStringList fnames = QDir(path).entryList(QStringList("*.command"), QDir::Files, QDir::Name);
    foreach(QString fname, fnames)
    {
        QString path0 = path + "/" + fname;
        qint64 elapsed_sec = QFileInfo(path0).lastModified().secsTo(QDateTime::currentDateTime());
        if (elapsed_sec > 20) {
            if (!QFile::remove(path0)) {
                qCritical() << "Unable to remove command file: " + path0;
            }
        } else {
            QString json = read_text_file(path0);
            QJsonParseError error;
            QJsonObject obj = QJsonDocument::fromJson(json.toLatin1(), &error).object();
            if (error.error != QJsonParseError::NoError) {
                qCritical() << "Error in slot_commands_directory_changed parsing json file";
            } else {
                d->process_command(obj);
            }
            if (!QFile::remove(path0)) {
                qCritical() << "Problem removing command file: " + path0;
            }
        }
    }
}

void MPDaemon::slot_pript_qprocess_finished()
{
    debug_log(__FUNCTION__, __FILE__, __LINE__);

    QProcess* P = qobject_cast<QProcess*>(sender());
    if (!P)
        return;
    QString pript_id = P->property("pript_id").toString();
    MPDaemonPript* S;
    if (d->m_pripts.contains(pript_id)) {
        S = &d->m_pripts[pript_id];
    } else {
        d->writeLogRecord("error", "message", "Unexpected problem in slot_pript_qprocess_finished. Unable to find script or process with id: " + pript_id);
        qCritical() << "Unexpected problem in slot_pript_qprocess_finished. Unable to find script or process with id: " + pript_id;
        return;
    }
    if (!S->output_fname.isEmpty()) {
        debug_log(__FUNCTION__, __FILE__, __LINE__);
        QString runtime_results_json = read_text_file(S->output_fname);
        if (runtime_results_json.isEmpty()) {
            S->success = false;
            S->error = "Could not read results file: " + S->output_fname;
        } else {
            QJsonParseError error;
            S->runtime_results = QJsonDocument::fromJson(runtime_results_json.toLatin1(), &error).object();
            if (error.error != QJsonParseError::NoError) {
                S->success = false;
                S->error = "Error parsing json of runtime results.";
            } else {
                S->success = S->runtime_results["success"].toBool();
                S->error = S->runtime_results["error"].toString();
            }
        }
    } else {
        S->success = true;
    }
    d->finish_and_finalize(*S);

    QJsonObject obj0;
    obj0["pript_id"] = pript_id;
    obj0["reason"] = "finished";
    obj0["success"] = S->success;
    obj0["error"] = S->error;
    if (S->prtype == ScriptType) {
        d->writeLogRecord("stop-script", obj0);
        printf("  Script %s finished ", pript_id.toLatin1().data());
    } else {
        d->writeLogRecord("stop-process", obj0);
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
            } else {
                d->writeLogRecord("error", "pript_id", pript_id, "message", "Process did not finish after waiting even though we are in the slot for finished!!");
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

void MPDaemon::slot_qprocess_output()
{

    QProcess* P = qobject_cast<QProcess*>(sender());
    if (!P)
        return;
    QByteArray str = P->readAll();
    QString pript_id = P->property("pript_id").toString();
    if ((d->m_pripts.contains(pript_id)) && (d->m_pripts[pript_id].stdout_file)) {
        if (d->m_pripts[pript_id].stdout_file->isOpen()) {
            d->m_pripts[pript_id].stdout_file->write(str);
            d->m_pripts[pript_id].stdout_file->flush();
        }
    } else {
        printf("%s", str.data());
    }
}

/*
void MPDaemonPrivate::write_daemon_state()
{
    /// Witold rather than starting at 100000, I'd like to format the num in the fname to be like 0000023. Could you please help?
    static long num = 100000;
    QString timestamp = MPDaemon::makeTimestamp();
    QString fname = QString("%1/daemon_state/%2.%3.json").arg(MPDaemon::daemonPath()).arg(timestamp).arg(num);
    num++;

    QJsonObject state;

    state["is_running"] = m_is_running;

    {
        QJsonObject scripts;
        QJsonObject processes;
        QStringList keys = m_pripts.keys();
        foreach (QString key, keys) {
            if (m_pripts[key].prtype == ScriptType)
                scripts[key] = pript_struct_to_obj(m_pripts[key]);
            else
                processes[key] = pript_struct_to_obj(m_pripts[key]);
        }
        state["scripts"] = scripts;
        state["processes"] = processes;
    }

    QString json = QJsonDocument(state).toJson();
    write_text_file(fname + ".tmp", json);
    /// Witold I don't think rename is an atomic operation. Is there a way to guarantee that I don't read the file halfway through the rename?
    QFile::rename(fname + ".tmp", fname);

    //remove the pripts that are finished
    {
        QStringList keys = m_pripts.keys();
        foreach (QString key, keys) {
            if (m_pripts[key].is_finished) {
                m_pripts.remove(key);
            }
        }
    }

    //finally, clean up
    QStringList list = QDir(MPDaemon::daemonPath() + "/daemin_state").entryList(QStringList("*.json"), QDir::Files, QDir::Name);
    foreach (QString fname, list) {
        QString path0 = MPDaemon::daemonPath() + "/daemon_state/" + fname;
        qint64 secs = QFileInfo(path0).lastModified().secsTo(QDateTime::currentDateTime());
        if ((secs <= -60) || (secs >= 60)) { //I feel a bit paranoid. That's why I allow some future stuff.
            QFile::remove(path0);
        }
    }
}
*/

void MPDaemonPrivate::process_command(QJsonObject obj)
{
    writeLogRecord("process-command", obj);
    QString command = obj.value("command").toString();
    if (command == "stop") {
        debug_log(__FUNCTION__, __FILE__, __LINE__);
        m_is_running = false;
    } else if (command == "queue-script") {
        debug_log(__FUNCTION__, __FILE__, __LINE__);
        MPDaemonPript S = pript_obj_to_struct(obj);
        S.prtype = ScriptType;
        S.timestamp_queued = QDateTime::currentDateTime();
        if (m_pripts.contains(S.id)) {
            writeLogRecord("error", "message", "Unable to queue script. Process or script with this id already exists: " + S.id);
            qCritical() << "Unable to queue script. Process or script with this id already exists: " + S.id;
            return;
        }
        writeLogRecord("queue-script", "pript_id", S.id);
        printf("QUEUING SCRIPT %s\n", S.id.toLatin1().data());
        m_pripts[S.id] = S;
        write_pript_file(S);
    } else if (command == "queue-process") {
        debug_log(__FUNCTION__, __FILE__, __LINE__);
        MPDaemonPript P = pript_obj_to_struct(obj);
        P.prtype = ProcessType;
        P.timestamp_queued = QDateTime::currentDateTime();
        if (m_pripts.contains(P.id)) {
            writeLogRecord("error", "message", "Unable to queue process. Process or script with this id already exists: " + P.id);
            qCritical() << "Unable to queue process. Process or script with this id already exists: " + P.id;
            return;
        }
        writeLogRecord("queue-process", P.id);
        printf("QUEUING PROCESS %s %s\n", P.processor_name.toLatin1().data(), P.id.toLatin1().data());
        m_pripts[P.id] = P;
        write_pript_file(P);
    } else if (command == "clear-processing") {
        q->clearProcessing();
    } else {
        qCritical() << "Unrecognized command: " + command;
        writeLogRecord("error", "message", "Unrecognized command: " + command);
    }
}

void MPDaemonPrivate::writeLogRecord(QString record_type, QString key1, QVariant val1, QString key2, QVariant val2, QString key3, QVariant val3)
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

void MPDaemonPrivate::writeLogRecord(QString record_type, const QJsonObject& obj)
{
    QJsonObject X;
    X["record_type"] = record_type;
    X["timestamp"] = QDateTime::currentDateTime().toString("yyyy-MM-dd|hh:mm:ss.zzz");
    X["data"] = obj;
    QString line = QJsonDocument(X).toJson(QJsonDocument::Compact);
    if (!m_log_path.isEmpty()) {
        append_line_to_file(m_log_path + "/mpdaemon.log", line);
    }
}

bool MPDaemonPrivate::handle_scripts()
{
    int max_simultaneous_scripts = 100;
    if (num_running_scripts() < max_simultaneous_scripts) {
        if (num_pending_scripts() > 0) {
            if (launch_next_script()) {
                printf("%d scripts running.\n", num_running_scripts());
            } else {
                qCritical() << "Unexpected problem. Failed to launch_next_script";
                return false;
            }
        }
    }
    return true;
}

bool MPDaemonPrivate::launch_next_script()
{
    debug_log(__FUNCTION__, __FILE__, __LINE__);

    QStringList keys = m_pripts.keys();
    foreach(QString key, keys)
    {
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

bool MPDaemonPrivate::launch_pript(QString pript_id)
{
    debug_log(__FUNCTION__, __FILE__, __LINE__);

    MPDaemonPript* S;
    if (!m_pripts.contains(pript_id))
        return false;
    S = &m_pripts[pript_id];
    if (S->is_running)
        return false;
    if (S->is_finished)
        return false;

    QString exe = qApp->applicationFilePath();
    QStringList args;

    if (S->prtype == ScriptType) {
        debug_log(__FUNCTION__, __FILE__, __LINE__);
        args << "run-script";
        if (!S->output_fname.isEmpty()) {
            args << "--~script_output=" + S->output_fname;
        }
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
            if (compute_checksum_of_file(fname) != S->script_path_checksums.value(ii)) {
                QString message = "Script checksums do not match. Script file has changed since queueing: " + fname + " Not launching process: " + pript_id;
                qWarning() << message;
                qWarning() << compute_checksum_of_file(fname) << "<>" << S->script_path_checksums.value(ii);
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
        write_text_file(par_fname, parameters_json);
        args << par_fname;
    } else if (S->prtype == ProcessType) {
        debug_log(__FUNCTION__, __FILE__, __LINE__);
        args << "run-process";
        args << S->processor_name;
        if (!S->output_fname.isEmpty())
            args << "--~process_output=" + S->output_fname;
        QStringList pkeys = S->parameters.keys();
        foreach(QString pkey, pkeys)
        {
            args << QString("--%1=%2").arg(pkey).arg(S->parameters[pkey].toString());
        }
        S->runtime_opts.num_threads_allotted = S->num_threads_requested;
        S->runtime_opts.memory_gb_allotted = S->memory_gb_requested;
    }
    debug_log(__FUNCTION__, __FILE__, __LINE__);
    QProcess* qprocess = new QProcess;
    qprocess->setProperty("pript_id", pript_id);
    qprocess->setProcessChannelMode(QProcess::MergedChannels);
    QObject::connect(qprocess, SIGNAL(readyRead()), q, SLOT(slot_qprocess_output()));
    if (S->prtype == ScriptType) {
        printf("   Launching script %s: ", pript_id.toLatin1().data());
        writeLogRecord("start-script", "pript_id", pript_id);
        foreach(QString fname, S->script_paths)
        {
            QString str = QFileInfo(fname).fileName();
            printf("%s ", str.toLatin1().data());
        }
        printf("\n");
    } else {
        printf("   Launching process %s %s: ", S->processor_name.toLatin1().data(), pript_id.toLatin1().data());
        writeLogRecord("start-process", "pript_id", pript_id);
        QString cmd = args.join(" ");
        printf("%s\n", cmd.toLatin1().data());
    }

    QObject::connect(qprocess, SIGNAL(finished(int)), q, SLOT(slot_pript_qprocess_finished()));

    debug_log(__FUNCTION__, __FILE__, __LINE__);

    qprocess->start(exe, args);
    if (qprocess->waitForStarted()) {
        S->qprocess = qprocess;
        if (!S->stdout_fname.isEmpty()) {
            S->stdout_file = new QFile(S->stdout_fname);
            if (!S->stdout_file->open(QFile::WriteOnly)) {
                qCritical() << "Unable to open stdout file for writing: " + S->stdout_fname;
                delete S->stdout_file;
                S->stdout_file = 0;
            }
        }
        S->is_running = true;
        S->timestamp_started = QDateTime::currentDateTime();
        write_pript_file(*S);
        return true;
    } else {
        debug_log(__FUNCTION__, __FILE__, __LINE__);
        if (S->prtype == ScriptType) {
            writeLogRecord("stop-script", "pript_id", pript_id, "reason", "Unable to start script.");
            qCritical() << "Unable to start script: " + S->id;
        } else {
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

ProcessResources MPDaemonPrivate::compute_process_resources_available()
{
    ProcessResources ret = m_total_resources_available;
    QStringList keys = m_pripts.keys();
    foreach(QString key, keys)
    {
        if (m_pripts[key].prtype == ProcessType) {
            if (m_pripts[key].is_running) {
                ProcessRuntimeOpts rtopts = m_pripts[key].runtime_opts;
                ret.num_threads -= rtopts.num_threads_allotted;
                ret.memory_gb -= rtopts.memory_gb_allotted;
            }
        }
    }
    return ret;
}

ProcessResources MPDaemonPrivate::compute_process_resources_needed(MPDaemonPript P)
{
    ProcessResources ret;
    ret.num_threads = P.num_threads_requested;
    ret.memory_gb = P.memory_gb_requested;
    return ret;
}

bool MPDaemonPrivate::handle_processes()
{
    ProcessResources pr_available = compute_process_resources_available();
    QStringList keys = m_pripts.keys();
    foreach(QString key, keys)
    {
        if (m_pripts[key].prtype == ProcessType) {
            if ((!m_pripts[key].is_running) && (!m_pripts[key].is_finished)) {
                ProcessResources pr_needed = compute_process_resources_needed(m_pripts[key]);
                if (is_at_most(pr_needed, pr_available)) {
                    if (process_parameters_are_okay(key)) {
                        if (okay_to_run_process(key)) { //check whether there are io file conflicts at the moment
                            if (launch_pript(key)) {
                                pr_available.num_threads -= m_pripts[key].runtime_opts.num_threads_allotted;
                                pr_available.memory_gb -= m_pripts[key].runtime_opts.memory_gb_allotted;
                            }
                        }
                    } else {
                        writeLogRecord("unqueue-process", "pript_id", key, "reason", "processor not found or parameters are incorrect.");
                        m_pripts.remove(key);
                    }
                }
            }
        }
    }
    return true;
}

bool MPDaemonPrivate::process_parameters_are_okay(const QString& key)
{
    //check that the processor is registered and that the parameters are okay
    ProcessManager* PM = ProcessManager::globalInstance();
    if (!m_pripts.contains(key))
        return false;
    MPDaemonPript P0 = m_pripts[key];
    MLProcessor MLP = PM->processor(P0.processor_name);
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

bool MPDaemonPrivate::okay_to_run_process(const QString& key)
{
    //next we check that there are no running processes where the input or output files conflict with the proposed input or output files
    //but note that it is okay for the same input file to be used in more than one simultaneous process
    QSet<QString> pending_input_paths;
    QSet<QString> pending_output_paths;
    QStringList pripts_keys = m_pripts.keys();
    foreach(QString key0, pripts_keys)
    {
        if (m_pripts[key0].prtype == ProcessType) {
            if (m_pripts[key0].is_running) {
                QStringList input_paths = get_input_paths(m_pripts[key0]);
                foreach(QString path0, input_paths)
                {
                    pending_input_paths.insert(path0);
                }
                QStringList output_paths = get_output_paths(m_pripts[key0]);
                foreach(QString path0, output_paths)
                {
                    pending_output_paths.insert(path0);
                }
            }
        }
    }
    {
        QStringList proposed_input_paths = get_input_paths(m_pripts[key]);
        QStringList proposed_output_paths = get_output_paths(m_pripts[key]);
        foreach(QString path0, proposed_input_paths)
        {
            if (pending_output_paths.contains(path0))
                return false;
        }
        foreach(QString path0, proposed_output_paths)
        {
            if (pending_input_paths.contains(path0))
                return false;
            if (pending_output_paths.contains(path0))
                return false;
        }
    }
    return true;
}

QStringList MPDaemonPrivate::get_input_paths(MPDaemonPript P)
{
    QStringList ret;
    if (P.prtype != ProcessType)
        return ret;
    ProcessManager* PM = ProcessManager::globalInstance();
    MLProcessor MLP = PM->processor(P.processor_name);
    QStringList pnames = MLP.inputs.keys();
    foreach(QString pname, pnames)
    {
        QString path0 = P.parameters.value(pname).toString();
        if (!path0.isEmpty()) {
            ret << path0;
        }
    }
    return ret;
}

QStringList MPDaemonPrivate::get_output_paths(MPDaemonPript P)
{
    QStringList ret;
    if (P.prtype != ProcessType)
        return ret;
    ProcessManager* PM = ProcessManager::globalInstance();
    MLProcessor MLP = PM->processor(P.processor_name);
    QStringList pnames = MLP.outputs.keys();
    foreach(QString pname, pnames)
    {
        QString path0 = P.parameters.value(pname).toString();
        if (!path0.isEmpty()) {
            ret << path0;
        }
    }
    return ret;
}

#include "signal.h"
bool MPDaemonPrivate::write_running_file()
{
    QString fname = cfp(mlTmpPath() + "/mpdaemon_running.pid");
    QString txt = QString("%1").arg(qApp->applicationPid());

    //check for another daemon running!
    if (QFileInfo(fname).lastModified().secsTo(QDateTime::currentDateTime()) <= 60) {
        QString txt0 = read_text_file(fname);
        if ((!txt0.isEmpty()) && (txt0 != txt)) {
            if (kill(txt0.toLong(), 0) == 0) {
                printf("Another daemon seems to be running. Closing.\n");
                return false;
            }
        }
    }
    write_text_file(fname, txt);
    //we will be forgiving if we cannot write the text file, since it is more important that the daemon stays up
    return true;
}

void MPDaemonPrivate::write_daemon_state()
{
    /// Witold rather than starting at 100000, I'd like to format the num in the fname to be like 0000023. Could you please help?
    static long num = 100000;
    QString timestamp = MPDaemon::makeTimestamp();
    QString fname = QString("%1/daemon_state/%2.%3.json").arg(MPDaemon::daemonPath()).arg(timestamp).arg(num);
    num++;

    QJsonObject state;

    state["is_running"] = m_is_running;

    {
        QJsonObject scripts;
        QJsonObject processes;
        QStringList keys = m_pripts.keys();
        foreach(QString key, keys)
        {
            if (m_pripts[key].prtype == ScriptType)
                scripts[key] = pript_struct_to_obj(m_pripts[key], AbbreviatedRecord);
            else
                processes[key] = pript_struct_to_obj(m_pripts[key], AbbreviatedRecord);
        }
        state["scripts"] = scripts;
        state["processes"] = processes;
    }

    QString json = QJsonDocument(state).toJson();
    write_text_file(fname + ".tmp", json);
    /// Witold I don't think rename is an atomic operation. Is there a way to guarantee that I don't read the file halfway through the rename?
    QFile::rename(fname + ".tmp", fname);

    //remove the pripts that have been finished for a while
    {
        QStringList keys = m_pripts.keys();
        foreach(QString key, keys)
        {
            if (m_pripts[key].is_finished) {
                double elapsed_sec = m_pripts[key].timestamp_finished.secsTo(QDateTime::currentDateTime());
                if (elapsed_sec > 20) {
                    m_pripts.remove(key);
                }
            }
        }
    }

    //finally, clean up
    QStringList list = QDir(MPDaemon::daemonPath() + "/daemon_state").entryList(QStringList("*.json"), QDir::Files, QDir::Name);
    foreach(QString fname, list)
    {
        QString path0 = MPDaemon::daemonPath() + "/daemon_state/" + fname;
        qint64 secs = QFileInfo(path0).lastModified().secsTo(QDateTime::currentDateTime());
        if ((secs <= -60) || (secs >= 60)) { //I feel a bit paranoid. That's why I allow some future stuff.
            QFile::remove(path0);
        }
    }
}

bool MPDaemonPrivate::stop_or_remove_pript(const QString& key)
{
    if (!m_pripts.contains(key))
        return false;
    MPDaemonPript* PP = &m_pripts[key];
    if ((PP->is_running)) {
        if (PP->qprocess) {
            qWarning() << "Terminating qprocess: " + key;
            PP->qprocess->terminate(); // I think it's okay to terminate a process. It won't cause this program to crash.
            delete PP->qprocess;
        }
        finish_and_finalize(*PP);
        m_pripts.remove(key);
        if (PP->prtype == ScriptType)
            writeLogRecord("stop-script", "pript_id", key, "reason", "requested");
        else
            writeLogRecord("stop-process", "pript_id", key, "reason", "requested");
    } else {
        m_pripts.remove(key);
        if (PP->prtype == ScriptType)
            writeLogRecord("unqueue-script", "pript_id", key, "reason", "requested");
        else
            writeLogRecord("unqueue-process", "pript_id", key, "reason", "requested");
    }
    return true;
}

void MPDaemonPrivate::write_pript_file(const MPDaemonPript& P)
{
    if (m_log_path.isEmpty())
        return;
    if (P.id.isEmpty())
        return;
    QString fname;
    if (P.prtype == ScriptType) {
        fname = QString("%1/scripts/%2.json").arg(m_log_path).arg(P.id);
    } else if (P.prtype == ProcessType) {
        fname = QString("%1/processes/%2.json").arg(m_log_path).arg(P.id);
    }
    QJsonObject obj = pript_struct_to_obj(P, RuntimeRecord);
    QString json = QJsonDocument(obj).toJson();
    write_text_file(fname, json);
}

void MPDaemonPrivate::finish_and_finalize(MPDaemonPript& P)
{
    P.is_finished = true;
    P.is_running = false;
    P.timestamp_finished = QDateTime::currentDateTime();
    if (!P.stdout_fname.isEmpty()) {
        P.runtime_results["stdout"] = read_text_file(P.stdout_fname);
    }
    write_pript_file(P);
}

void MPDaemonPrivate::stop_orphan_processes_and_scripts()
{
    QStringList keys = m_pripts.keys();
    foreach(QString key, keys)
    {
        if (!m_pripts[key].is_finished) {
            if ((m_pripts[key].parent_pid) && (!MPDaemon::pidExists(m_pripts[key].parent_pid))) {
                debug_log(__FUNCTION__, __FILE__, __LINE__);
                if (m_pripts[key].qprocess) {
                    if (m_pripts[key].prtype == ScriptType) {
                        writeLogRecord("stop-script", "pript_id", key, "reason", "orphan");
                        qWarning() << "Terminating orphan script qprocess: " + key;
                    } else {
                        writeLogRecord("stop-process", "pript_id", key, "reason", "orphan");
                        qWarning() << "Terminating orphan process qprocess: " + key;
                    }

                    m_pripts[key].qprocess->disconnect(); //so we don't go into the finished slot
                    m_pripts[key].qprocess->terminate();
                    delete m_pripts[key].qprocess;
                    finish_and_finalize(m_pripts[key]);
                    m_pripts.remove(key);
                } else {
                    if (m_pripts[key].prtype == ScriptType) {
                        writeLogRecord("unqueue-script", "pript_id", key, "reason", "orphan");
                        qWarning() << "Removing orphan script: " + key + " " + m_pripts[key].script_paths.value(0);
                    } else {
                        writeLogRecord("unqueue-process", "pript_id", key, "reason", "orphan");
                        qWarning() << "Removing orphan process: " + key + " " + m_pripts[key].processor_name;
                    }
                    finish_and_finalize(m_pripts[key]);
                    m_pripts.remove(key);
                }
            }
        }
    }
}

#include "signal.h"
bool MPDaemon::pidExists(qint64 pid)
{
    // check whether process exists (works on Linux)
    return (kill(pid, 0) == 0);
}

bool MPDaemon::waitForFinishedAndWriteOutput(QProcess* P)
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

int MPDaemonPrivate::num_running_pripts(PriptType prtype)
{
    int ret = 0;
    QStringList keys = m_pripts.keys();
    foreach(QString key, keys)
    {
        if (m_pripts[key].is_running) {
            if (m_pripts[key].prtype == prtype)
                ret++;
        }
    }
    return ret;
}

int MPDaemonPrivate::num_pending_pripts(PriptType prtype)
{
    int ret = 0;
    QStringList keys = m_pripts.keys();
    foreach(QString key, keys)
    {
        if ((!m_pripts[key].is_running) && (!m_pripts[key].is_finished)) {
            if (m_pripts[key].prtype == prtype)
                ret++;
        }
    }
    return ret;
}

QJsonArray stringlist_to_json_array(QStringList list)
{
    QJsonArray ret;
    foreach(QString str, list)
    {
        ret << QJsonValue(str);
    }
    return ret;
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
    QJsonObject ret;
    QStringList keys = map.keys();
    foreach(QString key, keys)
    {
        /// Witold I would like to map numbers to numbers here. Can you help?
        ret[key] = QJsonValue::fromVariant(map[key]);
    }
    return ret;
}

QVariantMap json_obj_to_variantmap(QJsonObject obj)
{
    QVariantMap ret;
    QStringList keys = obj.keys();
    foreach(QString key, keys)
    {
        ret[key] = obj[key].toVariant();
    }
    return ret;
}

QStringList paths_to_file_names(const QStringList& paths)
{
    QStringList ret;
    foreach(QString path, paths)
    {
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
        ret["stdout_fname"] = S.stdout_fname;
        ret["parent_pid"] = QString("%1").arg(S.parent_pid);
    }
    ret["id"] = S.id;
    ret["success"] = S.success;
    ret["error"] = S.error;
    ret["timestamp_queued"] = S.timestamp_queued.toString("yyyy-MM-dd|hh:mm:ss.zzz");
    ret["timestamp_started"] = S.timestamp_started.toString("yyyy-MM-dd|hh:mm:ss.zzz");
    ret["timestamp_finished"] = S.timestamp_finished.toString("yyyy-MM-dd|hh:mm:ss.zzz");
    if (S.prtype == ScriptType) {
        ret["prtype"] = "script";
        if (rt != AbbreviatedRecord) {
            ret["script_paths"] = stringlist_to_json_array(S.script_paths);
            ret["script_path_checksums"] = stringlist_to_json_array(S.script_path_checksums);
        } else {
            ret["script_names"] = stringlist_to_json_array(paths_to_file_names(S.script_paths));
        }
    } else {
        ret["prtype"] = "process";
        ret["processor_name"] = S.processor_name;
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
    ret.parameters = json_obj_to_variantmap(obj.value("parameters").toObject());
    ret.id = obj.value("id").toString();
    ret.output_fname = obj.value("output_fname").toString();
    ret.stdout_fname = obj.value("stdout_fname").toString();
    ret.success = obj.value("success").toBool();
    ret.error = obj.value("error").toString();
    ret.parent_pid = obj.value("parent_pid").toString().toLongLong();
    ret.timestamp_queued = QDateTime::fromString(obj.value("timestamp_queued").toString(), "yyyy-MM-dd|hh:mm:ss.zzz");
    ret.timestamp_started = QDateTime::fromString(obj.value("timestamp_started").toString(), "yyyy-MM-dd|hh:mm:ss.zzz");
    ret.timestamp_finished = QDateTime::fromString(obj.value("timestamp_finished").toString(), "yyyy-MM-dd|hh:mm:ss.zzz");
    if (obj.value("prtype").toString() == "script") {
        ret.prtype = ScriptType;
        ret.script_paths = json_array_to_stringlist(obj.value("script_paths").toArray());
        ret.script_path_checksums = json_array_to_stringlist(obj.value("script_path_checksums").toArray());
    } else {
        ret.prtype = ProcessType;
        ret.processor_name = obj.value("processor_name").toString();
    }
    return ret;
}

bool is_at_most(ProcessResources PR1, ProcessResources PR2)
{
    return ((PR1.num_threads <= PR2.num_threads) && (PR1.memory_gb <= PR2.memory_gb));
}

QJsonObject runtime_opts_struct_to_obj(ProcessRuntimeOpts opts)
{
    QJsonObject ret;
    ret["memory_gb_allotted"] = opts.memory_gb_allotted;
    ret["num_threads_allotted"] = opts.num_threads_allotted;
    return ret;
}
