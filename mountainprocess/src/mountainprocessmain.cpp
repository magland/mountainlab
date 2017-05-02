/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
** Created: 4/27/2016
*******************************************************/

/*
An example of what happens when running a processingn script (aka pipeline)

0. Daemon needs to be running in the background:
    > mountainprocess daemon-start

1. User queues the script
    > mountainprocess queue-script foo_script.pipeline [parameters...]
    A command is send to the daemon via: MPDaemonInterface::queueScript

2. The daemon receives the command and adds this job to its queue of scripts.
    When it becomes time to run the script it, the daemon launches (and tracks) a new QProcess
    > mountainprocess run-script foo_script.pipeline [parameters...]

3. The script is executed by QJSEngine (see run_script() below)
    During the course of execution, various processes will be queued. These involve system calls
    > mountainprocess queue-process [processor_name] [parameters...]

4. When "mountainprocess queue-process" is called, a command is sent to the daemon
    just as above MPDaemonInterface::queueProcess

5. The daemon receives the command (as above) and adds this job to its queue of processes.
    When it becomes time to run the process, the daemon launches (and tracks) a new QProcess
    > mountainprocess run-process [processor_name] [parameters...]

6. The process is actually run. This involves calling the appropriate processor library,
    for example mountainsort.mp

7. When the process QProcess ends, an output file with JSON info is written. That triggers things to stop:
    mountainprocess run-process stops ->
    -> mountainprocess queue-process stops
    once all of the queued processes for the script have stopped ->
    -> mountainprocess run-script stops ->
    -> mountainprocess queue-script stops

If anything crashes along the way, every involved QProcess is killed.

*/

#include "mpdaemon.h"
#include "mpdaemoninterface.h"

#include <QCoreApplication>
#include <QFile>
#include <QDebug>
#include <QSettings>
#include <QJSEngine>
#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonArray>
#include <QFileInfo>
#include <QDir>
#include "processmanager.h"

#include "cachemanager.h"
#include "mlcommon.h"
#include "scriptcontroller2.h"
#include <objectregistry.h>
#include <qprocessmanager.h>
#include <unistd.h>
#include "signal.h"

#ifndef Q_OS_LINUX
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#endif

//#define NO_FORK

/// TODO security in scripts that are able to be submitted
/// TODO title on mountainview from mountainbrowser
/// TODO on startup of mountainprocess daemon show all the loaded processors
/// TODO remove all references to datalaboratory.org and magland.org in the repository (don't just search .h/.cpp files)
/// TODO rigorously check mpserver for potential crashes, unhandled exceptions
/// TODO make MAX_SHORT_TERM_GB and MAX_LONG_TERM_GB configurable (but default should be zero - so we don't interrupt every run -- so like the daemon should be handling it)

#define MAX_SHORT_TERM_GB 100
#define MAX_LONG_TERM_GB 100
#define MAX_MDACHUNK_GB 20

struct run_script_opts;
void print_usage();
void print_daemon_instructions();
bool load_parameter_file(QVariantMap& params, const QString& fname);
bool run_script(const QStringList& script_fnames, const QVariantMap& params, const run_script_opts& opts, QString& error_message, QJsonObject& results);
bool initialize_process_manager();
void remove_system_parameters(QVariantMap& params);
bool queue_pript(PriptType prtype, const CLParams& CLP, QString working_path);
QString get_daemon_state_summary(const QJsonObject& state);

#define EXIT_ON_CRITICAL_ERROR
void mountainprocessMessageOutput(QtMsgType type, const QMessageLogContext& context, const QString& msg);

struct run_script_opts {
    run_script_opts()
    {
        nodaemon = false;
        force_run = false;
    }

    bool nodaemon;
    //QStringList server_urls;
    //QString server_base_path;
    bool force_run;
    QString working_path;
    bool preserve_tempdir = false;
};

QJsonArray monitor_stats_to_json_array(const QList<MonitorStats>& stats);
int compute_peak_mem_bytes(const QList<MonitorStats>& stats);
double compute_peak_cpu_pct(const QList<MonitorStats>& stats);
double compute_avg_cpu_pct(const QList<MonitorStats>& stats);
//void log_begin(int argc,char* argv[]);
//void //log_end();

bool get_missing_prvs(const QVariantMap& clparams, QJsonObject& missing_prvs, QString working_path);
QString get_default_daemon_id();
void set_default_daemon_id(QString daemon_id);
void print_default_daemon_id();
QJsonObject get_daemon_state(QString daemon_id);

void silent_message_output(QtMsgType type, const QMessageLogContext& context, const QString& msg)
{
    Q_UNUSED(type)
    Q_UNUSED(context)
    Q_UNUSED(msg)
    return;
}

void sig_handler(int signum)
{
    (void)signum;
    QProcessManager* manager = ObjectRegistry::getObject<QProcessManager>();
    if (manager) {
        manager->closeAll();
    }
    abort();
}

int main(int argc, char* argv[])
{
    QCoreApplication app(argc, argv);
    CLParams CLP(argc, argv);

    ObjectRegistry registry;

    //The process manager
    QProcessManager* processManager = new QProcessManager;
    registry.addAutoReleasedObject(processManager);
    signal(SIGINT, sig_handler);
    signal(SIGKILL, sig_handler);
    signal(SIGTERM, sig_handler);

    bool detach_mode = CLP.named_parameters.contains("_detach");

    // If _working_path is specified then we change the current directory
    QString working_path = CLP.named_parameters.value("_working_path").toString();
    if (working_path.isEmpty())
        working_path = QDir::currentPath();
    if (!working_path.isEmpty()) {
        if (!QDir::setCurrent(working_path)) {
            qWarning() << "Unable to set working path to: " << working_path;
        }
    }

    //log_begin(argc,argv);

    // Do not allow downloads or *processing* to resolve prv files
    // since this is now handled in a separate gui -- as it very much should!!!
    // The --_prvgui option can be used to launch the gui

    // Find .prv files that are missing, and handle that situation accordingly.
    QJsonObject missing_prvs;
    if (!get_missing_prvs(CLP.named_parameters, missing_prvs, working_path)) {
        if (missing_prvs.keys().isEmpty())
            return -1;
        if (CLP.named_parameters.contains("_prvgui")) {
            QString tmp_fname = CacheManager::globalInstance()->makeLocalFile() + ".prvs";
            QString json = QJsonDocument(missing_prvs).toJson();
            TextFile::write(tmp_fname, json);
            CacheManager::globalInstance()->setTemporaryFileDuration(tmp_fname, 600);
            QString cmd = QString("prv-gui %1").arg(tmp_fname);
            QProcess::startDetached(cmd);
            return -1;
        }
        else {
            qWarning() << missing_prvs.keys();
            QTextStream(stdout) << "One or more prvs could not be found. To resolve these using the GUI, run this command again with the --_prvgui option." << endl;
            return -1;
        }
    }

    // Make sure the MP_DAEMON_ID environment variable is set
    // This will be used in all calls to queue-process (and queue-script)
    // But note that the id argument needs to be used for daemon-start, daemon-stop, daemon-state, etc.
    {
        QString daemon_id = qgetenv("MP_DAEMON_ID");
        if (daemon_id.isEmpty()) {
            QSettings settings("Magland", "MountainLab");
            daemon_id = settings.value("default_daemon_id").toString();
        }
        if (CLP.named_parameters.contains("_daemon_id")) {
            daemon_id = CLP.named_parameters["_daemon_id"].toString();
        }
        qputenv("MP_DAEMON_ID", daemon_id.toUtf8().data());
    }

    QString arg1 = CLP.unnamed_parameters.value(0);
    QString arg2 = CLP.unnamed_parameters.value(1);

    QString log_path = MLUtil::mlLogPath() + "/mountainprocess";
    QString tmp_path = MLUtil::tempPath();
    CacheManager::globalInstance()->setLocalBasePath(tmp_path);

    QString iff;
    if (CLP.named_parameters.contains("_intermediate_file_folder"))
        iff = CLP.named_parameters["_intermediate_file_folder"].toString();
    if (CLP.named_parameters.contains("_iff"))
        iff = CLP.named_parameters["_iff"].toString();

    if (!iff.isEmpty()) {
        if (!detach_mode) {
            qDebug().noquote() << "### Setting intermediate file folder: " + iff;
        }
        CacheManager::globalInstance()->setIntermediateFileFolder(iff);
    }

    qInstallMessageHandler(mountainprocessMessageOutput);

    /// TODO don't need to always load the process manager?

    ProcessManager* PM = ProcessManager::globalInstance();

    setbuf(stdout, NULL);

    if (arg1 == "list-processors") { //Provide a human-readable list of the available processors
        if (!initialize_process_manager()) { //load the processor plugins etc
            //log_end();
            return -1;
        }
        QStringList pnames = PM->processorNames();
        qSort(pnames);
        foreach (QString pname, pnames) {
            printf("%s\n", pname.toLatin1().data());
        }
        //log_end();
        return 0;
    }
    else if (arg1 == "spec") { // Show the spec for a single processor
        qInstallMessageHandler(silent_message_output);
        if (!initialize_process_manager()) { // load the processor plugins etc
            //log_end();
            return -1;
        }
        if (!arg2.isEmpty()) {
            MLProcessor MLP = PM->processor(arg2);
            QString json = QJsonDocument(MLP.spec).toJson(QJsonDocument::Indented);
            printf("%s\n", json.toLatin1().data());
        }
        else {
            QJsonArray processors_array;
            QStringList pnames = PM->processorNames();
            foreach (QString pname, pnames) {
                MLProcessor MLP = PM->processor(pname);
                processors_array.push_back(MLP.spec);
            }
            QJsonObject obj;
            obj["processors"] = processors_array;
            QString json = QJsonDocument(obj).toJson(QJsonDocument::Indented);
            printf("%s\n", json.toLatin1().data());
        }
    }
    else if ((arg1 == "run-process") || (arg1 == "exec-process")) { //Run a process synchronously
        //clean up expired files
        CacheManager::globalInstance()->removeExpiredFiles();

        bool exec_mode = (arg1 == "exec-process");
        //printf("Initializing process manager...\n");
        if (!initialize_process_manager()) { // load the processor plugins etc
            //log_end();
            return -1;
        }
        QString output_fname = CLP.named_parameters.value("_process_output").toString(); //maybe the user/framework specified where output is to be saved
        if (!output_fname.isEmpty()) {
            output_fname = QDir::current().absoluteFilePath(output_fname); //make it absolute
        }
        QString processor_name = arg2; //name of the processor is the second user-supplied arg
        QVariantMap process_parameters = CLP.named_parameters;
        remove_system_parameters(process_parameters); //remove parameters starting with "_"
        int ret = 0;
        QString error_message;
        MLProcessInfo info;
        int parent_pid = CLP.named_parameters.value("_parent_pid", 0).toLongLong();

        bool preserve_tempdir = CLP.named_parameters.contains("_preserve_tempdir");

        bool force_run = CLP.named_parameters.contains("_force_run"); // do not check for already computed
        int request_num_threads = CLP.named_parameters.value("_request_num_threads", 0).toInt(); // the processor may or may not respect this request. But mountainsort/omp does.
        PM->setDefaultParameters(processor_name, process_parameters);
        if ((!force_run) && (!exec_mode) && (PM->processAlreadyCompleted(processor_name, process_parameters))) {
            // We have a record of this procesor already completed. If so, we save a lot of time by not re-running
            printf("Process already completed: %s\n", processor_name.toLatin1().data());
        }
        else {
            QString id;

            //check to see that all our parameters are in order for the given processor
            if (!PM->checkParameters(processor_name, process_parameters)) {
                error_message = "Problem checking process: " + processor_name;
                ret = -1;
            }
            else {
                RequestProcessResources RPR;
                RPR.request_num_threads = request_num_threads;
                id = PM->startProcess(processor_name, process_parameters, RPR, exec_mode, preserve_tempdir); //start the process and retrieve a unique id
                if (id.isEmpty()) {
                    error_message = "Problem starting process: " + processor_name;
                    ret = -1;
                }
                if (!PM->waitForFinished(id, parent_pid)) { //wait for the process to finish
                    error_message = "Problem waiting for process to finish: " + processor_name;
                    ret = -1;
                }
                info = PM->processInfo(id); //get the info about the process from the process manager
                if (info.exit_status == QProcess::CrashExit) {
                    error_message = "Process crashed: " + processor_name;
                    ret = -1;
                }
                if (info.exit_code != 0) {
                    error_message = "Exit code is non-zero: " + processor_name;
                    ret = -1;
                }
                PM->clearProcess(id); //clean up
            }

            printf("---------------------------------------------------------------\n");
            printf("PROCESS COMPLETED (exit code = %d): %s\n", info.exit_code, info.processor_name.toLatin1().data());
            if (!error_message.isEmpty())
                printf("ERROR: %s\n", error_message.toLatin1().data());
            int mb = compute_peak_mem_bytes(info.monitor_stats) / 1000000;
            double cpu = compute_peak_cpu_pct(info.monitor_stats);
            double cpu_avg = compute_avg_cpu_pct(info.monitor_stats);
            double sec = info.start_time.msecsTo(info.finish_time) * 1.0 / 1000;
            printf("Peak RAM: %d MB. Peak CPU: %g%%. Avg CPU: %g%%. Elapsed time: %g seconds.\n", mb, cpu, cpu_avg, sec);
            printf("---------------------------------------------------------------\n");
        }
        QJsonObject obj; //the output info to be saved
        obj["exe_command"] = info.exe_command;
        obj["exit_code"] = info.exit_code;
        if (info.exit_status == QProcess::CrashExit)
            obj["exit_status"] = "CrashExit";
        else
            obj["exit_status"] = "NormalExit";
        obj["finished"] = QJsonValue(info.finished);
        obj["parameters"] = variantmap_to_json_obj(info.parameters);
        obj["processor_name"] = info.processor_name;
        obj["standard_output"] = QString(info.standard_output);
        obj["standard_error"] = QString(info.standard_error);
        obj["success"] = error_message.isEmpty();
        obj["error"] = error_message;
        obj["peak_mem_bytes"] = (int)compute_peak_mem_bytes(info.monitor_stats);
        obj["peak_cpu_pct"] = compute_peak_cpu_pct(info.monitor_stats);
        obj["avg_cpu_pct"] = compute_avg_cpu_pct(info.monitor_stats);
        obj["start_time"] = info.start_time.toString("yyyy-MM-dd:hh-mm-ss.zzz");
        obj["finish_time"] = info.finish_time.toString("yyyy-MM-dd:hh-mm-ss.zzz");
        //obj["monitor_stats"]=monitor_stats_to_json_array(info.monitor_stats); -- at some point we can include this in the file. For now we only worry about the computed peak values
        if (!output_fname.isEmpty()) { //The user wants the results to go in this file
            QFile::remove(output_fname); //important -- added 9/9/16
            QString obj_json = QJsonDocument(obj).toJson();
            if (!TextFile::write(output_fname, obj_json)) {
                qCritical() << "Unable to write results to: " + output_fname;
                return -1;
            }
        }
        if (!error_message.isEmpty()) {
            qCritical() << "Error in mountainprocessmain" << error_message;
        }
        printf("\n");

        //log_end();
        return ret; //returns exit code 0 if okay
    }
    else if (arg1 == "run-script") { // run a script synchronously (although note that individual processes will be queued (unless --_nodaemon is specified), but the script will wait for them to complete before exiting)
        //clean up expired files
        CacheManager::globalInstance()->removeExpiredFiles();

        QString daemon_id = qgetenv("MP_DAEMON_ID");
        if (!CLP.named_parameters.contains("_nodaemon")) {
            if (daemon_id.isEmpty()) {
                printf("*****************************************************************\n");
                printf("Cannot run script because...\n");
                printf("\nNo daemon has been specified.\n");
                printf("Use the --_nodaemon flag or do the following:\n");
                print_daemon_instructions();
                return -1;
            }
            QJsonObject state = get_daemon_state(daemon_id);
            if (!state["is_running"].toBool()) {
                printf("*****************************************************************\n");
                printf("Cannot run script because...\n");
                printf("\nDaemon \"%s\" is not running.\n", daemon_id.toUtf8().data());
                printf("Use the --_nodaemon flag or do the following:\n");
                print_daemon_instructions();
                return -1;
            }
        }

        if (!initialize_process_manager()) { // load the processor plugins etc
            //log_end();
            return -1;
        }

        ProcessManager::globalInstance()->cleanUpCompletedProcessRecords(); //those .json files in completed_process tmp directory

        QString output_fname = CLP.named_parameters.value("_script_output").toString(); //maybe the user or framework specified where output is to be saved
        if (!output_fname.isEmpty()) {
            output_fname = QDir::current().absoluteFilePath(output_fname); //make it absolute
        }

        int ret = 0;
        QString error_message;

        QStringList script_fnames; //names of the javascript source files
        QVariantMap params; //parameters to be passed into the main() function of the javascript
        for (int i = 0; i < CLP.unnamed_parameters.count(); i++) {
            QString str = CLP.unnamed_parameters[i];
            if ((str.endsWith(".js")) || (str.endsWith(".pipeline"))) { // it is a javascript source file
                if (!QFile::exists(str)) {
                    QString str2 = MLUtil::mountainlabBasePath() + "/mountainprocess/scripts/" + str;
                    if (QFile::exists(str2)) {
                        str = str2;
                    }
                }
                if (QFile::exists(str)) {
                    script_fnames << str;
                }
                else {
                    qCritical() << "Unable to find script file: " + str;
                    //log_end();
                    return -1;
                }
            }
            if ((str.endsWith(".par")) || (str.endsWith(".json"))) { // note that we can have multiple parameter files! the later ones override the earlier ones.
                if (!load_parameter_file(params, str)) {
                    qWarning() << "Unable to load parameter file: " + str;
                    //log_end();
                    return -1;
                }
            }
        }
        // load the command-line parameters into params. These will override parameters set by .par or .json files
        QStringList keys0 = CLP.named_parameters.keys();
        foreach (QString key0, keys0) {
            params[key0] = CLP.named_parameters[key0];
        }
        remove_system_parameters(params);

        run_script_opts opts;
        opts.nodaemon = CLP.named_parameters.contains("_nodaemon");
        opts.force_run = CLP.named_parameters.contains("_force_run");
        opts.working_path = QDir::currentPath(); // this should get passed through to the processors
        opts.preserve_tempdir = CLP.named_parameters.contains("_preserve_tempdir");
        QJsonObject results;
        // actually run the script
        if (!run_script(script_fnames, params, opts, error_message, results)) {
            ret = -1;
        }

        QJsonObject obj; //the output info
        obj["script_fnames"] = QJsonArray::fromStringList(script_fnames);
        obj["parameters"] = variantmap_to_json_obj(params);
        obj["success"] = (ret == 0);
        obj["results"] = results;
        obj["error"] = error_message;
        if (!output_fname.isEmpty()) { //The user wants the results to go in this file
            QFile::remove(output_fname); //important -- added 9/9/16
            QString obj_json = QJsonDocument(obj).toJson();
            if (!TextFile::write(output_fname, obj_json)) {
                qCritical() << "Unable to write results to: " + output_fname;
                return -1;
            }
        }

        QJsonArray PP = results["processes"].toArray();

        if (ret == 0) { //success

            printf("\nPeak Memory (MB):\n");
            for (int i = 0; i < PP.count(); i++) {
                QJsonObject QQ = PP[i].toObject();
                QJsonObject RR = QQ["results"].toObject();
                QString processor_name = QQ["processor_name"].toString();
                double mem = RR["peak_mem_bytes"].toDouble() / 1e6;
                if (!processor_name.isEmpty()) {
                    QString tmp = QString("  %1 (%2)").arg(mem).arg(processor_name);
                    printf("%s\n", tmp.toUtf8().data());
                }
            }

            printf("\nPeak CPU percent:\n");
            for (int i = 0; i < PP.count(); i++) {
                QJsonObject QQ = PP[i].toObject();
                QJsonObject RR = QQ["results"].toObject();
                QString processor_name = QQ["processor_name"].toString();
                double cpu = RR["peak_cpu_pct"].toDouble();
                if (!processor_name.isEmpty()) {
                    QString tmp = QString("  %1 (%2)").arg(cpu).arg(processor_name);
                    printf("%s\n", tmp.toUtf8().data());
                }
            }

            printf("\nAvg CPU (pct):\n");
            for (int i = 0; i < PP.count(); i++) {
                QJsonObject QQ = PP[i].toObject();
                QJsonObject RR = QQ["results"].toObject();
                QString processor_name = QQ["processor_name"].toString();
                double cpu = RR["avg_cpu_pct"].toDouble();
                if (!processor_name.isEmpty()) {
                    QString tmp = QString("  %1 (%2)").arg(cpu).arg(processor_name);
                    printf("%s\n", tmp.toUtf8().data());
                }
            }

            printf("\nElapsed time (sec):\n");
            for (int i = 0; i < PP.count(); i++) {
                QJsonObject QQ = PP[i].toObject();
                QJsonObject RR = QQ["results"].toObject();
                QString processor_name = QQ["processor_name"].toString();
                QDateTime start_time = QDateTime::fromString(RR["start_time"].toString(), "yyyy-MM-dd:hh-mm-ss.zzz");
                QDateTime finish_time = QDateTime::fromString(RR["finish_time"].toString(), "yyyy-MM-dd:hh-mm-ss.zzz");
                double elapsed = start_time.msecsTo(finish_time) * 1.0 / 1000;
                if (!processor_name.isEmpty()) {
                    QString tmp = QString("  %1 (%2)").arg(elapsed).arg(processor_name);
                    printf("%s\n", tmp.toUtf8().data());
                }
            }
        }

        //log_end();
        return ret;
    }
    else if (arg1 == "set-default-daemon") {
        set_default_daemon_id(arg2);
        printf("default daemon id is set to: ");
        print_default_daemon_id();
        return 0;
    }
    else if (arg1 == "get-default-daemon") {
        print_default_daemon_id();
        return 0;
    }
    else if (arg1 == "queue-script") { //Queue a script -- to be executed when resources are available
        //clean up expired files
        CacheManager::globalInstance()->removeExpiredFiles();

        QString daemon_id = qgetenv("MP_DAEMON_ID");
        if (daemon_id.isEmpty()) {
            printf("\nCANNOT QUEUE SCRIPT: No daemon has been specified.\n");
            printf("Use the --_nodaemon flag or do the following:\n");
            print_daemon_instructions();
            return -1;
        }

        if (queue_pript(ScriptType, CLP, working_path)) {
            //log_end();
            return 0;
        }
        else {
            //log_end();
            return -1;
        }
    }
    else if (arg1 == "queue-process") {
        //clean up expired files
        CacheManager::globalInstance()->removeExpiredFiles();

        QString daemon_id = qgetenv("MP_DAEMON_ID");
        if (daemon_id.isEmpty()) {
            printf("\nCANNOT QUEUE PROCESS: No daemon has been specified.\n");
            printf("Use the --_nodaemon flag or do the following:\n");
            print_daemon_instructions();
            return -1;
        }

        if (!initialize_process_manager()) { // load the processor plugins etc... needed because we check to see if the processor exists, and attach the spec (for later checking), before it is queued
            //log_end();
            return -1;
        }

        if (queue_pript(ProcessType, CLP, working_path)) {
            //log_end();
            return 0;
        }
        else {
            //log_end();
            return -1;
        }
    }
    else if (arg1 == "list-daemons") {
        // hack by jfm to temporarily implement mp-list-daemons
        {
            QSettings settings("Magland", "MountainLab");
            QStringList candidates = settings.value("mp-list-daemons-candidates").toStringList();
            foreach (QString candidate, candidates) {
                qputenv("MP_DAEMON_ID", candidate.toUtf8().data());
                MPDaemonClient client;
                QJsonObject state = client.state();
                if (!state.isEmpty()) {
                    printf("%s\n", candidate.toUtf8().data());
                }
            }
        }
        return 0;
    }
    /*
    else if (arg1 == "daemon-start") {
        MPDaemonInterface MPDI;
        if (MPDI.daemonIsRunning()) {
            printf("Daemon is already running.\n");
            return 0;
        }
        QString cmd = qApp->applicationFilePath() + " daemon-start-internal " + QString("--_daemon_id=%1").arg(qApp->property("daemon_id").toString());
        if (!QProcess::startDetached(cmd)) {
            printf("Unexpected problem: Unable to start program for daemon-start-internal.\n");
            return -1;
        }
        QTime timer;
        timer.start();
        MPDaemon::wait(500);
        while (!MPDI.daemonIsRunning()) {
            printf(".");
            MPDaemon::wait(1000);
            if (timer.elapsed() > 10000) {
                printf("Daemon was not started after waiting.\n");
                return -1;
            }
        }
        printf("Daemon has been started ().\n");
        return 0;
    }
    else if (arg1 == "daemon-start-internal") {
    */
    else if (arg1 == "daemon-start") {
        QString daemon_id = arg2;

        qputenv("MP_DAEMON_ID", daemon_id.toUtf8().data());

        if (daemon_id.isEmpty()) {
            printf("You must specify a daemon id like this:\n");
            printf("daemon-start [some_id]\n");
            printf("For example you can use your user name\n");
            exit(-1);
        }
        /// Witold, I am changing the behavior to the following so that there is no interaction with stdin. I believe this is best.
        if (daemon_id != get_default_daemon_id()) {
            QString msg;
            if (get_default_daemon_id().isEmpty())
                msg = "Setting default daemon to: " + daemon_id;
            else
                msg = QString("Changing default daemon id from %1 to %2. To change it back use mp-set-default-daemon-id.").arg(get_default_daemon_id()).arg(daemon_id);
            printf("%s\n", msg.toUtf8().data());
            set_default_daemon_id(daemon_id);
        }
        //clean up expired files
        CacheManager::globalInstance()->removeExpiredFiles();

        /*
        bool user_wants_to_set_as_default = false;
        while (1) {
            if (daemon_id != get_default_daemon_id()) {
                printf("Would you like to set this as the default daemon (Y/N)? ");
                char response[1000];
                /// Witold, I suppose I should use something better than gets() to avoid the warning
                gets(response);
                QString resp = QString(response);
                if (resp.toLower() == "y") {
                    user_wants_to_set_as_default = true;
                    break;
                }
                else if (resp.toLower() == "n") {
                    user_wants_to_set_as_default = false;
                    break;
                }
            }
            else {
                break;
            }
        }
        qputenv("MP_DAEMON_ID", daemon_id.toUtf8().data());
        if (user_wants_to_set_as_default) {
            /// Witold, ideally this should be done *after* the daemon is started
            set_default_daemon_id(daemon_id);
        }
        */

        bool do_fork = (!CLP.named_parameters.contains("_nofork"));

        // Start the mountainprocess daemon
        if (do_fork) {
#ifndef NO_FORK
/*
         *  The following magic ensures we detach from the parent process
         *  and from the controlling terminal. This is to prevent process
         *  that spawned us to wait for our children to complete.
         */
#ifdef Q_OS_LINUX
            // fork, setsid(?), redirect stdout to /dev/null
            if (daemon(1, 0)) {
                exit(1);
            }
#else
            // fork, setsid, fork, redirect stdout to /dev/null
            if (fork() > 0) {
                exit(0);
            }
            Q_UNUSED(setsid());
            if (fork() > 0) {
                exit(0);
            }
            int devnull = open("/dev/null", O_WRONLY);
            Q_UNUSED(dup2(devnull, STDOUT_FILENO));
            Q_UNUSED(dup2(devnull, STDERR_FILENO));
#endif
#endif
        }
        qDebug().noquote() << "Initializing process manager...";
        if (!initialize_process_manager()) { // load the processor plugins etc
            //log_end();
            return -1;
        }
        //QString mdaserver_base_path = MLUtil::configResolvedPath("mountainprocess", "mdaserver_base_path");
        //QString mdachunk_data_path = MLUtil::configResolvedPath("mountainprocess", "mdachunk_data_path");
        //figure out what to do about this cleaner business
        /*
        TempFileCleaner cleaner;
        cleaner.addPath(MLUtil::tempPath() + "/tmp_short_term", MAX_SHORT_TERM_GB);
        cleaner.addPath(MLUtil::tempPath() + "/tmp_long_term", MAX_LONG_TERM_GB);
        cleaner.addPath(mdaserver_base_path + "/tmp_short_term", MAX_SHORT_TERM_GB);
        cleaner.addPath(mdaserver_base_path + "/tmp_long_term", MAX_LONG_TERM_GB);
        cleaner.addPath(mdachunk_data_path + "/tmp_short_term", MAX_MDACHUNK_GB);
        cleaner.addPath(mdachunk_data_path + "/tmp_long_term", MAX_MDACHUNK_GB);
        */

        MountainProcessServer server;
        server.setLogPath(log_path);

        ProcessResources RR; // these are the rules for determining how many processes to run simultaneously
        //RR.num_threads = qMax(0.0, MLUtil::configValue("mountainprocess", "max_num_simultaneous_threads").toDouble());
        //RR.memory_gb = qMax(0.0, MLUtil::configValue("mountainprocess", "max_total_memory_gb").toDouble());
        RR.num_threads = 0;
        RR.memory_gb = 0;
        RR.num_processes = MLUtil::configValue("mountainprocess", "max_num_simultaneous_processes").toDouble();
        server.setTotalResourcesAvailable(RR);
        qDebug().noquote() << "Starting server...";
        if (!server.start())
            return -1;
        return 0;
    }
    else if (arg1 == "daemon-stop") { //Stop the daemon
        QString daemon_id = arg2;
        if (daemon_id.isEmpty()) {
            daemon_id = get_default_daemon_id();
        }
        if (daemon_id.isEmpty()) {
            printf("You must specify a daemon id like this:\n");
            printf("daemon-stop [some_id]\n");
            exit(-1);
        }
        qputenv("MP_DAEMON_ID", daemon_id.toUtf8().data());

        MPDaemonClient client;
        return client.stop() ? 0 : -1;
    }
    else if (arg1 == "daemon-restart") { //Restart the daemon
        QString daemon_id = arg2;
        if (daemon_id.isEmpty()) {
            daemon_id = get_default_daemon_id();
        }
        if (daemon_id.isEmpty()) {
            printf("You must specify a daemon id like this:\n");
            printf("daemon-restart [some_id]\n");
            exit(-1);
        }
        qputenv("MP_DAEMON_ID", daemon_id.toUtf8().data());

        //clean up expired files
        CacheManager::globalInstance()->removeExpiredFiles();

        MPDaemonClient client;
        if (!client.stop() || !client.start())
            return -1;
        printf("Daemon has been restarted.\n");
        return 0;
    }
    /*
    else if (arg1 == "-internal-daemon-start") { //This is called internaly to start the daemon (which is the central program running in the background)
        MPDaemon X;
        if (!X.run())
            return -1;
        return 0;
    }
    */
    else if (arg1 == "print-log") {
        QString daemon_id = arg2;
        if (daemon_id.isEmpty()) {
            printf("You must specify a daemon id like this:\n");
            printf("print-log [some_id]\n");
            exit(-1);
        }
        qputenv("MP_DAEMON_ID", daemon_id.toUtf8().data());

        MPDaemonClient client;
        QJsonArray arr = client.log();
        QTextStream qout(stdout);
        foreach (QJsonValue v, arr) {
            QJsonObject logRecord = v.toObject();
            qout << logRecord["timestamp"].toString() << '\t'
                 << logRecord["record_type"].toString() << '\t'
                 << QJsonDocument(logRecord["data"].toObject()).toJson() << endl;
        }
        return 0;
    }
    else if (arg1 == "log") {
        QString daemon_id = arg2;
        if (daemon_id.isEmpty()) {
            printf("You must specify a daemon id like this:\n");
            printf("log [some_id]\n");
            exit(-1);
        }
        qputenv("MP_DAEMON_ID", daemon_id.toUtf8().data());

        MPDaemonClient client;
        client.contignousLog();
        return 0;
    }
    else if (arg1 == "daemon-state") { //Print some information on the state of the daemon
        QString daemon_id = arg2;
        if (daemon_id.isEmpty()) {
            daemon_id = get_default_daemon_id();
        }
        if (daemon_id.isEmpty()) {
            printf("You must specify a daemon id like this:\n");
            printf("daemon-state [some_id]\n");
            exit(-1);
        }
        QJsonObject state = get_daemon_state(daemon_id);
        if (CLP.named_parameters.contains("script_id")) {
            QString script_id = CLP.named_parameters["script_id"].toString();
            state = state.value("scripts").toObject().value(script_id).toObject();
        }
        QString json = QJsonDocument(state).toJson();
        printf("%s", json.toLatin1().constData());
        //log_end();
        return 0;
    }
    else if (arg1 == "daemon-state-summary") { //Print some information on the state of the daemon
        QString daemon_id = arg2;
        if (daemon_id.isEmpty()) {
            daemon_id = get_default_daemon_id();
        }
        if (daemon_id.isEmpty()) {
            printf("You must specify a daemon id like this:\n");
            printf("daemon-state-summary [some_id]\n");
            exit(-1);
        }
        qputenv("MP_DAEMON_ID", daemon_id.toUtf8().data());

        MPDaemonClient client;
        QString txt = get_daemon_state_summary(client.state());
        printf("%s", txt.toLatin1().constData());
        //log_end();
        return 0;
    }
    else if (arg1 == "clear-processing") { // stop all active processing including running or queued pripts
        QString daemon_id = arg2;
        if (daemon_id.isEmpty()) {
            daemon_id = get_default_daemon_id();
        }
        if (daemon_id.isEmpty()) {
            printf("You must specify a daemon id like this:\n");
            printf("clear-processing [some_id]\n");
            exit(-1);
        }
        qputenv("MP_DAEMON_ID", daemon_id.toUtf8().data());

        MPDaemonClient client;
        client.clearProcessing();
        //log_end();
        return 0;
    }
    else if (arg1 == "get-script") { // print information about a script with given id
        if (!log_path.isEmpty()) {
            QString str = TextFile::read(log_path + "/scripts/" + CLP.named_parameters["id"].toString() + ".json");
            printf("%s", str.toLatin1().data());
        }
    }
    else if (arg1 == "get-process") { // print information about a process with given id
        if (!log_path.isEmpty()) {
            QString str = TextFile::read(log_path + "/processes/" + CLP.named_parameters["id"].toString() + ".json");
            printf("%s", str.toLatin1().data());
        }
    }
    else if (arg1 == "cleanup-cache") { // remove old files in the temporary directory
        CacheManager::globalInstance()->cleanUp();
    }
    else if (arg1 == "temp") {
        QString tmp_path = CacheManager::globalInstance()->localTempPath();
        printf("%s\n", tmp_path.toUtf8().data());
        return 0;
    }
    else {
        print_usage(); //print usage information
        //log_end();
        return -1;
    }

    //log_end();
    return 0;
}

// load the processor plugins etc
bool initialize_process_manager()
{

    // Load the processor paths
    QStringList processor_paths = MLUtil::configResolvedPathList("mountainprocess", "processor_paths");
    if (processor_paths.isEmpty()) {
        qCritical() << "No processor paths found.";
        return false;
    }
    ProcessManager* PM = ProcessManager::globalInstance();

    PM->setProcessorPaths(processor_paths);
    return true;

    //everything below has been moved to ProcessManager

    /*
    // Load the processors

    foreach (QString processor_path, processor_paths) {
        //printf("Searching for processors in %s\n", p0.toLatin1().data());
        //int previous_num_processors = PM->processorNames().count();
        PM->loadProcessors(processor_path);
        //int num_processors = PM->processorNames().count() - previous_num_processors;
        //qDebug().noquote() << QString("Loaded %1 processors in %2").arg(num_processors).arg(processor_path);
    }

    return true;
    */
}

QString remove_comments_in_line(QString line)
{
    int ind = line.indexOf("//");
    if (ind >= 0) {
        return line.mid(0, ind);
    }
    else {
        return line;
    }
}

// thie allows us to have comments in the .json parameter file. Not necessarily a good idea overall
QString remove_comments(QString json)
{
    QStringList lines = json.split("\n");
    for (int i = 0; i < lines.count(); i++) {
        lines[i] = remove_comments_in_line(lines[i]);
    }
    return lines.join("\n");
}

QJsonValue evaluate_js(QString js)
{
    QJSEngine engine;
    return QJsonValue::fromVariant(engine.evaluate(js.toUtf8()).toVariant());
}

bool load_parameter_file(QVariantMap& params, const QString& fname)
{
    QString json = TextFile::read(fname);
    if (json.isEmpty()) {
        qCritical() << "Non-existent or empty parameter file: " + fname;
        return false;
    }
    json = remove_comments(json);
    QJsonParseError error;
    QJsonObject obj = QJsonDocument::fromJson(json.toUtf8(), &error).object();
    if (error.error != QJsonParseError::NoError) {
        qCritical() << "Error parsing json file: " + fname + " : " + error.errorString();
        return false;
    }
    QStringList keys = obj.keys();
    foreach (QString key, keys) {
        QString str = obj[key].toString();
        if ((str.startsWith("eval")) && (str.endsWith(")"))) { // in case we want to evaluate a numerical expression. eg, eval(6*1000)
            obj[key] = evaluate_js(str);
        }
        params[key] = obj[key].toVariant();
    }
    return true;
}

void display_error(QJSValue result)
{
    qDebug().noquote() << result.property("name").toString();
    qDebug().noquote() << result.property("message").toString();
    qDebug().noquote() << QString("%1 line %2").arg(result.property("fileName").toString()).arg(result.property("lineNumber").toInt());
}

// Run a processing script synchronously
bool run_script(const QStringList& script_fnames, const QVariantMap& params, const run_script_opts& opts, QString& error_message, QJsonObject& results)
{
    QJsonObject parameters = variantmap_to_json_obj(params);

    QJSEngine engine;

    ScriptController2 Controller2;
    Controller2.setNoDaemon(opts.nodaemon);
    //Controller2.setServerUrls(opts.server_urls);
    //Controller2.setServerBasePath(opts.server_base_path);
    Controller2.setForceRun(opts.force_run);
    Controller2.setWorkingPath(opts.working_path);
    Controller2.setPreserveTempdir(opts.preserve_tempdir);
    QJSValue MP2 = engine.newQObject(&Controller2);
    engine.globalObject().setProperty("_MP2", MP2);

    foreach (QString fname, script_fnames) {
        QJSValue result = engine.evaluate(TextFile::read(fname), fname);
        if (result.isError()) {
            display_error(result);
            error_message = "Error running script";
            qCritical() << "Error running script.";
            return false;
        }
    }
    {
        QString parameters_json = QJsonDocument(parameters).toJson(QJsonDocument::Compact);
        //Think about passing the parameters object directly in to the main() rather than converting to json string
        QString code = QString("main(JSON.parse('%1'));").arg(parameters_json);
        QJSValue result = engine.evaluate(code);
        if (result.isError()) {
            display_error(result);
            error_message = "Error running script: " + QString("%1 line %2: %3").arg(result.property("fileName").toString()).arg(result.property("lineNumber").toInt()).arg(result.property("message").toString());
            qCritical().noquote() << error_message;
            return false;
        }
    }

    results = Controller2.getResults();

    return true;
}

void print_usage()
{
    printf("Usage:\n");
    printf("mp-run-process [processor_name] --[param1]=[val1] --[param2]=[val2] ... [--_force_run] [--_request_num_threads=4]\n");
    printf("mp-run-script [script1].js [script2.js] ... [file1].par [file2].par ... [--_force_run] [--_request_num_threads=4]\n");
    printf("mp-list-daemons\n");
    printf("mp-daemon-start [some daemon id]\n");
    printf("mp-daemon-stop [some daemon id]\n");
    printf("mp-daemon-restart [some daemon id]\n");
    printf("mp-daemon-state [some daemon id]\n");
    printf("mp-daemon-state [some daemon id] [--script_id=id]\n");
    printf("mp-daemon-state-summary [some daemon id]\n");
    printf("mp-get-default-daemon\n");
    printf("mp-set-default-daemon [some daemon id]\n");
    printf("mountainprocess clear-processing [some daemon id]\n");
    printf("mp-queue-process [processor_name] --_process_output=[optional_output_fname] --[param1]=[val1] --[param2]=[val2] ... [--_force_run]\n");
    printf("mp-list-processors\n");
    printf("mp-spec [processor_name]\n");
    printf("mp-cleanup-cache\n");
}

void print_daemon_instructions()
{
    printf("\nTo connect to a currently running daemon:\n");
    printf("Set the MP_DAEMON_ID environment variable, or\n");
    printf("> mp-set-default-daemon [some daemon id]\n\n");
    printf("To start a new daemon:\n");
    printf("> mp-daemon-start [some daemon id, for example your user name]\n\n");
    printf("To list all running daemons:\n");
    printf("> mp-list-daemons\n\n");
    return;
}

void remove_system_parameters(QVariantMap& params)
{
    QStringList keys = params.keys();
    foreach (QString key, keys) {
        if (key.startsWith("_")) {
            params.remove(key);
        }
    }
}

// Queue a script or process
bool queue_pript(PriptType prtype, const CLParams& CLP, QString working_path)
{
    MPDaemonPript PP;

    if (CLP.named_parameters.contains("_pript_id")) {
        PP.id = CLP.named_parameters["_pript_id"].toString();
    }

    //bool detach = CLP.named_parameters.value("_detach", false).toBool();
    bool detach = CLP.named_parameters.contains("_detach");

    QVariantMap params;
    if (prtype == ScriptType) {
        // It is a script
        for (int i = 0; i < CLP.unnamed_parameters.count(); i++) {
            QString str = CLP.unnamed_parameters[i];
            if ((str.endsWith(".js")) || (str.endsWith(".pipeline"))) {
                PP.script_paths << str;
                PP.script_path_checksums << MLUtil::computeSha1SumOfFile(str);
            }
            if ((str.endsWith(".par")) || (str.endsWith(".json"))) { // note that we can have multiple parameter files! the later ones override the earlier ones.
                if (!load_parameter_file(params, str)) {
                    qWarning() << "Error loading parameter file" << str;
                    return false;
                }
            }
        }
        PP.prtype = ScriptType;
    }
    else {
        // It is a process
        int request_num_threads = CLP.named_parameters.value("_request_num_threads", 0).toInt();
        PP.prtype = ProcessType;
        PP.RPR.request_num_threads = request_num_threads;
        PP.processor_name = CLP.unnamed_parameters.value(1); //arg2 -- the name of the processor
        ProcessManager* PM = ProcessManager::globalInstance();
        if (!PM->processorNames().contains(PP.processor_name)) {
            qWarning() << QString("Unable to find processor %1. (%2 processors loaded)").arg(PP.processor_name).arg(PM->processorNames().count());
            return false;
        }
        PP.processor_spec = PM->processor(PP.processor_name).spec; //this will be checked later for consistency (make sure processor has not changed)
    }

    // Command-line parameters
    QStringList pkeys = CLP.named_parameters.keys();
    foreach (QString pkey, pkeys) {
        params[pkey] = CLP.named_parameters[pkey];
    }
    remove_system_parameters(params);
    PP.parameters = params;

    if (PP.id.isEmpty())
        PP.id = MLUtil::makeRandomId(20); // Random id for the process or script

    if (prtype == ScriptType) {
        // It is a script
        PP.output_fname = CLP.named_parameters["_script_output"].toString();
        if (!PP.output_fname.isEmpty()) {
            PP.output_fname = QDir::current().absoluteFilePath(PP.output_fname); //make it absolute
            QFile::remove(PP.output_fname); //important, added 9/9/16
        }
    }
    else {
        // It is a process
        PP.output_fname = CLP.named_parameters["_process_output"].toString();
        PP.preserve_tempdir = CLP.named_parameters.contains("_preserve_tempdir");
        if (!PP.output_fname.isEmpty()) {
            PP.output_fname = QDir::current().absoluteFilePath(PP.output_fname); //make it absolute
            QFile::remove(PP.output_fname); //important, added 9/9/16
        }
    }
    if (PP.output_fname.isEmpty()) {
        // If the user has not specified an output file, we create one with a name depending on the random id
        if (prtype == ScriptType)
            PP.output_fname = CacheManager::globalInstance()->makeLocalFile("script_output." + PP.id + ".json", CacheManager::ShortTerm);
        else
            PP.output_fname = CacheManager::globalInstance()->makeLocalFile("process_output." + PP.id + ".json", CacheManager::ShortTerm);
    }
    // We write the standard output to a file -- which will be read at the end, I believe
    if (prtype == ScriptType)
        PP.stdout_fname = CacheManager::globalInstance()->makeLocalFile("script_stdout." + PP.id + ".txt", CacheManager::ShortTerm);
    else
        PP.stdout_fname = CacheManager::globalInstance()->makeLocalFile("process_stdout." + PP.id + ".txt", CacheManager::ShortTerm);

    if (!detach) {
        // If we don't plan to detach, then we keep track of the parent pid, so if the parent is terminated, the pript will terminate itself
        PP.parent_pid = QCoreApplication::applicationPid();
    }

    PP.force_run = CLP.named_parameters.contains("_force_run"); // do not check if processes have already run
    PP.working_path = working_path; // all processes and scripts should be run with this working path

    MPDaemonClient client;
    if (!client.isConnected()) {
        QString daemon_id = qgetenv("MP_DAEMON_ID");
        printf("\nCould not connect to daemon: %s\n", daemon_id.toUtf8().data());
        printf("Use the --_nodaemon option or do the following:\n");
        print_daemon_instructions();
        return false;
    }

    if (prtype == ScriptType) {
        // It is a script
        if (!client.queueScript(PP)) { //queue the script
            qWarning() << "Error queueing script";
            return false;
        }
    }
    else {
        // It is a process
        if (!client.queueProcess(PP)) { //queue the process
            qWarning() << "Error queueing process";
            return false;
        }
    }
    if (!detach) {
        // If we are not detaching, then we are going to wait for the output file to appear
        qint64 parent_pid = CLP.named_parameters.value("_parent_pid", 0).toLongLong();
        MPDaemon::waitForFileToAppear(PP.output_fname, -1, false, parent_pid, PP.stdout_fname);
        QJsonParseError error;
        QJsonObject results_obj = QJsonDocument::fromJson(TextFile::read(PP.output_fname).toLatin1(), &error).object();
        if (error.error != QJsonParseError::NoError) {
            qCritical().noquote() << "Error in queue_pript in parsing output json file.";
            return false;
        }
        bool success = results_obj["success"].toBool();
        if (!success) {
            if (prtype == ScriptType) {
                qCritical().noquote() << "Error in script " + PP.id + ": " + results_obj["error"].toString();
            }
            else {
                qCritical().noquote() << "Error in process " + PP.processor_name + " " + PP.id + ": " + results_obj["error"].toString();
            }
            return false;
        }
    }
    return true;
}

void mountainprocessMessageOutput(QtMsgType type, const QMessageLogContext& context, const QString& msg)
{
    (void)context;
    QByteArray localMsg = msg.toLocal8Bit();
    switch (type) {
    case QtDebugMsg:
        //fprintf(stderr, "Debug: %s (%s:%u, %s)\n", localMsg.constData(), context.file, context.line, context.function);
        fprintf(stderr, "Debug: %s\n", localMsg.constData());
        break;
    case QtInfoMsg:
        fprintf(stderr, "%s\n", localMsg.constData());
        break;
    case QtWarningMsg:
        fprintf(stderr, "Warning: %s\n", localMsg.constData());
        break;
    case QtCriticalMsg:
        //write to error log here
        fprintf(stderr, "Critical: %s\n", localMsg.constData());
#ifdef EXIT_ON_CRITICAL_ERROR
        exit(-1);
#endif
        break;
    case QtFatalMsg:
        fprintf(stderr, "Fatal: %s\n", localMsg.constData());
        exit(-1);
    }
}

QJsonArray monitor_stats_to_json_array(const QList<MonitorStats>& stats)
{
    QJsonArray ret;
    for (int i = 0; i < stats.count(); i++) {
        MonitorStats X = stats[i];
        QJsonObject obj;
        obj["timestamp"] = X.timestamp.toMSecsSinceEpoch();
        obj["mem_bytes"] = (int)X.mem_bytes;
        obj["cpu_pct"] = X.cpu_pct;
        ret << obj;
    }
    return ret;
}
int compute_peak_mem_bytes(const QList<MonitorStats>& stats)
{
    int ret = 0;
    for (int i = 0; i < stats.count(); i++) {
        ret = qMax(ret, stats[i].mem_bytes);
    }
    return ret;
}

double compute_peak_cpu_pct(const QList<MonitorStats>& stats)
{
    double ret = 0;
    for (int i = 0; i < stats.count(); i++) {
        ret = qMax(ret, stats[i].cpu_pct);
    }
    return ret;
}

double compute_avg_cpu_pct(const QList<MonitorStats>& stats)
{
    double ret = 0;
    for (int i = 0; i < stats.count(); i++) {
        ret += stats[i].cpu_pct;
    }
    if (stats.count())
        ret /= stats.count();
    return ret;
}

struct ProcessorCount {
    int queued = 0;
    int running = 0;
    int finished = 0;
    int errors = 0;
};

QString get_daemon_state_summary(const QJsonObject& state)
{
    QString ret;
    if (state["is_running"].toBool()) {
        ret += "Daemon is running.\n";
    }
    else {
        ret += "Daemon is NOT running.\n";
        return ret;
    }

    //scripts
    {
        int num_running = 0;
        int num_finished = 0;
        int num_queued = 0;
        int num_errors = 0;
        QJsonObject scripts = state["scripts"].toObject();
        QStringList ids = scripts.keys();
        foreach (QString id, ids) {
            QJsonObject script = scripts[id].toObject();
            if (script["is_running"].toBool()) {
                num_running++;
            }
            else if (script["is_finished"].toBool()) {
                if (!script["error"].toString().isEmpty()) {
                    num_errors++;
                }
                else {
                    num_finished++;
                }
            }
            else {
                num_queued++;
            }
        }
        ret += QString("%1 scripts: %2 queued, %3 running, %4 finished, %5 errors\n").arg(ids.count()).arg(num_queued).arg(num_running).arg(num_finished).arg(num_errors);
    }

    //processes
    {
        int num_running = 0;
        int num_finished = 0;
        int num_queued = 0;
        int num_errors = 0;
        QMap<QString, ProcessorCount> processor_counts;
        QJsonObject processes = state["processes"].toObject();
        QStringList ids = processes.keys();
        foreach (QString id, ids) {
            QJsonObject process = processes[id].toObject();
            if (process["is_running"].toBool()) {
                num_running++;
                processor_counts[process["processor_name"].toString()].running++;
            }
            else if (process["is_finished"].toBool()) {
                if (!process["error"].toString().isEmpty()) {
                    num_errors++;
                    processor_counts[process["processor_name"].toString()].errors++;
                }
                else {
                    num_finished++;
                    processor_counts[process["processor_name"].toString()].finished++;
                }
            }
            else {
                num_queued++;
                processor_counts[process["processor_name"].toString()].queued++;
            }
        }
        ret += QString("%1 processes: %2 queued, %3 running, %4 finished, %5 errors\n").arg(ids.count()).arg(num_queued).arg(num_running).arg(num_finished).arg(num_errors);
        QStringList keys = processor_counts.keys();
        foreach (QString processor_name, keys) {
            ret += QString("  %1: %2 queued, %3 running, %4 finished, %5 errors\n").arg(processor_name).arg(processor_counts[processor_name].queued).arg(processor_counts[processor_name].running).arg(processor_counts[processor_name].finished).arg(processor_counts[processor_name].errors);
        }
    }

    return ret;
}

QJsonObject read_prv_file(QString path)
{
    QString txt = TextFile::read(path);
    if (txt.isEmpty()) {
        qWarning() << "Working directory: " + QDir::currentPath();
        qWarning() << QString("Unable to read .prv file: ") + "[" + path + "]";
        return QJsonObject();
    }
    QJsonParseError err;
    QJsonObject obj = QJsonDocument::fromJson(txt.toUtf8(), &err).object();
    if (err.error != QJsonParseError::NoError) {
        qWarning() << "Error parsing .prv file: " + path;
        return QJsonObject();
    }
    return obj;
}

void get_missing_prvs(QString key, QVariant clparam, QJsonObject& missing_prvs, QString working_path)
{
    if (clparam.type() == QVariant::List) {
        //handle the case of multiple values
        QVariantList list = clparam.toList();
        for (int i = 0; i < list.count(); i++) {
            QVariant val = list[i];
            get_missing_prvs(QString("%1-%2").arg(key).arg(i), val, missing_prvs, working_path);
        }
        return;
    }
    else {
        QString val = clparam.toString();
        if (val.endsWith(".prv")) {
            val = MLUtil::resolvePath(working_path, val);
            QJsonObject obj = read_prv_file(val);
            QString path0 = locate_prv(obj);
            if (path0.isEmpty()) {
                qWarning() << QString("prv is missing: %1 (parameter: %2, working_path:%3)").arg(MLUtil::resolvePath(working_path, val)).arg(key).arg(working_path);
                missing_prvs[key] = obj;
            }
        }
    }
}

bool get_missing_prvs(const QVariantMap& clparams, QJsonObject& missing_prvs, QString working_path)
{
    QStringList keys = clparams.keys();
    foreach (QString key, keys) {
        get_missing_prvs(key, clparams[key], missing_prvs, working_path);
    }
    if (!missing_prvs.isEmpty()) {
        return false;
    }
    return true; //something is missing
}

QString get_default_daemon_id()
{
    QSettings settings("Magland", "MountainLab");
    QString ret = settings.value("default_daemon_id").toString();
    return ret;
}

void set_default_daemon_id(QString daemon_id)
{
    QSettings settings("Magland", "MountainLab");
    settings.setValue("default_daemon_id", daemon_id);
}

void print_default_daemon_id()
{
    QSettings settings("Magland", "MountainLab");
    QString val = settings.value("default_daemon_id").toString();
    printf("%s\n", val.toUtf8().data());
}

/*
static QString s_log_file="";
void log_begin(int argc,char* argv[]) {
    s_log_file=MLUtil::tempPath()+QString("/mountainprocess_log.%1.%2.txt").arg(QCoreApplication::applicationPid()).arg(MLUtil::makeRandomId(4));
    QString txt;
    for (int i=0; i<argc; i++) txt+=QString(argv[i])+" ";
    TextFile::write(s_log_file,txt);
}

void //log_end() {
    QFile::remove(s_log_file);
}
*/

QJsonObject get_daemon_state(QString daemon_id)
{
    qputenv("MP_DAEMON_ID", daemon_id.toUtf8().data());

    MPDaemonClient client;
    QJsonObject state = client.state();
    if (state.isEmpty()) {
        qWarning() << "Can't connect to daemon";
        return QJsonObject();
    }
    return state;
}
