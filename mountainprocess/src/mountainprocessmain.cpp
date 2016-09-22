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
#include "scriptcontroller.h"

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
#include "tempfilecleaner.h"
#include "mlcommon.h"
#include "scriptcontroller2.h"
#include <unistd.h>

#ifndef Q_OS_LINUX
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#endif

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
bool load_parameter_file(QVariantMap& params, const QString& fname);
bool run_script(const QStringList& script_fnames, const QVariantMap& params, const run_script_opts& opts, QString& error_message, QJsonObject& results);
bool initialize_process_manager();
void remove_system_parameters(QVariantMap& params);
bool queue_pript(PriptType prtype, const CLParams& CLP);
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
    QStringList server_urls;
    QString server_base_path;
    bool force_run;
    QString working_path;
};

QJsonArray monitor_stats_to_json_array(const QList<MonitorStats>& stats);
long compute_peak_mem_bytes(const QList<MonitorStats>& stats);
double compute_peak_cpu_pct(const QList<MonitorStats>& stats);
//void log_begin(int argc,char* argv[]);
//void //log_end();

class MountainProcess {
public:
    MountainProcess() {

    }

    bool list_processors() {
        ProcessManager* PM = ProcessManager::globalInstance();
        if (!initialize_process_manager()) { //load the processor plugins etc
            //log_end();
            return false;
        }
        QStringList pnames = PM->processorNames();
        qSort(pnames);
        foreach (QString pname, pnames) {
            printf("%s\n", pname.toLatin1().data());
        }
        //log_end();
        return true;
    }
    bool spec(const QString &name) {
        ProcessManager* PM = ProcessManager::globalInstance();
        if (!initialize_process_manager()) { //load the processor plugins etc
            //log_end();
            return false;
        }
        MLProcessor MLP = PM->processor(name);
        QString json = QJsonDocument(MLP.spec).toJson(QJsonDocument::Indented);
        printf("%s\n", json.toLatin1().data());
        return true;
    }

    bool daemon_start() {
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
        if (!initialize_process_manager()) {
            //log_end();
            return false;
        }
        QString log_path = MLUtil::mlLogPath() + "/mountainprocess";
        QString mdaserver_base_path = MLUtil::configResolvedPath("mountainprocess", "mdaserver_base_path");
        QString mdachunk_data_path = MLUtil::configResolvedPath("server", "mdachunk_data_path");
        //figure out what to do about this cleaner business
        TempFileCleaner cleaner;
        cleaner.addPath(MLUtil::tempPath() + "/tmp_short_term", MAX_SHORT_TERM_GB);
        cleaner.addPath(MLUtil::tempPath() + "/tmp_long_term", MAX_LONG_TERM_GB);
        cleaner.addPath(mdaserver_base_path + "/tmp_short_term", MAX_SHORT_TERM_GB);
        cleaner.addPath(mdaserver_base_path + "/tmp_long_term", MAX_LONG_TERM_GB);
        cleaner.addPath(mdachunk_data_path + "/tmp_short_term", MAX_MDACHUNK_GB);
        cleaner.addPath(mdachunk_data_path + "/tmp_long_term", MAX_MDACHUNK_GB);
        MPDaemon X;
        X.setLogPath(log_path);
        ProcessResources RR;
        RR.num_threads = qMax(1.0, MLUtil::configValue("mountainprocess", "num_threads").toDouble());
        RR.memory_gb = qMax(1.0, MLUtil::configValue("mountainprocess", "memory_gb").toDouble());
        X.setTotalResourcesAvailable(RR);
        return X.run();
    }

    bool daemon_stop() {
        MPDaemonInterface X;
        return X.stop();
    }

    bool daemon_restart() {
        MPDaemonInterface X;
        if (!X.stop() || !X.start())
            return false;
        printf("Daemon has been restarted.\n");
        return true;
    }
    bool daemon_state() {
        MPDaemonInterface X;
        QJsonObject state = X.getDaemonState();
        QString json = QJsonDocument(state).toJson();
        printf("%s", json.toLatin1().data());
        //log_end();
        return true;
    }

    bool daemon_state_summary() {
        MPDaemonInterface X;
        QString txt = get_daemon_state_summary(X.getDaemonState());
        printf("%s", qPrintable(txt));
        //log_end();
        return true;
    }
    bool clear_processing() {
        MPDaemonInterface X;
        X.clearProcessing();
        //log_end();
        return true;
    }

private:

};

int retCode(bool v) { return v ? 0 : -1; }

int main(int argc, char* argv[])
{
    QCoreApplication app(argc, argv);
    CLParams CLP(argc, argv);

    //log_begin(argc,argv);

    if (!resolve_prv_files(CLP.named_parameters)) {
        qWarning() << "Could not resolve .prv file. Try adjusting the bigfile_paths in mountainlab.ini.";
        //log_end();
        return -1;
    }
    MountainProcess prog;

    QString arg1 = CLP.unnamed_parameters.value(0);
    QString arg2 = CLP.unnamed_parameters.value(1);

    QString log_path = MLUtil::mlLogPath() + "/mountainprocess";
    QString tmp_path = MLUtil::tempPath();
    CacheManager::globalInstance()->setLocalBasePath(tmp_path);

    if (arg1.endsWith(".prv")) {
        QString path0 = resolve_prv_file(arg1);
        printf("FILE: %s\n", path0.toLatin1().data());
        //log_end();
        return 0;
    }

    qInstallMessageHandler(mountainprocessMessageOutput);

    /// TODO don't need to always load the process manager?

    ProcessManager* PM = ProcessManager::globalInstance();
    QStringList server_urls = MLUtil::configResolvedPathList("mountainprocess", "server_urls");
    PM->setServerUrls(server_urls);
    QString server_base_path = MLUtil::configResolvedPath("server", "mdaserver_base_path");
    PM->setServerBasePath(server_base_path);

    setbuf(stdout, NULL);

    QString working_path = CLP.named_parameters.value("_working_path").toString();
    if (!working_path.isEmpty()) {
        if (!QDir::setCurrent(working_path)) {
            qWarning() << "Unable to set working path to: " << working_path;
        }
    }

    if (arg1 == "list-processors") { //Provide a human-readable list of the available processors
        return retCode(prog.list_processors());
    }
    else if (arg1 == "spec") {
        return retCode(prog.spec(arg2));
    }
    else if (arg1 == "run-process") { //Run a process synchronously
        if (!initialize_process_manager()) {
            //log_end();
            return -1; //load the processor plugins etc
        }
        QString output_fname = CLP.named_parameters.value("_process_output").toString(); //maybe the user specified where output is to be reported
        if (!output_fname.isEmpty()) {
            output_fname = QDir::current().absoluteFilePath(output_fname); //make it absolute
        }
        QString processor_name = arg2; //name of the processor is the second user-supplied arg
        QVariantMap process_parameters = CLP.named_parameters;
        remove_system_parameters(process_parameters); //remove parameters starting with "_"
        int ret = 0;
        QString error_message;
        MLProcessInfo info;

        bool force_run = CLP.named_parameters.contains("_force_run");
        if ((!force_run) && (PM->processAlreadyCompleted(processor_name, process_parameters))) { //do we have a record of this procesor already completing? If so, we save a lot of time by not re-running
            printf("Process already completed: %s\n", processor_name.toLatin1().data());
        }
        else {
            QString id;

            if (!PM->checkParameters(processor_name, process_parameters)) { //check to see that all our parameters are in order for the given processor
                error_message = "Problem checking process: " + processor_name;
                ret = -1;
            }
            else {
                id = PM->startProcess(processor_name, process_parameters); //start the process and retrieve a unique id
                if (id.isEmpty()) {
                    error_message = "Problem starting process: " + processor_name;
                    ret = -1;
                }
                if (!PM->waitForFinished(id)) { //wait for the process to finish
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
            printf("PROCESS COMPLETED: %s\n", info.processor_name.toLatin1().data());
            if (!error_message.isEmpty())
                printf("ERROR: %s\n", error_message.toLatin1().data());
            long mb = compute_peak_mem_bytes(info.monitor_stats) / 1000000;
            double cpu = compute_peak_cpu_pct(info.monitor_stats);
            if (cpu) {
                printf("Peak usage: %ld MB RAM / %g%% CPU\n", mb, cpu);
            }
            printf("---------------------------------------------------------------\n");
        }
        QJsonObject obj; //the output info
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
        obj["peak_mem_bytes"] = (long long)compute_peak_mem_bytes(info.monitor_stats);
        obj["peak_cpu_pct"] = compute_peak_cpu_pct(info.monitor_stats);
        //obj["monitor_stats"]=monitor_stats_to_json_array(info.monitor_stats); -- at some point we can include this in the file. For now we only worry about the computed peak values
        if (!output_fname.isEmpty()) { //The user wants the results to go in this file
            QFile::remove(output_fname); //important -- added 9/9/16
            QString obj_json = QJsonDocument(obj).toJson();
            if (!TextFile::write(output_fname, obj_json)) {
                qCritical() << "Unable to write results to: " + output_fname;
            }
        }
        if (!error_message.isEmpty()) {
            qCritical() << "Error in mountainprocessmain" << error_message;
        }
        printf("\n");

        //log_end();
        return ret; //returns exit code 0 if okay
    }
    else if (arg1 == "run-script") { //run a script synchronously (although note that individual processes will be queued (unless --_nodaemon is specified), but the script will wait for them to complete)
        if (!initialize_process_manager()) {
            //log_end();
            return -1;
        }

        QString output_fname = CLP.named_parameters.value("_script_output").toString(); //maybe the user specified where output is to be reported
        if (!output_fname.isEmpty()) {
            output_fname = QDir::current().absoluteFilePath(output_fname); //make it absolute
        }

        int ret = 0;
        QString error_message;

        QStringList script_fnames; //names of the javascript source files
        QVariantMap params; //parameters to be passed into the main() function of the javascript
        for (int i = 0; i < CLP.unnamed_parameters.count(); i++) {
            QString str = CLP.unnamed_parameters[i];
            if ((str.endsWith(".js")) || (str.endsWith(".pipeline"))) { //must be a javascript source file
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
                    qWarning() << "Unable to find script file: " + str;
                    //log_end();
                    return -1;
                }
            }
            if ((str.endsWith(".par")) || (str.endsWith(".json"))) { // note that we can have multiple parameter files! the later ones override the earlier ones.
                qDebug() << "Loading parameter file: " + str;
                if (!load_parameter_file(params, str)) {
                    qWarning() << "Unable to load parameter file: " + str;
                    //log_end();
                    return -1;
                }
            }
        }
        QStringList keys0 = CLP.named_parameters.keys();
        foreach (QString key0, keys0) {
            params[key0] = CLP.named_parameters[key0];
        }
        remove_system_parameters(params);

        run_script_opts opts;
        opts.nodaemon = CLP.named_parameters.contains("_nodaemon");
        opts.server_urls = server_urls;
        opts.server_base_path = server_base_path;
        opts.force_run = CLP.named_parameters.contains("_force_run");
        opts.working_path = QDir::currentPath();
        QJsonObject results;
        if (!run_script(script_fnames, params, opts, error_message, results)) { //actually run the script
            ret = -1;
        }

        QJsonObject obj; //the output info
        obj["script_fnames"] = stringlist_to_json_array(script_fnames);
        obj["parameters"] = variantmap_to_json_obj(params);
        obj["success"] = (ret == 0);
        obj["results"] = results;
        obj["error"] = error_message;
        if (!output_fname.isEmpty()) { //The user wants the results to go in this file
            QFile::remove(output_fname); //important -- added 9/9/16
            QString obj_json = QJsonDocument(obj).toJson();
            if (!TextFile::write(output_fname, obj_json)) {
                qCritical() << "Unable to write results to: " + output_fname;
            }
        }

        //log_end();
        return ret;
    }
    else if (arg1 == "daemon-start") {
        return retCode(prog.daemon_start());
    } else if (arg1 == "daemon-stop") { //Stop the daemon
        return retCode(prog.daemon_stop());
    } else if (arg1 == "daemon-restart") { //Restart the daemon
        return retCode(prog.daemon_restart());
    }
    /*
    else if (arg1 == "-internal-daemon-start") { //This is called internaly to start the daemon (which is the central program running in the background)
        MPDaemon X;
        if (!X.run())
            return -1;
        return 0;
    }
    */
    else if (arg1 == "daemon-state") { //Print some information on the state of the daemon
        prog.daemon_state();
        return 0;
    }
    else if (arg1 == "daemon-state-summary") { //Print some information on the state of the daemon
        prog.daemon_state_summary();
        return 0;
    }
    else if (arg1 == "clear-processing") {
        prog.clear_processing();
        return 0;
    }
    else if (arg1 == "queue-script") { //Queue a script -- to be executed when resources are available
        if (queue_pript(ScriptType, CLP)) {
            //log_end();
            return 0;
        }
        else {
            //log_end();
            return -1;
        }
    }
    else if (arg1 == "queue-process") {
        if (queue_pript(ProcessType, CLP)) {
            //log_end();
            return 0;
        }
        else {
            //log_end();
            return -1;
        }
    }
    else if (arg1 == "get-script") {
        if (!log_path.isEmpty()) {
            QString str = TextFile::read(log_path + "/scripts/" + CLP.named_parameters["id"].toString() + ".json");
            printf("%s", str.toLatin1().data());
        }
    }
    else if (arg1 == "get-process") {
        if (!log_path.isEmpty()) {
            QString str = TextFile::read(log_path + "/processes/" + CLP.named_parameters["id"].toString() + ".json");
            printf("%s", str.toLatin1().data());
        }
    }
    else if (arg1 == "create-prv") {

        if (arg2.isEmpty()) {
            print_usage();
            //log_end();
            return -1;
        }
        QString arg3 = CLP.unnamed_parameters.value(2);
        if (arg3.isEmpty())
            arg3 = arg2 + ".prv";

        if (!QFile::exists(arg2)) {
            qWarning() << "File does not exist: " + arg2;
            //log_end();
            return -1;
        }
        QJsonObject obj = make_prv_object(arg2);
        QString json = QJsonDocument(obj).toJson(QJsonDocument::Indented);
        if (!TextFile::write(arg3, json)) {
            qWarning() << "Error writing file: " + arg3;
            //log_end();
            return -1;
        }
    }
    else {
        print_usage(); //print usage information
        //log_end();
        return -1;
    }

    //log_end();
    return 0;
}

bool initialize_process_manager()
{
    /*
     * Load the processor paths
     */
    QStringList processor_paths = MLUtil::configResolvedPathList("mountainprocess", "processor_paths");
    if (processor_paths.isEmpty()) {
        qCritical() << "No processor paths found.";
        return false;
    }

    /*
     * Load the processors
     */
    ProcessManager* PM = ProcessManager::globalInstance();
    foreach (QString processor_path, processor_paths) {
        //printf("Searching for processors in %s\n", p0.toLatin1().data());
        PM->loadProcessors(processor_path);
        int num_processors = PM->processorNames().count();
        printf("Loaded %d processors in %s\n", num_processors, processor_path.toLatin1().data());
    }

    return true;
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

QString remove_comments(QString json)
{
    QStringList lines = json.split("\n");
    for (int i = 0; i < lines.count(); i++) {
        lines[i] = remove_comments_in_line(lines[i]);
    }
    return lines.join("\n");
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
    /// Witold I use toLatin1() everywhere. Is this the appropriate way to convert to byte array?
    /// Jeremy: toUtf8() or toLocal8Bit() might be better
    QJsonObject obj = QJsonDocument::fromJson(json.toLatin1(), &error).object();
    if (error.error != QJsonParseError::NoError) {
        qCritical() << "Error parsing json file: " + fname + " : " + error.errorString();
        return false;
    }
    QStringList keys = obj.keys();
    foreach (QString key, keys) {
        params[key] = obj[key].toVariant();
    }
    return true;
}

void display_error(QJSValue result)
{
    /// Witold there must be a better way to print the QJSValue error message out to the console.
    /// Jeremy: No, this is the proper way. You can use "stack" property to get a full trace, if you want
    qDebug().noquote() << result.property("name").toString(); //okay
    qDebug().noquote() << result.property("message").toString(); //okay
    qDebug().noquote() << QString("%1 line %2").arg(result.property("fileName").toString()).arg(result.property("lineNumber").toInt()); //okay
}

bool run_script(const QStringList& script_fnames, const QVariantMap& params, const run_script_opts& opts, QString& error_message, QJsonObject& results)
{
    QJsonObject parameters = variantmap_to_json_obj(params);

    QJSEngine engine;

    ScriptController Controller;
    Controller.setNoDaemon(opts.nodaemon);
    Controller.setServerUrls(opts.server_urls);
    Controller.setServerBasePath(opts.server_base_path);
    Controller.setForceRun(opts.force_run);
    QJSValue MP = engine.newQObject(&Controller);
    engine.globalObject().setProperty("MP", MP);

    ScriptController2 Controller2;
    Controller2.setNoDaemon(opts.nodaemon);
    Controller2.setServerUrls(opts.server_urls);
    Controller2.setServerBasePath(opts.server_base_path);
    Controller2.setForceRun(opts.force_run);
    Controller2.setWorkingPath(opts.working_path);
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
            qCritical() << error_message;
            return false;
        }
    }

    results = Controller2.getResults();

    return true;
}

void print_usage()
{
    printf("Usage:\n");
    printf("mountainprocess run-process [processor_name] --[param1]=[val1] --[param2]=[val2] ... [--_force_run]\n");
    printf("mountainprocess run-script [script1].js [script2.js] ... [file1].par [file2].par ... [--_force_run] \n");
    printf("mountainprocess daemon-start\n");
    //printf("mountainprocess daemon-stop\n");
    printf("mountainprocess daemon-restart\n");
    printf("mountainprocess daemon-state\n");
    printf("mountainprocess daemon-state-summary\n");
    printf("mountainprocess queue-script --_script_output=[optional_output_fname] [script1].js [script2.js] ... [file1].par [file2].par ...  [--_force_run]\n");
    printf("mountainprocess queue-process [processor_name] --_process_output=[optional_output_fname] --[param1]=[val1] --[param2]=[val2] ... [--_force_run]\n");
    printf("mountainprocess create-prv [filename]\n");
    printf("mountainprocess create-prv [filename] [output_filename].prv\n");
    printf("mountainprocess set-big-file-search-path\n");
    printf("mountainprocess set-big-file-search-path [path]\n");
    printf("mountainprocess list-processors\n");
    printf("mountainprocess spec [processor_name]\n");
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

bool queue_pript(PriptType prtype, const CLParams& CLP)
{
    MPDaemonPript PP;

    bool detach = CLP.named_parameters.value("_detach", false).toBool();

    if (prtype == ScriptType) {
        QVariantMap params;
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
        QStringList pkeys = CLP.named_parameters.keys();
        foreach (QString pkey, pkeys) {
            params[pkey] = CLP.named_parameters[pkey];
        }
        remove_system_parameters(params);
        PP.parameters = params;
        PP.prtype = ScriptType;
    }
    else {
        PP.parameters = CLP.named_parameters;
        PP.prtype = ProcessType;
        remove_system_parameters(PP.parameters);
        PP.processor_name = CLP.unnamed_parameters.value(1); //arg2
    }

    PP.id = MLUtil::makeRandomId(20);

    if (prtype == ScriptType) {
        PP.output_fname = CLP.named_parameters["_script_output"].toString();
        if (!PP.output_fname.isEmpty()) {
            PP.output_fname = QDir::current().absoluteFilePath(PP.output_fname); //make it absolute
            QFile::remove(PP.output_fname); //important, added 9/9/16
        }
    }
    else {
        PP.output_fname = CLP.named_parameters["_process_output"].toString();
        if (!PP.output_fname.isEmpty()) {
            PP.output_fname = QDir::current().absoluteFilePath(PP.output_fname); //make it absolute
            QFile::remove(PP.output_fname); //important, added 9/9/16
        }
    }
    if (PP.output_fname.isEmpty()) {
        if (prtype == ScriptType)
            PP.output_fname = CacheManager::globalInstance()->makeLocalFile("script_output." + PP.id + ".json", CacheManager::ShortTerm);
        else
            PP.output_fname = CacheManager::globalInstance()->makeLocalFile("process_output." + PP.id + ".json", CacheManager::ShortTerm);
    }
    if (prtype == ScriptType)
        PP.stdout_fname = CacheManager::globalInstance()->makeLocalFile("script_stdout." + PP.id + ".txt", CacheManager::ShortTerm);
    else
        PP.stdout_fname = CacheManager::globalInstance()->makeLocalFile("process_stdout." + PP.id + ".txt", CacheManager::ShortTerm);

    if (!detach) {
        PP.parent_pid = QCoreApplication::applicationPid();
    }

    PP.force_run = CLP.named_parameters.contains("_force_run");
    PP.working_path = QDir::currentPath();

    qDebug() << ":::::::::::::::::::::::::::::::::::::::::"
             << "queue_pript" << PP.force_run;

    MPDaemonInterface X;
    // ensure daemon is running

    if (prtype == ScriptType) {
        if (!X.queueScript(PP)) { //queue the script
            qWarning() << "Error queueing script";
            return false;
        }
    }
    else {
        if (!X.queueProcess(PP)) { //queue the process
            qWarning() << "Error queueing process";
            return false;
        }
    }
    if (!detach) {
        qint64 parent_pid = CLP.named_parameters.value("_parent_pid", 0).toLongLong();
        MPDaemon::waitForFileToAppear(PP.output_fname, -1, false, parent_pid, PP.stdout_fname);
        QJsonParseError error;
        QJsonObject results_obj = QJsonDocument::fromJson(TextFile::read(PP.output_fname).toLatin1(), &error).object();
        if (error.error != QJsonParseError::NoError) {
            qWarning() << "Error in queue_pript in parsing output json file.";
            return false;
        }
        bool success = results_obj["success"].toBool();
        if (!success) {
            if (prtype == ScriptType) {
                qWarning() << "Error in script " + PP.id + ": " + results_obj["error"].toString();
            }
            else {
                qWarning() << "Error in process " + PP.processor_name + " " + PP.id + ": " + results_obj["error"].toString();
            }
            return false;
        }
    }
    return true;
}

void mountainprocessMessageOutput(QtMsgType type, const QMessageLogContext& context, const QString& msg)
{
    QByteArray localMsg = msg.toLocal8Bit();
    switch (type) {
    case QtDebugMsg:
        fprintf(stderr, "Debug: %s (%s:%u, %s)\n", localMsg.constData(), context.file, context.line, context.function);
        break;
    case QtInfoMsg:
        fprintf(stderr, "Info: %s (%s:%u, %s)\n", localMsg.constData(), context.file, context.line, context.function);
        break;
    case QtWarningMsg:
        fprintf(stderr, "Warning: %s (%s:%u, %s)\n", localMsg.constData(), context.file, context.line, context.function);
        break;
    case QtCriticalMsg:
        //write to error log here
        fprintf(stderr, "Critical: %s (%s:%u, %s)\n", localMsg.constData(), context.file, context.line, context.function);
#ifdef EXIT_ON_CRITICAL_ERROR
        exit(-1);
#endif
        break;
    case QtFatalMsg:
        fprintf(stderr, "Fatal: %s (%s:%u, %s)\n", localMsg.constData(), context.file, context.line, context.function);
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
        obj["mem_bytes"] = (long long)X.mem_bytes;
        obj["cpu_pct"] = X.cpu_pct;
        ret << obj;
    }
    return ret;
}
long compute_peak_mem_bytes(const QList<MonitorStats>& stats)
{
    long ret = 0;
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
        long num_running = 0;
        long num_finished = 0;
        long num_queued = 0;
        long num_errors = 0;
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
        long num_running = 0;
        long num_finished = 0;
        long num_queued = 0;
        long num_errors = 0;
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
