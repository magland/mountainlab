/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
** Created: 4/27/2016
*******************************************************/

#include "commandlineparams.h"
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
#include "processmanager.h"
#include "textfile.h"
#include "cachemanager.h"
#include "mlutils.h"

/// TODO change all .ini config files to json
/// TODO convert python web servers to nodejs
/// TODO make central list of all servers/services to be running
/// TODO web page that checks status of all services
/// TODO add field linking process to parent script
/// TODO improve stdout view
/// TODO put mountainbrowser in html -- no qtwebkit
/// TODO consolidate all temporary and data directories in mountainlab/data mountainlab/tmp
/// TODO consolidate all configuration files in mountainlab/config
/// TODO security in scripts that are able to be submitted
/// TODO remove mscmdserver -- replace by mpserver

void print_usage();
bool load_parameter_file(QVariantMap& params, const QString& fname);
bool run_script(const QStringList& script_fnames, const QVariantMap& params, QString& error_message);
bool initialize_process_manager();
void remove_system_parameters(QVariantMap& params);
int queue_pript(PriptType prtype, const CLParams& CLP);

int main(int argc, char* argv[])
{
    QCoreApplication app(argc, argv);
    CLParams CLP = commandlineparams(argc, argv);

    CacheManager::globalInstance()->setLocalBasePath(cfp(qApp->applicationDirPath() + "/../tmp"));

    QString arg1 = CLP.unnamed_parameters.value(0);
    QString arg2 = CLP.unnamed_parameters.value(1);

    setbuf(stdout, NULL);

    QString config_fname = cfp(qApp->applicationDirPath() + "/mountainprocess.ini");
    QSettings config(config_fname, QSettings::IniFormat);
    QString log_path = config.value("log_path").toString();
    if (!log_path.isEmpty()) {
        if (QFileInfo(log_path).isRelative()) {
            log_path = cfp(qApp->applicationDirPath() + "/" + log_path);
        }
    }

    if (arg1 == "list-processors") { //Provide a human-readable list of the available processors
        if (!initialize_process_manager())
            return -1; //load the processor plugins etc
        ProcessManager* PM = ProcessManager::globalInstance();
        QStringList pnames = PM->processorNames();
        qSort(pnames);
        foreach (QString pname, pnames) {
            printf("%s\n", pname.toLatin1().data());
        }
        return 0;
    }
    else if (arg1 == "run-process") { //Run a process synchronously
        if (!initialize_process_manager())
            return -1; //load the processor plugins etc
        ProcessManager* PM = ProcessManager::globalInstance();

        QString output_fname = CLP.named_parameters.value("~process_output").toString(); //maybe the user specified where output is to be reported
        QString processor_name = arg2; //name of the processor is the second user-supplied arg
        QVariantMap process_parameters = CLP.named_parameters;
        remove_system_parameters(process_parameters); //remove parameters starting with "~"
        int ret = 0;
        QString error_message;
        MLProcessInfo info;

        if (PM->processAlreadyCompleted(processor_name, process_parameters)) { //do we have a record of this procesor already completing? If so, we save a lot of time by not re-running
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
                if (!PM->waitForFinished(id, -1)) { //wait for the process to finish
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
        if (!output_fname.isEmpty()) { //The user wants the results to go in this file
            QString obj_json = QJsonDocument(obj).toJson();
            if (!write_text_file(output_fname, obj_json)) {
                qWarning() << "Unable to write results to: " + output_fname;
            }
        }
        if (!error_message.isEmpty()) {
            qWarning() << error_message;
        }

        return ret; //returns exit code 0 if okay
    }
    else if (arg1 == "run-script") { //run a script synchronously (although note that individual processes will be queued, but the script will wait for them to complete)
        if (!initialize_process_manager())
            return -1;

        QString output_fname = CLP.named_parameters.value("~script_output").toString(); //maybe the user specified where output is to be reported

        int ret = 0;
        QString error_message;

        QStringList script_fnames; //names of the javascript source files
        QVariantMap params; //parameters to be passed into the main() function of the javascript
        for (int i = 0; i < CLP.unnamed_parameters.count(); i++) {
            QString str = CLP.unnamed_parameters[i];
            if (str.endsWith(".js")) { //must be a javascript source file
                script_fnames << str;
            }
            if (str.endsWith(".par")) { // note that we can have multiple parameter files! the later ones override the earlier ones.
                if (!load_parameter_file(params, str)) {
                    return -1;
                }
            }
        }
        QStringList keys0 = CLP.named_parameters.keys();
        foreach (QString key0, keys0) {
            params[key0] = CLP.named_parameters[key0];
        }
        remove_system_parameters(params);

        if (!run_script(script_fnames, params, error_message)) { //actually run the script
            ret = -1;
        }

        QJsonObject obj; //the output info
        obj["script_fnames"] = stringlist_to_json_array(script_fnames);
        obj["parameters"] = variantmap_to_json_obj(params);
        obj["success"] = (ret == 0);
        obj["error"] = error_message;
        if (!output_fname.isEmpty()) { //The user wants the results to go in this file
            QString obj_json = QJsonDocument(obj).toJson();
            if (!write_text_file(output_fname, obj_json)) {
                qWarning() << "Unable to write results to: " + output_fname;
            }
        }

        return ret;
    }
    else if (arg1 == "daemon-start") {
        if (!initialize_process_manager())
            return -1;
        MPDaemon X;
        X.setLogPath(log_path);
        ProcessResources RR;
        RR.num_threads = config.value("num_threads", 8).toDouble();
        RR.memory_gb = config.value("memory_gb", 8).toDouble();
        X.setTotalResourcesAvailable(RR);
        if (!X.run())
            return -1;
        return 0;
    }
    /*
    else if (arg1 == "-internal-daemon-start") { //This is called internaly to start the daemon (which is the central program running in the background)
        MPDaemon X;
        if (!X.run())
            return -1;
        return 0;
    } else if (arg1 == "daemon-start") { //Start the daemon
        MPDaemonInterface X;
        if (X.start()) {
            printf("Started daemon.\n");
            return 0;
        } else {
            printf("Failed to start daemon.\n");
            return -1;
        }
    } else if (arg1 == "daemon-stop") { //Stop the daemon
        MPDaemonInterface X;
        if (X.stop())
            return 0;
        else
            return -1;
    } else if (arg1 == "daemon-restart") { //Restart the daemon
        MPDaemonInterface X;
        if (!X.stop())
            return -1;
        if (!X.start())
            return -1;
        printf("Daemon has been restarted.\n");
        return 0;
    }
    */
    else if (arg1 == "daemon-state") { //Print some information on the state of the daemon
        MPDaemonInterface X;
        QJsonObject state = X.getDaemonState();
        QString json = QJsonDocument(state).toJson();
        printf("%s", json.toLatin1().data());
        return 0;
    }
    else if (arg1 == "clear-processing") {
        MPDaemonInterface X;
        X.clearProcessing();
        return 0;
    }
    else if (arg1 == "queue-script") { //Queue a script -- to be executed when resources are available
        return queue_pript(ScriptType, CLP);
    }
    else if (arg1 == "queue-process") {
        return queue_pript(ProcessType, CLP);
    }
    else if (arg1 == "get-script") {
        if (!log_path.isEmpty()) {
            QString str = read_text_file(log_path + "/scripts/" + CLP.named_parameters["id"].toString() + ".json");
            printf("%s", str.toLatin1().data());
        }
    }
    else if (arg1 == "get-process") {
        if (!log_path.isEmpty()) {
            QString str = read_text_file(log_path + "/processes/" + CLP.named_parameters["id"].toString() + ".json");
            printf("%s", str.toLatin1().data());
        }
    }
    else {
        print_usage(); //print usage information
        return -1;
    }

    return 0;
}

bool initialize_process_manager()
{
    /*
     * Load configuration file. If it doesn't exist, copy example configuration file.
     */
    QString config_fname = cfp(qApp->applicationDirPath() + "/mountainprocess.ini");
    if (!QFile::exists(config_fname)) {
        if (!QFile::copy(config_fname + ".example", config_fname)) {
            qWarning() << "Unable to copy example configuration file to " + config_fname;
            return false;
        }
    }
    QSettings config(config_fname, QSettings::IniFormat);
    /*
     * Load the processor paths
     */
    QStringList processor_paths = config.value("processor_paths").toStringList();
    if (processor_paths.isEmpty()) {
        qWarning() << "No processor paths found in " + config_fname;
    }

    /*
     * Load the processors
     */
    ProcessManager* PM = ProcessManager::globalInstance();
    foreach (QString processor_path, processor_paths) {
        QString p0 = processor_path;
        if (QFileInfo(p0).isRelative()) {
            p0 = cfp(qApp->applicationDirPath() + "/" + p0);
        }
        printf("Searching for processors in %s\n", p0.toLatin1().data());
        PM->loadProcessors(p0);
    }

    return true;
}

bool load_parameter_file(QVariantMap& params, const QString& fname)
{
    QString json = read_text_file(fname);
    if (json.isEmpty()) {
        qWarning() << "Non-existent or empty parameter file: " + fname;
        return false;
    }
    QJsonParseError error;
    /// Witold I use toLatin1() everywhere. Is this the appropriate way to convert to byte array?
    QJsonObject obj = QJsonDocument::fromJson(json.toLatin1(), &error).object();
    if (error.error != QJsonParseError::NoError) {
        qWarning() << "Error parsing json file: " + fname + " : " + error.errorString();
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
    /// Witold In general is it possible to not display quotes around strings for qDebug?
    qDebug() << result.property("name").toString();
    qDebug() << result.property("message").toString();
    qDebug() << QString("%1 line %2").arg(result.property("fileName").toString()).arg(result.property("lineNumber").toInt());
}

bool run_script(const QStringList& script_fnames, const QVariantMap& params, QString& error_message)
{
    QJsonObject parameters = variantmap_to_json_obj(params);

    QJSEngine engine;
    ScriptController Controller;
    QJSValue MP = engine.newQObject(&Controller);
    engine.globalObject().setProperty("MP", MP);
    foreach (QString fname, script_fnames) {
        QJSValue result = engine.evaluate(read_text_file(fname), fname);
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

    return true;
}

void print_usage()
{
    printf("Usage:\n");
    printf("mountainprocess runProcess [processor_name] --[param1]=[val1] --[param2]=[val2] ...\n");
    printf("mountainprocess runScript [script1].js [script2.js] ... [file1].par [file2].par ... \n");
    printf("mountainprocess daemon-start\n");
    printf("mountainprocess daemon-stop\n");
    printf("mountainprocess daemon-restart\n");
    printf("mountainprocess daemon-state\n");
    printf("mountainprocess queue-script --~script_output=[optional_output_fname] [script1].js [script2.js] ... [file1].par [file2].par ... \n");
    printf("mountainprocess queue-process [processor_name] --~process_output=[optional_output_fname] --[param1]=[val1] --[param2]=[val2] ...\n");
}

void remove_system_parameters(QVariantMap& params)
{
    QStringList keys = params.keys();
    foreach (QString key, keys) {
        if (key.startsWith("~")) {
            params.remove(key);
        }
    }
}

int queue_pript(PriptType prtype, const CLParams& CLP)
{
    MPDaemonPript PP;

    bool detach = CLP.named_parameters.value("~detach", false).toBool();

    if (prtype == ScriptType) {
        QVariantMap params;
        for (int i = 0; i < CLP.unnamed_parameters.count(); i++) {
            QString str = CLP.unnamed_parameters[i];
            if (str.endsWith(".js")) {
                PP.script_paths << str;
                PP.script_path_checksums << compute_checksum_of_file(str);
            }
            if (str.endsWith(".par")) { // note that we can have multiple parameter files! the later ones override the earlier ones.
                if (!load_parameter_file(params, str)) {
                    return -1;
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

    PP.id = make_random_id();

    if (prtype == ScriptType) {
        PP.output_fname = CLP.named_parameters["~script_output"].toString();
    }
    else {
        PP.output_fname = CLP.named_parameters["~process_output"].toString();
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

    MPDaemonInterface X;
    if (prtype == ScriptType) {
        if (!X.queueScript(PP)) { //queue the script
            return -1;
        }
    }
    else {
        if (!X.queueProcess(PP)) { //queue the process
            return -1;
        }
    }
    if (!detach) {
        qint64 parent_pid = CLP.named_parameters.value("~parent_pid", 0).toLongLong();
        MPDaemon::waitForFileToAppear(PP.output_fname, -1, false, parent_pid, PP.stdout_fname);
        QJsonObject results_obj = QJsonDocument::fromJson(read_text_file(PP.output_fname).toLatin1()).object();
        bool success = results_obj["success"].toBool();
        if (!success) {
            if (prtype == ScriptType) {
                qWarning() << "Error in script " + PP.id + ": " + results_obj["error"].toString();
            }
            else {
                qWarning() << "Error in process " + PP.processor_name + " " + PP.id + ": " + results_obj["error"].toString();
            }
            return -1;
        }
    }
    return 0;
}
