/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
** Created: 4/19/2016
*******************************************************/

#include "scriptcontroller2.h"

#include <QCryptographicHash>
#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include "cachemanager.h"
#include "processmanager.h"
#include <QTime>
#include <QCoreApplication>
#include <QDebug>
#include <unistd.h> //for usleep
#include "mpdaemon.h"
#include "mlcommon.h"

struct PipelineNode2 {
    // A node in the processing pipeline -- representing a single process
    PipelineNode2()
    {
        completed = false;
        running = false;
        qprocess = 0;
    }
    ~PipelineNode2()
    {
        if (qprocess) {
            delete qprocess;
        }
    }
    QMap<QString, QDateTime> timestamps;
    QString processor_name;
    QVariantMap inputs;
    QVariantMap parameters;
    QVariantMap outputs;
    bool create_prv = false;
    bool copy_file = false;
    bool remove_intermediate = false;
    bool completed;
    bool running;
    QString process_output_fname; //internal
    QProcess* qprocess;

    QStringList input_paths()
    {
        QStringList input_pnames = inputs.keys();
        QStringList ret;
        foreach (QString pname, input_pnames) {
            ret.append(MLUtil::toStringList(inputs[pname]));
        }
        return ret;
    }
    QStringList output_paths()
    {
        QStringList output_pnames = outputs.keys();
        QStringList ret;
        foreach (QString pname, output_pnames) {
            ret.append(MLUtil::toStringList(outputs[pname]));
        }
        return ret;
    }
};

// prv processes are special processing nodes that create a .prv file from an output file
QJsonArray get_prv_processes_2(const QList<PipelineNode2>& nodes, const QMap<QString, int>& node_indices_for_outputs, QString path, QSet<int>& node_indices_already_used, bool* ok);

QString create_temporary_path_for_output(QString processor_name, QVariantMap inputs, QVariantMap parameters, QString output_pname, int output_index = -1);

class ScriptController2Private {
public:
    ScriptController2Private()
    {
        m_nodaemon = false;
    }

    ScriptController2* q;
    bool m_nodaemon = false;
    //QStringList m_server_urls;
    //QString m_server_base_path;
    bool m_force_run = false;
    QString m_working_path;
    bool m_preserve_tempdir = false;
    QJsonObject m_results;
    int m_num_threads = 0;

    QList<PipelineNode2> m_pipeline_nodes;

    QProcess* queue_process(QString processor_name, const QVariantMap& parameters, bool use_run, bool force_run, bool preserve_tempdir, QString process_output_fname, int request_num_threads);
    QProcess* run_process(QString processor_name, const QVariantMap& parameters, bool force_run, bool preserve_tempdir, QString process_output_fname, int request_num_threads);

    void make_absolute_paths(QVariantMap& fnames);
    QString make_absolute_path(QString fname);
    QVariant make_absolute_path_variant(QVariant val);

    bool run_or_queue_node(PipelineNode2* node, const QMap<QString, int>& node_indices_for_outputs);
    PipelineNode2* find_node_ready_to_run();
    bool handle_running_processes();
    bool get_node_indices_for_outputs(QMap<QString, int>& node_indices_for_outputs);
    bool okay_to_remove_intermediate_file(const QString& path);
    bool create_rprv(const QString& path);
};

ScriptController2::ScriptController2()
{
    d = new ScriptController2Private;
    d->q = this;
}

ScriptController2::~ScriptController2()
{
    delete d;
}

void ScriptController2::setNoDaemon(bool val)
{
    d->m_nodaemon = val;
}

/*
void ScriptController2::setServerUrls(const QStringList& urls)
{
    d->m_server_urls = urls;
}

void ScriptController2::setServerBasePath(const QString& path)
{
    d->m_server_base_path = path;
}
*/

void ScriptController2::setForceRun(bool force_run)
{
    d->m_force_run = force_run;
}

void ScriptController2::setWorkingPath(QString working_path)
{
    d->m_working_path = working_path;
}

void ScriptController2::setPreserveTempdir(bool tempdir)
{
    d->m_preserve_tempdir = tempdir;
}

QJsonObject ScriptController2::getResults()
{
    return d->m_results;
}

QVariant filter_process_output(QVariant X, QString processor_name, QVariantMap inputs, QVariantMap parameters, QString pname)
{
    // Create temporary files for any output that is an empty string
    if (X.type() == QVariant::List) {
        QVariantList list = X.toList();
        for (int i = 0; i < list.count(); i++) {
            list[i] = filter_process_output(list[i], processor_name, inputs, parameters, QString("%1-%2").arg(pname).arg(i));
        }
        return list;
    }
    else {
        if (X.toString().isEmpty()) {
            X = create_temporary_path_for_output(processor_name, inputs, parameters, pname);
        }
        return X;
    }
}

// Add a process to the processing pipeline
QString ScriptController2::addProcess(QString processor_name, QString inputs_json, QString parameters_json, QString outputs_json)
{
    // Find the processor in the process manager
    MLProcessor PP = ProcessManager::globalInstance()->processor(processor_name);
    if (PP.name != processor_name) { //rather use PP.isNull()
        qWarning() << "Unable to find processor **: " + processor_name;
        return "{}";
    }
    QMap<QString, MLParameter> PP_outputs = PP.outputs; // The spec of the output parameters

    PipelineNode2 node;
    /// TODO: check for json parse errors here
    node.processor_name = processor_name;

    {
        QVariantMap inputs0 = QJsonDocument::fromJson(inputs_json.toLatin1()).object().toVariantMap();
        QStringList input_keys = PP.inputs.keys();
        foreach (QString key, input_keys) {
            node.inputs[key] = inputs0[key];
        }
    }
    {
        QVariantMap parameters0 = QJsonDocument::fromJson(parameters_json.toLatin1()).object().toVariantMap();
        QStringList parameter_keys = PP.parameters.keys();
        foreach (QString key, parameter_keys) {
            node.parameters[key] = parameters0[key];
        }
    }

    //node.inputs = QJsonDocument::fromJson(inputs_json.toLatin1()).object().toVariantMap();
    //node.parameters = QJsonDocument::fromJson(parameters_json.toLatin1()).object().toVariantMap();
    if (!outputs_json.isEmpty()) {
        node.outputs = QJsonDocument::fromJson(outputs_json.toLatin1()).object().toVariantMap();
        // check to see if outputs are consistent with PP_outputs
        QSet<QString> output_keys = PP.outputs.keys().toSet();
        QStringList node_output_keys = node.outputs.keys();
        foreach (QString key, node_output_keys) {
            if (!output_keys.contains(key)) {
                qWarning() << "Processor " + processor_name + " does not contain output: " + key;
                return "{}";
            }
        }
    }

    if ((outputs_json.isEmpty()) || (outputs_json == "\"\"")) {
        // any outputs that have not been specified will be marked as an empty string
        // later we will handle these cases by creating temporary files, and then returning this info from the function call
        QStringList output_pnames = PP_outputs.keys();
        foreach (QString pname, output_pnames) {
            if (!node.outputs.contains(pname)) {
                //node.outputs[pname] = d->create_temporary_path_for_output(node.processor_name, node.inputs, node.parameters, pname);
                node.outputs[pname] = "";
            }
        }
    }

    // Create temporary files for any output that is an empty string
    foreach (QString pname, node.outputs.keys()) {
        node.outputs[pname] = filter_process_output(node.outputs[pname], node.processor_name, node.inputs, node.parameters, pname);
    }

    d->make_absolute_paths(node.inputs);
    d->make_absolute_paths(node.outputs);
    d->m_pipeline_nodes << node;

    return QJsonDocument(QJsonObject::fromVariantMap(node.outputs)).toJson();
}

void ScriptController2::addPrv(QString input_path, QString output_path)
{
    PipelineNode2 node;
    node.inputs["input"] = d->make_absolute_path(input_path);
    node.outputs["output"] = d->make_absolute_path(output_path);
    node.create_prv = true;
    d->m_pipeline_nodes << node;
}

void ScriptController2::addCopyFile(QString input_path, QString output_path)
{
    PipelineNode2 node;
    node.inputs["input"] = d->make_absolute_path(input_path);
    node.outputs["output"] = d->make_absolute_path(output_path);
    node.copy_file = true;
    d->m_pipeline_nodes << node;
}

bool ScriptController2::runPipeline()
{
    QDateTime timestamp_start = QDateTime::currentDateTime();
    //check for empty output paths
    for (int i = 0; i < d->m_pipeline_nodes.count(); i++) {
        PipelineNode2* node = &d->m_pipeline_nodes[i];
        QStringList output_paths = node->output_paths();
        foreach (QString path, output_paths) {
            if (path.isEmpty()) {
                qWarning() << "Output path is empty" << node->processor_name << node->input_paths() << node->output_paths();
                return false;
            }
        }
    }

    //check for situation where input file matches output file
    for (int i = 0; i < d->m_pipeline_nodes.count(); i++) {
        PipelineNode2* node = &d->m_pipeline_nodes[i];
        QSet<QString> input_paths_set = node->input_paths().toSet();
        QStringList output_paths = node->output_paths();
        foreach (QString output_path, output_paths) {
            if (input_paths_set.contains(output_path)) {
                qWarning() << node->processor_name << output_path;
                qWarning() << "An input path is the same as an output path. This can happen sometimes when using .prv files (checksum lookups) in the case where a process creates an output file that matches an input file.";
                return false;
            }
        }
    }

    //record which outputs get created by which nodes (by index)
    QMap<QString, int> node_indices_for_outputs;
    if (!d->get_node_indices_for_outputs(node_indices_for_outputs)) {
        return false;
    }

    bool done = false;
    while (!done) {
        bool found = true;
        while (found) {
            found = false;
            PipelineNode2* node = d->find_node_ready_to_run();
            if (node) {
                if (!d->run_or_queue_node(node, node_indices_for_outputs)) {
                    return false;
                }
                found = true;
            }
        }

        if (!d->handle_running_processes()) {
            return false;
        }

        //check to see if we are done
        if (!d->find_node_ready_to_run()) {
            int num_running = 0;
            for (int i = 0; i < d->m_pipeline_nodes.count(); i++) {
                PipelineNode2* node = &d->m_pipeline_nodes[i];
                if (node->running)
                    num_running++;
            }
            if (num_running == 0)
                done = true; //if nothing is ready to run, and nothing is running, then we are done
        }

        if (!done) {
            //            MPDaemon::wait(100);
            usleep(100 * 1000);
            qApp->processEvents(); //important I think for detecting when processes end.
        }
    }

    //check whether everything got run
    for (int i = 0; i < d->m_pipeline_nodes.count(); i++) {
        PipelineNode2* node = &d->m_pipeline_nodes[i];
        if (!node->completed) {
            qWarning() << "Not every process in the pipeline was run. For example: " + node->processor_name + ". This could be due to a cyclic dependency.";
            return false;
        }
    }

    QDateTime timestamp_finish = QDateTime::currentDateTime();
    d->m_results["total_time_sec"] = timestamp_start.msecsTo(timestamp_finish) * 1.0 / 1000;
    return true;
}

void ScriptController2::setNumThreads(int num_threads)
{
    d->m_num_threads = num_threads;
}

QJsonObject make_prv_object_2(QString path)
{
    if (!QFile::exists(path)) {
        if (QFile::exists(path + ".rprv")) {
            QString json = TextFile::read(path + ".rprv");
            return QJsonDocument::fromJson(json.toUtf8()).object();
        }
        else {
            qWarning() << "Unable to find file (in make_prv_object_2):" << path;
            return QJsonObject();
        }
    }
    if (path.endsWith(".prv")) {
        //this important section added on 10/12/16 by jfm
        QString json = TextFile::read(path);
        return QJsonDocument::fromJson(json.toUtf8()).object();
    }
    return MLUtil::createPrvObject(path);
}

void ScriptController2::log(const QString& message)
{
    printf("SCRIPT: %s\n", message.toLatin1().data());
}

QString ScriptController2::dir(const QString& path)
{
    QStringList files = QDir(path).entryList(QStringList("*"), QDir::Files, QDir::Name);
    QStringList folders = QDir(path).entryList(QStringList("*"), QDir::Dirs | QDir::NoDotAndDotDot, QDir::Name);
    QJsonObject ret;
    {
        QJsonArray tmp;
        foreach (QString file, files) {
            QJsonObject obj;
            obj["name"] = file;
            tmp.push_back(obj);
        }
        ret["files"] = tmp;
    }
    {
        QJsonArray tmp;
        foreach (QString folder, folders) {
            QJsonObject obj;
            obj["name"] = folder;
            tmp.push_back(obj);
        }
        ret["folders"] = tmp;
    }
    return QJsonDocument(ret).toJson();
}

void ScriptController2::mkdir(const QString& path)
{
    QDir(QFileInfo(path).path()).mkpath(QFileInfo(path).fileName());
}

QString ScriptController2::createPrvObject(const QString& path)
{
    QMap<QString, int> node_indices_for_outputs;
    if (!d->get_node_indices_for_outputs(node_indices_for_outputs))
        return "";
    QJsonObject obj = make_prv_object_2(path);
    if (obj.isEmpty())
        return "";
    /*
    bool ok;
    QSet<int> node_indices_already_used; //to avoid infinite cycles, which can happen, for example, when an input file is the same as an output file
    obj["processes"] = get_prv_processes_2(d->m_pipeline_nodes, node_indices_for_outputs, path, node_indices_already_used, &ok);
    if (!ok) {
        qWarning() << "Error in get_prv_processes";
        return "";
    }
    */
    QString obj_json = QJsonDocument(obj).toJson(QJsonDocument::Indented);
    return obj_json;
}

bool ScriptController2::writeTextFile(const QString& path, const QString& text)
{
    return TextFile::write(path, text);
}

bool ScriptController2::fileExists(const QString& path)
{
    return QFile::exists(path);
}

void ScriptController2::removeFile(const QString& path)
{
    QFile::remove(path);
}

QProcess* ScriptController2Private::queue_process(QString processor_name, const QVariantMap& parameters, bool use_run, bool force_run, bool preserve_tempdir, QString process_output_fname, int request_num_threads)
{
    QString exe = qApp->applicationFilePath();
    QStringList args;
    if (use_run) {
        args << "run-process";
    }
    else {
        args << "queue-process";
    }
    args << processor_name;
    QStringList pkeys = parameters.keys();
    foreach (QString pkey, pkeys) {
        QStringList list = MLUtil::toStringList(parameters[pkey]);
        foreach (QString val, list) {
            args << QString("--%1=%2").arg(pkey).arg(val);
        }
    }
    //args << QString("--_parent_pid=%1").arg(QCoreApplication::applicationPid());
    if (force_run) {
        args << "--_force_run";
    }
    if (preserve_tempdir) {
        args << "--_preserve_tempdir";
    }
    if (!m_working_path.isEmpty()) {
        args << "--_working_path=" + m_working_path;
    }
    args << "--_process_output=" + process_output_fname;
    if (request_num_threads)
        args << QString("--_request_num_threads=%1").arg(request_num_threads);
    QProcess* P1 = new QProcess;
    P1->setReadChannelMode(QProcess::MergedChannels);
    //P1->start(exe, args);
    qDebug().noquote() << exe + " " + args.join(" ");
    MPDaemon::start_bash_command_and_kill_when_pid_is_gone(P1, exe, args, QCoreApplication::applicationPid());
    if (!P1->waitForStarted()) {
        qWarning() << "Error waiting for process to start: " + processor_name;
        delete P1;
        return 0;
    }
    return P1;
}

QProcess* ScriptController2Private::run_process(QString processor_name, const QVariantMap& parameters, bool force_run, bool preserve_tempdir, QString process_output_fname, int request_num_threads)
{
    return ScriptController2Private::queue_process(processor_name, parameters, true, force_run, preserve_tempdir, process_output_fname, request_num_threads);
}

void ScriptController2Private::make_absolute_paths(QVariantMap& fnames)
{
    QStringList pnames = fnames.keys();
    foreach (QString pname, pnames) {
        fnames[pname] = make_absolute_path_variant(fnames[pname]);
    }
}

QString ScriptController2Private::make_absolute_path(QString fname)
{
    if (fname.isEmpty())
        return "";

    if ((QDir::isRelativePath(fname)) && (!m_working_path.isEmpty())) {
        fname = m_working_path + "/" + fname;
    }
    return fname;
}

QVariant ScriptController2Private::make_absolute_path_variant(QVariant val)
{
    if (val.type() == QVariant::List) {
        QVariantList list = val.toList();
        for (int i = 0; i < list.count(); i++) {
            list[i] = make_absolute_path_variant(list[i]);
        }
        return list;
    }
    else {
        return make_absolute_path(val.toString());
    }
}

QVariant resolve_prv_files(QVariant f)
{
    if (f.type() == QVariant::List) {
        QVariantList ret;
        QVariantList list = f.toList();
        for (int i = 0; i < list.count(); i++) {
            ret << resolve_prv_files(list[i]);
        }
        return ret;
    }
    else {
        QString fname = f.toString();
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
    }
}

bool ScriptController2Private::run_or_queue_node(PipelineNode2* node, const QMap<QString, int>& node_indices_for_outputs)
{
    QVariantMap parameters0;
    {
        QStringList pnames = node->inputs.keys();
        foreach (QString pname, pnames) {
            parameters0[pname] = resolve_prv_files(node->inputs[pname]);
        }
    }
    {
        QStringList pnames = node->parameters.keys();
        foreach (QString pname, pnames) {
            parameters0[pname] = node->parameters[pname];
        }
    }
    {
        QStringList pnames = node->outputs.keys();
        foreach (QString pname, pnames) {
            if (!node->create_prv)
                parameters0[pname] = resolve_prv_files(node->outputs[pname]);
            else
                parameters0[pname] = node->outputs[pname];
        }
    }

    ProcessManager* PM = ProcessManager::globalInstance();
    if (node->create_prv) {
        QString input_path = node->inputs["input"].toString();
        QString output_path = node->outputs["output"].toString();
        if (input_path.isEmpty()) {
            qWarning() << "Input path for create_prv is empty.";
            return false;
        }
        if (!QFile::exists(input_path)) {
            qWarning() << "Input file for create_prv does not exist: " + input_path;
            return false;
        }
        if (!output_path.endsWith(".prv")) {
            qWarning() << "Problem creating .prv file. The file path must end with .prv.";
            return false;
        }

        if (!output_path.endsWith(".prv")) {
            qWarning() << ".prv file must end with .prv extension";
            return false;
        }
        QJsonObject obj = make_prv_object_2(input_path);
        if (obj.isEmpty())
            return false;
        bool ok;
        QSet<int> node_indices_already_used; //to avoid infinite cycles, which can happen, for example, when an input file is the same as an output file
        obj["processes"] = get_prv_processes_2(m_pipeline_nodes, node_indices_for_outputs, input_path, node_indices_already_used, &ok);
        if (!ok) {
            qWarning() << "Error in get_prv_processes";
            return false;
        }
        QString obj_json = QJsonDocument(obj).toJson(QJsonDocument::Indented);
        if (TextFile::write(output_path, obj_json)) {
            node->completed = true;
            return true;
        }
        else {
            qWarning() << "Unable to write prv file: " << output_path;
            return false;
        }
    }
    if (node->copy_file) {
        QString input_path = node->inputs["input"].toString();
        QString output_path = node->outputs["output"].toString();
        printf("COPYING FILE: %s -> %s", input_path.toUtf8().data(), output_path.toUtf8().data());
        if (input_path.isEmpty()) {
            qWarning() << "Input path for copy_file is empty.";
            return false;
        }
        if (!QFile::exists(input_path)) {
            qWarning() << "Input file for copy_file does not exist: " + input_path;
            return false;
        }

        if (QFile::exists(output_path)) {
            if (!QFile::remove(output_path)) {
                printf("Unable to remove output file: %s\n", output_path.toUtf8().data());
            }
        }
        if (!QFile::copy(input_path, output_path)) {
            printf("Unable to copy file: %s %s\n", input_path.toUtf8().data(), output_path.toUtf8().data());
            return false;
        }

        node->completed = true;
        return true;
    }
    if (node->remove_intermediate) {
        QString input_path = node->inputs["input"].toString();
        if (QFile::exists(input_path)) { //note that the file may not exist if the .rprv from a previous run already exists
            qDebug().noquote() << QString("Removing intermediate file: %1").arg(input_path);
            if (!create_rprv(input_path)) {
                qWarning() << "Unable to create .rprv file";
                return false;
            }
            if (!QFile::remove(input_path)) {
                qWarning() << "Unable to remove intermediate file: " + input_path;
                return false;
            }
        }
        node->completed = true;
        return true;
    }
    if (!PM->checkParameters(node->processor_name, parameters0)) {
        qWarning() << "Error checking parameters for processor: " + node->processor_name;
        return false;
    }
    if ((!m_force_run) && (PM->processAlreadyCompleted(node->processor_name, parameters0, true, true))) {
        q->log(QString("Process already completed: %1").arg(node->processor_name));
        node->completed = true;
        return true;
    }
    else {
        QProcess* P1;
        node->process_output_fname = CacheManager::globalInstance()->makeLocalFile() + ".process_output";
        if (m_nodaemon) {
            printf("Launching process %s\n", node->processor_name.toLatin1().data());
            P1 = run_process(node->processor_name, parameters0, m_force_run, m_preserve_tempdir, node->process_output_fname, m_num_threads);
            if (!P1) {
                qWarning() << "Unable to launch process: " + node->processor_name;
                return false;
            }
        }
        else {
            printf("Queuing process from script controller: %s\n", node->processor_name.toLatin1().data());
            P1 = queue_process(node->processor_name, parameters0, false, m_force_run, m_preserve_tempdir, node->process_output_fname, m_num_threads);
            if (!P1) {
                qWarning() << "Unable to queue process: " + node->processor_name;
                return false;
            }
        }

        node->running = true;
        node->qprocess = P1;
        return true;
    }
}

PipelineNode2* ScriptController2Private::find_node_ready_to_run()
{
    QSet<QString> file_paths_waiting_to_be_created;
    for (int i = 0; i < m_pipeline_nodes.count(); i++) {
        PipelineNode2* node = &m_pipeline_nodes[i];
        if (!node->completed) {
            QStringList pnames = node->outputs.keys();
            foreach (QString pname, pnames) {
                QStringList paths0 = MLUtil::toStringList(node->outputs[pname]);
                foreach (QString path0, paths0) {
                    file_paths_waiting_to_be_created.insert(path0);
                }
            }
        }
    }
    for (int i = 0; i < m_pipeline_nodes.count(); i++) {
        PipelineNode2* node = &m_pipeline_nodes[i];
        if ((!node->completed) && (!node->running)) {
            bool ready_to_run = true;
            QStringList pnames = node->inputs.keys();
            foreach (QString pname, pnames) {
                QStringList input_paths = MLUtil::toStringList(node->inputs[pname]);
                foreach (QString input_path, input_paths) {
                    if (file_paths_waiting_to_be_created.contains(input_path)) {
                        ready_to_run = false;
                    }
                }
            }
            if (ready_to_run) {
                if (node->remove_intermediate) {
                    if (okay_to_remove_intermediate_file(node->inputs["input"].toString())) {
                        return node;
                    }
                }
                else {
                    return node;
                }
            }
        }
    }
    return 0;
}

QString create_temporary_path_for_output(QString processor_name, QVariantMap inputs, QVariantMap parameters, QString output_pname, int output_index)
{
    QJsonObject obj;
    obj["processor_name"] = processor_name;
    obj["inputs"] = QJsonObject::fromVariantMap(inputs);
    obj["parameters"] = QJsonObject::fromVariantMap(parameters);
    obj["output_pname"] = output_pname;
    if (output_index >= 0)
        obj["output_index"] = output_index;
    QString json = QJsonDocument(obj).toJson();
    QString code = MLUtil::computeSha1SumOfString(json);
    QString str = "";
    if (output_index >= 0)
        str = QString("-%1").arg(output_index);
    return CacheManager::globalInstance()->makeIntermediateFile(code + "-" + processor_name + "-" + output_pname + str + ".tmp");
}

bool ScriptController2Private::handle_running_processes()
{
    for (int i = 0; i < m_pipeline_nodes.count(); i++) {
        PipelineNode2* node = &m_pipeline_nodes[i];
        if (node->running) {
            if (!node->qprocess) {
                qCritical() << "Unexpected problem" << __FILE__ << __LINE__;
                return false;
            }
            {
                QByteArray str = node->qprocess->readAll();
                if (!str.isEmpty()) {
                    printf("%s:: %s", node->processor_name.toLatin1().data(), str.data());
                }
            }
            if (node->qprocess->state() == QProcess::NotRunning) {
                printf("Process finished: %s\n", node->processor_name.toLatin1().data());
                node->qprocess->waitForReadyRead();
                QByteArray str = node->qprocess->readAll();
                if (!str.isEmpty()) {
                    printf("%s", str.data());
                }
                if (node->qprocess->exitStatus() == QProcess::CrashExit) {
                    qWarning() << "Process crashed: " + node->processor_name;
                    return false;
                }
                if (node->qprocess->exitCode() != 0) {
                    qWarning() << "Process returned with non-zero exit code: " + node->processor_name;
                    return false;
                }
                node->completed = true;
                node->running = false;

                if (!node->process_output_fname.isEmpty()) {
                    QString tmp_json = TextFile::read(node->process_output_fname);
                    if (tmp_json.isEmpty()) {
                        qWarning() << "process output file is empty or does not exist for processor: " + node->processor_name;
                    }
                    CacheManager::globalInstance()->setTemporaryFileDuration(node->process_output_fname, 600);
                    //QFile::remove(node->process_output_fname);
                    QJsonArray PP = m_results["processes"].toArray();
                    while (i >= PP.count())
                        PP.append(QJsonObject());
                    QJsonObject X;
                    X["processor_name"] = node->processor_name;
                    X["inputs"] = QJsonObject::fromVariantMap(node->inputs);
                    X["outputs"] = QJsonObject::fromVariantMap(node->outputs);
                    X["parameters"] = QJsonObject::fromVariantMap(node->parameters);
                    if (!tmp_json.isEmpty()) {
                        X["results"] = QJsonDocument::fromJson(tmp_json.toUtf8()).object();
                    }
                    PP[i] = X;
                    m_results["processes"] = PP;
                }

                delete node->qprocess;
                node->qprocess = 0;
            }
        }
    }
    return true;
}

bool ScriptController2Private::get_node_indices_for_outputs(QMap<QString, int>& node_indices_for_outputs)
{
    for (int i = 0; i < m_pipeline_nodes.count(); i++) {
        QStringList output_paths = m_pipeline_nodes[i].output_paths();
        foreach (QString path, output_paths) {
            if (!path.isEmpty()) {
                if (node_indices_for_outputs.contains(path)) {
                    qWarning() << "Same output is created twice in pipeline." << path;
                    return false;
                }
            }
            node_indices_for_outputs[path] = i;
        }
    }
    return true;
}

bool ScriptController2Private::okay_to_remove_intermediate_file(const QString& path)
{
    for (int i = 0; i < m_pipeline_nodes.count(); i++) {
        PipelineNode2* node = &m_pipeline_nodes[i];
        if (!node->completed) {
            if ((!node->remove_intermediate) && (node->input_paths().contains(path))) {
                return false;
            }
            if (node->output_paths().contains(path)) {
                return false;
            }
        }
    }
    return true;
}

bool ScriptController2Private::create_rprv(const QString& path)
{
    if (!QFile::exists(path)) {
        qWarning() << "In create_prv: file does not exist: " + path;
        return false;
    }
    QJsonObject obj = MLUtil::createPrvObject(path);

    //important for provenance tracking
    obj["original_last_modified"] = QFileInfo(path).lastModified().toString("yyyy-MM-dd-hh-mm-ss-zzz");

    QString json = QJsonDocument(obj).toJson();
    if (!TextFile::write(path + ".rprv", json)) {
        qWarning() << "Unable to write .rprv file: " + path + ".rprv";
        return false;
    }
    return true;
}

/*
QString resolve_file_name_2(QStringList server_urls, QString server_base_path, QString fname_in)
{
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

    return fname;
}
*/

QJsonArray get_prv_processes_2(const QList<PipelineNode2>& nodes, const QMap<QString, int>& node_indices_for_outputs, QString path, QSet<int>& node_indices_already_used, bool* ok)
{
    QJsonArray processes;
    int ind0 = node_indices_for_outputs.value(path, -1);
    if ((ind0 >= 0) && (!node_indices_already_used.contains(ind0))) { //avoid infinite cycles, which can happen, for example, if an output file is the same as an input file
        const PipelineNode2* node = &nodes[ind0];
        if (!node->processor_name.isEmpty()) {
            QJsonObject process;
            process["processor_name"] = node->processor_name;
            {
                QStringList input_pnames = node->inputs.keys();
                QJsonObject inputs;
                foreach (QString pname, input_pnames) {
                    QStringList list = MLUtil::toStringList(node->inputs[pname]);
                    if (list.count() == 1) {
                        bool is_list = (node->inputs[pname].type() == QVariant::List);
                        QJsonObject tmp;
                        tmp = make_prv_object_2(list[0]);
                        if (tmp.isEmpty()) {
                            *ok = false;
                            return processes;
                        }
                        inputs[pname] = tmp;
                        if (is_list) {
                            QJsonArray list2;
                            list2 << inputs[pname].toObject();
                            inputs[pname] = list2;
                        }
                    }
                    else {
                        QJsonArray list2;
                        foreach (QString str, list) {
                            QJsonObject tmp;
                            tmp = make_prv_object_2(str);
                            if (tmp.isEmpty()) {
                                *ok = false;
                                return processes;
                            }
                            list2.push_back(tmp);
                        }
                        inputs[pname] = list2;
                    }
                }
                process["inputs"] = inputs;
            }
            {
                QStringList output_pnames = node->outputs.keys();
                QJsonObject outputs;
                foreach (QString pname, output_pnames) {
                    QStringList list = MLUtil::toStringList(node->outputs[pname]);
                    if (list.count() == 1) {
                        bool is_list = (node->outputs[pname].type() == QVariant::List);
                        QJsonObject tmp;
                        tmp = make_prv_object_2(list[0]);
                        if (tmp.isEmpty()) {
                            *ok = false;
                            return processes;
                        }
                        outputs[pname] = tmp;
                        if (is_list) {
                            QJsonArray list2;
                            list2 << outputs[pname].toObject();
                            outputs[pname] = list2;
                        }
                    }
                    else {
                        QJsonArray list2;
                        foreach (QString str, list) {
                            QJsonObject tmp;
                            tmp = make_prv_object_2(str);
                            if (tmp.isEmpty()) {
                                *ok = false;
                                return processes;
                            }
                            list2.push_back(tmp);
                        }
                        outputs[pname] = list2;
                    }
                }
                process["outputs"] = outputs;
            }
            /*
            {

                QJsonObject outputs;
                QStringList output_pnames = node->outputs.keys();
                foreach (QString pname, output_pnames) {
                    QJsonObject tmp;
                    tmp = make_prv_object_2(node->outputs[pname].toString());
                    if (tmp.isEmpty()) {
                        *ok = false;
                        return processes;
                    }
                    outputs[pname] = tmp;
                }
                process["outputs"] = outputs;
            }
            */
            {
                QJsonObject parameters;
                QStringList pnames = node->parameters.keys();
                foreach (QString pname, pnames) {
                    parameters[pname] = node->parameters[pname].toString();
                }
                process["parameters"] = parameters;
            }
            node_indices_already_used.insert(ind0);
            processes.append(process);
            {
                QStringList input_pnames = node->inputs.keys();
                foreach (QString pname, input_pnames) {
                    QJsonArray X = get_prv_processes_2(nodes, node_indices_for_outputs, node->inputs[pname].toString(), node_indices_already_used, ok);
                    if (!ok)
                        return processes;
                    for (int a = 0; a < X.count(); a++) {
                        processes.append(X[a]);
                    }
                }
            }
        }
    }
    *ok = true;
    return processes;
}

void ScriptController2::removeIntermediate(const QString& path)
{
    PipelineNode2 node;
    node.inputs["input"] = d->make_absolute_path(path);
    node.remove_intermediate = true;
    d->m_pipeline_nodes << node;
}
