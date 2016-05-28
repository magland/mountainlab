/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
** Created: 4/19/2016
*******************************************************/

#include "scriptcontroller.h"

#include <QCryptographicHash>
#include <QFile>
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

class ScriptControllerPrivate {
public:
    ScriptControllerPrivate()
    {
        m_nodaemon = false;
    }

    ScriptController* q;
    bool m_nodaemon;
    QStringList m_server_urls;
    QString m_server_base_path;

    QProcess* queue_process(QString processor_name, const QVariantMap& parameters, bool use_run = false);
    QProcess* run_process(QString processor_name, const QVariantMap& parameters);
    //static bool queue_process_and_wait_for_finished(QString processor_name, const QVariantMap& parameters);
    static void wait(qint64 msec);

    QString resolve_file_name_p(QString fname);
};

ScriptController::ScriptController()
{
    d = new ScriptControllerPrivate;
    d->q = this;
}

ScriptController::~ScriptController()
{
    delete d;
}

void ScriptController::setNoDaemon(bool val)
{
    d->m_nodaemon = val;
}

void ScriptController::setServerUrls(const QStringList& urls)
{
    d->m_server_urls = urls;
}

void ScriptController::setServerBasePath(const QString& path)
{
    d->m_server_base_path = path;
}

QString ScriptController::fileChecksum(const QString& fname_in)
{
    QString fname = d->resolve_file_name_p(fname_in);
    QTime timer;
    timer.start();
    printf("Computing checksum for file %s\n", fname.toLatin1().data());
    QFile file(fname);
    if (!file.open(QIODevice::ReadOnly))
        return "";
    QCryptographicHash hash(QCryptographicHash::Sha1);
    hash.addData(&file);
    file.close();
    QString ret = QString(hash.result().toHex());
    printf("%s -- Elapsed: %g sec\n", ret.toLatin1().data(), timer.elapsed() * 1.0 / 1000);
    return ret;
}

QString ScriptController::stringChecksum(const QString& str)
{
    QByteArray X = str.toLatin1();
    QCryptographicHash hash(QCryptographicHash::Sha1);
    hash.addData(X);
    return QString(hash.result().toHex());
}

QString ScriptController::createTemporaryFileName(const QString& code)
{
    return CacheManager::globalInstance()->makeLocalFile(code, CacheManager::LongTerm);
}

/*
bool ScriptController::runProcess(const QString& processor_name, const QString& parameters_json)
{
    QJsonObject params = QJsonDocument::fromJson(parameters_json.toLatin1()).object();
    QStringList keys = params.keys();
    QMap<QString, QVariant> parameters;
    foreach (QString key, keys) {
        parameters[key] = params[key].toVariant();
    }
    ProcessManager* PM = ProcessManager::globalInstance();
    if (!PM->checkParameters(processor_name, parameters)) {
        return false;
    }
    if (PM->processAlreadyCompleted(processor_name, parameters)) {
        this->log(QString("Process already completed: %1").arg(processor_name));
        return true;
    }

    if (d->queue_process_and_wait_for_finished(processor_name, parameters)) {
        return true;
    }
    else {
        return false;
    }
}
*/

struct PipelineNode {
    PipelineNode()
    {
        completed = false;
        running = false;
        qprocess = 0;
    }
    ~PipelineNode()
    {
        if (qprocess) {
            delete qprocess;
        }
    }
    QString processor_name;
    QVariantMap parameters;
    QStringList input_paths;
    QStringList output_paths;
    bool completed;
    bool running;
    QProcess* qprocess;
};

bool ScriptController::runPipeline(const QString& json)
{
    ProcessManager* PM = ProcessManager::globalInstance();

    //parse the json
    QJsonParseError parse_err;
    QJsonObject pipe_obj = QJsonDocument::fromJson(json.toLatin1(), &parse_err).object();
    if (parse_err.error != QJsonParseError::NoError) {
        qWarning() << "runPipeline: error parsing json: " + parse_err.errorString();
        return false;
    }
    QJsonArray processes = pipe_obj["processes"].toArray();

    //set up the node tree
    QList<PipelineNode> nodes;
    for (int i = 0; i < processes.count(); i++) {
        QJsonObject process = processes[i].toObject();
        PipelineNode node;
        node.parameters = json_obj_to_variantmap(process["parameters"].toObject());
        node.processor_name = process["processor_name"].toString();
        MLProcessor PP = PM->processor(node.processor_name);
        if (PP.name != node.processor_name) {
            qWarning() << "Unable to find processor: " + node.processor_name;
            return false;
        }
        QStringList input_pnames = PP.inputs.keys();
        QStringList output_pnames = PP.outputs.keys();
        foreach(QString pname, input_pnames)
        {
            QString path0 = node.parameters.value(pname).toString();
            if (!path0.isEmpty())
                node.input_paths << path0;
        }
        foreach(QString pname, output_pnames)
        {
            QString path0 = node.parameters.value(pname).toString();
            if (!path0.isEmpty())
                node.output_paths << path0;
        }
        nodes << node;
    }

    //record which outputs get created by which nodes (by index)
    QMap<QString, int> node_indices_for_outputs;
    for (int i = 0; i < nodes.count(); i++) {
        foreach(QString path, nodes[i].output_paths)
        {
            if (!path.isEmpty()) {
                if (node_indices_for_outputs.contains(path)) {
                    qWarning() << "Same output is created twice in pipeline.";
                    return false;
                }
            }
            node_indices_for_outputs[path] = i;
        }
    }

    //run the processing in parallel
    bool done = false;
    while (!done) {
        QList<int> node_indices_ready_to_be_run;
        QList<int> node_indices_not_ready_to_be_run;
        QList<int> node_indices_running;
        QList<int> node_indices_completed;
        for (int i = 0; i < nodes.count(); i++) {
            if (nodes[i].completed) {
                node_indices_completed << i;
            } else if (nodes[i].running) {
                node_indices_running << i;
            } else {
                bool ready_to_go = true;
                foreach(QString path, nodes[i].input_paths)
                {
                    if (node_indices_for_outputs.contains(path)) {
                        int ii = node_indices_for_outputs[path];
                        if (!nodes[ii].completed) {
                            ready_to_go = false;
                        }
                    }
                }
                if (ready_to_go)
                    node_indices_ready_to_be_run << i;
                else
                    node_indices_not_ready_to_be_run << i;
            }
        }
        if ((node_indices_running.isEmpty()) && (node_indices_ready_to_be_run.isEmpty()) && (node_indices_not_ready_to_be_run.isEmpty())) {
            //we are done!
            done = true;
        }
        if ((node_indices_running.isEmpty()) && (node_indices_ready_to_be_run.isEmpty()) && (!node_indices_not_ready_to_be_run.isEmpty())) {
            //Somehow we are not done, but nothing is ready and nothing is running.
            qWarning() << "Unable to run all processes in pipeline.";
            return false;
        }
        if (!node_indices_ready_to_be_run.isEmpty()) {
            //something is ready to run
            int ii = node_indices_ready_to_be_run[0];
            PipelineNode* node = &nodes[ii];
            if (!PM->checkParameters(node->processor_name, node->parameters)) {
                qWarning() << "Error checking parameters for processor: " + node->processor_name;
                return false;
            }
            if (PM->processAlreadyCompleted(node->processor_name, node->parameters)) {
                this->log(QString("Process already completed: %1").arg(node->processor_name));
                node->completed = true;
            } else {

                QProcess* P1;
                if (d->m_nodaemon) {
                    printf("Launching process %s\n", node->processor_name.toLatin1().data());
                    P1 = d->run_process(node->processor_name, node->parameters);
                    if (!P1) {
                        qWarning() << "Unable to launch process: " + node->processor_name;
                        return false;
                    }
                } else {
                    printf("Queuing process %s\n", node->processor_name.toLatin1().data());
                    P1 = d->queue_process(node->processor_name, node->parameters);
                    if (!P1) {
                        qWarning() << "Unable to queue process: " + node->processor_name;
                        return false;
                    }
                }

                node->running = true;
                node->qprocess = P1;
            }
        }
        //check for completed processes
        for (int i = 0; i < node_indices_running.count(); i++) {
            PipelineNode* node = &nodes[node_indices_running[i]];
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
                printf("Done.\n");
                node->completed = true;
                node->running = false;
                delete node->qprocess;
                node->qprocess = 0;
            }
        }
        if (!done) {
            MPDaemon::wait(100);
            qApp->processEvents(); //important I think for detecting when processes end.
        }
    }

    return true;
}

void ScriptController::log(const QString& message)
{
    printf("SCRIPT: %s\n", message.toLatin1().data());
}

QProcess* ScriptControllerPrivate::queue_process(QString processor_name, const QVariantMap& parameters, bool use_run)
{
    QString exe = qApp->applicationFilePath();
    QStringList args;
    if (use_run) {
        args << "run-process";
    } else {
        args << "queue-process";
    }
    args << processor_name;
    QStringList pkeys = parameters.keys();
    foreach(QString pkey, pkeys)
    {
        args << QString("--%1=%2").arg(pkey).arg(parameters[pkey].toString());
    }
    /// TODO switch all --~* parameters to --_* (more readable)
    args << QString("--~parent_pid=%1").arg(QCoreApplication::applicationPid());
    QProcess* P1 = new QProcess;
    P1->setReadChannelMode(QProcess::MergedChannels);
    P1->start(exe, args);
    if (!P1->waitForStarted()) {
        qWarning() << "Error waiting for process to start: " + processor_name;
        delete P1;
        return 0;
    }
    return P1;
}

QProcess* ScriptControllerPrivate::run_process(QString processor_name, const QVariantMap& parameters)
{
    return ScriptControllerPrivate::queue_process(processor_name, parameters, true);
}

/*
bool ScriptControllerPrivate::queue_process_and_wait_for_finished(QString processor_name, const QVariantMap& parameters)
{
    QString exe = qApp->applicationFilePath();
    QStringList args;
    args << "queue-process";
    args << processor_name;
    QStringList pkeys = parameters.keys();
    foreach (QString pkey, pkeys) {
        args << QString("--%1=%2").arg(pkey).arg(parameters[pkey].toString());
    }
    args << QString("--~parent_pid=%1").arg(QCoreApplication::applicationPid());
    QProcess P1;
    P1.setReadChannelMode(QProcess::MergedChannels);
    P1.start(exe, args);
    if (!MPDaemon::waitForFinishedAndWriteOutput(&P1)) {
        //if (!P1.waitForFinished(-1)) {
        printf("Error waiting for queue-process to finish: %s\n", processor_name.toLatin1().data());
        return false;
    }
    if (P1.exitStatus() == QProcess::CrashExit) {
        printf("Error -- queue-process crashed: %s\n", processor_name.toLatin1().data());
        return false;
    }
    if (P1.exitCode() != 0) {
        printf("Error -- queue-process returned non-zero exit code: %s\n", processor_name.toLatin1().data());
        return false;
    }
    return true;
}
*/

void ScriptControllerPrivate::wait(qint64 msec)
{
    usleep(msec * 1000);
}

QString ScriptControllerPrivate::resolve_file_name_p(QString fname_in)
{
    return resolve_file_name(m_server_urls, m_server_base_path, fname_in);
}

QString resolve_file_name(QStringList server_urls, QString server_base_path, QString fname_in)
{
    QString fname = fname_in;
    foreach(QString str, server_urls)
    {
        if (fname.startsWith(str + "/mdaserver"))
            fname = server_base_path + "/" + fname.mid((str + "/mdaserver").count());
    }
    if (!fname.startsWith(server_base_path)) {
        qWarning() << "Path does not start with " + server_base_path + ": " + fname;
        fname = "";
    }
    if (fname.mid(server_base_path.count()).contains("..")) {
        qWarning() << "Illegal .. in file path: " + fname;
        fname = "";
    }
    return fname;
}
