/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
** Created: 4/27/2016
*******************************************************/

#ifndef PROCESSMANAGER_H
#define PROCESSMANAGER_H

#include <QString>
#include <QVariant>
#include <QProcess>
#include <QDateTime>
#include <QJsonObject>

struct RequestProcessResources {
    int request_num_threads = 0;
};

struct MLParameter {
    QString name;
    QString ptype;
    QString description;
    bool optional;
    QVariant default_value;
};

struct MLProcessor {
    QString name;
    QString version;
    QString description;
    QMap<QString, MLParameter> inputs;
    QMap<QString, MLParameter> outputs;
    QMap<QString, MLParameter> parameters;
    QString exe_command;
    QJsonObject spec;

    QString basepath;
};

struct MonitorStats {
    QDateTime timestamp;
    int mem_bytes = 0;
    double cpu_pct = 0;
};

struct MLProcessInfo {
    QDateTime start_time;
    QDateTime finish_time;
    QString processor_name;
    QVariantMap parameters;
    QString exe_command;
    bool finished;
    int exit_code;
    QList<MonitorStats> monitor_stats;
    QProcess::ExitStatus exit_status;
    QByteArray standard_output;
    QByteArray standard_error;
};

class ProcessManagerPrivate;
class ProcessManager : public QObject {
    Q_OBJECT
public:
    friend class ProcessManagerPrivate;
    ProcessManager();
    virtual ~ProcessManager();

    //void setServerUrls(const QStringList& urls);
    //void setServerBasePath(const QString& path);

    void setProcessorPaths(const QStringList& paths);
    void reloadProcessors();

    QStringList processorNames() const;
    MLProcessor processor(const QString& name);

    bool checkParameters(const QString& processor_name, const QVariantMap& parameters);
    void setDefaultParameters(const QString& processor_name, QVariantMap& parameters);
    bool processAlreadyCompleted(const QString& processor_name, const QVariantMap& parameters, bool allow_rprv_inputs = true, bool allow_rprv_outputs = false);
    QString startProcess(const QString& processor_name, const QVariantMap& parameters, const RequestProcessResources& RPR, bool exec_mode, bool preserve_tempdir); //returns the process id/handle (a random string)
    bool waitForFinished(const QString& process_id, int parent_pid);
    MLProcessInfo processInfo(const QString& id);
    void clearProcess(const QString& id);
    void clearAllProcesses();

    QStringList allProcessIds() const;

    bool isFinished(const QString& id);

    void cleanUpCompletedProcessRecords();

    static ProcessManager* globalInstance();

private:
    bool loadProcessors(const QString& path, bool recursive = true);
    bool loadProcessorFile(const QString& path);

signals:
    void processFinished(QString id);
private slots:
    void slot_process_finished();
    void slot_qprocess_output();
    void slot_monitor();

private:
    ProcessManagerPrivate* d;
};

#endif // PROCESSMANAGER_H
