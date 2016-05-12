/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
** Created: 5/2/2016
*******************************************************/

#ifndef MPDAEMON_H
#define MPDAEMON_H

#include <QObject>
#include <QString>
#include <QDateTime>
#include <QVariantMap>
#include <QJsonObject>
#include <QProcess>
#include <QFile>

struct ProcessResources {
    ProcessResources()
    {
        num_threads = 0;
        memory_gb = 0;
    }
    double num_threads;
    double memory_gb;
};

class MPDaemonPrivate;
class MPDaemon : public QObject {
    Q_OBJECT
public:
    friend class MPDaemonPrivate;
    MPDaemon();
    virtual ~MPDaemon();
    void setTotalResourcesAvailable(ProcessResources PR);
    void setLogPath(const QString& path);
    bool run();
    void clearProcessing();

    static QString daemonPath();
    static QString makeTimestamp(const QDateTime& dt = QDateTime::currentDateTime());
    static QDateTime parseTimestamp(const QString& timestamp);
    static bool waitForFileToAppear(QString fname, qint64 timeout_ms = -1, bool remove_on_appear = false, qint64 parent_pid = 0, QString stdout_fname = "");
    static void wait(qint64 msec);
    static bool pidExists(qint64 pid);
    static bool waitForFinishedAndWriteOutput(QProcess* P);

private
slots:
    void slot_commands_directory_changed();
    void slot_pript_qprocess_finished();
    void slot_qprocess_output();

private:
    MPDaemonPrivate* d;
};

struct ProcessRuntimeOpts {
    ProcessRuntimeOpts()
    {
        num_threads_allotted = 1;
        memory_gb_allotted = 1;
    }
    double num_threads_allotted;
    double memory_gb_allotted;
};

bool is_at_most(ProcessResources PR1, ProcessResources PR2);

enum PriptType {
    ScriptType,
    ProcessType
};

struct MPDaemonPript {
    //Represents a process or a script
    MPDaemonPript()
    {
        is_running = false;
        is_finished = false;
        success = false;
        parent_pid = 0;
        qprocess = 0;
        stdout_file = 0;
        prtype = ScriptType;
        num_threads_requested = 1;
        memory_gb_requested = 1;
    }
    PriptType prtype;
    QString id;
    QString output_fname;
    QString stdout_fname;
    QVariantMap parameters;
    bool is_running;
    bool is_finished;
    bool success;
    QString error;
    QJsonObject runtime_results;
    qint64 parent_pid;
    QDateTime timestamp_queued;
    QDateTime timestamp_started;
    QDateTime timestamp_finished;
    QProcess* qprocess;
    QFile* stdout_file;

    //For a script:
    QStringList script_paths;
    QStringList script_path_checksums; //to ensure that scripts have not changed at time of running

    //For a process:
    QString processor_name;
    double num_threads_requested;
    double memory_gb_requested;
    ProcessRuntimeOpts runtime_opts; //defined at run time
};

enum RecordType {
    AbbreviatedRecord,
    FullRecord,
    RuntimeRecord
};

QJsonObject pript_struct_to_obj(MPDaemonPript S, RecordType rt);
MPDaemonPript pript_obj_to_struct(QJsonObject obj);
QJsonArray stringlist_to_json_array(QStringList list);
QStringList json_array_to_stringlist(QJsonArray X);
QJsonObject variantmap_to_json_obj(QVariantMap map);
QVariantMap json_obj_to_variantmap(QJsonObject obj);
void mkdir_if_doesnt_exist(const QString& path);
QJsonObject runtime_opts_struct_to_obj(ProcessRuntimeOpts opts);

#endif // MPDAEMON_H
