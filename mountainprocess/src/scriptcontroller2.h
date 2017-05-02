/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
** Created: 4/27/2016
*******************************************************/

#ifndef SCRIPTCONTROLLER2_H
#define SCRIPTCONTROLLER2_H

#include <QObject>

class ScriptController2Private;
class ScriptController2 : public QObject {
    Q_OBJECT
public:
    friend class ScriptController2Private;
    ScriptController2();
    virtual ~ScriptController2();
    void setNoDaemon(bool val);
    void setServerUrls(const QStringList& urls);
    void setServerBasePath(const QString& path);
    void setForceRun(bool force_run);
    void setWorkingPath(QString working_path);
    void setPreserveTempdir(bool tempdir);
    QJsonObject getResults();

    Q_INVOKABLE QString addProcess(QString processor_name, QString inputs_json, QString parameters_json, QString outputs_json); //returns json
    Q_INVOKABLE void addPrv(QString input_path, QString output_path);
    Q_INVOKABLE void addCopyFile(QString input_path, QString output_path);
    Q_INVOKABLE void removeIntermediate(const QString& path);
    Q_INVOKABLE bool runPipeline();
    Q_INVOKABLE void setNumThreads(int num_threads);
    Q_INVOKABLE void log(const QString& message);
    Q_INVOKABLE QString dir(const QString& path);
    Q_INVOKABLE void mkdir(const QString& path);
    Q_INVOKABLE QString createPrvObject(const QString& path);
    Q_INVOKABLE bool writeTextFile(const QString& path, const QString& text);
    Q_INVOKABLE bool fileExists(const QString& path);
    Q_INVOKABLE void removeFile(const QString& path);

private:
    ScriptController2Private* d;
};

QString resolve_file_name_2(QStringList server_urls, QString server_base_path, QString fname_in);
QJsonObject make_prv_object_2(QString path);

#endif // SCRIPTCONTROLLER2_H
