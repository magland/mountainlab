#include "msprocessmanager.h"
#include "mountainsort_version.h"
#include "msprocessor.h"
#include <QList>
#include <QTime>
#include <QCoreApplication>
#include <QDir>
#include <QCryptographicHash>

#include "example_processor.h"
#include "bandpass_filter_processor.h"
#include "whiten_processor.h"
#include "detect_processor.h"
#include "adjust_times_processor.h"
#include "branch_cluster_v2_processor.h"
#include "remove_duplicate_clusters_processor.h"
#include "compute_outlier_scores_processor.h"
#include "compute_detectability_scores_processor.h"
#include "copy_processor.h"
#include "mda2txt_processor.h"
#include "mask_out_artifacts_processor.h"
#include "fit_stage_processor.h"
#include "compute_templates_processor.h"
#include "msmisc.h"
#include "mv_firings_filter_processor.h"
#include "mv_subfirings_processor.h"
#include "mv_compute_templates_processor.h"
#include "extract_clips_processor.h"
#include "extract_clips_features_processor.h"
#include "merge_labels_processor.h"
#include "filter_events_processor.h"
#include "confusion_matrix_processor.h"
#include "extract_raw_processor.h"
#include "merge_across_channels_processor.h"
#include "geom2adj_processor.h"
#include "mlutils.h"

#include "textfile.h"
#include <sys/stat.h>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>

class MSProcessManagerPrivate {
public:
    MSProcessManager* q;
    QMap<QString, MSProcessor*> m_processors;
    example_Processor m_dummy_processor;

    MSProcessor* find_processor(const QString& name);
    QVariantMap compute_process_info(const QString& processor_name, const QVariantMap& parameters);
};

MSProcessManager::MSProcessManager()
{
    d = new MSProcessManagerPrivate;
    d->q = this;
}

Q_GLOBAL_STATIC(MSProcessManager, theInstance)
MSProcessManager* MSProcessManager::globalInstance()
{
    return theInstance;
}

MSProcessManager::~MSProcessManager()
{
    qDeleteAll(d->m_processors);
    delete d;
}

void MSProcessManager::loadDefaultProcessors()
{
    loadProcessor(new example_Processor);
    loadProcessor(new bandpass_filter_Processor);
    loadProcessor(new whiten_Processor);
    loadProcessor(new detect_Processor);
    loadProcessor(new adjust_times_Processor);
    loadProcessor(new branch_cluster_v2_Processor);
    loadProcessor(new remove_duplicate_clusters_Processor);
    loadProcessor(new compute_outlier_scores_Processor);
    loadProcessor(new compute_detectability_scores_Processor);
    loadProcessor(new copy_Processor);
    loadProcessor(new mda2txt_Processor);
    loadProcessor(new mask_out_artifacts_Processor);
    loadProcessor(new fit_stage_Processor);
    loadProcessor(new compute_templates_Processor);
    loadProcessor(new mv_compute_templates_Processor);
    loadProcessor(new mv_firings_filter_Processor);
    loadProcessor(new mv_subfirings_Processor);
    loadProcessor(new extract_clips_Processor);
    loadProcessor(new extract_clips_features_Processor);
    loadProcessor(new merge_labels_Processor);
    loadProcessor(new filter_events_Processor);
    loadProcessor(new confusion_matrix_Processor);
    loadProcessor(new extract_raw_Processor);
    loadProcessor(new merge_across_channels_Processor);
    loadProcessor(new geom2adj_Processor);
}

bool MSProcessManager::containsProcessor(const QString& processor_name) const
{
    return d->m_processors.contains(processor_name);
}

bool MSProcessManager::checkProcess(const QString& processor_name, const QVariantMap& parameters) const
{
    return d->find_processor(processor_name)->check(parameters);
}

QString compute_process_code(const QString& processor_name, const QVariantMap& parameters)
{
    QJsonObject X;
    X["processor_name"] = processor_name;
    X["parameters"] = QJsonObject::fromVariantMap(parameters);
    QString json = QJsonDocument(X).toJson();
    /// Witold I need a canonical json here so that the hash is always the same. (relatively important)
    return compute_hash(json);
}

QString compute_file_code(const QString& path)
{
    //the code comprises the device,inode,size, and modification time (in seconds)
    //note that it is not dependent on the file name
    struct stat SS;
    stat(path.toLatin1().data(), &SS);
    QString id_string = QString("%1:%2:%3:%4").arg(SS.st_dev).arg(SS.st_ino).arg(SS.st_size).arg(SS.st_mtim.tv_sec);
    return id_string;
}

QVariantMap MSProcessManagerPrivate::compute_process_info(const QString& processor_name, const QVariantMap& parameters)
{
    QVariantMap ret;
    QString version;
    QStringList input_file_parameters;
    QStringList output_file_parameters;
    if (MSProcessor* PP = m_processors.value(processor_name)) {
        version = PP->version();
        input_file_parameters = PP->inputFileParameters();
        output_file_parameters = PP->outputFileParameters();
    } else
        return ret; //can't even find the processor (not registered)

    QVariantMap input_file_codes;
    foreach(const QString & pp, input_file_parameters)
    {
        QString path0 = parameters[pp].toString();
        if (!path0.isEmpty()) {
            QString code0 = compute_file_code(path0);
            input_file_codes[path0] = code0;
        }
    }
    QVariantMap output_file_codes;
    foreach(const QString & pp, output_file_parameters)
    {
        QString path0 = parameters[pp].toString();
        if (!path0.isEmpty()) {
            QString code0 = compute_file_code(path0);
            output_file_codes[path0] = code0;
        }
    }
    ret["processor_name"] = processor_name;
    ret["version"] = version;
    ret["parameters"] = parameters;
    ret["input_file_codes"] = input_file_codes;
    ret["output_file_codes"] = output_file_codes;
    return ret;
}

bool MSProcessManager::runProcess(const QString& processor_name, const QVariantMap& parameters_in)
{
    QVariantMap parameters = parameters_in;

    printf("RUNNING %s\n", processor_name.toLatin1().data());
    QTime timer;
    timer.start();
    MSProcessor* processor = d->find_processor(processor_name);
    if (!processor) {
        qWarning() << "Unable to find processor [202]" << processor_name;
        return false;
    }
    bool ret = processor->run(parameters);
    if (ret) {
        printf("Elapsed time for processor %s: %g sec\n", processor_name.toLatin1().data(), timer.elapsed() * 1.0 / 1000);
    } else {
        qWarning() << "Error in processor->run" << processor_name;
    }

    return ret;
}

bool MSProcessManager::checkAndRunProcess(const QString& processor_name, const QVariantMap& parameters)
{
    if (!this->containsProcessor(processor_name)) {
        printf("Unable to find processor: %s\n", processor_name.toLatin1().data());
        return false;
    }
    if (!this->checkProcess(processor_name, parameters)) {
        printf("Problem checking processor: %s\n", processor_name.toLatin1().data());
        return false;
    }
    {
        if (!this->runProcess(processor_name, parameters)) {
            printf("Problem running processor: %s\n", processor_name.toLatin1().data());
            return 0;
        } else
            return true;
    }
}

MSProcessor* MSProcessManager::processor(const QString& processor_name)
{
    return d->find_processor(processor_name);
}

void MSProcessManager::loadProcessor(MSProcessor* P)
{
    if (d->m_processors.contains(P->name())) {
        qWarning() << "Processor with this name has already been loaded:" << P->name();
        return;
    }
    d->m_processors[P->name()] = P;
}

QStringList MSProcessManager::allProcessorNames() const
{
    return d->m_processors.keys();
}

QString MSProcessManager::usageString() const
{
    QString str;
    str += QString("MountainSort version %1\n").arg(mountainsort_version());

    str += QString("Available processors:\n");
    QStringList names = allProcessorNames();
    str += names.join(", ");
    str += "\n";

    return str;
}

QString tostr(const QStringList& list)
{
    return list.join(" ");
}

void MSProcessManager::printDetails() const
{
    QStringList names = allProcessorNames();
    for (int i = 0; i < names.count(); i++) {
        QString name = names[i];
        MSProcessor* P = d->m_processors[name];
        QString str;
        str += QString("***** Processor %1 *****").arg(P->name()) + "\n";
        str += QString("    Input files: %1").arg(tostr(P->inputFileParameters())) + "\n";
        str += QString("    Output files: %1").arg(tostr(P->outputFileParameters())) + "\n";
        str += QString("    Required params: %1").arg(tostr(P->requiredParameters())) + "\n";
        str += QString("    Optional params: %1").arg(tostr(P->optionalParameters())) + "\n";
        printf("%s\n", str.toLatin1().data());
    }
}

void MSProcessManager::printJsonSpec() const
{
    QJsonArray processors;
    QStringList names = allProcessorNames();
    for (int i = 0; i < names.count(); i++) {
        QString name = names[i];
        MSProcessor* P = d->m_processors[name];

        QJsonArray inputs;
        {
            QStringList input_file_parameters = P->inputFileParameters();
            foreach(QString pname, input_file_parameters)
            {
                QJsonObject P;
                P["name"] = pname;
                inputs.append(P);
            }
        }

        QJsonArray outputs;
        {
            QStringList output_file_parameters = P->outputFileParameters();
            foreach(QString pname, output_file_parameters)
            {
                QJsonObject P;
                P["name"] = pname;
                outputs.append(P);
            }
        }

        QJsonArray parameters;
        {
            QStringList parameters0 = P->requiredParameters();
            foreach(QString pname, parameters0)
            {
                QJsonObject P;
                P["name"] = pname;
                P["optional"] = false;
                parameters.append(P);
            }
        }
        {
            QStringList parameters0 = P->optionalParameters();
            foreach(QString pname, parameters0)
            {
                QJsonObject P;
                P["name"] = pname;
                P["optional"] = true;
                parameters.append(P);
            }
        }

        QJsonObject obj;
        obj["name"] = P->name();
        obj["version"] = P->version();
        obj["description"] = P->description();
        obj["inputs"] = inputs;
        obj["outputs"] = outputs;
        obj["parameters"] = parameters;
        QString exe_command = QString("%1 %2 $(arguments)").arg(qApp->applicationFilePath()).arg(P->name());
        obj["exe_command"] = exe_command;

        processors.append(obj);
    }
    QJsonObject X;
    X["processors"] = processors;
    QString json = QJsonDocument(X).toJson();
    printf("%s", json.toLatin1().data());
}

MSProcessor* MSProcessManagerPrivate::find_processor(const QString& name)
{
    if (!m_processors.contains(name)) {
        qWarning() << "Unable to find processor: " << name;
        return &m_dummy_processor;
    }
    return m_processors[name];
}
