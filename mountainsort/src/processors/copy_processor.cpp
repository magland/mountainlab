#include "copy_processor.h"

#include <QFile>
#include <QFileInfo>

class copy_ProcessorPrivate {
public:
    copy_Processor* q;
};

copy_Processor::copy_Processor()
{
    d = new copy_ProcessorPrivate;
    d->q = this;

    this->setName("copy");
    this->setVersion("0.1");
    this->setInputFileParameters("input");
    this->setOutputFileParameters("output");
}

copy_Processor::~copy_Processor()
{
    delete d;
}

bool copy_Processor::check(const QMap<QString, QVariant>& params)
{
    if (!this->checkParameters(params))
        return false;
    return true;
}

bool copy_Processor::run(const QMap<QString, QVariant>& params)
{
    QString input_path = params["input"].toString();
    QString output_path = params["output"].toString();
    if (!QFile::exists(input_path)) {
        qWarning() << "File does not exist:" << input_path;
        return false;
    }
    if (QFileInfo(input_path).canonicalFilePath() == QFileInfo(output_path).canonicalFilePath()) {
        return true;
    }
    if (QFile::exists(output_path)) {
        if (!QFile::remove(output_path)) {
            qWarning() << "Unable to remove file:" << output_path;
            return false;
        }
    }
    if (!QFile::copy(input_path, output_path)) {
        qWarning() << "Error copying file" << input_path << output_path;
        return false;
    }
    return true;
}
