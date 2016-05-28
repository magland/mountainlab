#include "mv_compute_templates_processor.h"

#include "mv_compute_templates_processor.h"
#include "mv_compute_templates.h"

class mv_compute_templates_ProcessorPrivate {
public:
    mv_compute_templates_Processor* q;
};

mv_compute_templates_Processor::mv_compute_templates_Processor()
{
    d = new mv_compute_templates_ProcessorPrivate;
    d->q = this;

    this->setName("mv_compute_templates");
    this->setVersion("0.12");
    this->setInputFileParameters("timeseries", "firings");
    this->setOutputFileParameters("templates", "stdevs");
    this->setRequiredParameters("clip_size");
}

mv_compute_templates_Processor::~mv_compute_templates_Processor()
{
    delete d;
}

bool mv_compute_templates_Processor::check(const QMap<QString, QVariant>& params)
{
    if (!this->checkParameters(params))
        return false;
    return true;
}

bool mv_compute_templates_Processor::run(const QMap<QString, QVariant>& params)
{
    QString timeseries = params["timeseries"].toString();
    QString firings = params["firings"].toString();
    QString templates = params["templates"].toString();
    QString stdevs = params["stdevs"].toString();
    int clip_size = params["clip_size"].toInt();
    return mv_compute_templates(timeseries, firings, templates, stdevs, clip_size);
}
