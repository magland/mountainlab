#ifndef MV_COMPUTE_TEMPLATES_PROCESSOR_H
#define MV_COMPUTE_TEMPLATES_PROCESSOR_H

#include "msprocessor.h"

class mv_compute_templates_ProcessorPrivate;
class mv_compute_templates_Processor : public MSProcessor
{
public:
    friend class mv_compute_templates_ProcessorPrivate;
    mv_compute_templates_Processor();
    virtual ~mv_compute_templates_Processor();

    bool check(const QMap<QString,QVariant> &params);
    bool run(const QMap<QString,QVariant> &params);
private:
    mv_compute_templates_ProcessorPrivate *d;
};

#endif // MV_COMPUTE_TEMPLATES_PROCESSOR_H

