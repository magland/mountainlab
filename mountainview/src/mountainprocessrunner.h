/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
** Created: 4/4/2016
*******************************************************/

#ifndef MOUNTAINSORTTHREAD_H
#define MOUNTAINSORTTHREAD_H

#include "computationthread.h"
#include "haltagent.h"

class MountainProcessRunnerPrivate;
class MountainProcessRunner {
public:
    friend class MountainProcessRunnerPrivate;
    MountainProcessRunner();
    virtual ~MountainProcessRunner();

    void setProcessorName(const QString& pname);
    void setInputParameters(const QMap<QString, QVariant>& parameters);
    void setMLProxyUrl(const QString& url);
    QString makeOutputFilePath(const QString& pname);
    void runProcess(HaltAgent* halt_agent);

private:
    MountainProcessRunnerPrivate* d;
};

#endif // MOUNTAINSORTTHREAD_H
