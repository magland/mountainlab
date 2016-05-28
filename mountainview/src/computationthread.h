/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
** Created: 3/31/2016
*******************************************************/

#ifndef COMPUTATIONTHREAD
#define COMPUTATIONTHREAD

#include <QThread>
#include <QString>
#include <QMutex>

#include "haltagent.h"

class ComputationThreadPrivate;
class ComputationThread : public QThread, public HaltAgent {
    Q_OBJECT
public:
    friend class ComputationThreadPrivate;
    ComputationThread();
    virtual ~ComputationThread();

    virtual void compute() = 0;

    void setDeleteOnComplete(bool val);
    void startComputation(); //will stop existing computation
    bool stopComputation(int timeout = 0); //will wait for stop before returning, returns true if successfully stopped
    bool isComputing();
    bool isFinished();
    bool hasError();
    QString errorMessage();

    bool stopRequested();

signals:
    void computationFinished();

protected:
    void setErrorMessage(const QString& error);

private:
    void run();
private
slots:
    void slot_start();

private:
    ComputationThreadPrivate* d;
};

#endif // COMPUTATIONTHREAD
