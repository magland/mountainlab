/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
** Created: 3/25/2016
*******************************************************/

#include "closemehandler.h"
#include <QDebug>
#include <QTimer>
#include <QDateTime>
#include <QDir>
#include <QApplication>
#include <QFileInfo>
#include "mlutils.h"

class CloseMeHandlerPrivate {
public:
    CloseMeHandler* q;
    QDateTime m_start_time;
    void do_start();
};

CloseMeHandler::CloseMeHandler()
{
    d = new CloseMeHandlerPrivate;
    d->q = this;
}

CloseMeHandler::~CloseMeHandler()
{
    delete d;
}

void CloseMeHandler::start()
{
    CloseMeHandler* X = new CloseMeHandler();
    X->d->do_start();
}

void CloseMeHandler::slot_timer()
{
    QString fname = cfp(qApp->applicationDirPath() + "/closeme.tmp");
    if (QFile::exists(fname)) {
        QDateTime dt = QFileInfo(fname).created();
        if (dt > d->m_start_time) {
            exit(0);
        }
    }
    QTimer::singleShot(1000, this, SLOT(slot_timer()));
}

void CloseMeHandlerPrivate::do_start()
{
    m_start_time = QDateTime::currentDateTime();
    QTimer::singleShot(2000, q, SLOT(slot_timer()));
}
