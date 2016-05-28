#include "mvstatusbar.h"
#include "taskprogress.h"

#include <QHBoxLayout>
#include <QLabel>

class MVStatusBarPrivate {
public:
    MVStatusBar* q;
    TaskProgressAgent* m_tp_agent;
    QLabel m_bytes_downloaded_label;
    QLabel m_tasks_running_label;
    QLabel m_remote_processing_time_label;
    QLabel m_bytes_allocated_label;

    void update_font();
};

MVStatusBar::MVStatusBar()
{
    d = new MVStatusBarPrivate;
    d->q = this;

    QHBoxLayout* hlayout1 = new QHBoxLayout;
    hlayout1->setSpacing(10);
    hlayout1->setMargin(0);
    hlayout1->addWidget(&d->m_bytes_downloaded_label);
    hlayout1->addWidget(&d->m_tasks_running_label);
    hlayout1->addWidget(&d->m_remote_processing_time_label);
    hlayout1->addStretch();

    QHBoxLayout* hlayout2 = new QHBoxLayout;
    hlayout2->setSpacing(10);
    hlayout2->setMargin(0);
    hlayout2->addWidget(&d->m_bytes_allocated_label);
    hlayout2->addStretch();

    QVBoxLayout *layout=new QVBoxLayout;
    layout->setSpacing(0);
    layout->setMargin(0);
    layout->addLayout(hlayout1,1);
    layout->addLayout(hlayout2,1);
    layout->setAlignment(Qt::AlignLeft);
    setLayout(layout);

    d->m_tp_agent = TaskProgressAgent::globalInstance();
    connect(d->m_tp_agent, SIGNAL(quantitiesChanged()), this, SLOT(slot_update_quantities()));
    connect(d->m_tp_agent, SIGNAL(tasksChanged()), this, SLOT(slot_update_tasks()));
}

MVStatusBar::~MVStatusBar()
{
    delete d;
}

QString format_num_bytes(double num_bytes)
{
    if (num_bytes < 100)
        return QString("%1 bytes").arg(num_bytes);
    if (num_bytes < 1000 * 1e3)
        return QString("%1 KB").arg(num_bytes / 1e3, 0, 'f', 2);
    if (num_bytes < 1000 * 1e6)
        return QString("%1 MB").arg(num_bytes / 1e6, 0, 'f', 2);
    return QString("%1 GB").arg(num_bytes / 1e9, 0, 'f', 2);
}

QString format_duration(double msec)
{
    if (msec < 1e3)
        return QString("%1 msec").arg(msec);
    return QString("%1 sec").arg(msec / 1000, 0, 'f', 1);
}

void MVStatusBar::slot_update_quantities()
{
    {
        QString txt = QString("%1 downloaded |").arg(format_num_bytes(d->m_tp_agent->getQuantity("bytes_downloaded")));
        d->m_bytes_downloaded_label.setText(txt);
    }
    {
        QString txt = QString("%1 remote processing").arg(format_duration(d->m_tp_agent->getQuantity("remote_processing_time")));
        d->m_remote_processing_time_label.setText(txt);
    }
    {
        double using_bytes=d->m_tp_agent->getQuantity("bytes_allocated")-d->m_tp_agent->getQuantity("bytes_freed");
        double bytes_read=d->m_tp_agent->getQuantity("bytes_read");
        QString txt = QString("%1 RAM | %2 Read").arg(format_num_bytes(using_bytes)).arg(format_num_bytes(bytes_read));
        d->m_bytes_allocated_label.setText(txt);
    }
}

void MVStatusBar::slot_update_tasks()
{
    d->update_font();
    {
        QString txt = QString("%1 tasks running |").arg(d->m_tp_agent->activeTasks().count());
        d->m_tasks_running_label.setText(txt);
    }
}

void MVStatusBarPrivate::update_font()
{
    QFont fnt = q->font();
    fnt.setPixelSize(qMax(9, q->height()/2 - 4));
    q->setFont(fnt);
}
