/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
** Created: 4/4/2016
*******************************************************/

#include "mvclipswidget.h"
#include "mvclipsview.h"
#include "taskprogress.h"
#include <QHBoxLayout>
#include "computationthread.h"
#include "mountainprocessrunner.h"
#include "msmisc.h"

class MVClipsWidgetComputer : public ComputationThread {
public:
    //input
    DiskReadMda firings;
    DiskReadMda timeseries;
    //QString mscmdserver_url;
    QString mlproxy_url;
    int clip_size;
    QList<int> labels_to_use;

    //output
    Mda clips;
    QList<double> times;

    void compute();
};

class MVClipsWidgetPrivate {
public:
    MVClipsWidget* q;
    //QString m_mscmdserver_url;
    QString m_mlproxy_url;
    DiskReadMda m_timeseries;
    DiskReadMda m_firings;
    QList<int> m_labels_to_use;
    int m_clip_size;
    MVClipsView* m_view;
    MVClipsWidgetComputer m_computer;

    void start_computation();
};

MVClipsWidget::MVClipsWidget()
{
    d = new MVClipsWidgetPrivate;
    d->q = this;

    d->m_view = MVClipsView::newInstance();

    QHBoxLayout* hlayout = new QHBoxLayout;
    hlayout->addWidget(d->m_view);
    this->setLayout(hlayout);

    connect(&d->m_computer, SIGNAL(computationFinished()), this, SLOT(slot_computation_finished()));
    connect(d->m_view, SIGNAL(currentEventChanged()), this, SIGNAL(currentEventChanged()));
}

MVClipsWidget::~MVClipsWidget()
{
    d->m_computer.stopComputation(); // important do take care of this before things start getting destructed!
    delete d;
}

/*
void MVClipsWidget::setMscmdServerUrl(const QString& url)
{
    d->m_mscmdserver_url = url;
}
*/

void MVClipsWidget::setMLProxyUrl(const QString& url)
{
    d->m_mlproxy_url = url;
}

void MVClipsWidget::setTimeseries(DiskReadMda& X)
{
    d->m_timeseries = X;
    d->start_computation();
}

void MVClipsWidget::setFirings(DiskReadMda& F)
{
    d->m_firings = F;
    d->start_computation();
}

void MVClipsWidget::setLabelsToUse(const QList<int>& labels)
{
    d->m_labels_to_use = labels;
    d->start_computation();
}

void MVClipsWidget::setClipSize(int clip_size)
{
    d->m_clip_size = clip_size;
    d->start_computation();
}

int MVClipsWidget::currentClipIndex()
{
    return d->m_view->currentClipIndex();
}

MVEvent MVClipsWidget::currentEvent()
{
    return d->m_view->currentEvent();
}

void MVClipsWidget::setCurrentEvent(MVEvent evt)
{
    d->m_view->setCurrentEvent(evt);
}

void MVClipsWidget::slot_computation_finished()
{
    d->m_computer.stopComputation(); //because I'm paranoid
    d->m_view->setClips(d->m_computer.clips);
    d->m_view->setTimes(d->m_computer.times);
}

void MVClipsWidgetComputer::compute()
{
    TaskProgress task("Clips Widget Computer");
    QString firings_out_path;
    {
        QString labels_str;
        foreach(int x, labels_to_use)
        {
            if (!labels_str.isEmpty())
                labels_str += ",";
            labels_str += QString("%1").arg(x);
        }

        MountainProcessRunner MT;
        QString processor_name = "mv_subfirings";
        MT.setProcessorName(processor_name);

        QMap<QString, QVariant> params;
        params["firings"] = firings.path();
        params["labels"] = labels_str;
        MT.setInputParameters(params);
        //MT.setMscmdServerUrl(mscmdserver_url);
        MT.setMLProxyUrl(mlproxy_url);

        firings_out_path = MT.makeOutputFilePath("firings_out");

        MT.runProcess(this);
        if (this->stopRequested()) {
            task.error(QString("Halted while running process: " + processor_name));
            return;
        }
    }

    QString clips_path;
    {
        MountainProcessRunner MT;
        QString processor_name = "extract_clips";
        MT.setProcessorName(processor_name);

        QMap<QString, QVariant> params;
        params["timeseries"] = timeseries.path();
        params["firings"] = firings_out_path;
        params["clip_size"] = clip_size;
        MT.setInputParameters(params);
        //MT.setMscmdServerUrl(mscmdserver_url);
        MT.setMLProxyUrl(mlproxy_url);

        clips_path = MT.makeOutputFilePath("clips");

        MT.runProcess(this);
        if (this->stopRequested()) {
            task.error(QString("Halted while running process: " + processor_name));
            return;
        }
    }
    task.log("Reading: " + firings_out_path);
    DiskReadMda firings_out(firings_out_path);
    task.log(QString("firings_out: %1 x %2").arg(firings_out.N1()).arg(firings_out.N2()));
    times.clear();
    for (long j = 0; j < firings_out.N2(); j++) {
        times << firings_out.value(1, j);
    }
    task.log("Reading: " + clips_path);
    DiskReadMda CC(clips_path);
    task.log(QString("CC: %1 x %2 x %3").arg(CC.N1()).arg(CC.N2()).arg(CC.N3()));
    CC.setHaltAgent(this);
    CC.readChunk(clips, 0, 0, 0, CC.N1(), CC.N2(), CC.N3());
    if (this->stopRequested()) {
        task.error(QString("Halted while reading chunk from: " + clips_path));
        return;
    }
}

void MVClipsWidgetPrivate::start_computation()
{
    m_computer.stopComputation();
    //m_computer.mscmdserver_url = m_mscmdserver_url;
    m_computer.mlproxy_url = m_mlproxy_url;
    m_computer.firings = m_firings;
    m_computer.timeseries = m_timeseries;
    m_computer.labels_to_use = m_labels_to_use;
    m_computer.clip_size = m_clip_size;
    m_computer.startComputation();
}
