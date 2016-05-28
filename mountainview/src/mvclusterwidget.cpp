#include "mvclusterwidget.h"
#include "mvclusterview.h"
#include "computationthread.h"
#include <QCheckBox>
#include <QHBoxLayout>
#include <QLabel>
#include <QList>
#include "mvclipsview.h"
#include "msmisc.h"
#include <math.h>
#include "extract_clips.h"
#include "mountainprocessrunner.h"

class MVClusterWidgetComputer : public ComputationThread {
public:
    //input
    //QString mscmdserver_url;
    QString mlproxy_url;
    DiskReadMda timeseries;
    DiskReadMda firings;
    int clip_size;
    QList<int> labels_to_use;

    //output
    Mda data;
    DiskReadMda firings_subset;

    void compute();
};

class MVClusterWidgetPrivate {
public:
    MVClusterWidget* q;
    QList<MVClusterView*> m_views;
    MVClipsView* m_clips_view;
    QLabel* m_info_bar;
    Mda m_data;
    //QString m_mscmdserver_url;
    QString m_mlproxy_url;
    DiskReadMda m_timeseries;
    DiskReadMda m_firings;
    QList<int> m_labels_to_use;
    int m_clip_size;
    QList<double> m_outlier_scores;
    MVClusterWidgetComputer m_computer;

    void connect_view(MVClusterView* V);
    void update_clips_view();
    int current_event_index();
    void set_data_on_visible_views_that_need_it();
    void start_computation();
};

MVClusterWidget::MVClusterWidget()
{
    d = new MVClusterWidgetPrivate;
    d->q = this;

    d->m_clip_size = 200;

    d->m_clips_view = MVClipsView::newInstance();

    {
        MVClusterView* X = new MVClusterView;
        X->setMode(MVCV_MODE_HEAT_DENSITY);
        d->m_views << X;
    }
    {
        MVClusterView* X = new MVClusterView;
        X->setMode(MVCV_MODE_LABEL_COLORS);
        d->m_views << X;
    }
    {
        MVClusterView* X = new MVClusterView;
        X->setMode(MVCV_MODE_TIME_COLORS);
        X->setVisible(false);
        d->m_views << X;
    }
    {
        MVClusterView* X = new MVClusterView;
        X->setMode(MVCV_MODE_AMPLITUDE_COLORS);
        X->setVisible(false);
        d->m_views << X;
    }

    QVBoxLayout* mainlayout = new QVBoxLayout;
    this->setLayout(mainlayout);

    QHBoxLayout* bottom_panel = new QHBoxLayout;
    {
        QCheckBox* CB = new QCheckBox;
        CB->setText("Clip View");
        connect(CB, SIGNAL(toggled(bool)), this, SLOT(slot_show_clip_view_toggled(bool)));
        CB->setChecked(true);
        bottom_panel->addWidget(CB);
    }
    {
        QCheckBox* CB = new QCheckBox;
        CB->setText("Density Plot");
        CB->setProperty("view_index", 0);
        connect(CB, SIGNAL(toggled(bool)), this, SLOT(slot_show_view_toggled(bool)));
        CB->setChecked(true);
        bottom_panel->addWidget(CB);
    }
    {
        QCheckBox* CB = new QCheckBox;
        CB->setText("Label Colors");
        CB->setProperty("view_index", 1);
        connect(CB, SIGNAL(toggled(bool)), this, SLOT(slot_show_view_toggled(bool)));
        CB->setChecked(true);
        bottom_panel->addWidget(CB);
    }
    {
        QCheckBox* CB = new QCheckBox;
        CB->setText("Time Colors");
        CB->setProperty("view_index", 2);
        connect(CB, SIGNAL(toggled(bool)), this, SLOT(slot_show_view_toggled(bool)));
        CB->setChecked(false);
        bottom_panel->addWidget(CB);
    }
    {
        QCheckBox* CB = new QCheckBox;
        CB->setText("Amplitude Colors");
        CB->setProperty("view_index", 3);
        connect(CB, SIGNAL(toggled(bool)), this, SLOT(slot_show_view_toggled(bool)));
        CB->setChecked(false);
        bottom_panel->addWidget(CB);
    }
    bottom_panel->addStretch();

    QHBoxLayout* hlayout = new QHBoxLayout;
    mainlayout->addLayout(hlayout);
    mainlayout->addLayout(bottom_panel);

    QVBoxLayout* vlayout = new QVBoxLayout;
    vlayout->addWidget(d->m_clips_view);
    d->m_info_bar = new QLabel;
    d->m_info_bar->setFixedHeight(20);
    vlayout->addWidget(d->m_info_bar);
    //d->m_clips_view->setFixedWidth(250);

    QSizePolicy clips_size_policy(QSizePolicy::Expanding, QSizePolicy::Preferred);
    clips_size_policy.setHorizontalStretch(1);

    QWidget* clips_widget = new QWidget;
    clips_widget->setLayout(vlayout);
    clips_widget->setSizePolicy(clips_size_policy);
    hlayout->addWidget(clips_widget);

    QSizePolicy view_size_policy(QSizePolicy::Expanding, QSizePolicy::Preferred);
    view_size_policy.setHorizontalStretch(1);

    foreach(MVClusterView * V, d->m_views)
    {
        V->setSizePolicy(view_size_policy);
        hlayout->addWidget(V);
    }

    foreach(MVClusterView * V, d->m_views)
    {
        d->connect_view(V);
    }

    connect(&d->m_computer, SIGNAL(computationFinished()), this, SLOT(slot_computation_finished()));
}

MVClusterWidget::~MVClusterWidget()
{
    d->m_computer.stopComputation(); // important do take care of this before things start getting destructed!
    delete d;
}

/*
void MVClusterWidget::setMscmdServerUrl(const QString& url)
{
    d->m_mscmdserver_url = url;
}
*/

void MVClusterWidget::setMLProxyUrl(const QString& url)
{
    d->m_mlproxy_url = url;
}

void MVClusterWidget::setData(const Mda& X)
{
    d->m_data = X;
    foreach(MVClusterView * V, d->m_views)
    {
        V->setData(Mda());
    }
    double max_abs_val = 0;
    int NN = X.totalSize();
    for (int i = 0; i < NN; i++) {
        if (fabs(X.get(i)) > max_abs_val)
            max_abs_val = fabs(X.get(i));
    }
    if (max_abs_val) {
        AffineTransformation T;
        T.setIdentity();
        double factor = 1 / 1.2;
        T.scale(1.0 / max_abs_val * factor, 1.0 / max_abs_val * factor, 1.0 / max_abs_val * factor);
        this->setTransformation(T);
    }
    d->set_data_on_visible_views_that_need_it();
}

void MVClusterWidget::setTimes(const QList<double>& times)
{
    foreach(MVClusterView * V, d->m_views)
    {
        V->setTimes(times);
    }
}

void MVClusterWidget::setLabels(const QList<int>& labels)
{
    foreach(MVClusterView * V, d->m_views)
    {
        V->setLabels(labels);
    }
}

void MVClusterWidget::setAmplitudes(const QList<double>& amps)
{
    foreach(MVClusterView * V, d->m_views)
    {
        V->setAmplitudes(amps);
    }
}

void MVClusterWidget::setOutlierScores(const QList<double>& outlier_scores)
{
    d->m_outlier_scores = outlier_scores;
}

void MVClusterWidget::setCurrentEvent(const MVEvent& evt)
{
    foreach(MVClusterView * V, d->m_views)
    {
        V->setCurrentEvent(evt);
    }
    d->update_clips_view();
}

void MVClusterWidget::setClipSize(int clip_size)
{
    d->m_clip_size = clip_size;
    d->start_computation();
}

void MVClusterWidget::setTimeseries(const DiskReadMda& X)
{
    d->m_timeseries = X;
    d->start_computation();
}

void MVClusterWidget::setFirings(const DiskReadMda& F)
{
    d->m_firings = F;

    d->start_computation();
}

void MVClusterWidget::setLabelsToUse(const QList<int>& labels)
{
    d->m_labels_to_use = labels;
    d->start_computation();
}

void MVClusterWidget::setTransformation(const AffineTransformation& T)
{
    foreach(MVClusterView * V, d->m_views)
    {
        V->setTransformation(T);
    }
}

MVEvent MVClusterWidget::currentEvent()
{
    return d->m_views[0]->currentEvent();
}

void MVClusterWidget::slot_view_current_event_changed()
{
    MVClusterView* V0 = (MVClusterView*)sender();
    this->setCurrentEvent(V0->currentEvent());
    emit currentEventChanged();
}

void MVClusterWidget::slot_view_transformation_changed()
{
    MVClusterView* V0 = (MVClusterView*)sender();
    AffineTransformation T = V0->transformation();
    foreach(MVClusterView * V, d->m_views)
    {
        V->setTransformation(T);
    }
}

void MVClusterWidget::slot_show_clip_view_toggled(bool val)
{
    d->m_clips_view->setVisible(val);
}

void MVClusterWidget::slot_show_view_toggled(bool val)
{
    int index = sender()->property("view_index").toInt();
    if ((index >= 0) && (index < d->m_views.count())) {
        d->m_views[index]->setVisible(val);
    }
    d->set_data_on_visible_views_that_need_it();
}

void MVClusterWidget::slot_computation_finished()
{
    d->m_computer.stopComputation(); //because I'm paranoid

    DiskReadMda F = d->m_computer.firings_subset;

    QList<double> times;
    QList<int> labels;
    QList<double> amplitudes;
    QList<double> outlier_scores;
    for (long j = 0; j < F.N2(); j++) {
        times << F.value(1, j);
        labels << (int)F.value(2, j);
        amplitudes << F.value(3, j);
        outlier_scores << F.value(4, j);
    }

    int K = compute_max(labels);
    QList<int> labels_map;
    labels_map << 0;
    int aa = 1;
    for (int k = 1; k <= K; k++) {
        if (d->m_labels_to_use.indexOf(k) >= 0) {
            labels_map << aa;
            aa++;
        } else {
            labels_map << 0;
        }
    }
    for (long j = 0; j < labels.count(); j++) {
        labels[j] = labels_map[labels[j]];
    }

    this->setTimes(times);
    this->setLabels(labels);
    this->setAmplitudes(amplitudes);
    this->setOutlierScores(outlier_scores);

    this->setData(d->m_computer.data);
}

void MVClusterWidgetPrivate::connect_view(MVClusterView* V)
{
    QObject::connect(V, SIGNAL(currentEventChanged()), q, SLOT(slot_view_current_event_changed()));
    QObject::connect(V, SIGNAL(transformationChanged()), q, SLOT(slot_view_transformation_changed()));
}

void MVClusterWidgetPrivate::update_clips_view()
{
    MVEvent evt = q->currentEvent();
    QString info_txt;
    if (evt.time >= 0) {
        QList<double> times;
        times << evt.time;
        Mda clip0 = extract_clips(m_timeseries, times, m_clip_size);
        double ppp = m_outlier_scores.value(current_event_index());
        if (ppp) {
            info_txt = QString("Outlier score: %1").arg(ppp);
        }
        m_clips_view->setClips(clip0);
    } else {
        m_clips_view->setClips(Mda());
    }
    m_info_bar->setText(info_txt);
}

int MVClusterWidgetPrivate::current_event_index()
{
    return m_views[0]->currentEventIndex();
}

void MVClusterWidgetPrivate::set_data_on_visible_views_that_need_it()
{
    foreach(MVClusterView * V, m_views)
    {
        if (V->isVisible()) {
            if (!V->hasData()) {
                V->setData(m_data);
            }
        }
    }
}

void MVClusterWidgetPrivate::start_computation()
{
    m_computer.stopComputation();
    //m_computer.mscmdserver_url = m_mscmdserver_url;
    m_computer.mlproxy_url = m_mlproxy_url;
    m_computer.timeseries = m_timeseries;
    m_computer.firings = m_firings;
    m_computer.clip_size = m_clip_size;
    m_computer.labels_to_use = m_labels_to_use;
    m_computer.startComputation();
}

void MVClusterWidgetComputer::compute()
{
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
    }

    QString features_path;
    {
        MountainProcessRunner MT;
        QString processor_name = "extract_clips_features";
        MT.setProcessorName(processor_name);

        QMap<QString, QVariant> params;
        params["timeseries"] = timeseries.path();
        params["firings"] = firings_out_path;
        params["clip_size"] = clip_size;
        params["num_features"] = 3;
        MT.setInputParameters(params);
        //MT.setMscmdServerUrl(mscmdserver_url);
        MT.setMLProxyUrl(mlproxy_url);

        features_path = MT.makeOutputFilePath("features");

        MT.runProcess(this);
    }
    firings_subset = DiskReadMda(firings_out_path);
    firings_subset.setHaltAgent(this);

    DiskReadMda features(features_path);
    features.setHaltAgent(this);
    features.readChunk(data, 0, 0, features.N1(), features.N2());
}
