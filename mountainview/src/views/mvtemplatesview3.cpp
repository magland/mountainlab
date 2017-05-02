/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
** Created: 7/29/2016
*******************************************************/

#include "mvtemplatesview3.h"
#include "mvtemplatesview2panel.h"

#include <QLabel>
#include <QSpinBox>
#include <QVBoxLayout>
#include <mountainprocessrunner.h>
#include <taskprogress.h>
#include "actionfactory.h"

struct ClusterData2 {
    int k = 0;
    int channel = 0;
    Mda template0;
    Mda stdev0;
    int num_events = 0;

    QJsonObject toJsonObject();
    void fromJsonObject(const QJsonObject& X);
};

class MVTemplatesView3Calculator {
public:
    //input
    QString mlproxy_url;
    DiskReadMda32 timeseries;
    DiskReadMda firings;
    int clip_size;

    //output
    QList<ClusterData2> cluster_data;

    virtual void compute();

    bool loaded_from_static_output = false;
    QJsonObject exportStaticOutput();
    void loadStaticOutput(const QJsonObject& X);
    static void mv_compute_templates_stdevs(DiskReadMda& templates_out, DiskReadMda& stdevs_out, const QString& mlproxy_url, const QString& timeseries, const QString& firings, int clip_size);
};

class MVTemplatesView3Private {
public:
    MVTemplatesView3* q;

    MVTemplatesView3Calculator m_calculator;
    QList<ClusterData2> m_cluster_data;
    double m_max_absval = 1;
    QList<MVTemplatesView2Panel*> m_panels;
    //int m_num_rows = 1;

    double m_total_time_sec = 0;
    bool m_zoomed_out_once = false;
    double m_vscale_factor = 4;
    double m_hscale_factor = 2;

    void compute_total_time();
    void update_panels();
    void update_scale_factors();
    static double get_disksize_for_firing_rate(double firing_rate);
};

MVTemplatesView3::MVTemplatesView3(MVAbstractContext* mvcontext)
    : MVGridView(mvcontext)
{
    d = new MVTemplatesView3Private;
    d->q = this;

    ActionFactory::addToToolbar(ActionFactory::ActionType::ZoomInVertical, this, SLOT(slot_vertical_zoom_in()));
    ActionFactory::addToToolbar(ActionFactory::ActionType::ZoomOutVertical, this, SLOT(slot_vertical_zoom_out()));

    MVContext* c = qobject_cast<MVContext*>(mvcontext);
    Q_ASSERT(c);

    this->recalculateOn(c, SIGNAL(firingsChanged()), false);
    this->recalculateOn(c, SIGNAL(currentTimeseriesChanged()));
    this->recalculateOnOptionChanged("clip_size");

    QObject::connect(c, SIGNAL(clusterMergeChanged()), this, SLOT(slot_update_highlighting()));
    QObject::connect(c, SIGNAL(currentClusterChanged()), this, SLOT(slot_update_highlighting()));
    QObject::connect(c, SIGNAL(selectedClustersChanged()), this, SLOT(slot_update_highlighting()));
    connect(c, SIGNAL(clusterVisibilityChanged()), this, SLOT(slot_update_panels()));
    connect(c, SIGNAL(viewMergedChanged()), this, SLOT(slot_update_panels()));

    connect(this, SIGNAL(signalViewClicked(int, Qt::KeyboardModifiers)), this, SLOT(slot_panel_clicked(int, Qt::KeyboardModifiers)));

    /*
    {
        QSpinBox* SB = new QSpinBox;
        SB->setRange(1, 10);
        SB->setValue(1);
        QObject::connect(SB, SIGNAL(valueChanged(int)), this, SLOT(slot_set_num_rows(int)));
        this->addToolbarControl(new QLabel("# Rows:"));
        this->addToolbarControl(SB);
    }
    */

    this->recalculate();
}

MVTemplatesView3::~MVTemplatesView3()
{
    this->stopCalculation();
    delete d;
}

void MVTemplatesView3::prepareCalculation()
{
    if (!d->m_calculator.loaded_from_static_output) {
        d->compute_total_time();
    }

    MVContext* c = qobject_cast<MVContext*>(mvContext());
    Q_ASSERT(c);

    d->m_calculator.mlproxy_url = c->mlProxyUrl();
    d->m_calculator.timeseries = c->currentTimeseries();
    d->m_calculator.firings = c->firings();
    d->m_calculator.clip_size = c->option("clip_size", 100).toInt();
    d->m_calculator.cluster_data.clear();
    update();
}

void MVTemplatesView3::runCalculation()
{
    d->m_calculator.compute();
}

void MVTemplatesView3::onCalculationFinished()
{
    d->m_cluster_data = d->m_calculator.cluster_data;
    if (!d->m_zoomed_out_once) {
        this->zoomAllTheWayOut();
        d->m_zoomed_out_once = true;
    }
    d->update_panels();
}

void MVTemplatesView3::zoomAllTheWayOut()
{
    //d->m_panel_widget->zoomAllTheWayOut();
}

void MVTemplatesView3::keyPressEvent(QKeyEvent* evt)
{
    if (evt->key() == Qt::Key_Up) {
        slot_vertical_zoom_in();
    }
    if (evt->key() == Qt::Key_Down) {
        slot_vertical_zoom_out();
    }
}

void MVTemplatesView3::prepareMimeData(QMimeData& mimeData, const QPoint& pos)
{
    MVContext* c = qobject_cast<MVContext*>(mvContext());
    Q_ASSERT(c);

    /*
    int view_index = d->find_view_index_at(pos);
    if (view_index >= 0) {
        int k = d->m_views[view_index]->k();
        if (!c->selectedClusters().contains(k)) {
            c->clickCluster(k, Qt::NoModifier);
        }
    }
    */

    QByteArray ba;
    QDataStream ds(&ba, QIODevice::WriteOnly);
    ds << c->selectedClusters();
    mimeData.setData("application/x-msv-clusters", ba); // selected cluster data

    MVAbstractView::prepareMimeData(mimeData, pos); // call base class implementation
}

void MVTemplatesView3::slot_update_panels()
{
    d->update_panels();
}

void MVTemplatesView3::slot_update_highlighting()
{
    MVContext* c = qobject_cast<MVContext*>(mvContext());
    Q_ASSERT(c);

    int current_k = c->currentCluster();
    QSet<int> selected_ks = c->selectedClusters().toSet();
    for (int i = 0; i < d->m_panels.count(); i++) {
        int indx = d->m_panels[i]->property("cluster_data_index").toInt();
        int k0 = d->m_cluster_data.value(indx).k;
        d->m_panels[i]->setCurrent(current_k == k0);
        d->m_panels[i]->setSelected(selected_ks.contains(k0));
        if (current_k == k0) {
            //d->m_panel_widget->setCurrentPanelIndex(i); //for purposes of zooming
        }
    }
    updateViews();
}

void MVTemplatesView3::slot_panel_clicked(int index, Qt::KeyboardModifiers modifiers)
{
    MVContext* c = qobject_cast<MVContext*>(mvContext());
    Q_ASSERT(c);

    if (modifiers & Qt::ShiftModifier) {
        int i1 = MVGridView::currentViewIndex();
        int i2 = index;
        int j1 = qMin(i1, i2);
        int j2 = qMax(i1, i2);
        if ((j1 >= 0) && (j2 >= 0)) {
            QSet<int> set = c->selectedClusters().toSet();
            for (int j = j1; j <= j2; j++) {
                MVTemplatesView2Panel* panel = d->m_panels.value(j);
                if (panel) {
                    int indx = panel->property("cluster_data_index").toInt();
                    int k = d->m_cluster_data.value(indx).k;
                    if (k > 0)
                        set.insert(k);
                }
            }
            c->setSelectedClusters(set.toList());
        }
    }
    else {
        MVTemplatesView2Panel* panel = d->m_panels.value(index);
        int indx = panel->property("cluster_data_index").toInt();
        ClusterData2 CD = d->m_cluster_data.value(indx);
        if (CD.k > 0) {
            c->clickCluster(CD.k, modifiers);
        }
    }
}

void MVTemplatesView3::slot_vertical_zoom_in()
{
    d->m_vscale_factor *= 1.1;
    d->update_scale_factors();
}

void MVTemplatesView3::slot_vertical_zoom_out()
{
    d->m_vscale_factor /= 1.1;
    d->update_scale_factors();
}

/*
void MVTemplatesView3::slot_set_num_rows(int num_rows)
{
    d->m_num_rows = num_rows;
    d->update_panels();
}
*/

void MVTemplatesView3Calculator::mv_compute_templates_stdevs(DiskReadMda& templates_out, DiskReadMda& stdevs_out, const QString& mlproxy_url, const QString& timeseries, const QString& firings, int clip_size)
{
    TaskProgress task(TaskProgress::Calculate, "mv_compute_templates_stdevs");
    task.log("mlproxy_url: " + mlproxy_url);
    MountainProcessRunner X;
    QString processor_name = "mv_compute_templates";
    X.setProcessorName(processor_name);

    QMap<QString, QVariant> params;
    params["timeseries"] = timeseries;
    params["firings"] = firings;
    params["clip_size"] = clip_size;
    X.setInputParameters(params);
    X.setMLProxyUrl(mlproxy_url);

    QString templates_fname = X.makeOutputFilePath("templates");
    QString stdevs_fname = X.makeOutputFilePath("stdevs");

    task.log("X.compute()");
    X.runProcess();
    task.log("Returning DiskReadMda: " + templates_fname + " " + stdevs_fname);
    templates_out.setPath(templates_fname);
    stdevs_out.setPath(stdevs_fname);

    //templates_out.setRemoteDataType("float32_q8");
    //stdevs_out.setRemoteDataType("float32_q8");
}

void MVTemplatesView3Calculator::compute()
{
    TaskProgress task(TaskProgress::Calculate, "Cluster Details");
    if (this->loaded_from_static_output) {
        task.log("Loaded from static output");
        return;
    }

    QTime timer;
    timer.start();
    task.setProgress(0.1);

    int M = timeseries.N1();
    //int N = timeseries.N2();
    int L = firings.N2();
    int T = clip_size;

    QVector<int> labels;

    task.log("Setting up labels");
    task.setProgress(0.2);
    for (int i = 0; i < L; i++) {
        labels << (int)firings.value(2, i);
    }

    if (MLUtil::threadInterruptRequested()) {
        task.error("Halted *");
        return;
    }
    task.log("Clearing data");
    cluster_data.clear();

    task.setLabel("Computing templates");
    task.setProgress(0.4);
    int K = MLCompute::max(labels);

    QString timeseries_path = timeseries.makePath();
    QString firings_path = firings.makePath();

    task.log("mp_compute_templates_stdevs: " + mlproxy_url + " timeseries_path=" + timeseries_path + " firings_path=" + firings_path);
    task.setProgress(0.6);
    DiskReadMda templates0, stdevs0;
    mv_compute_templates_stdevs(templates0, stdevs0, mlproxy_url, timeseries_path, firings_path, T);
    if (MLUtil::threadInterruptRequested()) {
        task.error("Halted **");
        return;
    }

    task.setLabel("Setting cluster data");
    task.setProgress(0.75);
    for (int k = 1; k <= K; k++) {
        if (MLUtil::threadInterruptRequested()) {
            task.error("Halted ***");
            return;
        }
        ClusterData2 CD;
        CD.k = k;
        CD.channel = 0;
        for (int i = 0; i < L; i++) {
            if (labels[i] == k) {
                CD.num_events++;
            }
        }
        if (MLUtil::threadInterruptRequested()) {
            task.error("Halted ****");
            return;
        }
        if (!templates0.readChunk(CD.template0, 0, 0, k - 1, M, T, 1)) {
            qWarning() << "Unable to read chunk of templates in templates view 3";
            return;
        }
        if (!stdevs0.readChunk(CD.stdev0, 0, 0, k - 1, M, T, 1)) {
            qWarning() << "Unable to read chunk of stdevs in templates view 3";
        }
        if (!MLUtil::threadInterruptRequested()) {
            if (CD.num_events > 0) {
                cluster_data << CD;
            }
        }
    }
}

void MVTemplatesView3Private::compute_total_time()
{
    MVContext* c = qobject_cast<MVContext*>(q->mvContext());
    Q_ASSERT(c);

    m_total_time_sec = c->currentTimeseries().N2() / c->sampleRate();
}

double MVTemplatesView3Private::get_disksize_for_firing_rate(double firing_rate)
{
    return qMin(1.0, sqrt(firing_rate / 5));
}

void MVTemplatesView3Private::update_panels()
{
    MVContext* c = qobject_cast<MVContext*>(q->mvContext());
    Q_ASSERT(c);

    m_panels.clear();
    q->clearViews();
    QList<QColor> channel_colors;
    {
        int M = m_cluster_data.value(0).template0.N1();
        for (int m = 0; m < M; m++) {
            channel_colors << c->channelColor(m);
        }
    }
    m_max_absval = 0;
    for (int i = 0; i < m_cluster_data.count(); i++) {
        m_max_absval = qMax(qMax(m_max_absval, qAbs(m_cluster_data[i].template0.minimum())), qAbs(m_cluster_data[i].template0.maximum()));
    }

    for (int i = 0; i < m_cluster_data.count(); i++) {
        ClusterData2 CD = m_cluster_data[i];
        if (c->clusterIsVisible(CD.k)) {
            MVTemplatesView2Panel* panel = new MVTemplatesView2Panel;
            panel->setElectrodeGeometry(c->electrodeGeometry());
            panel->setTemplate(CD.template0);
            panel->setChannelColors(channel_colors);
            panel->setColors(c->colors());
            panel->setTitle(QString::number(CD.k));
            panel->setProperty("cluster_data_index", i);
            if (c->sampleRate()) {
                double total_time_sec = c->currentTimeseries().N2() / c->sampleRate();
                if (total_time_sec) {
                    double firing_rate = CD.num_events / total_time_sec;
                    panel->setFiringRateDiskDiameter(get_disksize_for_firing_rate(firing_rate));
                }
            }
            m_panels << panel;
        }
    }
    //int num_rows = m_num_rows;
    //int num_cols = ceil(m_panels.count() * 1.0 / num_rows);
    for (int i = 0; i < m_panels.count(); i++) {
        q->addView(new PaintLayerWidget(m_panels[i]));
    }
    update_scale_factors();
    q->slot_update_highlighting();

    q->updateLayout();
}

void MVTemplatesView3Private::update_scale_factors()
{
    double vfactor = m_vscale_factor;
    if (m_max_absval)
        vfactor /= m_max_absval;
    for (int i = 0; i < m_panels.count(); i++) {
        m_panels[i]->setVerticalScaleFactor(vfactor);
    }
    q->updateViews();
}

MVTemplatesView3Factory::MVTemplatesView3Factory(MVMainWindow* mw, QObject* parent)
    : MVAbstractViewFactory(mw, parent)
{
}

QString MVTemplatesView3Factory::id() const
{
    return QStringLiteral("open-templates");
}

QString MVTemplatesView3Factory::name() const
{
    return tr("Templates");
}

QString MVTemplatesView3Factory::title() const
{
    return tr("Templates");
}

MVAbstractView* MVTemplatesView3Factory::createView(MVAbstractContext* context)
{
    MVTemplatesView3* X = new MVTemplatesView3(context);
    return X;
}
