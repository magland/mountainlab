/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
** Created: 6/8/2016
*******************************************************/

#include "mvspikesprayview.h"
#include "extract_clips.h"
#include "mountainprocessrunner.h"
#include "mvmainwindow.h"
#include "mvspikespraypanel.h"

#include <QHBoxLayout>
#include <QJsonDocument>
#include <QMessageBox>
#include <QPainter>
#include <QScrollArea>
#include <QSettings>
#include <taskprogress.h>
#include "mlcommon.h"
#include "actionfactory.h"
#include <QFileDialog>
#include <QLabel>
#include <QSpinBox>

/// TODO: (MEDIUM) spike spray should respond to mouse wheel and show current position with marker
/// TODO: (MEDIUM) much more responsive rendering of spike spray
/// TODO: (HIGH) Implement the panels as layers, make a layer stack, and handle zoom/pan without using a scroll widget
/// TODO: (HIGH) GUI option to set max. number of spikespray spikes to render

class MVSpikeSprayComputer {
public:
    //input
    DiskReadMda32 timeseries;
    DiskReadMda firings;
    QString mlproxy_url;
    QSet<int> labels_to_use;
    int clip_size;
    int max_per_label;

    //output
    Mda clips_to_render;
    QVector<int> labels_to_render;

    void compute();

    bool loaded_from_static_output = false;
    QJsonObject exportStaticOutput();
    void loadStaticOutput(const QJsonObject& X);
};

QJsonObject MVSpikeSprayComputer::exportStaticOutput()
{
    QJsonObject ret;
    ret["version"] = "MVSpikeSprayComputer-0.1";
    {
        QByteArray ba = clips_to_render.toByteArray32();
        ret["clips_to_render"] = MLUtil::toJsonValue(ba);
    }
    ret["labels_to_render"] = MLUtil::toJsonValue(labels_to_render);
    return ret;
}

void MVSpikeSprayComputer::loadStaticOutput(const QJsonObject& X)
{
    {
        QByteArray ba;
        MLUtil::fromJsonValue(ba, X["clips_to_render"]);
        clips_to_render.fromByteArray(ba);
    }
    MLUtil::fromJsonValue(labels_to_render, X["labels_to_render"]);
    loaded_from_static_output = true;
}

class MVSpikeSprayViewPrivate {
public:
    MVSpikeSprayView* q;
    MVContext* m_context;
    QSet<int> m_labels_to_use;
    int m_max_spikes_per_label = 500;

    double m_amplitude_factor = 0; //zero triggers auto-calculation
    double m_brightness_factor = 1;
    double m_weight_factor = 1;
    int m_panel_width = 0;

    Mda m_clips_to_render;
    Mda m_clips_to_render_subchannels;
    QVector<int> m_labels_to_render;
    MVSpikeSprayComputer m_computer;

    QHBoxLayout* m_panel_layout;
    QList<MVSpikeSprayPanel*> m_panels;

    MVSpikeSprayPanel* add_panel();
    void set_amplitude_factor(double val);
    void set_panel_width(double w);
    double actual_panel_width();
    void update_panel_factors();
};

MVSpikeSprayView::MVSpikeSprayView(MVAbstractContext* context)
    : MVAbstractView(context)
{
    d = new MVSpikeSprayViewPrivate;
    d->q = this;

    MVContext* c = qobject_cast<MVContext*>(mvContext());
    Q_ASSERT(c);

    d->m_context = c;

    recalculateOnOptionChanged("clip_size");
    //recalculateOnOptionChanged("timeseries_for_spikespray");
    recalculateOn(c, SIGNAL(clusterMergeChanged()));
    recalculateOn(c, SIGNAL(timeseriesNamesChanged()));
    recalculateOn(c, SIGNAL(visibleChannelsChanged()));
    this->recalculateOn(context, SIGNAL(firingsChanged()), false);
    onOptionChanged("cluster_color_index_shift", this, SLOT(onCalculationFinished()));

    QWidget* panel_widget = new QWidget;
    d->m_panel_layout = new QHBoxLayout;
    panel_widget->setLayout(d->m_panel_layout);

    QScrollArea* SA = new QScrollArea;
    SA->setWidget(panel_widget);
    SA->setWidgetResizable(true);

    QHBoxLayout* layout = new QHBoxLayout;
    layout->setSpacing(0);
    layout->setMargin(0);
    layout->addWidget(SA);
    this->setLayout(layout);

    ActionFactory::addToToolbar(ActionFactory::ActionType::ZoomIn, this, SLOT(slot_zoom_in()));
    ActionFactory::addToToolbar(ActionFactory::ActionType::ZoomOut, this, SLOT(slot_zoom_out()));
    ActionFactory::addToToolbar(ActionFactory::ActionType::ZoomInVertical, this, SLOT(slot_vertical_zoom_in()));
    ActionFactory::addToToolbar(ActionFactory::ActionType::ZoomOutVertical, this, SLOT(slot_vertical_zoom_out()));

    /*
    {
        QAction* A = new QAction("<-col", this);
        A->setProperty("action_type", "toolbar");
        QObject::connect(A, SIGNAL(triggered(bool)), this, SLOT(slot_shift_colors_left()));
        this->addAction(A);
    }
    {
        QAction* A = new QAction("col->", this);
        A->setProperty("action_type", "toolbar");
        QObject::connect(A, SIGNAL(triggered(bool)), this, SLOT(slot_shift_colors_right()));
        this->addAction(A);
    }
    */

    {
        QAction* A = new QAction("Export static view", this);
        A->setProperty("action_type", "");
        QObject::connect(A, SIGNAL(triggered(bool)), this, SLOT(slot_export_static_view()));
        this->addAction(A);
    }

    {
        QSlider* S = new QSlider;
        S->setOrientation(Qt::Horizontal);
        S->setRange(0, 300);
        S->setValue(100);
        S->setMaximumWidth(150);
        this->addToolbarControl(new QLabel("brightness:"));
        this->addToolbarControl(S);
        QObject::connect(S, SIGNAL(valueChanged(int)), this, SLOT(slot_brightness_slider_changed(int)));
    }

    {
        QSlider* S = new QSlider;
        S->setOrientation(Qt::Horizontal);
        S->setRange(0, 200);
        S->setValue(100);
        S->setMaximumWidth(150);
        this->addToolbarControl(new QLabel("weight:"));
        this->addToolbarControl(S);
        QObject::connect(S, SIGNAL(valueChanged(int)), this, SLOT(slot_weight_slider_changed(int)));
    }

    {
        QSpinBox* SB = new QSpinBox;
        SB->setRange(1, 10000);
        SB->setValue(d->m_max_spikes_per_label);
        QObject::connect(SB, SIGNAL(valueChanged(int)), this, SLOT(slot_set_max_spikes_per_label(int)));
        this->addToolbarControl(new QLabel("Max. #spikes:"));
        this->addToolbarControl(SB);
    }
}

MVSpikeSprayView::~MVSpikeSprayView()
{
    this->stopCalculation();
    delete d;
}

void MVSpikeSprayView::setLabelsToUse(const QSet<int>& labels_to_use)
{
    d->m_labels_to_use = labels_to_use;
    qDeleteAll(d->m_panels);
    d->m_panels.clear();

    QList<int> list = labels_to_use.toList();
    qSort(list);

    {
        MVSpikeSprayPanel* P = d->add_panel();
        P->setLabelsToUse(labels_to_use);
        P->setLegendVisible(false);
    }
    if (list.count() > 1) {
        for (int i = 0; i < list.count(); i++) {
            MVSpikeSprayPanel* P = d->add_panel();
            QSet<int> tmp;
            tmp.insert(list[i]);
            P->setLabelsToUse(tmp);
        }
    }
    d->update_panel_factors();
    recalculate();
}

void MVSpikeSprayView::prepareCalculation()
{
    //QString timeseries_name = d->m_context->option("timeseries_for_spikespray").toString();
    //if (timeseries_name.isEmpty())
    //    timeseries_name = d->m_context->currentTimeseriesName();
    QString timeseries_name = d->m_context->currentTimeseriesName();
    this->setCalculatingMessage(QString("Calculating using %1...").arg(timeseries_name));
    d->m_labels_to_render.clear();
    d->m_computer.mlproxy_url = d->m_context->mlProxyUrl();
    d->m_computer.timeseries = d->m_context->timeseries(timeseries_name);
    d->m_computer.firings = d->m_context->firings();
    d->m_computer.labels_to_use = d->m_labels_to_use;
    d->m_computer.max_per_label = d->m_max_spikes_per_label;
    d->m_computer.clip_size = d->m_context->option("clip_size").toInt();

    d->m_amplitude_factor = 0;
}

void MVSpikeSprayView::runCalculation()
{
    d->m_computer.compute();
}

Mda extract_channels_from_clips(const Mda& clips, QList<int> channels)
{
    Mda ret(channels.count(), clips.N2(), clips.N3());
    for (int i = 0; i < clips.N3(); i++) {
        for (int t = 0; t < clips.N2(); t++) {
            for (int j = 0; j < channels.count(); j++) {
                ret.setValue(clips.value(channels[j] - 1, t, i), j, t, i);
            }
        }
    }
    return ret;
}

void MVSpikeSprayView::onCalculationFinished()
{
    MVContext* c = qobject_cast<MVContext*>(mvContext());
    Q_ASSERT(c);

    d->m_clips_to_render = d->m_computer.clips_to_render;
    d->m_labels_to_render = d->m_computer.labels_to_render;
    if (c->viewMerged()) {
        d->m_labels_to_render = d->m_context->clusterMerge().mapLabels(d->m_labels_to_render);
    }

    if (!d->m_amplitude_factor) {
        double maxval = qMax(qAbs(d->m_clips_to_render.minimum()), qAbs(d->m_clips_to_render.maximum()));
        if (maxval) {
            d->set_amplitude_factor(1.5 / maxval);
        }
    }

    if (!c->visibleChannels().isEmpty())
        d->m_clips_to_render_subchannels = extract_channels_from_clips(d->m_clips_to_render, c->visibleChannels());
    else
        d->m_clips_to_render_subchannels = d->m_clips_to_render;
    for (int i = 0; i < d->m_panels.count(); i++) {
        d->m_panels[i]->setClipsToRender(&d->m_clips_to_render_subchannels);
        d->m_panels[i]->setLabelsToRender(d->m_labels_to_render);
    }
}

QJsonObject MVSpikeSprayView::exportStaticView()
{
    QJsonObject ret = MVAbstractView::exportStaticView();
    ret["view-type"] = "MVSpikeSprayView";
    ret["computer-output"] = d->m_computer.exportStaticOutput();
    ret["labels-to-use"] = MLUtil::toJsonValue(d->m_labels_to_use.toList().toVector());
    return ret;
}

void MVSpikeSprayView::loadStaticView(const QJsonObject& X)
{
    MVAbstractView::loadStaticView(X);
    QJsonObject computer_output = X["computer-output"].toObject();
    d->m_computer.loadStaticOutput(computer_output);
    QVector<int> labels_to_use;
    MLUtil::fromJsonValue(labels_to_use, X["labels-to-use"]);
    this->setLabelsToUse(labels_to_use.toList().toSet());
    this->recalculate();
}

void MVSpikeSprayView::paintEvent(QPaintEvent* evt)
{
    Q_UNUSED(evt)

    QPainter painter(this);
}

void MVSpikeSprayView::keyPressEvent(QKeyEvent* evt)
{
    if (evt->key() == Qt::Key_Up) {
        this->slot_vertical_zoom_in();
    }
    else if (evt->key() == Qt::Key_Down) {
        this->slot_vertical_zoom_out();
    }
    else {
        QWidget::keyPressEvent(evt);
    }
}

void MVSpikeSprayView::wheelEvent(QWheelEvent* evt)
{
    if (evt->delta() > 0) {
        slot_zoom_in();
    }
    else if (evt->delta() < 0) {
        slot_zoom_out();
    }
}

void MVSpikeSprayView::slot_zoom_in()
{
    d->set_panel_width(d->actual_panel_width() + 10);
}

void MVSpikeSprayView::slot_zoom_out()
{
    d->set_panel_width(qMax(d->actual_panel_width() - 10, 30.0));
}

void MVSpikeSprayView::slot_vertical_zoom_in()
{
    d->set_amplitude_factor(d->m_amplitude_factor * 1.2);
}

void MVSpikeSprayView::slot_vertical_zoom_out()
{
    d->set_amplitude_factor(d->m_amplitude_factor / 1.2);
}

void MVSpikeSprayView::slot_shift_colors_left(int step)
{
    MVContext* c = qobject_cast<MVContext*>(mvContext());
    Q_ASSERT(c);

    int shift = c->option("cluster_color_index_shift", 0).toInt();
    shift += step;
    c->setOption("cluster_color_index_shift", shift);
    this->onCalculationFinished();
}

void MVSpikeSprayView::slot_shift_colors_right()
{
    slot_shift_colors_left(-1);
}

void MVSpikeSprayView::slot_export_static_view()
{
    //QSettings settings("SCDA", "MountainView");
    //QString default_dir = settings.value("default_export_dir", "").toString();
    QString default_dir = QDir::currentPath();
    QString fname = QFileDialog::getSaveFileName(this, "Export static spikespray view", default_dir, "*.smv");
    if (fname.isEmpty())
        return;
    //settings.setValue("default_export_dir", QFileInfo(fname).path());
    if (QFileInfo(fname).suffix() != "smv")
        fname = fname + ".smv";
    QJsonObject obj = exportStaticView();
    if (!TextFile::write(fname, QJsonDocument(obj).toJson())) {
        qWarning() << "Unable to write file: " + fname;
    }
}

void MVSpikeSprayView::slot_brightness_slider_changed(int val)
{
    d->m_brightness_factor = val * 1.0 / 100;
    d->update_panel_factors();
}

void MVSpikeSprayView::slot_weight_slider_changed(int val)
{
    //d->m_weight_factor=val*1.0/100;
    d->m_weight_factor = exp((val * 1.0 - 100) / 20);
    if (val == 200) {
        d->m_weight_factor = 0;
    }
    d->update_panel_factors();
}

void MVSpikeSprayView::slot_set_max_spikes_per_label(int val)
{
    d->m_max_spikes_per_label = val;
    this->suggestRecalculate();
}

MVSpikeSprayPanel* MVSpikeSprayViewPrivate::add_panel()
{
    MVSpikeSprayPanel* P = new MVSpikeSprayPanel(m_context);
    m_panel_layout->addWidget(P);
    m_panels << P;
    P->setMinimumWidth(m_panel_width);
    return P;
}

void MVSpikeSprayViewPrivate::set_amplitude_factor(double val)
{
    m_amplitude_factor = val;
    update_panel_factors();
}

void MVSpikeSprayViewPrivate::set_panel_width(double w)
{
    m_panel_width = w;
    for (int i = 0; i < m_panels.count(); i++) {
        m_panels[i]->setMinimumWidth(w);
    }
}

double MVSpikeSprayViewPrivate::actual_panel_width()
{
    if (m_panels.isEmpty())
        return m_panel_width;
    else
        return m_panels[0]->width();
}

void MVSpikeSprayViewPrivate::update_panel_factors()
{
    for (int i = 0; i < m_panels.count(); i++) {
        MVSpikeSprayPanel* P = m_panels[i];
        P->setAmplitudeFactor(m_amplitude_factor);
        P->setBrightnessFactor(m_brightness_factor);
        P->setWeightFactor(m_weight_factor);
    }
}

void MVSpikeSprayComputer::compute()
{
    TaskProgress task("Spike spray computer");
    if (loaded_from_static_output) {
        task.log("Loaded from static output");
        return;
    }

    labels_to_render.clear();

    QString firings_out_path;
    {
        QList<int> list = labels_to_use.toList();
        qSort(list);
        QString labels_str;
        foreach (int x, list) {
            if (!labels_str.isEmpty())
                labels_str += ",";
            labels_str += QString("%1").arg(x);
        }

        MountainProcessRunner MT;
        QString processor_name = "mv_subfirings";
        MT.setProcessorName(processor_name);

        QMap<QString, QVariant> params;
        params["firings"] = firings.makePath();
        params["labels"] = labels_str;
        params["max_per_label"] = max_per_label;
        MT.setInputParameters(params);
        MT.setMLProxyUrl(mlproxy_url);

        firings_out_path = MT.makeOutputFilePath("firings_out");

        MT.runProcess();

        if (MLUtil::threadInterruptRequested()) {
            task.error(QString("Halted while running process: " + processor_name));
            return;
        }
    }

    task.setProgress(0.25);
    task.log("firings_out_path: " + firings_out_path);

    QString clips_path;
    {
        MountainProcessRunner MT;
        QString processor_name = "extract_clips";
        MT.setProcessorName(processor_name);

        QMap<QString, QVariant> params;
        params["timeseries"] = timeseries.makePath();
        params["firings"] = firings_out_path;
        params["clip_size"] = clip_size;
        MT.setInputParameters(params);
        MT.setMLProxyUrl(mlproxy_url);

        clips_path = MT.makeOutputFilePath("clips");

        MT.runProcess();
        if (MLUtil::threadInterruptRequested()) {
            task.error(QString("Halted while running process: " + processor_name));
            return;
        }
    }

    task.setProgress(0.5);
    task.log("clips_path: " + clips_path);

    DiskReadMda clips0(clips_path);
    if (!clips0.readChunk(clips_to_render, 0, 0, 0, clips0.N1(), clips0.N2(), clips0.N3())) {
        qWarning() << "Unable to read chunk of clips in spikespray view";
        return;
    }

    task.setProgress(0.75);

    DiskReadMda firings2(firings_out_path);
    task.log(QString("%1x%2 from %3x%4 (%5x%6x%7) (%8)").arg(firings2.N1()).arg(firings2.N2()).arg(firings.N1()).arg(firings.N2()).arg(clips0.N1()).arg(clips0.N2()).arg(clips0.N3()).arg(clips_to_render.N3()));
    Mda firings0;
    if (!firings2.readChunk(firings0, 0, 0, firings2.N1(), firings2.N2())) {
        qWarning() << "Unable to read chunk of firings in spikespray view";
        return;
    }
    task.setProgress(0.9);
    labels_to_render.clear();
    for (int i = 0; i < firings0.N2(); i++) {
        int label0 = (int)firings0.value(2, i);
        labels_to_render << label0;
    }
}

MVSpikeSprayFactory::MVSpikeSprayFactory(MVMainWindow* mw, QObject* parent)
    : MVAbstractViewFactory(mw, parent)
{
}

QString MVSpikeSprayFactory::id() const
{
    return QStringLiteral("open-spike-spray");
}

QString MVSpikeSprayFactory::name() const
{
    return tr("Spike Spray");
}

QString MVSpikeSprayFactory::title() const
{
    return tr("Spike Spray");
}

MVAbstractView* MVSpikeSprayFactory::createView(MVAbstractContext* context)
{
    MVContext* c = qobject_cast<MVContext*>(context);
    Q_ASSERT(c);

    QList<int> ks = c->selectedClusters();
    if (ks.isEmpty())
        ks = c->clusterVisibilityRule().subset.toList();
    qSort(ks);
    if (ks.isEmpty()) {
        QMessageBox::warning(0, "Unable to open spike spray", "You must select at least one cluster.");
        return Q_NULLPTR;
    }
    MVSpikeSprayView* X = new MVSpikeSprayView(context);
    X->setLabelsToUse(ks.toSet());
    return X;
}

/*
QList<QAction*> MVSpikeSprayFactory::actions(const QMimeData& md)
{
    QSet<int> clusters;
    QDataStream ds(md.data("application/x-msv-clusters"));
    ds >> clusters;
    QList<QAction*> actions;
    if (!clusters.isEmpty()) {
        QAction* action = new QAction("Spike spray", 0);
        MVMainWindow* mw = this->mainWindow();
        /// Witold, do I need to worry about mw object being deleted before this action is triggered?
        connect(action, &QAction::triggered, [mw]() {
            mw->openView("open-spike-spray");
        });
        actions << action;
    }
    return actions;
}
*/

void MVSpikeSprayFactory::updateEnabled(MVContext* context)
{
    Q_UNUSED(context)
    //setEnabled(!c->selectedClusters().isEmpty());
}
