//#include "mvclustercontextmenu.h"
#include "clusterdetailview.h"
#include "mvmainwindow.h"
#include "mvcontext.h"
#include "tabber.h"
#include "taskprogress.h"
#include "actionfactory.h"
#include "clusterdetailviewpropertiesdialog.h"

#include <QPainter>
#include "mvutils.h"
#include "mlcommon.h"
#include <QProgressDialog>
#include <QTime>
#include <QMap>
#include <QDebug>
#include <QMouseEvent>
#include <QSet>
#include <QSettings>
#include <QImageWriter>
#include "extract_clips.h"
#include <QApplication>
#include <QFileDialog>
#include <QJsonArray>
#include <QJsonDocument>
#include <QMenu>
#include <QMessageBox>
#include <QTimer>
#include "compute_templates_0.h"
#include "mountainprocessrunner.h"
#include <math.h>
#include "mlcommon.h"
#include "mvmisc.h"
#include "get_sort_indices.h"

struct ClusterData {
    ClusterData()
    {
        k = 0;
        channel = 0;
    }

    int k;
    int channel;
    Mda32 template0;
    Mda32 stdev0;
    int num_events = 0;

    QJsonObject toJsonObject();
    void fromJsonObject(const QJsonObject& X);
};

struct ChannelSpacingInfo {
    QVector<double> channel_locations;
    double channel_location_spacing;
    double vert_scaling_factor;
};

class ClusterDetailViewCalculator {
public:
    //input
    //QString mscmdserver_url;
    QString mlproxy_url;
    DiskReadMda32 timeseries;
    DiskReadMda firings;
    int clip_size;
    QSet<int> clusters_to_force_show;

    //output
    QList<ClusterData> cluster_data;

    virtual void compute();

    bool loaded_from_static_output = false;
    QJsonObject exportStaticOutput();
    void loadStaticOutput(const QJsonObject& X);
};

class ClusterDetailView;
class ClusterDetailViewPrivate;
class ClusterView {
public:
    friend class ClusterDetailViewPrivate;
    friend class ClusterDetailView;
    ClusterView(ClusterDetailView* q0, ClusterDetailViewPrivate* d0);
    void setClusterData(const ClusterData& CD);
    void setAttributes(QJsonObject aa);

    int k();
    void setChannelSpacingInfo(const ChannelSpacingInfo& csi);
    void setHighlighted(bool val);
    void setHovered(bool val);
    void setSelected(bool val);
    void setStdevShading(bool val);
    bool stdevShading();

    void paint(QPainter* painter, QRectF rect, bool render_image_mode = false);
    double spaceNeeded();
    ClusterData* clusterData();
    QRectF rect();
    QPointF template_coord2pix(int m, double t, double val);

    ClusterDetailView* q;
    ClusterDetailViewPrivate* d;
    double x_position_before_scaling;

private:
    ClusterData m_CD;
    ChannelSpacingInfo m_csi;
    int m_T;
    QRectF m_rect;
    QRectF m_top_rect;
    QRectF m_template_rect;
    QRectF m_bottom_rect;
    bool m_highlighted;
    bool m_hovered;
    bool m_selected;
    QJsonObject m_attributes;
    bool m_stdev_shading;

    QColor get_firing_rate_text_color(double rate);
};

class ClusterDetailViewPrivate {
public:
    ClusterDetailView* q;

    QList<ClusterData> m_cluster_data;

    double m_vscale_factor;
    double m_space_ratio;
    double m_scroll_x;

    double m_total_time_sec;
    int m_hovered_k;
    double m_anchor_x;
    double m_anchor_scroll_x;
    int m_anchor_view_index;
    ClusterDetailViewCalculator m_calculator;
    bool m_stdev_shading;
    bool m_zoomed_out_once;
    ClusterDetailViewProperties m_properties;

    QComboBox* m_sort_order_chooser;

    QList<ClusterView*> m_views;

    void compute_total_time();
    void set_hovered_k(int k);
    int find_view_index_at(QPoint pos);
    ClusterView* find_view_for_k(int k);
    int find_view_index_for_k(int k);
    void ensure_view_visible(ClusterView* V);
    void zoom(double factor);
    QString group_label_for_k(int k);
    int get_current_view_index();
    void do_paint(QPainter& painter, int W, int H, bool render_image_mode = false);
    void export_image();
    void toggle_stdev_shading();
    void shift_select_clusters_between(int k1, int k2);
    void sort_cluster_data(QList<ClusterData>& cluster_data);
    double compute_sort_score(const ClusterData& CD);

    static QList<ClusterData> merge_cluster_data(const ClusterMerge& CM, const QList<ClusterData>& CD);
    Mda32 get_template_waveforms_for_export();
    Mda32 get_template_waveform_stdevs_for_export();
};

ClusterDetailView::ClusterDetailView(MVAbstractContext* context)
    : MVAbstractView(context)
{
    d = new ClusterDetailViewPrivate;
    d->q = this;
    d->m_vscale_factor = 2;
    d->m_total_time_sec = 1;
    d->m_hovered_k = -1;
    d->m_space_ratio = 50;
    d->m_scroll_x = 0;
    d->m_anchor_x = -1;
    d->m_anchor_scroll_x = -1;
    d->m_anchor_view_index = -1;
    d->m_stdev_shading = false;
    d->m_zoomed_out_once = false;

    MVContext* c = qobject_cast<MVContext*>(context);
    Q_ASSERT(c);

    QObject::connect(c, SIGNAL(clusterMergeChanged()), this, SLOT(update()));
    QObject::connect(c, SIGNAL(currentClusterChanged()), this, SLOT(update()));
    QObject::connect(c, SIGNAL(selectedClustersChanged()), this, SLOT(update()));
    connect(c, SIGNAL(clusterVisibilityChanged()), this, SLOT(update()));
    connect(c, SIGNAL(viewMergedChanged()), this, SLOT(update()));

    this->recalculateOn(c, SIGNAL(firingsChanged()), false);
    recalculateOn(c, SIGNAL(currentTimeseriesChanged()));
    recalculateOn(c, SIGNAL(visibleChannelsChanged()));
    recalculateOnOptionChanged("clip_size");

    this->setMouseTracking(true);

    /*
    {
        QAction* A = new QAction("Export waveforms (.mda)", this);
        A->setProperty("action_type", "");
        QObject::connect(A, SIGNAL(triggered(bool)), this, SLOT(slot_export_waveforms()));
        this->addAction(A);
    }
    */
    {
        QAction* a = new QAction("Export image", this);
        this->addAction(a);
        connect(a, SIGNAL(triggered(bool)), this, SLOT(slot_export_image()));
    }
    {
        QAction* a = new QAction("Toggle std. dev. shading", this);
        this->addAction(a);
        connect(a, SIGNAL(triggered(bool)), this, SLOT(slot_toggle_stdev_shading()));
    }
    ActionFactory::addToToolbar(ActionFactory::ActionType::ZoomIn, this, SLOT(slot_zoom_in()));
    ActionFactory::addToToolbar(ActionFactory::ActionType::ZoomOut, this, SLOT(slot_zoom_out()));
    ActionFactory::addToToolbar(ActionFactory::ActionType::ZoomInVertical, this, SLOT(slot_vertical_zoom_in()));
    ActionFactory::addToToolbar(ActionFactory::ActionType::ZoomOutVertical, this, SLOT(slot_vertical_zoom_out()));

    /*
    {
        QAction* A = new QAction("Export static view", this);
        A->setProperty("action_type", "");
        QObject::connect(A, SIGNAL(triggered(bool)), this, SLOT(slot_export_static_view()));
        this->addAction(A);
    }
    */

    {
        QAction* A = new QAction("Export template waveforms...", this);
        A->setProperty("action_type", "");
        QObject::connect(A, SIGNAL(triggered(bool)), this, SLOT(slot_export_template_waveforms()));
        this->addAction(A);
    }

    {
        QAction* A = new QAction("Export template waveform stdevs...", this);
        A->setProperty("action_type", "");
        QObject::connect(A, SIGNAL(triggered(bool)), this, SLOT(slot_export_template_waveform_stdevs()));
        this->addAction(A);
    }

    {
        QAction* A = new QAction(QString("View properties"), this);
        A->setProperty("action_type", "toolbar");
        QObject::connect(A, SIGNAL(triggered(bool)), this, SLOT(slot_view_properties()));
        this->addAction(A);
    }

    /*
    {
        QComboBox* CB = new QComboBox;
        CB->addItem("Default order", "default");
        CB->addItem("Quality order", "quality");
        d->m_sort_order_chooser = CB;
        this->addToolbarControl(CB);
        QObject::connect(CB, SIGNAL(currentIndexChanged(int)), this, SLOT(slot_update_sort_order()));
    }
    */

    recalculate();
}

ClusterDetailView::~ClusterDetailView()
{
    this->stopCalculation();
    qDeleteAll(d->m_views);
    delete d;
}

void ClusterDetailView::prepareCalculation()
{
    MVContext* c = qobject_cast<MVContext*>(mvContext());
    Q_ASSERT(c);

    if (!d->m_calculator.loaded_from_static_output) {
        d->compute_total_time();
    }
    d->m_calculator.mlproxy_url = c->mlProxyUrl();
    d->m_calculator.timeseries = c->currentTimeseries();
    d->m_calculator.firings = c->firings();
    d->m_calculator.clip_size = c->option("clip_size", 100).toInt();
    d->m_calculator.clusters_to_force_show = c->clustersToForceShow().toSet();
    update();
}

void ClusterDetailView::runCalculation()
{
    d->m_calculator.compute();
}

Mda32 extract_channels_from_template(const Mda32& X, QList<int> channels)
{
    Mda32 ret(channels.count(), X.N2(), X.N3());
    for (int i = 0; i < X.N3(); i++) { //handle case of more than 1 template, just because
        for (int t = 0; t < X.N2(); t++) {
            for (int j = 0; j < channels.count(); j++) {
                ret.setValue(X.value(channels[j] - 1, t, i), j, t, i);
            }
        }
    }
    return ret;
}

void extract_channels_from_templates(QList<ClusterData>& CD, QList<int> channels)
{
    for (int i = 0; i < CD.count(); i++) {
        CD[i].template0 = extract_channels_from_template(CD[i].template0, channels);
    }
}

void ClusterDetailView::onCalculationFinished()
{
    MVContext* c = qobject_cast<MVContext*>(mvContext());
    Q_ASSERT(c);

    d->m_cluster_data = d->m_calculator.cluster_data;
    if (!c->visibleChannels().isEmpty()) {
        extract_channels_from_templates(d->m_cluster_data, c->visibleChannels());
    }
    if (!d->m_zoomed_out_once) {
        this->zoomAllTheWayOut();
        d->m_zoomed_out_once = true;
    }
    this->update();
}

void ClusterDetailView::zoomAllTheWayOut()
{
    d->m_space_ratio = 0;
    update();
}

QImage ClusterDetailView::renderImage(int W, int H)
{
    MVContext* c = qobject_cast<MVContext*>(mvContext());
    Q_ASSERT(c);

    if (!W)
        W = 1600;
    if (!H)
        H = 800;
    QImage ret = QImage(W, H, QImage::Format_RGB32);
    QPainter painter(&ret);

    int current_k = c->currentCluster();
    QList<int> selected_ks = c->selectedClusters();
    c->setCurrentCluster(-1);
    c->setSelectedClusters(QList<int>());
    d->do_paint(painter, W, H, true);
    c->setCurrentCluster(current_k);
    c->setSelectedClusters(selected_ks);
    this->update(); //make sure we update, because some internal stuff has changed!

    return ret;
}

QJsonObject ClusterDetailView::exportStaticView()
{
    QJsonObject ret = MVAbstractView::exportStaticView();
    ret["view-type"] = "ClusterDetailView";
    ret["version"] = "0.1";
    ret["calculator-output"] = d->m_calculator.exportStaticOutput();

    QJsonObject info;
    info["vscale_factor"] = d->m_vscale_factor;
    info["space_ratio"] = d->m_space_ratio;
    info["total_time_sec"] = d->m_total_time_sec;
    info["stdev_shading"] = d->m_stdev_shading;
    ret["info"] = info;

    return ret;
}

void ClusterDetailView::loadStaticView(const QJsonObject& X)
{
    MVAbstractView::loadStaticView(X);
    d->m_calculator.loadStaticOutput(X["calculator-output"].toObject());
    QJsonObject info = X["info"].toObject();
    d->m_vscale_factor = info["vscale_factor"].toDouble();
    d->m_space_ratio = info["space_ratio"].toDouble();
    d->m_total_time_sec = info["total_time_sec"].toDouble();
    d->m_stdev_shading = info["stdev_shading"].toBool();
    this->recalculate();
}

void ClusterDetailView::leaveEvent(QEvent*)
{
    d->set_hovered_k(-1);
}

ChannelSpacingInfo compute_channel_spacing_info(QList<ClusterData>& cdata, double vscale_factor)
{
    ChannelSpacingInfo info;
    info.vert_scaling_factor = 1;
    if (cdata.count() == 0)
        return info;
    int M = cdata[0].template0.N1();
    int T = cdata[0].template0.N2();
    double minval = 0, maxval = 0;
    for (int i = 0; i < cdata.count(); i++) {
        for (int t = 0; t < T; t++) {
            for (int m = 0; m < M; m++) {
                double val = cdata[i].template0.value(m, t);
                if (val < minval)
                    minval = val;
                if (val > maxval)
                    maxval = val;
            }
        }
    }
    info.channel_location_spacing = 1.0 / M;
    double y0 = 0.5 * info.channel_location_spacing;
    for (int m = 0; m < M; m++) {
        info.channel_locations << y0;
        y0 += info.channel_location_spacing;
    }
    double maxabsval = qMax(maxval, -minval);
    //info.vert_scaling_factor = 0.5 / M / maxabsval * vscale_factor;
    info.vert_scaling_factor = 1.0 / maxabsval * vscale_factor;
    return info;
}

void ClusterDetailView::paintEvent(QPaintEvent* evt)
{
    Q_UNUSED(evt)

    QPainter painter(this);

    d->do_paint(painter, width(), height());
}

void ClusterDetailView::keyPressEvent(QKeyEvent* evt)
{
    MVContext* c = qobject_cast<MVContext*>(mvContext());
    Q_ASSERT(c);

    double factor = 1.15;
    if (evt->key() == Qt::Key_Up) {
        d->m_vscale_factor *= factor;
        update();
    }
    else if (evt->key() == Qt::Key_Down) {
        d->m_vscale_factor /= factor;
        update();
    }
    else if ((evt->key() == Qt::Key_Plus) || (evt->key() == Qt::Key_Equal)) {
        d->zoom(1.1);
    }
    else if (evt->key() == Qt::Key_Minus) {
        d->zoom(1 / 1.1);
    }
    else if ((evt->key() == Qt::Key_A) && (evt->modifiers() & Qt::ControlModifier)) {
        QList<int> ks;
        for (int i = 0; i < d->m_views.count(); i++) {
            ks << d->m_views[i]->k();
        }
        c->setSelectedClusters(ks);
    }
    else if (evt->key() == Qt::Key_Left) {
        int view_index = d->get_current_view_index();
        if (view_index > 0) {
            int k = d->m_views[view_index - 1]->k();
            QList<int> ks;
            if (evt->modifiers() & Qt::ShiftModifier) {
                ks = c->selectedClusters();
                ks << k;
            }
            c->setSelectedClusters(ks);
            c->setCurrentCluster(k);
        }
    }
    else if (evt->key() == Qt::Key_Right) {
        int view_index = d->get_current_view_index();
        if ((view_index >= 0) && (view_index + 1 < d->m_views.count())) {
            int k = d->m_views[view_index + 1]->k();
            QList<int> ks;
            if (evt->modifiers() & Qt::ShiftModifier) {
                ks = c->selectedClusters();
                ks << k;
            }
            c->setSelectedClusters(ks);
            c->setCurrentCluster(k);
        }
    }
    else if (evt->matches(QKeySequence::SelectAll)) {
        QList<int> all_ks;
        for (int i = 0; i < d->m_views.count(); i++) {
            all_ks << d->m_views[i]->k();
        }
        c->setSelectedClusters(all_ks);
    }
    else
        evt->ignore();
}

void ClusterDetailView::mousePressEvent(QMouseEvent* evt)
{
    if (evt->button() == Qt::LeftButton) {
        QPoint pt = evt->pos();
        d->m_anchor_x = pt.x();
        d->m_anchor_scroll_x = d->m_scroll_x;
    }
}

void ClusterDetailView::mouseReleaseEvent(QMouseEvent* evt)
{
    MVContext* c = qobject_cast<MVContext*>(mvContext());
    Q_ASSERT(c);

    QPoint pt = evt->pos();

    int view_index = d->find_view_index_at(pt);

    if (evt->button() == Qt::LeftButton) {
        if ((d->m_anchor_x >= 0) && (qAbs(pt.x() - d->m_anchor_x) > 5)) {
            d->m_scroll_x = d->m_anchor_scroll_x - (pt.x() - d->m_anchor_x);
            d->m_anchor_x = -1;
            update();
            return;
        }

        d->m_anchor_x = -1;

        if (view_index >= 0) {
            int k = d->m_views[view_index]->k();
            if (evt->modifiers() & Qt::ShiftModifier) {
                int k0 = c->currentCluster();
                d->shift_select_clusters_between(k0, k);
            }
            else {
                c->clickCluster(k, evt->modifiers());
            }
        }
    }
}

void ClusterDetailView::mouseMoveEvent(QMouseEvent* evt)
{
    QPoint pt = evt->pos();
    if ((d->m_anchor_x >= 0) && (qAbs(pt.x() - d->m_anchor_x) > 5)) {
        d->m_scroll_x = d->m_anchor_scroll_x - (pt.x() - d->m_anchor_x);
        update();
        return;
    }

    int view_index = d->find_view_index_at(pt);
    if (view_index >= 0) {
        d->set_hovered_k(d->m_views[view_index]->k());
    }
    else {
        d->set_hovered_k(-1);
    }
}

void ClusterDetailView::mouseDoubleClickEvent(QMouseEvent* evt)
{
    // send yourself a context menu request
    QContextMenuEvent e(QContextMenuEvent::Mouse, evt->pos());
    QCoreApplication::sendEvent(this, &e);
}

void ClusterDetailView::wheelEvent(QWheelEvent* evt)
{
    int delta = evt->delta();
    double factor = 1;
    if (delta > 0)
        factor = 1.1;
    else
        factor = 1 / 1.1;
    d->zoom(factor);
}

void ClusterDetailView::prepareMimeData(QMimeData& mimeData, const QPoint& pos)
{
    MVContext* c = qobject_cast<MVContext*>(mvContext());
    Q_ASSERT(c);

    int view_index = d->find_view_index_at(pos);
    if (view_index >= 0) {
        int k = d->m_views[view_index]->k();
        if (!c->selectedClusters().contains(k)) {
            c->clickCluster(k, Qt::NoModifier);
        }
    }

    QByteArray ba;
    QDataStream ds(&ba, QIODevice::WriteOnly);
    ds << c->selectedClusters();
    mimeData.setData("application/x-msv-clusters", ba); // selected cluster data

    MVAbstractView::prepareMimeData(mimeData, pos); // call base class implementation
}

/*
void ClusterDetailView::slot_export_waveforms()
{
    int K = d->m_views.count();
    if (!K) {
        QMessageBox::warning(0, "Unable to export waveforms", "Unable to export waveforms. No waveforms are visible");
        return;
    }
    int M, T;
    {
        Mda32 template0 = d->m_views[0]->clusterData()->template0;
        M = template0.N1();
        T = template0.N2();
    }
    Mda32 templates(M, T, K);
    for (int k = 0; k < K; k++) {
        Mda32 template0 = d->m_views[k]->clusterData()->template0;
        templates.setChunk(template0, 0, 0, k);
    }
    QString default_dir = QDir::currentPath();
    QString fname = QFileDialog::getSaveFileName(this, "Export waveforms", default_dir, "*.mda");
    if (fname.isEmpty())
        return;
    if (QFileInfo(fname).suffix() != "mda")
        fname = fname + ".mda";
    TaskProgress task("Export waveforms");
    if (templates.write32(fname)) {
        task.log(QString("Exported waveforms %1x%2x%3 to %4").arg(M).arg(T).arg(K).arg(fname));
    }
    else {
        task.error(QString("Unable to export waveforms. Unable to write file: ") + fname);
        QMessageBox::warning(0, "Unable to export waveforms", "Unable to write file: " + fname);
    }
}
*/

/*
void ClusterDetailView::slot_context_menu(const QPoint& pos)
{
    QMenu M;
    QAction* export_image = M.addAction("Export Image");
    QAction* toggle_stdev_shading = M.addAction("Toggle std. dev. shading");
    QAction* selected = M.exec(this->mapToGlobal(pos));
    if (selected == export_image) {
        d->export_image();
    } else if (selected == toggle_stdev_shading) {
        d->toggle_stdev_shading();
    }
}
*/

void ClusterDetailView::slot_export_image()
{
    d->export_image();
}

void ClusterDetailView::slot_toggle_stdev_shading()
{
    d->toggle_stdev_shading();
}

void ClusterDetailView::slot_zoom_in()
{
    d->zoom(1.2);
}

void ClusterDetailView::slot_zoom_out()
{
    d->zoom(1 / 1.2);
}

void ClusterDetailView::slot_vertical_zoom_in()
{
    d->m_vscale_factor *= 1.15;
    update();
}

void ClusterDetailView::slot_vertical_zoom_out()
{
    d->m_vscale_factor /= 1.15;
    update();
}

void ClusterDetailView::slot_export_static_view()
{
    //QSettings settings("SCDA", "MountainView");
    //QString default_dir = settings.value("default_export_dir", "").toString();
    QString default_dir = QDir::currentPath();
    QString fname = QFileDialog::getSaveFileName(this, "Export static cluster details view", default_dir, "*.smv");
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

void ClusterDetailView::slot_update_sort_order()
{
    this->update();
}

void ClusterDetailView::slot_view_properties()
{
    ClusterDetailViewPropertiesDialog dlg;
    dlg.setProperties(d->m_properties);
    if (dlg.exec() == QDialog::Accepted) {
        d->m_properties = dlg.properties();
        this->update();
    }
}

void ClusterDetailView::slot_export_template_waveforms()
{
    Mda32 X = d->get_template_waveforms_for_export();
    QString default_dir = QDir::currentPath();
    QString fname = QFileDialog::getSaveFileName(this, "Export template waveforms...", default_dir, "*.mda");
    if (fname.isEmpty())
        return;
    if (!X.write32(fname)) {
        qWarning() << "Problem writing file: " + fname;
    }
}

void ClusterDetailView::slot_export_template_waveform_stdevs()
{
    Mda32 X = d->get_template_waveform_stdevs_for_export();
    QString default_dir = QDir::currentPath();
    QString fname = QFileDialog::getSaveFileName(this, "Export template waveform stdevs...", default_dir, "*.mda");
    if (fname.isEmpty())
        return;
    if (!X.write32(fname)) {
        qWarning() << "Problem writing file: " + fname;
    }
}

void ClusterDetailViewPrivate::compute_total_time()
{
    MVContext* c = qobject_cast<MVContext*>(q->mvContext());
    Q_ASSERT(c);

    m_total_time_sec = c->currentTimeseries().N2() / c->sampleRate();
}

void ClusterDetailViewPrivate::set_hovered_k(int k)
{
    if (k == m_hovered_k)
        return;
    m_hovered_k = k;
    q->update();
}

int ClusterDetailViewPrivate::find_view_index_at(QPoint pos)
{
    for (int i = 0; i < m_views.count(); i++) {
        if (m_views[i]->rect().contains(pos))
            return i;
    }
    return -1;
}

ClusterView* ClusterDetailViewPrivate::find_view_for_k(int k)
{
    for (int i = 0; i < m_views.count(); i++) {
        if (m_views[i]->k() == k)
            return m_views[i];
    }
    return 0;
}

int ClusterDetailViewPrivate::find_view_index_for_k(int k)
{
    for (int i = 0; i < m_views.count(); i++) {
        if (m_views[i]->k() == k)
            return i;
    }
    return -1;
}

void ClusterDetailViewPrivate::ensure_view_visible(ClusterView* V)
{
    double x0 = V->x_position_before_scaling * m_space_ratio;
    if (x0 < m_scroll_x) {
        m_scroll_x = x0 - 100;
        if (m_scroll_x < 0)
            m_scroll_x = 0;
    }
    else if (x0 > m_scroll_x + q->width()) {
        m_scroll_x = x0 - q->width() + 100;
    }
}

void ClusterDetailViewPrivate::zoom(double factor)
{
    MVContext* c = qobject_cast<MVContext*>(q->mvContext());
    Q_ASSERT(c);

    int current_k = c->currentCluster();
    if ((current_k >= 0) && (find_view_for_k(current_k))) {
        ClusterView* view = find_view_for_k(current_k);
        double current_screen_x = view->x_position_before_scaling * m_space_ratio - m_scroll_x;
        m_space_ratio *= factor;
        m_scroll_x = view->x_position_before_scaling * m_space_ratio - current_screen_x;
        if (m_scroll_x < 0)
            m_scroll_x = 0;
    }
    else {
        m_space_ratio *= factor;
    }
    q->update();
}

QString ClusterDetailViewPrivate::group_label_for_k(int k)
{
    MVContext* c = qobject_cast<MVContext*>(q->mvContext());
    Q_ASSERT(c);

    if (c->viewMerged()) {
        return QString("%1").arg(c->clusterMerge().clusterLabelText(k));
    }
    else {
        return QString("%1").arg(k);
    }
}

int ClusterDetailViewPrivate::get_current_view_index()
{
    MVContext* c = qobject_cast<MVContext*>(q->mvContext());
    Q_ASSERT(c);

    int k = c->currentCluster();
    if (k < 0)
        return -1;
    return find_view_index_for_k(k);
}

QColor lighten(QColor col, float val)
{
    int r = col.red() * val;
    if (r > 255)
        r = 255;
    int g = col.green() * val;
    if (g > 255)
        g = 255;
    int b = col.blue() * val;
    if (b > 255)
        b = 255;
    return QColor(r, g, b);
}

QString truncate_based_on_font_and_width(QString txt, QFont font, double width)
{
    QFontMetrics fm(font);
    if (fm.width(txt) > width) {
        while ((txt.count() > 3) && (fm.width(txt + "...") > width)) {
            txt = txt.mid(0, txt.count() - 1);
        }
        return txt + "...";
    }
    return txt;
}

void ClusterView::paint(QPainter* painter, QRectF rect, bool render_image_mode)
{
    MVContext* c = qobject_cast<MVContext*>(q->mvContext());
    Q_ASSERT(c);

    int xmargin = 1;
    int ymargin = 8;
    QRectF rect2(rect.x() + xmargin, rect.y() + ymargin, rect.width() - xmargin * 2, rect.height() - ymargin * 2);
    painter->setClipRect(rect, Qt::IntersectClip);

    QColor background_color = c->color("view_background");
    if (m_highlighted)
        background_color = c->color("view_background_highlighted");
    else if (m_selected)
        background_color = c->color("view_background_selected");
    else if (m_hovered)
        background_color = c->color("view_background_hovered");
    if (render_image_mode)
        background_color = Qt::white;
    painter->fillRect(rect, QColor(220, 220, 225));
    painter->fillRect(rect2, background_color);

    QPen pen_frame;
    pen_frame.setWidth(1);
    if (render_image_mode)
        pen_frame.setWidth(2);
    if (m_selected)
        pen_frame.setColor(c->color("view_frame_selected"));
    else
        pen_frame.setColor(c->color("view_frame"));
    if (render_image_mode)
        pen_frame.setColor(Qt::white);
    painter->setPen(pen_frame);
    painter->drawRect(rect2);

    Mda32 template0 = m_CD.template0;
    Mda32 stdev0 = m_CD.stdev0;
    int M = template0.N1();
    int T = template0.N2();
    int Tmid = (int)((T + 1) / 2) - 1;
    m_T = T;

    int top_height = 10 + d->m_properties.cluster_number_font_size, bottom_height = 60;
    if (render_image_mode)
        bottom_height = 20;
    m_rect = rect;
    m_top_rect = QRectF(rect2.x(), rect2.y(), rect2.width(), top_height);
    m_template_rect = QRectF(rect2.x(), rect2.y() + top_height, rect2.width(), rect2.height() - bottom_height - top_height);
    m_bottom_rect = QRectF(rect2.x(), rect2.y() + rect2.height() - bottom_height, rect2.width(), bottom_height);

    {
        //the midline
        QColor midline_color = lighten(background_color, 0.9);
        QPointF pt0 = template_coord2pix(0, Tmid, 0);
        QPen pen;
        pen.setWidth(1);
        pen.setColor(midline_color);
        painter->setPen(pen);
        painter->drawLine(pt0.x(), rect2.bottom() - bottom_height, pt0.x(), rect2.top() + top_height);
    }

    for (int m = 0; m < M; m++) {
        QColor col = c->channelColor(m);
        QPen pen;
        pen.setWidth(1);
        pen.setColor(col);
        if (render_image_mode)
            pen.setWidth(3);
        painter->setPen(pen);
        if (d->m_stdev_shading) {
            QColor quite_light_gray(200, 200, 205);
            QPainterPath path;
            for (int t = 0; t < T; t++) {
                QPointF pt = template_coord2pix(m, t, template0.value(m, t) - stdev0.value(m, t));
                if (t == 0)
                    path.moveTo(pt);
                else
                    path.lineTo(pt);
            }
            for (int t = T - 1; t >= 0; t--) {
                QPointF pt = template_coord2pix(m, t, template0.value(m, t) + stdev0.value(m, t));
                path.lineTo(pt);
            }
            for (int t = 0; t <= 0; t++) {
                QPointF pt = template_coord2pix(m, t, template0.value(m, t) - stdev0.value(m, t));
                path.lineTo(pt);
            }
            painter->fillPath(path, QBrush(quite_light_gray));
        }
        { // the template
            QPainterPath path;
            for (int t = 0; t < T; t++) {
                QPointF pt = template_coord2pix(m, t, template0.value(m, t));
                if (t == 0)
                    path.moveTo(pt);
                else
                    path.lineTo(pt);
            }
            painter->drawPath(path);
        }
    }

    QFont font = painter->font();
    QString txt;
    QRectF RR;

    bool compressed_info = false;
    if (rect2.width() < 60)
        compressed_info = true;

    // Group (or cluster) label
    QString group_label = d->group_label_for_k(m_CD.k);
    {
        txt = QString("%1").arg(group_label);
        //font.setPixelSize(16);
        //if (compressed_info)
        //font.setPixelSize(12);
        font.setPixelSize(d->m_properties.cluster_number_font_size);

        txt = truncate_based_on_font_and_width(txt, font, m_top_rect.width());

        QPen pen;
        pen.setWidth(1);
        pen.setColor(c->color("cluster_label"));
        painter->setFont(font);
        painter->setPen(pen);
        painter->drawText(m_top_rect, Qt::AlignCenter | Qt::AlignBottom, txt);
    }

    font.setPixelSize(11);
    int text_height = 13;

    if (!render_image_mode) {
        // Stats at and tags at bottom
        if (!compressed_info) {
            RR = QRectF(m_bottom_rect.x(), m_bottom_rect.y() + m_bottom_rect.height() - text_height, m_bottom_rect.width(), text_height);
            txt = QString("%1 spikes").arg(m_CD.num_events);
            QPen pen;
            pen.setWidth(1);
            pen.setColor(c->color("info_text"));
            painter->setFont(font);
            painter->setPen(pen);
            painter->drawText(RR, Qt::AlignCenter | Qt::AlignBottom, txt);
        }

        {
            QPen pen;
            pen.setWidth(1);
            RR = QRectF(m_bottom_rect.x(), m_bottom_rect.y() + m_bottom_rect.height() - text_height * 2, m_bottom_rect.width(), text_height);
            double rate = m_CD.num_events * 1.0 / d->m_total_time_sec;
            pen.setColor(get_firing_rate_text_color(rate));
            if (!compressed_info)
                txt = QString("%1 sp/sec").arg(QString::number(rate, 'g', 2));
            else
                txt = QString("%1").arg(QString::number(rate, 'g', 2));
            painter->setFont(font);
            painter->setPen(pen);
            painter->drawText(RR, Qt::AlignCenter | Qt::AlignBottom, txt);
        }

        {
            QPen pen;
            pen.setWidth(1);
            RR = QRectF(m_bottom_rect.x(), m_bottom_rect.y() + m_bottom_rect.height() - text_height * 3, m_bottom_rect.width(), text_height);
            QJsonArray aa = m_attributes["tags"].toArray();
            QStringList aa_strlist = jsonarray2stringlist(aa);
            txt = aa_strlist.join(" ");
            painter->setFont(font);
            painter->setPen(pen);
            painter->drawText(RR, Qt::AlignCenter | Qt::AlignBottom, txt);
        }
    }
}

double ClusterView::spaceNeeded()
{
    return 1;
}

void ClusterDetailViewPrivate::do_paint(QPainter& painter, int W_in, int H_in, bool render_image_mode)
{
    MVContext* c = qobject_cast<MVContext*>(q->mvContext());
    Q_ASSERT(c);

    painter.setRenderHint(QPainter::Antialiasing);

    QColor background_color = c->color("background");
    if (render_image_mode)
        background_color = Qt::white;
    painter.fillRect(0, 0, W_in, H_in, background_color);

    int right_margin = 10; //make some room for the icon
    int left_margin = 30; //make room for the axis
    int W = W_in - right_margin - left_margin;
    int H = H_in;

    painter.setClipRect(QRectF(left_margin, 0, W, H));

    QList<ClusterData> cluster_data_merged;
    if (c->viewMerged()) {
        cluster_data_merged = merge_cluster_data(c->clusterMerge(), m_cluster_data);
    }
    else {
        cluster_data_merged = m_cluster_data;
    }

    sort_cluster_data(cluster_data_merged);

    int old_view_count = m_views.count();

    qDeleteAll(m_views);
    m_views.clear();
    QList<int> selected_clusters = c->selectedClusters();
    for (int i = 0; i < cluster_data_merged.count(); i++) {
        ClusterData CD = cluster_data_merged[i];
        if (c->clusterIsVisible(CD.k)) {
            ClusterView* V = new ClusterView(q, this);
            V->setStdevShading(m_stdev_shading);
            V->setHighlighted(CD.k == c->currentCluster());
            V->setSelected(selected_clusters.contains(CD.k));
            V->setHovered(CD.k == m_hovered_k);
            V->setClusterData(CD);
            V->setAttributes(c->clusterAttributes(CD.k));
            m_views << V;
        }
    }

    if (old_view_count != m_views.count()) {
        m_space_ratio = 0; //trigger zoom out
    }

    double total_space_needed = 0;
    for (int i = 0; i < m_views.count(); i++) {
        total_space_needed += m_views[i]->spaceNeeded();
    }
    if (m_scroll_x < 0)
        m_scroll_x = 0;
    if (total_space_needed * m_space_ratio - m_scroll_x < W) {
        m_scroll_x = total_space_needed * m_space_ratio - W;
        if (m_scroll_x < 0)
            m_scroll_x = 0;
    }
    if ((m_scroll_x == 0) && (total_space_needed * m_space_ratio < W)) {
        m_space_ratio = W / total_space_needed;
        if (m_space_ratio > 300)
            m_space_ratio = 300;
    }

    double space_ratio = m_space_ratio;
    double scroll_x = m_scroll_x;

    if (render_image_mode) {
        space_ratio = W / total_space_needed;
        scroll_x = 0;
    }

    ChannelSpacingInfo csi = compute_channel_spacing_info(cluster_data_merged, m_vscale_factor);

    float x0_before_scaling = 0;
    ClusterView* first_view = 0;
    for (int i = 0; i < m_views.count(); i++) {
        ClusterView* V = m_views[i];
        QRectF rect(left_margin + x0_before_scaling * space_ratio - scroll_x, 0, V->spaceNeeded() * space_ratio, H);
        V->setChannelSpacingInfo(csi);
        if ((rect.x() + rect.width() >= left_margin) && (rect.x() <= left_margin + W)) {
            first_view = V;
            QRegion save_clip_region = painter.clipRegion();
            V->paint(&painter, rect, render_image_mode);
            painter.setClipRegion(save_clip_region);
        }
        V->x_position_before_scaling = x0_before_scaling;
        x0_before_scaling += V->spaceNeeded();
    }

    painter.setClipRect(0, 0, W_in, H_in);
    if ((first_view) && (!render_image_mode)) {
        double fac0 = 0.95; //to leave a bit of a gap
        ClusterView* V = first_view;
        int M = cluster_data_merged[0].template0.N1();
        Q_UNUSED(M)
        //for (int m = 0; m < M; m++) {
        for (int m = 0; m <= 0; m++) {
            QPointF pt1 = V->template_coord2pix(m, 0, -0.5 / csi.vert_scaling_factor * fac0);
            QPointF pt2 = V->template_coord2pix(m, 0, 0.5 / csi.vert_scaling_factor * fac0);
            pt1.setX(left_margin);
            pt2.setX(left_margin);
            pt1.setY(pt1.y());
            pt2.setY(pt2.y());
            draw_axis_opts opts;
            opts.pt1 = pt1;
            opts.pt2 = pt2;
            opts.draw_tick_labels = false;
            opts.tick_length = 0;
            opts.draw_range = true;
            opts.minval = -0.5 / csi.vert_scaling_factor * fac0;
            opts.maxval = 0.5 / csi.vert_scaling_factor * fac0;
            opts.orientation = Qt::Vertical;
            draw_axis(&painter, opts);
        }
    }

    if ((q->isCalculating()) && (!render_image_mode)) {
        QFont font = painter.font();
        font.setPointSize(20);
        painter.setFont(font);
        painter.fillRect(QRectF(0, 0, q->width(), q->height()), c->color("calculation-in-progress"));
        painter.drawText(QRectF(left_margin, 0, W, H), Qt::AlignCenter | Qt::AlignVCenter, "Calculating...");
    }
}

void ClusterDetailViewPrivate::export_image()
{
    int W = m_properties.export_image_width;
    int H = m_properties.export_image_height;
    QImage img = q->renderImage(W, H);
    user_save_image(img);
}

void ClusterDetailViewPrivate::toggle_stdev_shading()
{
    m_stdev_shading = !m_stdev_shading;
    q->update();
}

void ClusterDetailViewPrivate::shift_select_clusters_between(int k1, int k2)
{
    MVContext* c = qobject_cast<MVContext*>(q->mvContext());
    Q_ASSERT(c);

    QSet<int> selected_clusters = c->selectedClusters().toSet();
    int ind1 = find_view_index_for_k(k1);
    int ind2 = find_view_index_for_k(k2);
    if ((ind1 >= 0) && (ind2 >= 0)) {
        for (int ii = qMin(ind1, ind2); ii <= qMax(ind1, ind2); ii++) {
            selected_clusters.insert(m_views[ii]->k());
        }
    }
    else if (ind1 >= 0) {
        selected_clusters.insert(m_views[ind1]->k());
    }
    else if (ind2 >= 0) {
        selected_clusters.insert(m_views[ind2]->k());
    }
    c->setSelectedClusters(QList<int>::fromSet(selected_clusters));
}

ClusterData combine_cluster_data_group(const QList<ClusterData>& group, ClusterData main_CD)
{
    ClusterData ret;
    ret.k = main_CD.k;
    ret.channel = main_CD.channel;
    Mda sum0;
    Mda sumsqr0;
    if (group.count() > 0) {
        sum0.allocate(group[0].template0.N1(), group[0].template0.N2(), group[0].template0.N3());
        sumsqr0.allocate(group[0].template0.N1(), group[0].template0.N2(), group[0].template0.N3());
    }
    double total_weight = 0;
    ret.num_events = group.count();
    for (int i = 0; i < group.count(); i++) {
        ret.num_events += group[i].num_events;
        double weight = group[i].num_events;
        for (int i3 = 0; i3 < sum0.N3(); i3++) {
            for (int i2 = 0; i2 < sum0.N2(); i2++) {
                for (int i1 = 0; i1 < sum0.N1(); i1++) {
                    double val1 = group[i].template0.value(i1, i2, i3) * weight;
                    double val2 = group[i].stdev0.value(i1, i2, i3) * group[i].stdev0.value(i1, i2, i3) * weight;
                    sum0.setValue(sum0.value(i1, i2, i3) + val1, i1, i2, i3);
                    sumsqr0.setValue(sumsqr0.value(i1, i2, i3) + (val2 + val1 * val1 / weight), i1, i2, i3);
                }
            }
        }
        total_weight += weight;
    }
    ret.template0.allocate(sum0.N1(), sum0.N2(), sum0.N3());
    ret.stdev0.allocate(sum0.N1(), sum0.N2(), sum0.N3());
    if (total_weight) {
        for (int i = 0; i < ret.template0.totalSize(); i++) {
            double val_sum = sum0.get(i);
            double val_sumsqr = sumsqr0.get(i);
            ret.template0.set(val_sum / total_weight, i);
            ret.stdev0.set(sqrt((val_sumsqr - val_sum * val_sum / total_weight) / total_weight), i);
        }
    }
    return ret;
}

QList<ClusterData> ClusterDetailViewPrivate::merge_cluster_data(const ClusterMerge& CM, const QList<ClusterData>& CD)
{
    QList<ClusterData> ret;
    for (int i = 0; i < CD.count(); i++) {
        if (CM.representativeLabel(CD[i].k) == CD[i].k) {
            QList<ClusterData> group;
            for (int j = 0; j < CD.count(); j++) {
                if (CM.representativeLabel(CD[j].k) == CD[i].k) {
                    group << CD[j];
                }
            }
            ret << combine_cluster_data_group(group, CD[i]);
        }
        else {
            ClusterData CD0;
            CD0.k = CD[i].k;
            CD0.channel = CD[i].channel;
            ret << CD0;
        }
    }
    return ret;
}

Mda32 ClusterDetailViewPrivate::get_template_waveforms_for_export()
{
    MVContext* c = qobject_cast<MVContext*>(q->mvContext());
    Q_ASSERT(c);

    QList<Mda32> templates;
    for (int ii = 0; ii < m_cluster_data.count(); ii++) {
        int k = m_cluster_data[ii].k;
        if (c->clusterIsVisible(k)) {
            templates << m_cluster_data[ii].template0;
        }
    }
    int A = templates.count();
    int M = templates.value(0).N1();
    int T = templates.value(0).N2();
    Mda32 X(M, T, A);
    for (int a = 0; a < A; a++) {
        X.setChunk(templates[a], 0, 0, a);
    }
    return X;
}

Mda32 ClusterDetailViewPrivate::get_template_waveform_stdevs_for_export()
{
    MVContext* c = qobject_cast<MVContext*>(q->mvContext());
    Q_ASSERT(c);

    QList<Mda32> stdevs;
    for (int ii = 0; ii < m_cluster_data.count(); ii++) {
        int k = m_cluster_data[ii].k;
        if (c->clusterIsVisible(k)) {
            stdevs << m_cluster_data[ii].stdev0;
        }
    }
    int A = stdevs.count();
    int M = stdevs.value(0).N1();
    int T = stdevs.value(0).N2();
    Mda32 X(M, T, A);
    for (int a = 0; a < A; a++) {
        X.setChunk(stdevs[a], 0, 0, a);
    }
    return X;
}

QPointF ClusterView::template_coord2pix(int m, double t, double val)
{
    double pcty = m_csi.channel_locations.value(m) - m_csi.channel_location_spacing * val * m_csi.vert_scaling_factor; //negative because (0,0) is top-left, not bottom-right
    double pctx = 0;
    if (m_T)
        pctx = (t + 0.5) / m_T;
    int margx = 4;
    int margy = 4;
    float x0 = m_template_rect.x() + margx + pctx * (m_template_rect.width() - margx * 2);
    float y0 = m_template_rect.y() + margy + pcty * (m_template_rect.height() - margy * 2);
    return QPointF(x0, y0);
}

void ClusterDetailViewPrivate::sort_cluster_data(QList<ClusterData>& CD)
{
    MVContext* c = qobject_cast<MVContext*>(q->mvContext());
    Q_ASSERT(c);

    QList<double> cluster_order_scores = c->clusterOrderScores();
    QVector<double> sort_scores;
    for (int i = 0; i < CD.count(); i++) {
        sort_scores << cluster_order_scores.value(CD[i].k - 1, 0);
    }

    /*
    QVector<double> sort_scores;
    for (int i = 0; i < CD.count(); i++) {
        sort_scores << compute_sort_score(CD[i]);
    }

    */

    QList<ClusterData> ret;
    QList<int> inds = get_sort_indices(sort_scores);
    for (int i = 0; i < inds.count(); i++) {
        ret << CD[inds[i]];
    }
    CD = ret;
}

double compute_centeredness_score(const Mda32& template0)
{
    int M = template0.N1();
    int T = template0.N2();
    QVector<double> tvals;
    QVector<double> yvals;
    for (int t = 0; t < T; t++) {
        for (int m = 0; m < M; m++) {
            yvals << qAbs(template0.value(m, t));
            tvals << t;
        }
    }
    double sum_yt = 0, sum_y = 0;
    for (int i = 0; i < tvals.count(); i++) {
        sum_yt += yvals[i] * tvals[i];
        sum_y += yvals[i];
    }
    if (sum_y) {
        double center_of_mass = sum_yt / sum_y;
        double score = 1 - qAbs((center_of_mass - T / 2) / (T / 2));
        return score;
    }
    else
        return 0;
}

double ClusterDetailViewPrivate::compute_sort_score(const ClusterData& CD)
{
    double score = 0;

    QString score_mode = m_sort_order_chooser->currentData().toString();

    if (score_mode == "quality") {
        //norm of signal divided by norm of noise
        double norm1 = MLCompute::norm(CD.template0.totalSize(), CD.template0.constDataPtr());
        double norm2 = MLCompute::norm(CD.stdev0.totalSize(), CD.stdev0.constDataPtr());
        if (norm2)
            score += norm1 / norm2;

        score += compute_centeredness_score(CD.template0);

        return -score;
    }
    else {
        return 0;
    }
}

QColor ClusterView::get_firing_rate_text_color(double rate)
{
    if (rate <= 0.1)
        return QColor(220, 220, 220);
    if (rate <= 1)
        return QColor(150, 150, 150);
    if (rate <= 10)
        return QColor(0, 50, 0);
    return QColor(50, 0, 0);
}

DiskReadMda mp_compute_templates(const QString& mlproxy_url, const QString& timeseries, const QString& firings, int clip_size)
{
    TaskProgress task(TaskProgress::Calculate, "mp_compute_templates");
    task.log("mlproxy_url: " + mlproxy_url);
    MountainProcessRunner X;
    QString processor_name = "compute_templates";
    X.setProcessorName(processor_name);

    QMap<QString, QVariant> params;
    params["timeseries"] = timeseries;
    params["firings"] = firings;
    params["clip_size"] = clip_size;
    X.setInputParameters(params);
    X.setMLProxyUrl(mlproxy_url);

    QString templates_fname = X.makeOutputFilePath("templates");

    task.log("X.compute()");
    X.runProcess();
    task.log("Returning DiskReadMda: " + templates_fname);
    DiskReadMda ret(templates_fname);
    //ret.setRemoteDataType("float32");
    //ret.setRemoteDataType("float32_q8"); //to save download time!
    return ret;
}

void mp_compute_templates_stdevs(DiskReadMda32& templates_out, DiskReadMda32& stdevs_out, const QString& mlproxy_url, const QString& timeseries, const QString& firings, int clip_size)
{
    TaskProgress task(TaskProgress::Calculate, "mp_compute_templates_stdevs");
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

void ClusterDetailViewCalculator::compute()
{
    TaskProgress task(TaskProgress::Calculate, "Cluster Detail");
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

    QVector<double> times;
    QVector<int> channels, labels;
    QVector<double> peaks;

    task.log("Setting up times/channels/labels/peaks");
    task.setProgress(0.2);
    for (int i = 0; i < L; i++) {
        times << firings.value(1, i) - 1; //convert to 0-based indexing
        channels << (int)firings.value(0, i) - 1; //convert to 0-based indexing
        labels << (int)firings.value(2, i);
        peaks << firings.value(3, i);
    }

    if (MLUtil::threadInterruptRequested()) {
        task.error("Halted *");
        return;
    }
    task.log("Clearing data");
    cluster_data.clear();

    task.setLabel("Computing templates");
    task.setProgress(0.4);
    int K = 0;
    for (int i = 0; i < L; i++)
        if (labels[i] > K)
            K = labels[i];

    QString timeseries_path = timeseries.makePath();
    QString firings_path = firings.makePath();
    /*
    this->setStatus("", "mscmd_compute_templates: "+mscmdserver_url+" timeseries_path="+timeseries_path+" firings_path="+firings_path, 0.6);
    DiskReadMda templates0 = mscmd_compute_templates(mscmdserver_url, timeseries_path, firings_path, T);
    */

    task.log("mp_compute_templates_stdevs: " + mlproxy_url + " timeseries_path=" + timeseries_path + " firings_path=" + firings_path);
    task.setProgress(0.6);
    //DiskReadMda templates0 = mp_compute_templates(mlproxy_url, timeseries_path, firings_path, T);
    DiskReadMda32 templates0, stdevs0;
    mp_compute_templates_stdevs(templates0, stdevs0, mlproxy_url, timeseries_path, firings_path, T);
    if (MLUtil::threadInterruptRequested()) {
        task.error("Halted **");
        return;
    }

    task.setLabel("Setting cluster data");
    task.setProgress(0.75);
    int K2 = K;
    foreach (int kk, clusters_to_force_show) {
        if (kk > K2)
            K2 = kk;
    }
    for (int k = 1; k <= K2; k++) {
        if (MLUtil::threadInterruptRequested()) {
            task.error("Halted ***");
            return;
        }
        ClusterData CD;
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
            qWarning() << "Unable to read chunk of templates in cluster detail view";
            return;
        }
        if (!stdevs0.readChunk(CD.stdev0, 0, 0, k - 1, M, T, 1)) {
            qWarning() << "Unable to read chunk of stdevs in cluster detail view";
        }
        if (!MLUtil::threadInterruptRequested()) {
            if ((CD.num_events > 0) || (clusters_to_force_show.contains(k))) {
                cluster_data << CD;
            }
        }
    }
}

QJsonObject ClusterDetailViewCalculator::exportStaticOutput()
{
    QJsonObject ret;
    ret["version"] = "ClusterDetailViewCalculator-0.1";

    QJsonArray cd_list;
    for (int i = 0; i < cluster_data.count(); i++) {
        QJsonObject cd = cluster_data[i].toJsonObject();
        cd_list.push_back(cd);
    }
    ret["cluster_data"] = cd_list;
    return ret;
}

void ClusterDetailViewCalculator::loadStaticOutput(const QJsonObject& X)
{
    cluster_data.clear();
    QJsonArray cd_list = X["cluster_data"].toArray();
    for (int i = 0; i < cd_list.count(); i++) {
        ClusterData CD;
        CD.fromJsonObject(cd_list[i].toObject());
        cluster_data << CD;
    }
    loaded_from_static_output = true;
}

ClusterView::ClusterView(ClusterDetailView* q0, ClusterDetailViewPrivate* d0)
{
    q = q0;
    d = d0;
    m_T = 1;
    m_highlighted = false;
    m_hovered = false;
    x_position_before_scaling = 0;
    m_stdev_shading = false;
}
void ClusterView::setClusterData(const ClusterData& CD)
{
    m_CD = CD;
}
void ClusterView::setAttributes(QJsonObject aa)
{
    m_attributes = aa;
}

int ClusterView::k()
{
    return this->clusterData()->k;
}
void ClusterView::setChannelSpacingInfo(const ChannelSpacingInfo& csi)
{
    m_csi = csi;
    m_T = 0;
}
void ClusterView::setHighlighted(bool val)
{
    m_highlighted = val;
}
void ClusterView::setHovered(bool val)
{
    m_hovered = val;
}
void ClusterView::setSelected(bool val)
{
    m_selected = val;
}
void ClusterView::setStdevShading(bool val)
{
    m_stdev_shading = val;
}
bool ClusterView::stdevShading()
{
    return m_stdev_shading;
}

ClusterData* ClusterView::clusterData()
{
    return &m_CD;
}
QRectF ClusterView::rect()
{
    return m_rect;
}

QJsonObject ClusterData::toJsonObject()
{
    QJsonObject ret;
    ret["channel"] = this->channel;
    //ret["inds"]=MLUtil::toJsonValue(this->inds); //to big, and unnecessary.... same with peaks and times
    ret["num_events"] = this->num_events;
    ret["k"] = this->k;
    ret["stdev0"] = MLUtil::toJsonValue(this->stdev0.toByteArray32());
    ret["template0"] = MLUtil::toJsonValue(this->template0.toByteArray32());
    return ret;
}

void ClusterData::fromJsonObject(const QJsonObject& X)
{
    this->channel = X["channel"].toInt();
    this->num_events = X["num_events"].toInt();
    this->k = X["k"].toInt();
    QByteArray tmp;
    {
        MLUtil::fromJsonValue(tmp, X["stdev0"]);
        this->stdev0.fromByteArray(tmp);
    }
    {
        MLUtil::fromJsonValue(tmp, X["template0"]);
        this->template0.fromByteArray(tmp);
    }
}
