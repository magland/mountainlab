#include "mvclusterdetailwidget.h"
#include "taskprogress.h"

#include <QPainter>
#include "mvutils.h"
#include "msmisc.h"
#include <QProgressDialog>
#include <QTime>
#include <QMap>
#include <QDebug>
#include <QMouseEvent>
#include <QSet>
#include <QSettings>
#include <QImageWriter>
#include "extract_clips.h"
#include <QFileDialog>
#include <QMenu>
#include <QMessageBox>
#include <QTimer>
#include "compute_templates_0.h"
#include "computationthread.h"
#include "mountainprocessrunner.h"
#include "toolbuttonmenu.h"

struct ClusterData {
    ClusterData()
    {
        k = 0;
        channel = 0;
    }

    int k;
    int channel;
    Mda template0;
    Mda stdev0;
    QList<int> inds;
    QList<double> times;
    QList<double> peaks;
};

struct ChannelSpacingInfo {
    QList<double> channel_locations;
    double vert_scaling_factor;
};

class MVClusterDetailWidgetCalculator : public ComputationThread {
public:
    //input
    //QString mscmdserver_url;
    QString mlproxy_url;
    DiskReadMda timeseries;
    DiskReadMda firings;
    int clip_size;

    virtual void compute();

    //output
    QList<ClusterData> cluster_data;
};

class ClusterView {
public:
    friend class MVClusterDetailWidgetPrivate;
    friend class MVClusterDetailWidget;
    ClusterView(MVClusterDetailWidget* q0, MVClusterDetailWidgetPrivate* d0)
    {
        q = q0;
        d = d0;
        m_T = 1;
        m_highlighted = false;
        m_hovered = false;
        x_position_before_scaling = 0;
        m_stdev_shading = false;
    }
    void setClusterData(const ClusterData& CD)
    {
        m_CD = CD;
    }
    void setAttributes(QJsonObject aa)
    {
        m_attributes = aa;
    }

    int k()
    {
        return this->clusterData()->k;
    }
    void setChannelSpacingInfo(const ChannelSpacingInfo& csi)
    {
        m_csi = csi;
        m_T = 0;
    }
    void setHighlighted(bool val)
    {
        m_highlighted = val;
    }
    void setHovered(bool val)
    {
        m_hovered = val;
    }
    void setSelected(bool val)
    {
        m_selected = val;
    }
    void setStdevShading(bool val)
    {
        m_stdev_shading = val;
    }
    bool stdevShading()
    {
        return m_stdev_shading;
    }

    void paint(QPainter* painter, QRectF rect);
    double spaceNeeded();
    ClusterData* clusterData()
    {
        return &m_CD;
    }
    QRectF rect()
    {
        return m_rect;
    }

    MVClusterDetailWidget* q;
    MVClusterDetailWidgetPrivate* d;
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

    QPointF template_coord2pix(int m, double t, double val);
    QColor get_firing_rate_text_color(double rate);
    QColor get_cluster_assessment_text_color(QString aa);
};

class MVClusterDetailWidgetPrivate {
public:
    MVClusterDetailWidget* q;

    //QString m_mscmdserver_url;
    QString m_mlproxy_url;
    DiskReadMda m_timeseries;
    DiskReadMda m_firings;
    double m_samplerate;
    QList<int> m_group_numbers;

    int m_clip_size;
    QList<ClusterData> m_cluster_data;

    double m_vscale_factor;
    double m_space_ratio;
    double m_scroll_x;

    QProgressDialog* m_progress_dialog;
    QMap<QString, QColor> m_colors;
    QList<QColor> m_channel_colors;
    double m_total_time_sec;
    int m_current_k;
    QSet<int> m_selected_ks;
    int m_hovered_k;
    double m_anchor_x;
    double m_anchor_scroll_x;
    int m_anchor_view_index;
    MVClusterDetailWidgetCalculator m_calculator;
    bool m_stdev_shading;

    MVViewAgent* m_view_agent;

    QList<ClusterView*> m_views;

    void compute_total_time();
    void set_current_k(int k);
    void set_hovered_k(int k);
    int find_view_index_at(QPoint pos);
    ClusterView* find_view_for_k(int k);
    int find_view_index_for_k(int k);
    void ensure_view_visible(ClusterView* V);
    void zoom(double factor);
    QString group_label_for_k(int k);
    bool has_nontrivial_group_numbers();
    int get_current_view_index();
    void do_paint(QPainter& painter, int W, int H);
    void export_image();
    void start_calculation();
    void toggle_stdev_shading();

    static QList<ClusterData> merge_cluster_data(const ClusterMerge& CM, const QList<ClusterData>& CD);
};

MVClusterDetailWidget::MVClusterDetailWidget(QWidget* parent)
    : QWidget(parent)
{
    d = new MVClusterDetailWidgetPrivate;
    d->q = this;
    d->m_clip_size = 100;
    d->m_progress_dialog = 0;
    d->m_vscale_factor = 2;
    d->m_total_time_sec = 1;
    d->m_samplerate = 0;
    d->m_current_k = -1;
    d->m_hovered_k = -1;
    d->m_space_ratio = 50;
    d->m_scroll_x = 0;
    d->m_anchor_x = -1;
    d->m_anchor_scroll_x = -1;
    d->m_anchor_view_index = -1;
    d->m_stdev_shading = false;

    d->m_colors["background"] = QColor(240, 240, 240);
    d->m_colors["frame1"] = QColor(245, 245, 245);
    d->m_colors["info_text"] = QColor(80, 80, 80);
    d->m_colors["view_background"] = QColor(245, 245, 245);
    d->m_colors["view_background_highlighted"] = QColor(250, 220, 200);
    d->m_colors["view_background_selected"] = QColor(250, 240, 230);
    d->m_colors["view_background_hovered"] = QColor(240, 245, 240);
    d->m_channel_colors << Qt::black;

    d->m_view_agent = 0;

    this->setFocusPolicy(Qt::StrongFocus);
    this->setMouseTracking(true);
    //this->setContextMenuPolicy(Qt::CustomContextMenu);
    //connect(this, SIGNAL(customContextMenuRequested(QPoint)), this, SLOT(slot_context_menu(QPoint)));

    connect(&d->m_calculator, SIGNAL(computationFinished()), this, SLOT(slot_calculator_finished()));

    this->setContextMenuPolicy(Qt::ActionsContextMenu);
    ToolButtonMenu *MM=new ToolButtonMenu(this);
    MM->setOffset(QSize(4, 4));
    QToolButton *tb=MM->activateOn(this);
    tb->setIconSize(QSize(24,24));
    QIcon icon(":images/gear.png");
    tb->setIcon(icon);
    {
        QAction *a=new QAction("Export image",this);
        this->addAction(a);
        connect(a,SIGNAL(triggered(bool)),this,SLOT(slot_export_image()));
    }
    {
        QAction *a=new QAction("Toggle std. dev. shading",this);
        this->addAction(a);
        connect(a,SIGNAL(triggered(bool)),this,SLOT(slot_toggle_stdev_shading()));
    }
}

MVClusterDetailWidget::~MVClusterDetailWidget()
{
    d->m_calculator.stopComputation(); // important do take care of this before things start getting destructed!
    qDeleteAll(d->m_views);
    delete d;
}

/*
void MVClusterDetailWidget::setMscmdServerUrl(const QString& url)
{
    d->m_mscmdserver_url = url;
}
*/

void MVClusterDetailWidget::setMLProxyUrl(const QString& url)
{
    d->m_mlproxy_url = url;
}

void MVClusterDetailWidget::setTimeseries(DiskReadMda& X)
{
    d->m_timeseries = X;
    d->compute_total_time();
    d->start_calculation();
}

void MVClusterDetailWidget::setFirings(const DiskReadMda& X)
{
    d->m_firings = X;
    d->start_calculation();
}

void MVClusterDetailWidget::setClipSize(int T)
{
    if (d->m_clip_size == T)
        return;
    d->m_clip_size = T;
    d->start_calculation();
}

void MVClusterDetailWidget::setGroupNumbers(const QList<int>& group_numbers)
{
    d->m_group_numbers = group_numbers;
    this->update();
}

void MVClusterDetailWidget::setSampleRate(double freq)
{
    d->m_samplerate = freq;
    d->compute_total_time();
    this->update();
}

void MVClusterDetailWidget::setChannelColors(const QList<QColor>& colors)
{
    d->m_channel_colors = colors;
    update();
}

void MVClusterDetailWidget::setColors(const QMap<QString, QColor>& colors)
{
    d->m_colors = colors;
}

int MVClusterDetailWidget::currentK()
{
    return d->m_current_k;
}

QList<int> MVClusterDetailWidget::selectedKs()
{
    return QList<int>::fromSet(d->m_selected_ks);
}

void MVClusterDetailWidget::setCurrentK(int k)
{
    d->set_current_k(k);
}

bool sets_are_equal(const QSet<int>& S1, const QSet<int>& S2)
{
    foreach(int val, S1)
    {
        if (!S2.contains(val))
            return false;
    }
    foreach(int val, S2)
    {
        if (!S1.contains(val))
            return false;
    }
    return true;
}

void MVClusterDetailWidget::setSelectedKs(const QList<int>& ks_in)
{
    QList<int> ks = ks_in;
    ks.removeAll(0);
    if (sets_are_equal(d->m_selected_ks, ks.toSet()))
        return;
    d->m_selected_ks = ks.toSet();
    emit this->signalSelectedKsChanged();
    update();
}

void MVClusterDetailWidget::setViewAgent(MVViewAgent* agent)
{
    d->m_view_agent = agent;
    if (agent) {
        QObject::connect(agent, SIGNAL(clusterMergeChanged()), this, SLOT(update()));
    }
}

void MVClusterDetailWidget::zoomAllTheWayOut()
{
    d->m_space_ratio = 0;
    update();
}

QImage MVClusterDetailWidget::renderImage(int W, int H)
{
    if (!W)
        W = 1600;
    if (!H)
        H = 800;
    QImage ret = QImage(W, H, QImage::Format_RGB32);
    QPainter painter(&ret);

    int current_k = d->m_current_k;
    QSet<int> selected_ks = d->m_selected_ks;
    d->m_current_k = -1;
    d->m_selected_ks.clear();
    d->do_paint(painter, W, H);
    d->m_current_k = current_k;
    d->m_selected_ks = selected_ks;
    this->update(); //make sure we update, because some internal stuff has changed!

    return ret;
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
    double y0 = 0.5 / M;
    for (int m = 0; m < M; m++) {
        info.channel_locations << y0;
        y0 += 1.0 / M;
    }
    double maxabsval = qMax(maxval, -minval);
    info.vert_scaling_factor = 0.5 / M / maxabsval * vscale_factor;
    return info;
}

void MVClusterDetailWidget::paintEvent(QPaintEvent* evt)
{
    Q_UNUSED(evt)

    QPainter painter(this);

    d->do_paint(painter, width(), height());
}

void MVClusterDetailWidget::keyPressEvent(QKeyEvent* evt)
{
    double factor = 1.15;
    if (evt->key() == Qt::Key_Up) {
        d->m_vscale_factor *= factor;
        update();
    } else if (evt->key() == Qt::Key_Down) {
        d->m_vscale_factor /= factor;
        update();
    } else if ((evt->key() == Qt::Key_Plus) || (evt->key() == Qt::Key_Equal)) {
        d->zoom(1.1);
    } else if (evt->key() == Qt::Key_Minus) {
        d->zoom(1 / 1.1);
    } else if ((evt->key() == Qt::Key_A) && (evt->modifiers() & Qt::ControlModifier)) {
        QList<int> ks;
        for (int i = 0; i < d->m_views.count(); i++) {
            ks << d->m_views[i]->k();
        }
        this->setSelectedKs(ks);
    } else if (evt->key() == Qt::Key_Left) {
        int view_index = d->get_current_view_index();
        if (view_index > 0) {
            int k = d->m_views[view_index - 1]->k();
            QList<int> ks;
            if (evt->modifiers() & Qt::ShiftModifier) {
                ks = this->selectedKs();
                ks << k;
            }
            this->setSelectedKs(ks);
            this->setCurrentK(k);
        }
    } else if (evt->key() == Qt::Key_Right) {
        int view_index = d->get_current_view_index();
        if ((view_index >= 0) && (view_index + 1 < d->m_views.count())) {
            int k = d->m_views[view_index + 1]->k();
            QList<int> ks;
            if (evt->modifiers() & Qt::ShiftModifier) {
                ks = this->selectedKs();
                ks << k;
            }
            this->setSelectedKs(ks);
            this->setCurrentK(k);
        }
    } else
        evt->ignore();
}

void MVClusterDetailWidget::mousePressEvent(QMouseEvent* evt)
{
    QPoint pt = evt->pos();
    d->m_anchor_x = pt.x();
    d->m_anchor_scroll_x = d->m_scroll_x;
}

void MVClusterDetailWidget::mouseReleaseEvent(QMouseEvent* evt)
{
    QPoint pt = evt->pos();

    if ((d->m_anchor_x >= 0) && (qAbs(pt.x() - d->m_anchor_x) > 5)) {
        d->m_scroll_x = d->m_anchor_scroll_x - (pt.x() - d->m_anchor_x);
        d->m_anchor_x = -1;
        update();
        return;
    }
    d->m_anchor_x = -1;

    if (evt->modifiers() & Qt::ControlModifier) {
        int view_index = d->find_view_index_at(pt);
        d->m_anchor_view_index = -1;
        if (view_index >= 0) {
            int k = d->m_views[view_index]->k();
            if (d->m_current_k == k) {
                d->set_current_k(-1);
            }
            if (d->m_selected_ks.contains(k)) {
                d->m_selected_ks.remove(k);
                emit signalSelectedKsChanged();
                update();
            } else {
                d->m_anchor_view_index = view_index;
                if (k)
                    d->m_selected_ks.insert(k);
                emit signalSelectedKsChanged();
                update();
            }
        }
    } else if (evt->modifiers() & Qt::ShiftModifier) {
        int view_index = d->find_view_index_at(pt);
        if (view_index >= 0) {
            if (d->m_anchor_view_index >= 0) {
                int min_index = qMin(d->m_anchor_view_index, view_index);
                int max_index = qMax(d->m_anchor_view_index, view_index);
                for (int i = min_index; i <= max_index; i++) {
                    if (i < d->m_views.count()) {
                        int k = d->m_views[i]->k();
                        if (k)
                            d->m_selected_ks.insert(k);
                    }
                }
                emit signalSelectedKsChanged();
                update();
            }
        }
    } else {
        d->m_anchor_view_index = -1;
        int view_index = d->find_view_index_at(pt);
        if (view_index >= 0) {
            d->m_anchor_view_index = view_index;
            int k = d->m_views[view_index]->k();
            if (d->m_current_k == k) {
            } else {
                d->set_current_k(k);
                d->m_selected_ks.clear();
                if (k)
                    d->m_selected_ks.insert(k);
                emit signalSelectedKsChanged();
                update();
            }
        } else {
            d->set_current_k(-1);
            d->m_selected_ks.clear();
            emit signalSelectedKsChanged();
            update();
        }
    }
}

void MVClusterDetailWidget::mouseMoveEvent(QMouseEvent* evt)
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
    } else {
        d->set_hovered_k(-1);
    }
}

void MVClusterDetailWidget::mouseDoubleClickEvent(QMouseEvent* evt)
{
    Q_UNUSED(evt);
    emit this->signalTemplateActivated();
}

void MVClusterDetailWidget::wheelEvent(QWheelEvent* evt)
{
    int delta = evt->delta();
    double factor = 1;
    if (delta > 0)
        factor = 1.1;
    else
        factor = 1 / 1.1;
    d->zoom(factor);
}

/*
void MVClusterDetailWidget::slot_context_menu(const QPoint& pos)
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

void MVClusterDetailWidget::slot_calculator_finished()
{
    d->m_calculator.stopComputation(); //because I'm paranoid!
    d->m_cluster_data = d->m_calculator.cluster_data;
    this->update();
}

void MVClusterDetailWidget::slot_export_image()
{
    d->export_image();
}

void MVClusterDetailWidget::slot_toggle_stdev_shading()
{
    d->toggle_stdev_shading();
}

void MVClusterDetailWidgetPrivate::compute_total_time()
{
    m_total_time_sec = m_timeseries.N2() / m_samplerate;
}

void MVClusterDetailWidgetPrivate::set_current_k(int k)
{
    if (k == m_current_k)
        return;
    m_current_k = k;
    if (k)
        m_selected_ks.insert(k);
    ClusterView* V = find_view_for_k(k);
    if (V)
        ensure_view_visible(V);
    q->update();
    emit q->signalCurrentKChanged();
}

void MVClusterDetailWidgetPrivate::set_hovered_k(int k)
{
    if (k == m_hovered_k)
        return;
    m_hovered_k = k;
    q->update();
}

int MVClusterDetailWidgetPrivate::find_view_index_at(QPoint pos)
{
    for (int i = 0; i < m_views.count(); i++) {
        if (m_views[i]->rect().contains(pos))
            return i;
    }
    return -1;
}

ClusterView* MVClusterDetailWidgetPrivate::find_view_for_k(int k)
{
    for (int i = 0; i < m_views.count(); i++) {
        if (m_views[i]->k() == k)
            return m_views[i];
    }
    return 0;
}

int MVClusterDetailWidgetPrivate::find_view_index_for_k(int k)
{
    for (int i = 0; i < m_views.count(); i++) {
        if (m_views[i]->k() == k)
            return i;
    }
    return -1;
}

void MVClusterDetailWidgetPrivate::ensure_view_visible(ClusterView* V)
{
    double x0 = V->x_position_before_scaling * m_space_ratio;
    if (x0 < m_scroll_x) {
        m_scroll_x = x0 - 100;
        if (m_scroll_x < 0)
            m_scroll_x = 0;
    } else if (x0 > m_scroll_x + q->width()) {
        m_scroll_x = x0 - q->width() + 100;
    }
}

void MVClusterDetailWidgetPrivate::zoom(double factor)
{
    if ((m_current_k >= 0) && (find_view_for_k(m_current_k))) {
        ClusterView* view = find_view_for_k(m_current_k);
        double current_screen_x = view->x_position_before_scaling * m_space_ratio - m_scroll_x;
        m_space_ratio *= factor;
        m_scroll_x = view->x_position_before_scaling * m_space_ratio - current_screen_x;
        if (m_scroll_x < 0)
            m_scroll_x = 0;
    } else {
        m_space_ratio *= factor;
    }
    q->update();
}

QString MVClusterDetailWidgetPrivate::group_label_for_k(int k)
{
    if (m_group_numbers.isEmpty())
        return QString("%1").arg(k);
    int g = m_group_numbers.value(k);
    for (int i = 1; i < k; i++) {
        if (m_group_numbers[i] == g)
            return "";
    }
    return QString("%1").arg(g);
}

bool MVClusterDetailWidgetPrivate::has_nontrivial_group_numbers()
{
    if (m_group_numbers.isEmpty())
        return false;
    for (int i = 1; i < m_group_numbers.count(); i++) {
        if (m_group_numbers[i] != i)
            return true;
    }
    return false;
}

int MVClusterDetailWidgetPrivate::get_current_view_index()
{
    int k = m_current_k;
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

void ClusterView::paint(QPainter* painter, QRectF rect)
{
    int xmargin = 1;
    int ymargin = 8;
    QRectF rect2(rect.x() + xmargin, rect.y() + ymargin, rect.width() - xmargin * 2, rect.height() - ymargin * 2);
    painter->setClipRect(rect);

    QColor background_color = d->m_colors["view_background"];
    if (m_highlighted)
        background_color = d->m_colors["view_background_highlighted"];
    else if (m_selected)
        background_color = d->m_colors["view_background_selected"];
    else if (m_hovered)
        background_color = d->m_colors["view_background_hovered"];
    painter->fillRect(rect, QColor(220, 220, 225));
    painter->fillRect(rect2, background_color);

    QPen pen_frame;
    pen_frame.setWidth(1);
    pen_frame.setColor(d->m_colors["frame1"]);
    if (m_selected)
        pen_frame.setColor(d->m_colors["view_frame_selected"]);
    painter->setPen(pen_frame);
    painter->drawRect(rect2);

    Mda template0 = m_CD.template0;
    Mda stdev0 = m_CD.stdev0;
    int M = template0.N1();
    int T = template0.N2();
    int Tmid = (int)((T + 1) / 2) - 1;
    m_T = T;

    int top_height = 20, bottom_height = 60;
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
        QColor col = d->m_channel_colors.value(m % d->m_channel_colors.count());
        QPen pen;
        pen.setWidth(1);
        pen.setColor(col);
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

    QString group_label = d->group_label_for_k(m_CD.k);
    if ((!group_label.isEmpty()) && (d->has_nontrivial_group_numbers())) {
        QPen pen;
        pen.setWidth(1);
        pen.setColor(d->m_colors["divider_line"]);
        pen.setStyle(Qt::DashLine);
        painter->setPen(pen);
        painter->drawLine(rect2.x(), rect2.y(), rect2.x(), rect2.y() + rect2.height());
    }
    {
        txt = QString("%1").arg(group_label);
        font.setPixelSize(16);
        if (compressed_info)
            font.setPixelSize(12);
        QPen pen;
        pen.setWidth(1);
        pen.setColor(Qt::darkBlue);
        painter->setFont(font);
        painter->setPen(pen);
        painter->drawText(m_top_rect, Qt::AlignCenter | Qt::AlignBottom, txt);
    }

    font.setPixelSize(11);
    int text_height = 13;

    if (!compressed_info) {
        RR = QRectF(m_bottom_rect.x(), m_bottom_rect.y() + m_bottom_rect.height() - text_height, m_bottom_rect.width(), text_height);
        txt = QString("%1 spikes").arg(m_CD.inds.count());
        QPen pen;
        pen.setWidth(1);
        pen.setColor(d->m_colors["info_text"]);
        painter->setFont(font);
        painter->setPen(pen);
        painter->drawText(RR, Qt::AlignCenter | Qt::AlignBottom, txt);
    }

    {
        QPen pen;
        pen.setWidth(1);
        RR = QRectF(m_bottom_rect.x(), m_bottom_rect.y() + m_bottom_rect.height() - text_height * 2, m_bottom_rect.width(), text_height);
        double rate = m_CD.inds.count() * 1.0 / d->m_total_time_sec;
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
        QString aa = m_attributes["assessment"].toString();
        pen.setColor(get_cluster_assessment_text_color(aa));
        txt = aa;
        painter->setFont(font);
        painter->setPen(pen);
        painter->drawText(RR, Qt::AlignCenter | Qt::AlignBottom, txt);
    }
}

double ClusterView::spaceNeeded()
{
    return 1;
}

void MVClusterDetailWidgetPrivate::do_paint(QPainter& painter, int W, int H)
{
    painter.fillRect(0, 0, W, H, m_colors["background"]);

    int right_margin=10; //make some room for the icon
    W=W-right_margin;

    if (m_calculator.isComputing()) {
        QFont font = painter.font();
        font.setPointSize(30);
        painter.setFont(font);
        painter.drawText(QRectF(0, 0, W, H), Qt::AlignCenter | Qt::AlignVCenter, "Calculating...");
        return;
    }

    QList<ClusterData> cluster_data_merged;
    QMap<int, QJsonObject> cluster_attributes;
    if (m_view_agent) {
        cluster_data_merged = merge_cluster_data(m_view_agent->clusterMerge(), m_cluster_data);
        cluster_attributes = m_view_agent->clusterAttributes();
    } else {
        cluster_data_merged = m_cluster_data;
    }

    qDeleteAll(m_views);
    m_views.clear();
    for (int i = 0; i < cluster_data_merged.count(); i++) {
        ClusterData CD = cluster_data_merged[i];
        ClusterView* V = new ClusterView(q, this);
        V->setStdevShading(m_stdev_shading);
        V->setHighlighted(CD.k == m_current_k);
        V->setSelected(m_selected_ks.contains(CD.k));
        V->setHovered(CD.k == m_hovered_k);
        V->setClusterData(CD);
        V->setAttributes(cluster_attributes[CD.k]);
        m_views << V;
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

    ChannelSpacingInfo csi = compute_channel_spacing_info(cluster_data_merged, m_vscale_factor);

    float x0_before_scaling = 0;
    for (int i = 0; i < m_views.count(); i++) {
        ClusterView* V = m_views[i];
        QRectF rect(x0_before_scaling * m_space_ratio - m_scroll_x, 0, V->spaceNeeded() * m_space_ratio, H);
        V->setChannelSpacingInfo(csi);
        V->paint(&painter, rect);
        V->x_position_before_scaling = x0_before_scaling;
        x0_before_scaling += V->spaceNeeded();
    }
}

void MVClusterDetailWidgetPrivate::export_image()
{
    QImage img = q->renderImage();
    user_save_image(img);
}

void MVClusterDetailWidgetPrivate::toggle_stdev_shading()
{
    m_stdev_shading = !m_stdev_shading;
    q->update();
}

void MVClusterDetailWidgetPrivate::start_calculation()
{
    m_calculator.stopComputation();
    //m_calculator.mscmdserver_url = m_mscmdserver_url;
    m_calculator.mlproxy_url = m_mlproxy_url;
    m_calculator.timeseries = m_timeseries;
    m_calculator.firings = m_firings;
    m_calculator.clip_size = m_clip_size;
    m_calculator.startComputation();
    q->update();
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
    for (int i = 0; i < group.count(); i++) {
        ret.inds << group[i].inds;
        ret.peaks << group[i].peaks;
        ret.times << group[i].times;
        double weight = ret.inds.count();
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
        for (long i = 0; i < ret.template0.totalSize(); i++) {
            double val_sum = sum0.get(i);
            double val_sumsqr = sumsqr0.get(i);
            ret.template0.set(val_sum / total_weight, i);
            ret.stdev0.set(sqrt((val_sumsqr - val_sum * val_sum / total_weight) / total_weight), i);
        }
    }
    return ret;
}

QList<ClusterData> MVClusterDetailWidgetPrivate::merge_cluster_data(const ClusterMerge& CM, const QList<ClusterData>& CD)
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
        } else {
            ClusterData CD0;
            CD0.k = CD[i].k;
            CD0.channel = CD[i].channel;
            ret << CD0;
        }
    }
    return ret;
}

QPointF ClusterView::template_coord2pix(int m, double t, double val)
{
    double pcty = m_csi.channel_locations.value(m) - val * m_csi.vert_scaling_factor; //negative because (0,0) is top-left, not bottom-right
    double pctx = 0;
    if (m_T)
        pctx = (t + 0.5) / m_T;
    int margx = 4;
    int margy = 4;
    float x0 = m_template_rect.x() + margx + pctx * (m_template_rect.width() - margx * 2);
    float y0 = m_template_rect.y() + margy + pcty * (m_template_rect.height() - margy * 2);
    return QPointF(x0, y0);
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

QColor ClusterView::get_cluster_assessment_text_color(QString aa)
{
    return Qt::black;
}

DiskReadMda mp_compute_templates(const QString& mlproxy_url, const QString& timeseries, const QString& firings, int clip_size, HaltAgent* halt_agent)
{
    TaskProgress task("mp_compute_templates");
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
    X.runProcess(halt_agent);
    task.log("Returning DiskReadMda: " + templates_fname);
    DiskReadMda ret(templates_fname);
    return ret;
}

void mp_compute_templates_stdevs(DiskReadMda& templates_out, DiskReadMda& stdevs_out, const QString& mlproxy_url, const QString& timeseries, const QString& firings, int clip_size, HaltAgent* halt_agent)
{
    TaskProgress task("mp_compute_templates_stdevs");
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
    X.runProcess(halt_agent);
    task.log("Returning DiskReadMda: " + templates_fname + " " + stdevs_fname);
    templates_out.setPath(templates_fname);
    stdevs_out.setPath(stdevs_fname);
}

void MVClusterDetailWidgetCalculator::compute()
{
    TaskProgress task("Cluster Details");

    QTime timer;
    timer.start();
    task.setProgress(0.1);

    int M = timeseries.N1();
    //int N = timeseries.N2();
    int L = firings.N2();
    int T = clip_size;

    QList<double> times;
    QList<int> channels, labels;
    QList<double> peaks;

    task.log("Setting up times/channels/labels/peaks");
    task.setProgress(0.2);
    for (int i = 0; i < L; i++) {
        times << firings.value(1, i) - 1; //convert to 0-based indexing
        channels << (int)firings.value(0, i) - 1; //convert to 0-based indexing
        labels << (int)firings.value(2, i);
        peaks << firings.value(3, i);
    }

    if (this->stopRequested()) {
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
    //DiskReadMda templates0 = mp_compute_templates(mlproxy_url, timeseries_path, firings_path, T, this);
    DiskReadMda templates0, stdevs0;
    mp_compute_templates_stdevs(templates0, stdevs0, mlproxy_url, timeseries_path, firings_path, T, this);
    if (this->stopRequested()) {
        task.error("Halted **");
        return;
    }

    task.setLabel("Setting cluster data");
    task.setProgress(0.75);
    for (int k = 1; k <= K; k++) {
        if (this->stopRequested()) {
            task.error("Halted ***");
            return;
        }
        ClusterData CD;
        CD.k = k;
        CD.channel = 0;
        for (int i = 0; i < L; i++) {
            if (labels[i] == k) {
                CD.inds << i;
                CD.times << times[i];
                CD.channel = channels[i];
                CD.peaks << peaks[i];
            }
        }
        if (this->stopRequested()) {
            task.error("Halted ****");
            return;
        }
        templates0.readChunk(CD.template0, 0, 0, k - 1, M, T, 1);
        stdevs0.readChunk(CD.stdev0, 0, 0, k - 1, M, T, 1);
        cluster_data << CD;
    }
}
