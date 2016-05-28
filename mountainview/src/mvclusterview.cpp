#include "mvclusterview.h"
#include <QImage>
#include <QColor>
#include <QPainter>
#include <stdio.h>
#include <QMouseEvent>
#include "affinetransformation.h"
#include <QTimer>
#include <math.h>
#include "mvutils.h"
#include "msmisc.h"
#include <QMenu>

QList<QColor> generate_colors(const QColor& bg, const QColor& fg, int noColors);

class MVClusterViewPrivate : public QObject {
    Q_OBJECT
public:
    MVClusterView* q;
    Mda m_data;
    Mda m_data_proj;
    bool m_data_proj_needed;
    Mda m_data_trans;
    bool m_data_trans_needed;
    int m_current_event_index;
    int m_mode;

    QImage m_grid_image;
    QRectF m_image_target;
    Mda m_heat_map_grid;
    Mda m_point_grid;
    Mda m_time_grid;
    Mda m_amplitude_grid;
    int m_grid_N1, m_grid_N2;
    bool m_grid_update_needed;
    QPointF m_anchor_point;
    QPointF m_last_mouse_release_point;
    AffineTransformation m_anchor_transformation;
    bool m_moved_from_anchor;
    QList<double> m_times;
    double m_max_time;
    QList<double> m_amplitudes; //absolute amplitudes, actually
    double m_max_amplitude;
    QList<int> m_labels;
    QSet<int> m_closest_inds_to_exclude;
    QList<QColor> m_label_colors;
    bool m_emit_transformation_changed_scheduled;

    Mda proj_matrix; //3xnum_features
    AffineTransformation m_transformation; //3x4

    void compute_data_proj();
    void compute_data_trans();
    void update_grid();
    void coord2gridindex(double x0, double y0, double& i1, double& i2);
    QPointF pixel2coord(QPointF pix);
    QPointF coord2pixel(QPointF coord);
    QColor get_label_color(int label);
    QColor get_time_color(double pct);
    int find_closest_event_index(double x, double y, const QSet<int>& inds_to_exclude);
    void set_current_event_index(int ind, bool do_emit = true);
    void schedule_emit_transformation_changed();
    void do_paint(QPainter& painter, int W, int H);
    void export_image();
public
slots:
    void slot_emit_transformation_changed();
};
#include "mvclusterview.moc"

MVClusterView::MVClusterView(QWidget* parent)
    : QWidget(parent)
{
    d = new MVClusterViewPrivate;
    d->q = this;
    d->m_grid_N1 = d->m_grid_N2 = 300;
    d->m_grid_update_needed = false;
    d->m_data_proj_needed = true;
    d->m_data_trans_needed = true;
    d->m_anchor_point = QPointF(-1, -1);
    d->m_last_mouse_release_point = QPointF(-1, -1);
    d->m_transformation.setIdentity();
    d->m_moved_from_anchor = false;
    d->m_current_event_index = -1;
    d->m_mode = MVCV_MODE_LABEL_COLORS;
    d->m_emit_transformation_changed_scheduled = false;
    d->m_max_time = 1;
    d->m_max_amplitude = 1;
    this->setMouseTracking(true);

    //QList<QString> color_strings;
    //color_strings << "black" << "blue" << "red" << "green" << "white" << "magenta";
    /*
    color_strings << "#F7977A"
                  << "#FDC68A"
                  << "#C4DF9B"
                  << "#82CA9D"
                  << "#6ECFF6"
                  << "#8493CA"
                  << "#A187BE"
                  << "#F49AC2"
                  << "#F9AD81"
                  << "#FFF79A"
                  << "#A2D39C"
                  << "#7BCDC8"
                  << "#7EA7D8"
                  << "#8882BE"
                  << "#BC8DBF"
                  << "#F6989D";
                  */
    //for (int i = 0; i < color_strings.size(); i++) {
    //    d->m_label_colors << QColor(color_strings[i]);
    //}

    //these should be relatively prime
    //int num1=5;
    //int num2=22;
    int num1 = 3;
    int num2 = 10;
    QList<QColor> colors = generate_colors(Qt::gray, Qt::white, num2);
    for (int j = 0; j < colors.count(); j++) {
        d->m_label_colors << colors.value((j * num1) % num2);
    }

    this->setContextMenuPolicy(Qt::CustomContextMenu);
    connect(this, SIGNAL(customContextMenuRequested(QPoint)), this, SLOT(slot_context_menu(QPoint)));
}

MVClusterView::~MVClusterView()
{
    delete d;
}

void MVClusterView::setData(const Mda& X)
{
    d->m_data = X;
    d->m_data_proj_needed = true;
    d->m_data_trans_needed = true;
    d->m_grid_update_needed = true;

    update();
}

bool MVClusterView::hasData()
{
    return (d->m_data.totalSize() > 1);
}

void MVClusterView::setTimes(const QList<double>& times)
{
    d->m_times = times;
    d->m_max_time = compute_max(times);
}

void MVClusterView::setLabels(const QList<int>& labels)
{
    d->m_labels = labels;
    if (d->m_times.count() != labels.count()) {
        d->m_times.clear();
        for (int i = 0; i < labels.count(); i++) {
            d->m_times << i;
        }
    }
}

void MVClusterView::setAmplitudes(const QList<double>& amps)
{
    QList<double> tmp;
    for (int i = 0; i < amps.count(); i++)
        tmp << fabs(amps[i]);
    d->m_amplitudes = tmp;
    d->m_max_amplitude = compute_max(tmp);
}

void MVClusterView::setMode(int mode)
{
    d->m_mode = mode;
    update();
}

void MVClusterView::setCurrentEvent(MVEvent evt, bool do_emit)
{
    if ((d->m_times.isEmpty()) || (d->m_labels.isEmpty()))
        return;
    for (int i = 0; i < d->m_times.count(); i++) {
        if ((d->m_times[i] == evt.time) && (d->m_labels.value(i) == evt.label)) {
            d->set_current_event_index(i, do_emit);
            return;
        }
    }
    d->set_current_event_index(-1);
}

MVEvent MVClusterView::currentEvent()
{
    MVEvent ret;
    if (d->m_current_event_index < 0) {
        ret.time = -1;
        ret.label = -1;
        return ret;
    }
    ret.time = d->m_times.value(d->m_current_event_index);
    ret.label = d->m_labels.value(d->m_current_event_index);
    return ret;
}

int MVClusterView::currentEventIndex()
{
    return d->m_current_event_index;
}

AffineTransformation MVClusterView::transformation()
{
    return d->m_transformation;
}

void MVClusterView::setTransformation(const AffineTransformation& T)
{
    if (d->m_transformation.equals(T))
        return;
    d->m_transformation = T;
    d->m_data_trans_needed = true;
    d->m_grid_update_needed = true;
    update();
    //do not emit to avoid excessive signals
}

QImage MVClusterView::renderImage(int W, int H)
{
    QImage ret = QImage(W, H, QImage::Format_RGB32);
    QPainter painter(&ret);

    int current_event_index = d->m_current_event_index;
    d->m_current_event_index = -1;

    d->do_paint(painter, W, H);

    d->m_current_event_index = current_event_index;

    return ret;
}

void MVClusterView::slot_context_menu(const QPoint& pos)
{
    QMenu M;
    QAction* export_image = M.addAction("Export Image");
    QAction* selected = M.exec(this->mapToGlobal(pos));
    if (selected == export_image) {
        d->export_image();
    }
}

QRectF compute_centered_square(QRectF R)
{
    int margin = 15;
    int W0 = R.width() - margin * 2;
    int H0 = R.height() - margin * 2;
    if (W0 > H0) {
        return QRectF(margin + (W0 - H0) / 2, margin, H0, H0);
    } else {
        return QRectF(margin, margin + (H0 - W0) / 2, W0, W0);
    }
}

void MVClusterView::paintEvent(QPaintEvent* evt)
{
    Q_UNUSED(evt)
    QPainter painter(this);
    d->do_paint(painter, width(), height());
}

void MVClusterView::mouseMoveEvent(QMouseEvent* evt)
{
    QPointF pt = evt->pos();
    if (d->m_anchor_point.x() >= 0) {
        //We have the mouse button down!
        QPointF diff = pt - d->m_anchor_point;
        if ((qAbs(diff.x()) >= 5) || (qAbs(diff.y()) >= 5) || (d->m_moved_from_anchor)) {
            d->m_moved_from_anchor = true;
            double factor = 1.0 / 2;
            double deg_x = -diff.x() * factor;
            double deg_y = diff.y() * factor;
            d->m_transformation = d->m_anchor_transformation;
            d->m_transformation.rotateY(deg_x * M_PI / 180);
            d->m_transformation.rotateX(deg_y * M_PI / 180);
            d->schedule_emit_transformation_changed();
            d->m_data_trans_needed = true;
            d->m_grid_update_needed = true;
            update();
        }
    }
}

void MVClusterView::mousePressEvent(QMouseEvent* evt)
{
    if (evt->button() == Qt::LeftButton) {
        QPointF pt = evt->pos();
        d->m_anchor_point = pt;
        d->m_anchor_transformation = d->m_transformation;
        d->m_moved_from_anchor = false;
    }
}

void MVClusterView::mouseReleaseEvent(QMouseEvent* evt)
{
    Q_UNUSED(evt)
    if (evt->button() == Qt::LeftButton) {
        d->m_anchor_point = QPointF(-1, -1);
        if (evt->pos() != d->m_last_mouse_release_point)
            d->m_closest_inds_to_exclude.clear();
        if (!d->m_moved_from_anchor) {
            QPointF coord = d->pixel2coord(evt->pos());
            int ind = d->find_closest_event_index(coord.x(), coord.y(), d->m_closest_inds_to_exclude);
            if (ind >= 0)
                d->m_closest_inds_to_exclude.insert(ind);
            d->set_current_event_index(ind);
        }
        d->m_last_mouse_release_point = evt->pos();
    }
}

void MVClusterView::wheelEvent(QWheelEvent* evt)
{
    double delta = evt->delta();
    double factor = 1;
    if (delta > 0) {
        factor = 1.1;
    } else if (delta < 0) {
        factor = 1 / 1.1;
    }
    if (delta != 1) {
        d->m_transformation.scale(factor, factor, factor);
        d->m_data_trans_needed = true;
        d->m_grid_update_needed = true;
        d->schedule_emit_transformation_changed();
        update();
    }
}

void MVClusterViewPrivate::slot_emit_transformation_changed()
{
    m_emit_transformation_changed_scheduled = false;
    emit q->transformationChanged();
}

void MVClusterViewPrivate::compute_data_proj()
{
    if (m_data.N2() <= 1)
        return;
    m_data_proj.allocate(3, m_data.N2());
    for (int i = 0; i < m_data.N2(); i++) {
        m_data_proj.setValue(m_data.value(0, i), 0, i);
        m_data_proj.setValue(m_data.value(1, i), 1, i);
        m_data_proj.setValue(m_data.value(2, i), 2, i);
    }
}

void MVClusterViewPrivate::compute_data_trans()
{
    int N2 = m_data_proj.N2();
    m_data_trans.allocate(3, N2);
    if (m_data_proj.N1() != 3)
        return; //prevents a crash when data hasn't been set, which is always a good thing
    double* AA = m_data_proj.dataPtr();
    double* BB = m_data_trans.dataPtr();
    double MM[16];
    m_transformation.getMatrixData(MM);
    int aaa = 0;
    for (int i = 0; i < N2; i++) {
        BB[aaa + 0] = AA[aaa + 0] * MM[0] + AA[aaa + 1] * MM[1] + AA[aaa + 2] * MM[2] + MM[3];
        BB[aaa + 1] = AA[aaa + 0] * MM[4] + AA[aaa + 1] * MM[5] + AA[aaa + 2] * MM[6] + MM[7];
        BB[aaa + 2] = AA[aaa + 0] * MM[8] + AA[aaa + 1] * MM[9] + AA[aaa + 2] * MM[10] + MM[11];
        aaa += 3;
    }
}

double compute_max(int N, double* X)
{
    if (N == 0)
        return 0;
    double ret = X[0];
    for (int i = 0; i < N; i++) {
        if (X[i] > ret)
            ret = X[i];
    }
    return ret;
}

QColor make_color(double r, double g, double b)
{
    if (r < 0)
        r = 0;
    if (r > 1)
        r = 1;
    if (g < 0)
        g = 0;
    if (g > 1)
        g = 1;
    if (b < 0)
        b = 0;
    if (b > 1)
        b = 1;
    return QColor((int)(r * 255), (int)(g * 255), (int)(b * 255));
}

void MVClusterViewPrivate::update_grid()
{
    if (m_data_proj_needed) {
        compute_data_proj();
        m_data_proj_needed = false;
    }
    if (m_data_trans_needed) {
        compute_data_trans();
        m_data_trans_needed = false;
    }
    int kernel_rad = 10;
    double kernel_tau = 3;
    double kernel[(kernel_rad * 2 + 1) * (kernel_rad * 2 + 1)];
    {
        int aa = 0;
        for (int dy = -kernel_rad; dy <= kernel_rad; dy++) {
            for (int dx = -kernel_rad; dx <= kernel_rad; dx++) {
                kernel[aa] = exp(-0.5 * (dx * dx + dy * dy) / (kernel_tau * kernel_tau));
                aa++;
            }
        }
    }
    int N1 = m_grid_N1, N2 = m_grid_N2;

    m_point_grid.allocate(N1, N2);
    for (int i = 0; i < N1 * N2; i++)
        m_point_grid.set(-1, i);
    double* m_point_grid_ptr = m_point_grid.dataPtr();

    if (m_mode == MVCV_MODE_TIME_COLORS) {
        m_time_grid.allocate(N1, N2);
        for (int i = 0; i < N1 * N2; i++)
            m_time_grid.set(-1, i);
    }
    double* m_time_grid_ptr = m_time_grid.dataPtr();

    if (m_mode == MVCV_MODE_AMPLITUDE_COLORS) {
        m_amplitude_grid.allocate(N1, N2);
        for (int i = 0; i < N1 * N2; i++)
            m_amplitude_grid.set(0, i);
    }
    double* m_amplitude_grid_ptr = m_amplitude_grid.dataPtr();

    if (m_mode == MVCV_MODE_HEAT_DENSITY) {
        m_heat_map_grid.allocate(N1, N2);
    }
    double* m_heat_map_grid_ptr = m_heat_map_grid.dataPtr();

    Mda z_grid;
    if (m_mode != MVCV_MODE_HEAT_DENSITY) {
        z_grid.allocate(N1, N2);
    }
    double* z_grid_ptr = z_grid.dataPtr();

    double max_abs_val = 0;
    int NN = m_data.totalSize();
    for (int i = 0; i < NN; i++) {
        if (fabs(m_data.get(i)) > max_abs_val)
            max_abs_val = fabs(m_data.get(i));
    }

    QList<double> x0s, y0s, z0s;
    QList<int> label0s;
    QList<double> time0s;
    QList<double> amp0s;
    for (int j = 0; j < m_data_trans.N2(); j++) {
        double x0 = m_data_trans.value(0, j);
        double y0 = m_data_trans.value(1, j);
        double z0 = m_data_trans.value(2, j);
        x0s << x0;
        y0s << y0;
        z0s << z0;
        label0s << m_labels.value(j);
        time0s << m_times.value(j);
        amp0s << m_amplitudes.value(j);
    }

    if (m_mode != MVCV_MODE_HEAT_DENSITY) {
        //3 axes
        double factor = 1.2;
        if (max_abs_val) {
            for (double aa = -max_abs_val * factor; aa <= max_abs_val * factor; aa += max_abs_val * factor / 50) {
                {
                    Point3D pt1 = Point3D(aa, 0, 0);
                    Point3D pt2 = m_transformation.map(pt1);
                    x0s << pt2.x;
                    y0s << pt2.y;
                    z0s << pt2.z;
                    label0s << -2;
                    time0s << -1;
                    amp0s << 0;
                }
                {
                    Point3D pt1 = Point3D(0, aa, 0);
                    Point3D pt2 = m_transformation.map(pt1);
                    x0s << pt2.x;
                    y0s << pt2.y;
                    z0s << pt2.z;
                    label0s << -2;
                    time0s << -1;
                    amp0s << 0;
                }
                {
                    Point3D pt1 = Point3D(0, 0, aa);
                    Point3D pt2 = m_transformation.map(pt1);
                    x0s << pt2.x;
                    y0s << pt2.y;
                    z0s << pt2.z;
                    label0s << -2;
                    time0s << -1;
                    amp0s << 0;
                }
            }
        }
    }
    for (int i = 0; i < x0s.count(); i++) {
        double x0 = x0s[i], y0 = y0s[i], z0 = z0s[i];
        int label0 = label0s[i];
        double time0 = time0s[i];
        double amp0 = amp0s[i];
        double factor = 1;
        //if (m_data_trans.value(2,j)>0) {
        //factor=0.2;
        //}
        double i1, i2;
        coord2gridindex(x0, y0, i1, i2);
        int ii1 = (int)(i1 + 0.5);
        int ii2 = (int)(i2 + 0.5);
        if ((ii1 - kernel_rad >= 0) && (ii1 + kernel_rad < N1) && (ii2 - kernel_rad >= 0) && (ii2 + kernel_rad < N2)) {
            int iiii = ii1 + N1 * ii2;
            if (m_mode != MVCV_MODE_HEAT_DENSITY) {
                if ((m_point_grid_ptr[iiii] == -1) || (z_grid_ptr[iiii] > z0)) {
                    m_point_grid_ptr[iiii] = label0;
                    if (m_mode == MVCV_MODE_TIME_COLORS) {
                        m_time_grid_ptr[iiii] = time0;
                    }
                    if (m_mode == MVCV_MODE_AMPLITUDE_COLORS) {
                        m_amplitude_grid_ptr[iiii] = amp0;
                    }
                    z_grid_ptr[iiii] = z0;
                }
            } else {
                m_point_grid_ptr[iiii] = 1;
                if (m_mode == MVCV_MODE_TIME_COLORS) {
                    m_time_grid_ptr[iiii] = time0;
                }
                if (m_mode == MVCV_MODE_AMPLITUDE_COLORS) {
                    m_amplitude_grid_ptr[iiii] = amp0;
                }
            }

            if (m_mode == MVCV_MODE_HEAT_DENSITY) {
                int aa = 0;
                for (int dy = -kernel_rad; dy <= kernel_rad; dy++) {
                    int ii2b = ii2 + dy;
                    for (int dx = -kernel_rad; dx <= kernel_rad; dx++) {
                        int ii1b = ii1 + dx;
                        m_heat_map_grid_ptr[ii1b + N1 * ii2b] += kernel[aa] * factor;
                        aa++;
                    }
                }
            }
        }
    }

    m_grid_image = QImage(N1, N2, QImage::Format_ARGB32);

    QColor dark(50, 50, 50);
    QColor axes_color(150, 150, 150);

    m_grid_image.fill(dark);

    if (m_mode == MVCV_MODE_HEAT_DENSITY) {
        //3 Axes
        QPainter paintr(&m_grid_image);
        paintr.setRenderHint(QPainter::Antialiasing);
        paintr.setPen(QPen(axes_color));
        double factor = 1.2;
        for (int pass = 1; pass <= 3; pass++) {
            double a1 = 0, a2 = 0, a3 = 0;
            if (pass == 1)
                a1 = max_abs_val * factor;
            if (pass == 2)
                a2 = max_abs_val * factor;
            if (pass == 3)
                a3 = max_abs_val * factor;
            Point3D pt1 = Point3D(-a1, -a2, -a3);
            Point3D pt2 = Point3D(a1, a2, a3);
            Point3D pt1b = m_transformation.map(pt1);
            Point3D pt2b = m_transformation.map(pt2);
            double x1, y1, x2, y2;
            coord2gridindex(pt1b.x, pt1b.y, x1, y1);
            coord2gridindex(pt2b.x, pt2b.y, x2, y2);
            paintr.drawLine(QPointF(x1, y1), QPointF(x2, y2));
        }
    }

    if (m_mode == MVCV_MODE_HEAT_DENSITY) {
        double max_heat_map_grid_val = compute_max(N1 * N2, m_heat_map_grid.dataPtr());
        for (int i2 = 0; i2 < N2; i2++) {
            for (int i1 = 0; i1 < N1; i1++) {
                double val = m_point_grid.value(i1, i2);
                if (val > 0) {
                    double val2 = m_heat_map_grid.value(i1, i2) / max_heat_map_grid_val;
                    QColor CC = get_heat_map_color(val2);
                    //m_grid_image.setPixel(i1,i2,white.rgb());
                    m_grid_image.setPixel(i1, i2, CC.rgb());
                }
            }
        }
    } else {
        for (int i2 = 0; i2 < N2; i2++) {
            for (int i1 = 0; i1 < N1; i1++) {
                double val = m_point_grid.value(i1, i2);
                if (m_mode == MVCV_MODE_TIME_COLORS) {
                    if (val >= 0) {
                        double time0 = m_time_grid.value(i1, i2);
                        if (m_max_time) {
                            QColor CC = get_time_color(time0 / m_max_time);
                            m_grid_image.setPixel(i1, i2, CC.rgb());
                        }
                    } else if (val == -2) {
                        m_grid_image.setPixel(i1, i2, axes_color.rgb());
                    }
                } else if (m_mode == MVCV_MODE_AMPLITUDE_COLORS) {
                    if (val >= 0) {
                        double amp0 = m_amplitude_grid.value(i1, i2);
                        if (m_max_amplitude) {
                            QColor CC = get_time_color(amp0 / m_max_amplitude);
                            m_grid_image.setPixel(i1, i2, CC.rgb());
                        }
                    } else if (val == -2) {
                        m_grid_image.setPixel(i1, i2, axes_color.rgb());
                    }
                } else if (m_mode == MVCV_MODE_LABEL_COLORS) {
                    if (val > 0) {
                        QColor CC = get_label_color((int)val);
                        m_grid_image.setPixel(i1, i2, CC.rgb());
                    } else if (val == 0) {
                        QColor CC = Qt::white;
                        m_grid_image.setPixel(i1, i2, CC.rgb());
                    } else if (val == -2) {
                        m_grid_image.setPixel(i1, i2, axes_color.rgb());
                    }
                }
            }
        }
    }
}

void MVClusterViewPrivate::coord2gridindex(double x0, double y0, double& i1, double& i2)
{
    int N1 = m_grid_N1;
    int N2 = m_grid_N2;
    int N1mid = (N1 + 1) / 2 - 1;
    int N2mid = (N2 + 1) / 2 - 1;
    double delta1 = 2.0 / N1;
    double delta2 = 2.0 / N2;
    i1 = N1mid + x0 / delta1;
    i2 = N2mid + y0 / delta2;
}

QPointF MVClusterViewPrivate::pixel2coord(QPointF pt)
{
    int N1 = m_grid_N1;
    int N2 = m_grid_N2;
    double delta1 = 2.0 / N1;
    double delta2 = 2.0 / N2;
    int N1mid = (N1 + 1) / 2 - 1;
    int N2mid = (N2 + 1) / 2 - 1;
    double pctx = (pt.x() - m_image_target.x()) / (m_image_target.width());
    double pcty = (pt.y() - m_image_target.y()) / (m_image_target.height());
    double xx = (-N1mid + pctx * N1) * delta1;
    double yy = (-N2mid + pcty * N2) * delta2;
    return QPointF(xx, yy);
}

QPointF MVClusterViewPrivate::coord2pixel(QPointF coord)
{
    int N1 = m_grid_N1;
    int N2 = m_grid_N2;
    double delta1 = 2.0 / N1;
    double delta2 = 2.0 / N2;
    int N1mid = (N1 + 1) / 2 - 1;
    int N2mid = (N2 + 1) / 2 - 1;
    double xx = coord.x();
    double yy = coord.y();
    double pctx = (xx / delta1 + N1mid) / N1;
    double pcty = (yy / delta2 + N2mid) / N2;
    double pt_x = pctx * m_image_target.width() + m_image_target.x();
    double pt_y = pcty * m_image_target.height() + m_image_target.y();
    return QPointF(pt_x, pt_y);
}

QColor MVClusterViewPrivate::get_label_color(int label)
{
    return m_label_colors[label % m_label_colors.size()];
}

QColor MVClusterViewPrivate::get_time_color(double pct)
{
    if (pct < 0)
        pct = 0;
    if (pct > 1)
        pct = 1;
    int h = (int)(pct * 359);
    return QColor::fromHsv(h, 200, 255);
}

int MVClusterViewPrivate::find_closest_event_index(double x, double y, const QSet<int>& inds_to_exclude)
{
    double best_distsqr = 0;
    int best_ind = 0;
    for (int i = 0; i < m_data_trans.N2(); i++) {
        if (!inds_to_exclude.contains(i)) {
            double distx = m_data_trans.value(0, i) - x;
            double disty = m_data_trans.value(1, i) - y;
            double distsqr = distx * distx + disty * disty;
            if ((distsqr < best_distsqr) || (i == 0)) {
                best_distsqr = distsqr;
                best_ind = i;
            }
        }
    }
    return best_ind;
}

void MVClusterViewPrivate::set_current_event_index(int ind, bool do_emit)
{
    if (m_current_event_index == ind)
        return;
    m_current_event_index = ind;
    if (do_emit) {
        emit q->currentEventChanged();
    }
    q->update();
}

void MVClusterViewPrivate::schedule_emit_transformation_changed()
{
    if (m_emit_transformation_changed_scheduled)
        return;
    m_emit_transformation_changed_scheduled = true;
    QTimer::singleShot(100, this, SLOT(slot_emit_transformation_changed()));
}

void MVClusterViewPrivate::do_paint(QPainter& painter, int W, int H)
{
    if (m_grid_update_needed) {
        update_grid();
        m_grid_update_needed = false;
    }

    painter.fillRect(0, 0, W, H, QColor(40, 40, 40));
    QRectF target = compute_centered_square(QRectF(0, 0, W, H));
    painter.drawImage(target, m_grid_image.scaled(W, H, Qt::IgnoreAspectRatio, Qt::SmoothTransformation));
    m_image_target = target;
    //QPen pen; pen.setColor(Qt::yellow);
    //painter.setPen(pen);
    //painter.drawRect(target);

    if (m_current_event_index >= 0) {
        double x = m_data_trans.value(0, m_current_event_index);
        double y = m_data_trans.value(1, m_current_event_index);
        QPointF pix = coord2pixel(QPointF(x, y));
        painter.setBrush(QBrush(Qt::darkGreen));
        painter.drawEllipse(pix, 6, 6);
    }
}

void MVClusterViewPrivate::export_image()
{
    QImage img = q->renderImage(1000, 800);
    user_save_image(img);
}

// generate_colors() is adapted from code by...
/*
 * Copyright (c) 2008 Helder Correia
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
*/

QList<QColor> generate_colors(const QColor& bg, const QColor& fg, int noColors)
{
    QList<QColor> colors;
    const int HUE_BASE = (bg.hue() == -1) ? 90 : bg.hue();
    int h, s, v;

    for (int i = 0; i < noColors; i++) {
        h = int(HUE_BASE + (360.0 / noColors * i)) % 360;
        s = 240;
        v = int(qMax(bg.value(), fg.value()) * 0.85);

        // take care of corner cases
        const int M = 35;
        if ((h < bg.hue() + M && h > bg.hue() - M)
            || (h < fg.hue() + M && h > fg.hue() - M)) {
            h = ((bg.hue() + fg.hue()) / (i + 1)) % 360;
            s = ((bg.saturation() + fg.saturation() + 2 * i) / 2) % 256;
            v = ((bg.value() + fg.value() + 2 * i) / 2) % 256;
        }

        colors.append(QColor::fromHsv(h, s, v));
    }

    return colors;
}
