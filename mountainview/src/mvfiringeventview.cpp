#include "mvfiringeventview.h"

#include <QPainter>
#include <QDebug>
#include <math.h>
#include <QHBoxLayout>
#include <QTimer>
#include <QMouseEvent>
#include <QMenu>
#include "sstimeseriesview.h"
#include "mvutils.h"
#include "msmisc.h"

class MVFiringEventViewPrivate {
public:
    MVFiringEventView* q;

    //set by user
    QList<double> m_times;
    QList<int> m_labels;
    QList<double> m_amplitudes;
    double m_samplerate;
    QList<Epoch> m_epochs;

    //internal
    double m_max_timepoint;
    double m_min_amplitude, m_max_amplitude;
    bool m_update_scheduled;
    QImage m_image;
    QRectF m_target_rect;
    Mda m_event_index_grid;
    int m_current_event_index;

    //settings
    double m_hmargin, m_vmargin;

    QPointF coord2imagepix(const QPointF& p, int W, int H);
    QPointF windowpix2coord(const QPointF& p);
    QPointF imagepix2coord(const QPointF& p, int W, int H);
    void schedule_update();
    int find_closest_event_index(const QPointF& pt);
    int find_closest_event_index(int i1, int i2);
    void set_current_event_index(int ind);
    void export_image();
    void do_paint(QPainter& painter, int W, int H);
};

MVFiringEventView::MVFiringEventView()
{
    d = new MVFiringEventViewPrivate;
    d->q = this;
    d->m_samplerate = 30000; //Hz
    d->m_max_timepoint = 0;
    d->m_min_amplitude = 0;
    d->m_max_amplitude = 0;
    d->m_update_scheduled = false;
    d->m_current_event_index = -1;

    d->m_hmargin = 25;
    d->m_vmargin = 25;

    this->setContextMenuPolicy(Qt::CustomContextMenu);
    connect(this, SIGNAL(customContextMenuRequested(QPoint)), this, SLOT(slot_context_menu(QPoint)));
}

MVFiringEventView::~MVFiringEventView()
{
    delete d;
}

void MVFiringEventView::setFirings(const Mda& firings)
{
    d->m_times.clear();
    d->m_labels.clear();
    d->m_amplitudes.clear();
    for (int i = 0; i < firings.N2(); i++) {
        d->m_times << firings.value(1, i);
        d->m_labels << (int)firings.value(2, i);
        d->m_amplitudes << firings.value(3, i);
    }

    d->m_max_timepoint = compute_max(d->m_times);
    d->m_min_amplitude = compute_min(d->m_amplitudes);
    d->m_max_amplitude = compute_max(d->m_amplitudes);

    d->schedule_update();
}

void MVFiringEventView::setSampleRate(double ff)
{
    d->m_samplerate = ff;
    d->schedule_update();
}

void MVFiringEventView::setCurrentEvent(MVEvent evt)
{
    Q_UNUSED(evt)
    //finish this!
}

void MVFiringEventView::setEpochs(const QList<Epoch>& epochs)
{
    d->m_epochs = epochs;
    update();
}

QImage MVFiringEventView::renderImage(int W, int H)
{
    QImage ret = QImage(W, H, QImage::Format_RGB32);
    QPainter painter(&ret);
    d->do_paint(painter, W, H);
    return ret;
}

void apply_kernel(int N, double* X, int kernel_rad, double* kernel)
{
    double Y[N];
    for (int i = 0; i < N; i++) {
        double val = 0;
        for (int dd = -kernel_rad; dd <= kernel_rad; dd++) {
            int i2 = i + dd;
            if ((0 <= i2) && (i2 < N)) {
                val += X[i2] * kernel[kernel_rad + dd];
            }
        }
        Y[i] = val;
    }
    for (int i = 0; i < N; i++)
        X[i] = Y[i];
}

void smooth_grid(Mda& X, double kernel_tau)
{
    int kernel_rad = (int)kernel_tau * 3;
    double kernel[(kernel_rad * 2 + 1)];
    for (int dd = -kernel_rad; dd <= kernel_rad; dd++) {
        kernel[dd + kernel_rad] = exp(-0.5 * (dd * dd) / (kernel_tau * kernel_tau));
    }
    double* ptr = X.dataPtr();
    double bufx[X.N1()];
    int aa = 0;
    for (int y = 0; y < X.N2(); y++) {
        for (int x = 0; x < X.N1(); x++) {
            bufx[x] = ptr[aa];
            aa++;
        }
        aa -= X.N1();
        apply_kernel(X.N1(), bufx, kernel_rad, kernel);
        for (int x = 0; x < X.N1(); x++) {
            ptr[aa] = bufx[x];
            aa++;
        }
    }
    double bufy[X.N2()];
    for (int x = 0; x < X.N1(); x++) {
        int aa = x;
        for (int y = 0; y < X.N2(); y++) {
            bufy[y] = ptr[aa];
            aa += X.N1();
        }
        apply_kernel(X.N2(), bufy, kernel_rad, kernel);
        aa = x;
        for (int y = 0; y < X.N2(); y++) {
            ptr[aa] = bufy[y];
            aa += X.N1();
        }
    }
}

void MVFiringEventView::slot_update()
{
    d->m_update_scheduled = false;

    int W = 500;
    int H = (int)(this->height() * 1.0 / this->width() * W);

    d->m_event_index_grid.allocate(W, H);
    double* ptr_eig = d->m_event_index_grid.dataPtr();
    for (int i = 0; i < W * H; i++)
        ptr_eig[i] = -1;
    Mda grid;
    grid.allocate(W, H);
    double* ptr = grid.dataPtr();
    for (int i = 0; i < d->m_times.count(); i++) {
        double t0 = d->m_times[i];
        double a0 = d->m_amplitudes.value(i);
        QPointF pt = d->coord2imagepix(QPointF(t0, a0), W, H);
        int x = (int)pt.x();
        int y = (int)pt.y();
        if ((x >= 0) && (x < W) && (y >= 0) && (y < H)) {
            ptr[x + W * y]++;
            ptr_eig[x + W * y] = i;
        }
    }

    Mda smoothed_grid = grid;
    smooth_grid(smoothed_grid, 4);
    double max_density = compute_max(smoothed_grid.totalSize(), smoothed_grid.dataPtr());
    double* smoothed_ptr = smoothed_grid.dataPtr();

    d->m_image = QImage(W, H, QImage::Format_ARGB32);
    d->m_image.fill(QColor(0, 0, 0, 0).toRgb());

    QColor white(255, 255, 255);
    QColor red(255, 0, 0);
    for (int y = 0; y < H; y++) {
        for (int x = 0; x < W; x++) {
            double val = ptr[x + W * y];
            if (val) {
                double pct = smoothed_ptr[x + W * y] / max_density;
                pct = sqrt(pct);
                QColor CC = get_heat_map_color(pct);
                d->m_image.setPixel(x, y, CC.rgb());
            }
        }
    }

    this->update();
    /*d->m_update_scheduled=false;
    d->do_compute();
    int Kv=d->m_event_counts.N1();
    long num_bins=d->m_event_counts.N2();
    Mda data; data.allocate(Kv,num_bins);
    for (int k=0; k<Kv; k++) {
        QList<double> rates;
        for (int i=0; i<num_bins; i++) {
            rates << d->m_event_counts.value(k,i)/d->m_window_width_sec;
        }
        QList<double> rates_smoothed=do_smoothing(rates,d->m_smoothing_kernel_size/d->m_window_width_sec);
        for (int i=0; i<rates_smoothed.count(); i++) {
            double val=0;
            if (rates_smoothed[i]) val=log(rates_smoothed[i]);
            data.setValue(val,k,i);
        }
    }
    DiskArrayModel_New *DAM=new DiskArrayModel_New();
    DAM->setFromMda(data);
    d->m_view->setData(DAM,this);*/
}

void MVFiringEventView::slot_context_menu(const QPoint& pos)
{
    QMenu M;
    QAction* export_image = M.addAction("Export Image");
    QAction* selected = M.exec(this->mapToGlobal(pos));
    if (selected == export_image) {
        d->export_image();
    }
}

void MVFiringEventView::paintEvent(QPaintEvent* evt)
{
    Q_UNUSED(evt);

    QPainter painter(this);
    d->do_paint(painter, width(), height());
}

void MVFiringEventView::mouseReleaseEvent(QMouseEvent* evt)
{
    QPointF pt = evt->pos();
    int index = d->find_closest_event_index(pt);
}

QPointF MVFiringEventViewPrivate::coord2imagepix(const QPointF& p, int W, int H)
{
    if ((m_max_timepoint == 0) || (m_max_amplitude == m_min_amplitude))
        return QPointF(0, 0);
    double pctx = p.x() / m_max_timepoint;
    double pcty = (p.y() - m_min_amplitude) / (m_max_amplitude - m_min_amplitude);

    double x0 = pctx * W;
    double y0 = (1 - pcty) * H;

    return QPointF(x0, y0);
}

QPointF MVFiringEventViewPrivate::windowpix2coord(const QPointF& p)
{
    QPointF q = p - m_target_rect.topLeft();
    return imagepix2coord(q, m_target_rect.width(), m_target_rect.height());
}

QPointF MVFiringEventViewPrivate::imagepix2coord(const QPointF& p, int W, int H)
{
    double pctx = p.x() / W;
    double pcty = 1 - p.y() / H;
    return QPointF(pctx * m_max_timepoint, pcty * (m_max_amplitude - m_min_amplitude) + m_min_amplitude);
}

void MVFiringEventViewPrivate::schedule_update()
{
    if (m_update_scheduled)
        return;
    m_update_scheduled = true;
    QTimer::singleShot(500, q, SLOT(slot_update()));
}

int MVFiringEventViewPrivate::find_closest_event_index(const QPointF& pt)
{
    QPointF coord = windowpix2coord(pt);
    QPointF qq = coord2imagepix(coord, m_event_index_grid.N1(), m_event_index_grid.N2());
    return find_closest_event_index((int)(qq.x()), (int)(qq.y()));
}

int MVFiringEventViewPrivate::find_closest_event_index(int i1, int i2)
{
    int best_ind = -1;
    int best_x = -1, best_y = -1;
    double best_dist = 0;
    for (int y = 0; y < m_event_index_grid.N2(); y++) {
        for (int x = 0; x < m_event_index_grid.N1(); x++) {
            int val = (int)m_event_index_grid.value(x, y);
            if (val >= 0) {
                double dist = sqrt((x - i1) * (x - i1) + (y - i2) * (y - i2));
                if ((best_ind < 0) || (dist < best_dist)) {
                    best_dist = dist;
                    best_ind = val;
                    best_x = x;
                    best_y = y;
                }
            }
        }
    }
    return best_ind;
}

void MVFiringEventViewPrivate::set_current_event_index(int ind)
{
    if (ind == m_current_event_index)
        return;
    m_current_event_index = ind;
    //emit!!!!!
    q->update();
}

void MVFiringEventViewPrivate::export_image()
{
    QImage img = q->renderImage(1800, 900);
    user_save_image(img);
}

void MVFiringEventViewPrivate::do_paint(QPainter& painter, int W, int H)
{
    painter.fillRect(0, 0, W, H, QColor(160, 160, 160));

    QRectF target = QRectF(m_hmargin, m_vmargin, W - 2 * m_hmargin, H - 2 * m_vmargin);
    m_target_rect = target;

    QFont font = painter.font();
    font.setPixelSize(24);
    painter.setFont(font);
    QColor epoch_color = QColor(150, 150, 150);
    //double margin=(m_max_amplitude-m_min_amplitude)*0.3;
    double margin = 0;
    for (int i = 0; i < m_epochs.count(); i++) {
        Epoch epoch = m_epochs[i];
        QPointF pt1 = coord2imagepix(QPointF(epoch.t_begin, m_max_amplitude + margin), target.width(), target.height());
        QPointF pt2 = coord2imagepix(QPointF(epoch.t_end, m_min_amplitude - margin), target.width(), target.height());
        QRect RR(target.x() + pt1.x(), target.y() + pt1.y(), pt2.x() - pt1.x(), pt2.y() - pt1.y());
        painter.fillRect(RR, epoch_color);
        RR.adjust(0, m_vmargin, 0, 0);
        if (fabs(m_max_amplitude) > fabs(m_min_amplitude)) {
            painter.drawText(RR, Qt::AlignTop | Qt::AlignHCenter, epoch.name);
        }
        else {
            painter.drawText(RR, Qt::AlignBottom | Qt::AlignHCenter, epoch.name);
        }
    }

    painter.drawImage(target, m_image.scaled(target.width(), target.height(), Qt::IgnoreAspectRatio, Qt::SmoothTransformation));
}
