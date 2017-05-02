/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
** Created: 6/1/2016
*******************************************************/

#include "mvtimeseriesrendermanager.h"

#include <QImage>
#include <QPainter>
#include <QThread>
#include <QTimer>
#include <QCoreApplication>
#include <QImageWriter>
#include "taskprogress.h"

#define PANEL_NUM_POINTS 1200
#define PANEL_WIDTH PANEL_NUM_POINTS * 2
#define PANEL_HEIGHT_PER_CHANNEL 30
#define MIN_PANEL_HEIGHT 600
#define MAX_PANEL_HEIGHT 1800
#define PANEL_HEIGHT(M) (int) qMin(MAX_PANEL_HEIGHT * 1.0, qMax(MIN_PANEL_HEIGHT * 1.0, PANEL_HEIGHT_PER_CHANNEL * M * 1.0))

#define MAX_NUM_IMAGE_PIXELS 50 * 1e6

struct ImagePanel {
    int ds_factor;
    int panel_width;
    int panel_num_points;
    int index;
    double amp_factor;
    QImage image;
    Mda32 min_data, max_data;
    QString make_code();
};

class MVTimeSeriesRenderManagerPrivate {
public:
    MVTimeSeriesRenderManager* q;
    MultiScaleTimeSeries m_ts;
    QMap<QString, ImagePanel> m_image_panels;
    QSet<QString> m_queued_or_running_panel_codes;
    ThreadManager m_thread_manager;
    double m_total_num_image_pixels;
    QList<QColor> m_channel_colors;
    double m_visible_minimum, m_visible_maximum;

    ImagePanel render_panel(ImagePanel p);
    void start_compute_panel(ImagePanel p);
    void stop_compute_panel(const QString& code);
    ImagePanel* closest_ancestor_panel(ImagePanel p);
    void cleanup_images(double t1, double t2, double amp_factor);
};

MVTimeSeriesRenderManager::MVTimeSeriesRenderManager()
{
    d = new MVTimeSeriesRenderManagerPrivate;
    d->q = this;
    d->m_total_num_image_pixels = 0;
    d->m_visible_minimum = d->m_visible_maximum = 0;
}

MVTimeSeriesRenderManager::~MVTimeSeriesRenderManager()
{
    delete d;
}

void MVTimeSeriesRenderManager::clear()
{
    d->m_image_panels.clear();
    d->m_total_num_image_pixels = 0;
    d->m_queued_or_running_panel_codes.clear();
    d->m_thread_manager.clear();
}

void MVTimeSeriesRenderManager::setMultiScaleTimeSeries(MultiScaleTimeSeries ts)
{
    d->m_ts = ts;
    this->clear();
}

void MVTimeSeriesRenderManager::setChannelColors(const QList<QColor>& colors)
{
    d->m_channel_colors = colors;
}

double MVTimeSeriesRenderManager::visibleMinimum() const
{
    return d->m_visible_minimum;
}

double MVTimeSeriesRenderManager::visibleMaximum() const
{
    return d->m_visible_maximum;
}

QImage MVTimeSeriesRenderManager::getImage(double t1, double t2, double amp_factor, double W, double H)
{
    // this is called in the gui thread
    if (t1 >= t2) {
        qWarning() << "t1>=t2 in MVTimeSeriesRenderManager::getImage" << t1 << t2;
        QImage tmp(W, H, QImage::Format_ARGB32);
        tmp.fill(QColor(200, 255, 200));
        return tmp;
    }

    QImage ret(W, H, QImage::Format_ARGB32);
    QColor transparent(0, 0, 0, 6);
    ret.fill(transparent);
    QPainter painter(&ret);

    int ds_factor = 1;
    int panel_width = PANEL_WIDTH;
    int panel_num_points = PANEL_NUM_POINTS;
    //points per pixel should be around 1
    while (((t2 - t1) / ds_factor) / W > 3) {
        ds_factor *= 3;
    }
    double points_per_pixel = ((t2 - t1) / ds_factor) / W;
    if (points_per_pixel < 1.0 / 3) {
        panel_num_points /= 3;
    }
    if ((t2 - t1 < panel_num_points)) {
        panel_num_points /= 3;
    }

    QSet<QString> panel_codes_needed;
    QList<ImagePanel> panels_to_start;

    d->m_visible_minimum = d->m_visible_maximum = 0;

    int ind1 = (int)(t1 / (ds_factor * panel_num_points));
    int ind2 = (int)(t2 / (ds_factor * panel_num_points));
    for (int iii = ind1; iii <= ind2; iii++) {
        ImagePanel p;
        p.amp_factor = amp_factor;
        p.ds_factor = ds_factor;
        p.panel_width = panel_width;
        p.panel_num_points = panel_num_points;
        p.index = iii;
        panel_codes_needed.insert(p.make_code());
        if (!d->m_image_panels.contains(p.make_code())) {
            if (!d->m_queued_or_running_panel_codes.contains(p.make_code())) {
                panels_to_start << p;
            }
        }

        p = d->render_panel(p);
        d->m_visible_minimum = qMin(d->m_visible_minimum, 1.0 * p.min_data.minimum());
        d->m_visible_maximum = qMax(d->m_visible_maximum, 1.0 * p.max_data.maximum());
        if (p.image.width()) {
            double a1 = (iii * panel_num_points * ds_factor - t1) * 1.0 / (t2 - t1) * W;
            double a2 = ((iii + 1) * panel_num_points * ds_factor - t1) * 1.0 / (t2 - t1) * W;
            if (a2 - a1 < 4000) { //we avoid running out of memory
                painter.drawImage(a1, 0, p.image.scaled(a2 - a1, H, Qt::IgnoreAspectRatio, Qt::SmoothTransformation));
            }
            else {
                QRectF rect(0, 0, ret.width(), ret.height());
                painter.fillRect(rect, Qt::lightGray);
                QString txt = "Image stretched too far.";
                painter.setPen(Qt::black);
                painter.drawText(rect, Qt::AlignCenter | Qt::AlignVCenter, txt);
            }
        }
    }

    //stop the threads that aren't needed
    foreach (QString code, d->m_queued_or_running_panel_codes) {
        if (!panel_codes_needed.contains(code)) {
            d->stop_compute_panel(code);
        }
    }

    for (int i = 0; i < panels_to_start.count(); i++) {
        d->start_compute_panel(panels_to_start[i]);
    }

    if (d->m_total_num_image_pixels > MAX_NUM_IMAGE_PIXELS) {
        d->cleanup_images(t1, t2, amp_factor);
    }

    return ret;
}

void MVTimeSeriesRenderManager::slot_thread_finished()
{
    MVTimeSeriesRenderManagerThread* thread = qobject_cast<MVTimeSeriesRenderManagerThread*>(sender());
    if (!thread)
        return;
    ImagePanel p;
    p.amp_factor = thread->amp_factor;
    p.ds_factor = thread->ds_factor;
    p.panel_width = thread->panel_width;
    p.panel_num_points = thread->panel_num_points;
    p.index = thread->index;
    QString code = p.make_code();
    d->m_queued_or_running_panel_codes.remove(code);

    if (thread->image.width()) {
        d->m_image_panels[code] = p;
        d->m_image_panels[code].image = thread->image;
        d->m_image_panels[code].min_data = thread->min_data;
        d->m_image_panels[code].max_data = thread->max_data;
        d->m_total_num_image_pixels += thread->image.width() * thread->image.height();
        emit updated();
    }
    else {
    }

    thread->deleteLater();
}

QString ImagePanel::make_code()
{
    return QString("amp=%1.ds=%2.pw=%3.pnp=%4.ind=%5").arg(this->amp_factor).arg(this->ds_factor).arg(this->panel_width).arg(this->panel_num_points).arg(this->index);
}

void MVTimeSeriesRenderManagerPrivate::start_compute_panel(ImagePanel p)
{
    if (m_queued_or_running_panel_codes.contains(p.make_code()))
        return;
    MVTimeSeriesRenderManagerThread* thread = new MVTimeSeriesRenderManagerThread;
    QObject::connect(thread, SIGNAL(finished()), q, SLOT(slot_thread_finished()));
    thread->amp_factor = p.amp_factor;
    thread->ds_factor = p.ds_factor;
    thread->panel_width = p.panel_width;
    thread->panel_num_points = p.panel_num_points;
    thread->index = p.index;
    thread->ts = m_ts;
    thread->channel_colors = m_channel_colors;
    m_queued_or_running_panel_codes.insert(p.make_code());
    m_thread_manager.start(p.make_code(), thread);
}

void MVTimeSeriesRenderManagerPrivate::stop_compute_panel(const QString& code)
{
    m_queued_or_running_panel_codes.remove(code);
    m_thread_manager.stop(code);
}

ImagePanel* MVTimeSeriesRenderManagerPrivate::closest_ancestor_panel(ImagePanel p)
{
    QList<ImagePanel*> candidates;
    QStringList keys = m_image_panels.keys();
    foreach (QString key, keys) {
        ImagePanel* pp = &m_image_panels[key];
        if (pp->amp_factor == p.amp_factor) {
            double t1 = p.index * p.ds_factor * p.panel_num_points;
            double t2 = (p.index + 1) * p.ds_factor * p.panel_num_points;
            double s1 = pp->index * pp->ds_factor * p.panel_num_points;
            double s2 = (pp->index + 1) * pp->ds_factor * p.panel_num_points;
            if ((s1 <= t1) && (t2 <= s2)) {
                candidates << pp;
            }
        }
    }
    if (candidates.isEmpty())
        return 0;
    ImagePanel* ret = candidates[0];
    int best_ds_factor = ret->ds_factor;
    for (int i = 0; i < candidates.count(); i++) {
        if (candidates[i]->ds_factor < best_ds_factor) {
            ret = candidates[i];
            best_ds_factor = ret->ds_factor;
        }
    }
    return ret;
}

void MVTimeSeriesRenderManagerPrivate::cleanup_images(double t1, double t2, double amp_factor)
{
    QStringList keys = m_image_panels.keys();
    foreach (QString key, keys) {
        ImagePanel* P = &m_image_panels[key];
        double s1 = P->index * P->ds_factor * P->panel_num_points;
        double s2 = (P->index + 1) * P->ds_factor * P->panel_num_points;
        if (((s2 < t1) || (s1 > t2)) || (P->amp_factor != amp_factor)) {
            //does not intersect or different amplitude
            m_total_num_image_pixels -= P->image.width() * P->image.height();
            m_image_panels.remove(key);
        }
        if (m_total_num_image_pixels < MAX_NUM_IMAGE_PIXELS * 0.5)
            return;
    }
}

QColor MVTimeSeriesRenderManagerThread::get_channel_color(int m)
{
    if (channel_colors.isEmpty())
        return Qt::black;
    return channel_colors[m % channel_colors.count()];
}

void MVTimeSeriesRenderManagerThread::run()
{
    int M = ts.N1();
    if (!M)
        return;

    if (MLUtil::threadInterruptRequested())
        return;

    QImage image0 = QImage(panel_width, PANEL_HEIGHT(M), QImage::Format_ARGB32);
    QColor transparent(0, 0, 0, 0);
    image0.fill(transparent);

    if (MLUtil::threadInterruptRequested())
        return;

    QPainter painter(&image0);
    painter.setRenderHint(QPainter::Antialiasing);

    int t1 = index * panel_num_points;
    int t2 = (index + 1) * panel_num_points;

    Mda32 Xmin, Xmax;
    ts.getData(Xmin, Xmax, t1, t2, ds_factor);

    if (MLUtil::threadInterruptRequested())
        return;

    double space = 0;
    double channel_height = (PANEL_HEIGHT(M) - (M - 1) * space) / M;
    int y0 = 0;
    QPen pen = painter.pen();
    for (int m = 0; m < M; m++) {
        if (MLUtil::threadInterruptRequested())
            return;
        pen.setColor(get_channel_color(m));
        if (ds_factor == 1)
            pen.setWidth(3);
        painter.setPen(pen);
        QRectF geom(0, y0, panel_width, channel_height);
        /*
        QPainterPath path;
        for (int ii = t1; ii <= t2; ii++) {
            double val_min = Xmin.value(m, ii - t1);
            double val_max = Xmax.value(m, ii - t1);
            double pctx = (ii - t1) * 1.0 / (t2 - t1);
            double pcty_min = 1 - (val_min * amp_factor + 1) / 2;
            double pcty_max = 1 - (val_max * amp_factor + 1) / 2;
            QPointF pt_min = QPointF(geom.x() + pctx * geom.width(), geom.y() + pcty_min * geom.height());
            QPointF pt_max = QPointF(geom.x() + pctx * geom.width(), geom.y() + pcty_max * geom.height());
            if (ii == t1)
                path.moveTo(pt_min);
            else
                path.lineTo(pt_min);
            path.lineTo(pt_max);
        }
        painter.drawPath(path);
        */
        QPainterPath path;
        QPointF first;
        for (int ii = t1; ii <= t2; ii++) {
            double val_min = Xmin.value(m, ii - t1);
            double pctx = (ii - t1) * 1.0 / (t2 - t1);
            double pcty_min = 1 - (val_min * amp_factor + 1) / 2;
            QPointF pt_min = QPointF(geom.x() + pctx * geom.width(), geom.y() + pcty_min * geom.height());
            if (ii == t1) {
                first = pt_min;
                path.moveTo(pt_min);
            }
            else
                path.lineTo(pt_min);
        }
        if (ds_factor > 1) {
            for (int ii = t2; ii >= t1; ii--) {
                double val_max = Xmax.value(m, ii - t1);
                double pctx = (ii - t1) * 1.0 / (t2 - t1);
                double pcty_max = 1 - (val_max * amp_factor + 1) / 2;
                QPointF pt_max = QPointF(geom.x() + pctx * geom.width(), geom.y() + pcty_max * geom.height());
                path.lineTo(pt_max);
            }
            path.lineTo(first);
            painter.fillPath(path, painter.pen().color());
        }
        else {
            painter.drawPath(path);
        }

        if (MLUtil::threadInterruptRequested())
            return;

        painter.drawPath(path);

        y0 += channel_height + space;
    }

    if (MLUtil::threadInterruptRequested())
        return;
    min_data = Xmin;
    max_data = Xmax;
    image = image0; //only copy on successful exit
}

ThreadManager::ThreadManager()
{
    QTimer::singleShot(100, this, SLOT(slot_timer()));
}

ThreadManager::~ThreadManager()
{
    m_queued_threads.clear();
    QStringList keys = m_running_threads.keys();
    foreach (QString key, keys) {
        this->stop(key);
    }
    while (!m_running_threads.isEmpty()) {
        qApp->processEvents();
    }
}

void ThreadManager::start(QString id, QThread* thread)
{
    if (m_queued_threads.contains(id))
        return;
    m_queued_threads[id] = thread;
    thread->setProperty("threadmanager_id", id);
    QObject::connect(thread, SIGNAL(finished()), this, SLOT(slot_thread_finished()));
}

void ThreadManager::stop(QString id)
{
    if (m_queued_threads.contains(id))
        m_queued_threads.remove(id);
    if (m_running_threads.contains(id)) {
        m_running_threads[id]->requestInterruption();
    }
}

void ThreadManager::clear()
{
    QStringList keys = m_running_threads.keys();
    foreach (QString key, keys) {
        this->stop(key);
    }
    m_running_threads.clear();
    m_queued_threads.clear();
}

void ThreadManager::slot_timer()
{
    while ((m_running_threads.count() < 4) && (!m_queued_threads.isEmpty())) {
        QString key = m_queued_threads.keys().value(0);
        QThread* T = m_queued_threads[key];
        m_queued_threads.remove(key);
        m_running_threads[key] = T;
        T->start();
    }
    QTimer::singleShot(100, this, SLOT(slot_timer()));
}

void ThreadManager::slot_thread_finished()
{
    QThread* thread = qobject_cast<QThread*>(sender());
    if (!thread)
        return;
    m_running_threads.remove(thread->property("threadmanager_id").toString());
}

ImagePanel MVTimeSeriesRenderManagerPrivate::render_panel(ImagePanel p)
{
    QString code = p.make_code();
    if (m_image_panels.contains(code)) {
        return m_image_panels[code];
    }
    else {
        ImagePanel* p2 = closest_ancestor_panel(p);
        if (p2) {
            double s1 = p2->ds_factor * p2->index * p2->panel_num_points;
            double s2 = p2->ds_factor * (p2->index + 1) * p2->panel_num_points;
            double t1 = p.ds_factor * p.index * p.panel_num_points;
            double t2 = p.ds_factor * (p.index + 1) * p.panel_num_points;
            double a1 = (t1 - s1) / (s2 - s1) * p2->image.width();
            double a2 = (t2 - s1) / (s2 - s1) * p2->image.width();
            p.min_data = p2->min_data;
            p.max_data = p2->max_data;
            p.image = p2->image.copy(a1, 0, a2 - a1, p2->image.height());
            QPainter painter(&p.image);
            QColor col(255, 0, 0, 6);
            painter.fillRect(0, 0, p.image.width(), p.image.height(), col);
            return p;
        }
        else {
            p.image = QImage(100, 100, QImage::Format_ARGB32);
            QColor col(0, 0, 0, 10);
            p.image.fill(col);
            return p;
        }
    }
}
