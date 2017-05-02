#include "mvclipsview.h"
#include <QDebug>
#include <mvcontext.h>

struct mvclipsview_coord {
    mvclipsview_coord(int i, int m, double t0, double v0)
    {
        clip_index = i;
        channel = m;
        t = t0;
        val = v0;
    }
    int clip_index = 0;
    int channel = 0;
    double t = 0;
    double val = 0;
};

class MVClipsViewPrivate {
public:
    MVClipsView* q;

    DiskReadMda32 m_clips;
    //QVector<double> m_times;
    //QVector<int> m_labels;
    MVContext* m_context;
    double m_clip_index_offset = 0;
    double m_pct_space_per_clip = 1;
    double m_vert_scale_factor = 1;

    QPointF coord2pix(const mvclipsview_coord& C);
    void auto_set_vert_scale_factor();
    void auto_set_pct_space_per_clip();
};

MVClipsView::MVClipsView(MVAbstractContext* context)
{
    d = new MVClipsViewPrivate;
    d->q = this;

    MVContext* c = qobject_cast<MVContext*>(context);
    Q_ASSERT(c);

    d->m_context = c;
}

MVClipsView::~MVClipsView()
{
    delete d;
}

void MVClipsView::setClips(const DiskReadMda32& clips)
{
    d->m_clips = clips;
    d->m_clip_index_offset = 0;
    d->auto_set_vert_scale_factor();
    d->auto_set_pct_space_per_clip();
    this->update();
}

void MVClipsView::paintEvent(QPaintEvent* evt)
{
    Q_UNUSED(evt)
    QPainter painter(this);

    int M = d->m_clips.N1();
    int T = d->m_clips.N2();
    int L = d->m_clips.N3();

    for (int i = 0; i < L; i++) {
        for (int m = 0; m < M; m++) {
            QColor col = d->m_context->channelColor(m + 1);
            painter.setPen(QPen(QBrush(col), 3));
            QPainterPath path;
            for (int t = 0; t < T; t++) {
                double val = d->m_clips.value(m, t, i);
                mvclipsview_coord C(i, m, t, val);
                QPointF pp = d->coord2pix(C);
                if (t == 0)
                    path.moveTo(pp);
                else
                    path.lineTo(pp);
            }
            painter.drawPath(path);
        }
    }
}

/*
void MVClipsView::setTimes(const QVector<double>& times)
{
    d->m_times = times;
}

void MVClipsView::setLabels(const QVector<int>& labels)
{
    d->m_labels = labels;
}
*/

QPointF MVClipsViewPrivate::coord2pix(const mvclipsview_coord& C)
{
    int M = m_clips.N1();
    int T = m_clips.N2();
    double x = (C.clip_index - m_clip_index_offset + C.t / T) * m_pct_space_per_clip * q->width();
    double y = (C.channel + 0.5 - C.val * m_vert_scale_factor * 0.5) / M * q->height();
    return QPointF(x, y);
}

void MVClipsViewPrivate::auto_set_vert_scale_factor()
{
    Mda32 X;
    if (!m_clips.readChunk(X, 0, 0, 0, m_clips.N1(), m_clips.N2(), m_clips.N3())) {
        qWarning() << "Unable to read chunk of clips in auto_set_vert_scale_factor of clips view.";
        return;
        ;
    }
    double val = qMax(qAbs(X.minimum()), qAbs(X.maximum()));
    if (!val) {
        m_vert_scale_factor = 1;
    }
    else {
        m_vert_scale_factor = 1 / val;
    }
}

void MVClipsViewPrivate::auto_set_pct_space_per_clip()
{
    if (!m_clips.N3())
        return;
    m_pct_space_per_clip = 1 / m_clips.N3();
}
