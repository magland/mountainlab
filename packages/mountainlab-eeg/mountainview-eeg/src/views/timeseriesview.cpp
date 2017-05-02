/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
** Created: 5/27/2016
*******************************************************/

#include "timeseriesview.h"
#include "multiscaletimeseries.h"
#include "mvtimeseriesrendermanager.h"
#include <math.h>

#include <QAction>
#include <QIcon>
#include <QImageWriter>
#include <QMouseEvent>
#include <QPainter>
#include <mveegcontext.h>

struct mvtsv_channel {
    int channel;
    QString label;
    QRectF geometry;
};

class TimeseriesViewCalculator {
public:
    //input
    DiskReadMda32 timeseries;
    QString timeseries_freq_range;
    //QString mlproxy_url;

    //output
    MultiScaleTimeSeries msts;
    double minval, maxval;
    int num_channels;

    void compute();
};

class TimeseriesViewPrivate {
public:
    TimeseriesView* q;
    MultiScaleTimeSeries m_msts;

    TimeseriesViewCalculator m_calculator;

    double m_amplitude_factor;
    QList<mvtsv_channel> m_channels;
    bool m_layout_needed;
    int m_num_channels;

    MVTimeSeriesRenderManager m_render_manager;

    QList<mvtsv_channel> make_channel_layout(int M);
    void paint_channel_labels(QPainter* painter, double W, double H);

    double val2ypix(int m, double val);
    double ypix2val(int m, double ypix);
};

TimeseriesView::TimeseriesView(MVAbstractContext* context)
    : MVTimeSeriesViewBase(context)
{
    d = new TimeseriesViewPrivate;
    d->q = this;

    MVEEGContext* c = qobject_cast<MVEEGContext*>(context);
    Q_ASSERT(c);

    d->m_amplitude_factor = 1.0;
    d->m_layout_needed = true;
    d->m_num_channels = 1;

    {
        QAction* a = new QAction(QIcon(":/images/vertical-zoom-in.png"), "Vertical Zoom In", this);
        a->setProperty("action_type", "toolbar");
        a->setToolTip("Vertical zoom in. Alternatively, use the UP arrow.");
        this->addAction(a);
        connect(a, SIGNAL(triggered(bool)), this, SLOT(slot_vertical_zoom_in()));
    }
    {
        QAction* a = new QAction(QIcon(":/images/vertical-zoom-out.png"), "Vertical Zoom Out", this);
        a->setProperty("action_type", "toolbar");
        a->setToolTip("Vertical zoom out. Alternatively, use the DOWN arrow.");
        this->addAction(a);
        connect(a, SIGNAL(triggered(bool)), this, SLOT(slot_vertical_zoom_out()));
    }

    QObject::connect(&d->m_render_manager, SIGNAL(updated()), this, SLOT(update()));
    this->recalculateOn(context, SIGNAL(currentTimeseriesChanged()));

    this->recalculate();
}

TimeseriesView::~TimeseriesView()
{
    this->stopCalculation(); //this is needed because of deletion of the multi-scale timeseries
    delete d;
}

void TimeseriesView::prepareCalculation()
{
    MVEEGContext* c = qobject_cast<MVEEGContext*>(mvContext());
    Q_ASSERT(c);

    d->m_layout_needed = true;
    d->m_calculator.timeseries = c->currentTimeseries();
    d->m_calculator.timeseries_freq_range = c->option("timeseries_freq_range").toString();
    //d->m_calculator.mlproxy_url = c->mlProxyUrl();

    MVTimeSeriesViewBase::prepareCalculation();
}

void TimeseriesView::runCalculation()
{
    d->m_calculator.compute();

    MVTimeSeriesViewBase::runCalculation();
}

void TimeseriesView::onCalculationFinished()
{
    MVEEGContext* c = qobject_cast<MVEEGContext*>(mvContext());
    Q_ASSERT(c);

    d->m_msts = d->m_calculator.msts;
    //d->m_render_manager.setChannelColors(c->channelColors());
    d->m_render_manager.setMultiScaleTimeSeries(d->m_msts);
    d->m_num_channels = d->m_calculator.num_channels;

    double max_range = qMax(qAbs(d->m_msts.minimum()), qAbs(d->m_msts.maximum()));
    if (max_range) {
        this->setAmplitudeFactor(1.5 / max_range);
    }

    d->m_layout_needed = true;

    MVTimeSeriesViewBase::onCalculationFinished();
}

double TimeseriesView::amplitudeFactor() const
{
    return d->m_amplitude_factor;
}

void TimeseriesView::resizeEvent(QResizeEvent* evt)
{
    d->m_layout_needed = true;
    MVTimeSeriesViewBase::resizeEvent(evt);
}

void TimeseriesView::setAmplitudeFactor(double factor)
{
    d->m_amplitude_factor = factor;
    update();
}

void TimeseriesView::autoSetAmplitudeFactorWithinTimeRange()
{
    double min0 = d->m_render_manager.visibleMinimum();
    double max0 = d->m_render_manager.visibleMaximum();
    double factor = qMax(qAbs(min0), qAbs(max0));
    if (factor)
        this->setAmplitudeFactor(1.5 / factor);
}

void TimeseriesView::paintContent(QPainter* painter)
{
    MVEEGContext* c = qobject_cast<MVEEGContext*>(mvContext());
    Q_ASSERT(c);

    // Geometry of channels
    if (d->m_layout_needed) {
        int M = d->m_num_channels;
        d->m_layout_needed = false;
        d->m_channels = d->make_channel_layout(M);
        for (int m = 0; m < M; m++) {
            mvtsv_channel* CH = &d->m_channels[m];
            CH->channel = m;
            CH->label = QString("%1").arg(m + 1);
        }
    }

    double WW = this->contentGeometry().width();
    double HH = this->contentGeometry().height();
    QImage img;
    img = d->m_render_manager.getImage(c->currentTimeRange().min, c->currentTimeRange().max, d->m_amplitude_factor, WW, HH);
    painter->drawImage(this->contentGeometry().left(), this->contentGeometry().top(), img);

    // Channel labels
    d->paint_channel_labels(painter, this->width(), this->height());
}

void TimeseriesView::keyPressEvent(QKeyEvent* evt)
{
    if (evt->key() == Qt::Key_Up) {
        d->m_amplitude_factor *= 1.2;
        update();
    }
    else if (evt->key() == Qt::Key_Down) {
        d->m_amplitude_factor /= 1.2;
        update();
    }
    else {
        MVTimeSeriesViewBase::keyPressEvent(evt);
    }
}

void TimeseriesView::slot_vertical_zoom_in()
{
    d->m_amplitude_factor *= 1.2;
    update();
}

void TimeseriesView::slot_vertical_zoom_out()
{
    d->m_amplitude_factor /= 1.2;
    update();
}

QList<mvtsv_channel> TimeseriesViewPrivate::make_channel_layout(int M)
{
    QList<mvtsv_channel> channels;
    if (!M)
        return channels;
    QRectF RR = q->contentGeometry();
    double space = 0;
    double channel_height = (RR.height() - (M - 1) * space) / M;
    double y0 = RR.top();
    for (int m = 0; m < M; m++) {
        mvtsv_channel X;
        X.geometry = QRectF(RR.left(), y0, RR.width(), channel_height);
        channels << X;
        y0 += channel_height + space;
    }
    return channels;
}

void TimeseriesViewPrivate::paint_channel_labels(QPainter* painter, double W, double H)
{
    Q_UNUSED(W)
    Q_UNUSED(H)
    QPen pen = painter->pen();
    pen.setColor(Qt::black);
    painter->setPen(pen);

    QFont font = painter->font();
    font.setPixelSize(13);
    painter->setFont(font);

    int M = m_num_channels;
    for (int m = 0; m < M; m++) {
        double ypix = val2ypix(m, 0);
        QRectF rect(0, ypix - 30, q->contentGeometry().left() - 5, 60);
        QString str = QString("[%1]").arg(m + 1);
        painter->drawText(rect, Qt::AlignRight | Qt::AlignVCenter, str);
    }
}

double TimeseriesViewPrivate::val2ypix(int m, double val)
{
    if (m >= m_channels.count())
        return 0;

    mvtsv_channel* CH = &m_channels[m];

    double py = CH->geometry.y() + CH->geometry.height() / 2 + CH->geometry.height() / 2 * (val * m_amplitude_factor);
    return py;
}

double TimeseriesViewPrivate::ypix2val(int m, double ypix)
{
    if (m >= m_channels.count())
        return 0;

    mvtsv_channel* CH = &m_channels[m];

    if (m_amplitude_factor) {
        double val = (ypix - (CH->geometry.y() + CH->geometry.height() / 2)) / m_amplitude_factor / (CH->geometry.height() / 2);
        return val;
    }
    else
        return 0;
}

DiskReadMda apply_filter(DiskReadMda X, QString freq_range)
{
    DiskReadMda ret = X;
    //apply filter here
    return ret;
}

void TimeseriesViewCalculator::compute()
{
    DiskReadMda32 timeseries2 = apply_filter(timeseries, this->timeseries_freq_range);
    msts.setData(timeseries2);
    //msts.setMLProxyUrl(mlproxy_url);
    msts.initialize();
    minval = msts.minimum();
    maxval = msts.maximum();
    num_channels = timeseries.N1();
}

TimeseriesDataFactory::TimeseriesDataFactory(MVMainWindow* mw, QObject* parent)
    : MVAbstractViewFactory(mw, parent)
{
}

QString TimeseriesDataFactory::id() const
{
    return QStringLiteral("open-timeseries-data");
}

QString TimeseriesDataFactory::name() const
{
    return tr("Timeseries Data");
}

QString TimeseriesDataFactory::title() const
{
    return tr("Timeseries");
}

MVAbstractView* TimeseriesDataFactory::createView(MVAbstractContext* context)
{
    MVEEGContext* c = qobject_cast<MVEEGContext*>(context);
    Q_ASSERT(c);

    TimeseriesView* X = new TimeseriesView(context);
    //QList<int> ks = c->selectedClusters();
    //X->setLabelsToView(ks.toSet());
    return X;
}
