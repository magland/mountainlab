#include "sstimeseriesplot.h"
#include <QPainter>
#include <QList>
#include "plotarea.h"
#include <QDebug>
#include <QTime>
#include <QMap>
#include <QTimer>
#include <QMouseEvent>
#include "sscommon.h"
#include "computationthread.h"
#include "taskprogress.h"
#define MIN_NUM_PIXELS 4000

class DataLoaderThread : public ComputationThread {
public:
    //input
    int x1, x2;
    int msfactor;
    DiskArrayModel_New* array_model;

    //output
    Mda data;

    void compute();
};

struct Some_Info {
    int x1, x2;
    int msfactor;
};

class SSTimeSeriesPlotPrivate {
public:
    SSTimeSeriesPlot* q;
    DiskArrayModel_New* m_data;
    SSLabelsModel1* m_labels;
    bool m_labels_is_owner;
    SSLabelsModel1* m_compare_labels;
    bool m_compare_labels_is_owner;
    int m_max_timepoint;
    int m_num_channels;
    int m_dim3;
    int m_control_panel_height;
    int m_bar_width;

    QList<double> m_minvals;
    QList<double> m_maxvals;
    PlotArea m_plot_area;
    QList<double> m_yrange_minvals, m_yrange_maxvals;
    QList<double> m_plot_offsets;
    QList<double> m_plot_minvals;
    QList<double> m_plot_maxvals;
    QList<double> m_plot_y1;
    QList<double> m_plot_y2;
    QList<QColor> m_channel_colors;
    QStringList m_channel_labels;
    QList<QColor> m_label_colors;
    QPixmap m_image;
    bool m_image_needs_update;
    //QMap<int,Mda *> m_multiscale_min_arrays;
    //QMap<int,Mda *> m_multiscale_max_arrays;
    int m_margins[4];
    bool m_uniform_vertical_channel_spacing;
    bool m_use_fixed_vertical_channel_spacing;
    double m_fixed_vertical_channel_spacing;
    DataLoaderThread m_data_loader;

    QColor get_channel_color(int num);
    void setup_plot_area();
    void schedule_setup_plot_area();
    bool m_setup_plot_area_scheduled;

    void draw_control_panel(QPainter* P);
};

SSTimeSeriesPlot::SSTimeSeriesPlot(QWidget* parent)
    : SSAbstractPlot(parent)
{
    d = new SSTimeSeriesPlotPrivate();
    d->q = this;

    d->m_image_needs_update = false;

    d->m_data = 0;
    d->m_labels = 0;
    d->m_labels_is_owner = 0;
    d->m_compare_labels = 0;
    d->m_compare_labels_is_owner = 0;

    d->m_max_timepoint = 0;
    d->m_num_channels = 1;
    d->m_dim3 = 1;
    d->m_control_panel_height = 16;
    d->m_bar_width = 1;

    d->m_uniform_vertical_channel_spacing = true;
    d->m_use_fixed_vertical_channel_spacing = false;
    d->m_fixed_vertical_channel_spacing = 0;

    QList<QString> color_strings;

    //color_strings << "darkblue" << "darkgreen" << "darkred" << "darkcyan" << "darkmagenta" << "darkorange" << "black";
    color_strings
        << "#282828"
        << "#402020"
        << "#204020"
        << "#202070";
    //<< "#403010"
    //            << "#401050";
    //            << "#103050";

    for (int i = 0; i < color_strings.size(); i++) {
        d->m_channel_colors << QColor(color_strings[i]);
    }

    color_strings.clear();
    /*
	color_strings << "#F7977A" << "#FDC68A"
				<< "#C4DF9B" << "#82CA9D"
				<< "#6ECFF6" << "#8493CA"
				<< "#A187BE" << "#F49AC2"
				<< "#F9AD81" << "#FFF79A"
				<< "#A2D39C" << "#7BCDC8"
				<< "#7EA7D8" << "#8882BE"
				<< "#BC8DBF" << "#F6989D";
                */

    /*
    color_strings
            << "#D0D0D0"
            << "#80D0D0"
            << "#D080D0"
            << "#D0D080"
            << "#8080D0"
            << "#80D080"
            << "#D08080";
    */
    //color_strings << "#CCBBFF";
    //color_strings << "#AA0000";
    color_strings << "#CCAAAA";

    for (int i = 0; i < color_strings.size(); i++) {
        d->m_label_colors << QColor(color_strings[i]);
    }

    d->m_plot_area.setMarkerColors(d->m_label_colors);
    QList<QString> mlabels;
    for (int i = 0; i < 5000; i++)
        mlabels << QString("%1").arg(i);
    d->m_plot_area.setMarkerLabels(mlabels);

    d->m_margins[0] = d->m_margins[1] = d->m_margins[2] = d->m_margins[3] = 0;

    connect(this, SIGNAL(replotNeeded()), this, SLOT(slot_replot_needed()));
    connect(&d->m_data_loader, SIGNAL(computationFinished()), this, SLOT(slot_data_loaded()));
}

SSTimeSeriesPlot::~SSTimeSeriesPlot()
{
    /*
	foreach (Mda *X,d->m_multiscale_min_arrays) {
		delete X;
	}
	foreach (Mda *X,d->m_multiscale_max_arrays) {
		delete X;
	}*/
    d->m_data_loader.stopComputation(); // important do take care of this before things start getting destructed!
    if (d->m_labels) {
        if (d->m_labels_is_owner)
            delete d->m_labels;
    }
    if (d->m_compare_labels) {
        if (d->m_compare_labels_is_owner)
            delete d->m_compare_labels;
    }
    delete d;
}

void SSTimeSeriesPlot::updateSize(int W, int H)
{
    d->m_plot_area.setPosition(d->m_margins[0], d->m_margins[2]);
    d->m_plot_area.setSize(W - d->m_margins[0] - d->m_margins[1], H - d->m_margins[2] - d->m_margins[3]);
}

void SSTimeSeriesPlot::paintPlot(QPainter* painter, int W, int H)
{
    updateSize(W, H);

    if ((width() != W) || (height() != d->m_image.height())) {
        d->m_image_needs_update = true;
    }

    int control_panel_height = d->m_control_panel_height;

    int W0 = W;
    int H0 = H - control_panel_height;

    if (d->m_image_needs_update) {
        d->m_image = QPixmap(W0, H0);
        d->m_image.fill(QColor(0, 0, 0, 0));
        QPainter painter2(&d->m_image);
        d->m_plot_area.refresh(&painter2);
        d->m_image_needs_update = false;
    }

    painter->drawPixmap(0, control_panel_height, d->m_image);
    if (d->m_control_panel_height) {
        d->draw_control_panel(painter);
    }
}

void SSTimeSeriesPlot::mousePressEvent(QMouseEvent* evt)
{
    if (evt->pos().y() <= d->m_control_panel_height) {
    } else {
        SSAbstractPlot::mousePressEvent(evt);
    }
}

void SSTimeSeriesPlot::mouseReleaseEvent(QMouseEvent* evt)
{
    if (evt->pos().y() <= d->m_control_panel_height) {
        float frac0 = evt->pos().x() * 1.0 / d->m_bar_width;
        int t0 = (int)((d->m_max_timepoint + 1) * frac0);
        emit requestMoveToTimepoint(t0);
    } else {
        SSAbstractPlot::mouseReleaseEvent(evt);
    }
}

void SSTimeSeriesPlot::mouseMoveEvent(QMouseEvent* evt)
{
    if (evt->pos().y() <= d->m_control_panel_height) {
    } else {
        SSAbstractPlot::mouseMoveEvent(evt);
    }
}

void SSTimeSeriesPlot::setData(SSARRAY* data)
{
    if (!data) {
        qWarning() << "Unexpected problem in SSTimeSeriesPlot::setData. data is null.";
        return;
    }
    d->m_data = data;
    d->m_max_timepoint = d->m_data->N2() - 1;
    d->m_num_channels = d->m_data->N1();
    d->m_dim3 = d->m_data->dim3();

    int M = d->m_num_channels;
    int N = d->m_max_timepoint + 1;

    if ((xRange().y >= N - 1) || (xRange().y < 200)) {
        setXRange(vec2(0, N - 1));
    }

    int scale0 = 1;
    while (N / scale0 > 100)
        scale0 *= MULTISCALE_FACTOR;
    Mda data0 = d->m_data->loadData(scale0, 0, N / scale0 + 1);

    d->m_minvals.clear();
    d->m_maxvals.clear();
    for (int ch = 0; ch < M; ch++) {
        d->m_minvals << 0;
        d->m_maxvals << 0;
    }
    //fixed on 12/4/15
    if (data0.N2() > 1) {
        for (int ch = 0; ch < M; ch++) {
            d->m_minvals[ch] = data0.value(ch, 0L, 0L);
            d->m_maxvals[ch] = data0.value(ch, 0L, 0L);
            for (long i3 = 0; i3 < data0.N3(); i3++) {
                for (long i2 = 0; i2 < data0.N2(); i2++) {
                    double val = data0.value(ch, i2, i3);
                    if (val < d->m_minvals[ch]) {
                        d->m_minvals[ch] = val;
                    }
                    if (val > d->m_maxvals[ch]) {
                        d->m_maxvals[ch] = val;
                    }
                }
            }
        }
    } else {
        for (int ch = 0; ch < M; ch++) {
            d->m_minvals[ch] = 0;
            d->m_maxvals[ch] = 1;
        }
    }

    d->m_image_needs_update = true;
    emit this->replotNeeded();
    update();
}

void SSTimeSeriesPlot::setLabels(SSLabelsModel1* L, bool is_owner)
{

    d->m_labels = L;
    d->m_labels_is_owner = is_owner;
    d->m_image_needs_update = true;
    update();
}

void SSTimeSeriesPlot::setCompareLabels(SSLabelsModel1* L, bool is_owner)
{
    d->m_compare_labels = L;
    d->m_compare_labels_is_owner = is_owner;
    d->m_image_needs_update = true;
    update();
}

void SSTimeSeriesPlot::initialize()
{
    d->schedule_setup_plot_area();
}

DiskArrayModel_New* SSTimeSeriesPlot::data()
{
    return d->m_data;
}

void SSTimeSeriesPlot::setChannelLabels(const QStringList& labels)
{
    d->m_channel_labels = labels;
}

void SSTimeSeriesPlot::setUniformVerticalChannelSpacing(bool val)
{
    d->m_uniform_vertical_channel_spacing = val;
    this->slot_replot_needed();
}

void SSTimeSeriesPlot::setFixedVerticalChannelSpacing(double fixed_val)
{
    d->m_use_fixed_vertical_channel_spacing = true;
    d->m_fixed_vertical_channel_spacing = fixed_val;
    this->slot_replot_needed();
}

bool SSTimeSeriesPlot::uniformVerticalChannelSpacing()
{
    return d->m_uniform_vertical_channel_spacing;
}

void SSTimeSeriesPlot::setShowMarkerLines(bool val)
{
    d->m_plot_area.setShowMarkerLines(val);
}

void SSTimeSeriesPlot::setControlPanelVisible(bool val)
{
    if (val)
        d->m_control_panel_height = 16;
    else
        d->m_control_panel_height = 0;
}

QColor SSTimeSeriesPlotPrivate::get_channel_color(int ch)
{
    if (m_channel_colors.size() == 0)
        return QColor(0, 0, 0);
    return m_channel_colors[ch % m_channel_colors.size()];
}

void SSTimeSeriesPlotPrivate::setup_plot_area()
{

    if (!m_data)
        return;

    QTime timer;
    timer.start();

    int M = m_num_channels;

    m_plot_offsets.clear();
    m_plot_minvals.clear();
    m_plot_maxvals.clear();
    m_plot_y1.clear();
    m_plot_y2.clear();
    for (int ch = 0; ch < M; ch++) {
        m_plot_minvals << 0;
        m_plot_maxvals << 0;
        m_plot_offsets << 0;
        m_plot_y1 << 0;
        m_plot_y2 << 0;
    }

    if ((m_yrange_maxvals.size() > 0) && (m_yrange_maxvals.size() == M)) {
        for (int i = 0; i < M; i++) {
            m_plot_minvals[i] = m_yrange_minvals[i] / q->verticalZoomFactor();
            m_plot_maxvals[i] = m_yrange_maxvals[i] / q->verticalZoomFactor();
        }
    } else if ((M > 0) && (m_minvals.size() > 0)) {
        for (int i = 0; i < M; i++) {
            m_plot_minvals[i] = m_minvals.value(i) / q->verticalZoomFactor();
            m_plot_maxvals[i] = m_maxvals.value(i) / q->verticalZoomFactor();
        }
    }

    double max00 = 0;
    for (int ch = 0; ch < M; ch++) {
        if (qAbs(m_plot_minvals[ch]) > max00) {
            max00 = qAbs(m_plot_minvals[ch]);
        }
        if (qAbs(m_plot_maxvals[ch]) > max00) {
            max00 = qAbs(m_plot_maxvals[ch]);
        }
    }

    if (m_uniform_vertical_channel_spacing) {
        if (m_plot_minvals.count() > 0) {
            float global_plot_minval = m_plot_minvals[0];
            float global_plot_maxval = m_plot_maxvals[0];
            for (int i = 0; i < m_plot_minvals.count(); i++) {
                if (m_plot_minvals[i] < global_plot_minval)
                    global_plot_minval = m_plot_minvals[i];
                if (m_plot_maxvals[i] > global_plot_maxval)
                    global_plot_maxval = m_plot_maxvals[i];
            }
            for (int i = 0; i < m_plot_minvals.count(); i++) {
                m_plot_minvals[i] = global_plot_minval;
                m_plot_maxvals[i] = global_plot_maxval;
            }
        }
    }

    double offset = 0;
    if (!q->channelFlip()) { // ahb
        for (int ch = 0; ch < M; ch++) { // jfm's upwards ordering

            if (m_use_fixed_vertical_channel_spacing) {
                m_plot_y1[ch] = offset + m_plot_minvals[ch] - max00;
                offset += m_fixed_vertical_channel_spacing;
            } else {
                m_plot_y1[ch] = offset;
                offset += (-m_plot_minvals[ch]);
            }
            m_plot_offsets[ch] = offset;
            if (m_use_fixed_vertical_channel_spacing) {
                m_plot_y1[ch] = offset + m_plot_maxvals[ch] + max00;
            } else {
                offset += m_plot_maxvals[ch];
                offset += max00 / 20;
                m_plot_y2[ch] = offset;
            }
        }
    } else {
        for (int ch = M - 1; ch >= 0; ch--) { // downwards ordering

            if (m_use_fixed_vertical_channel_spacing) {
                m_plot_y1[ch] = offset + m_plot_minvals[ch] - max00;
                offset -= m_fixed_vertical_channel_spacing;
            } else {
                m_plot_y1[ch] = offset;
                offset += (-m_plot_minvals[ch]);
            }
            m_plot_offsets[ch] = offset;
            if (m_use_fixed_vertical_channel_spacing) {
                m_plot_y2[ch] = offset + m_plot_maxvals[ch] + max00;
            } else {
                offset += m_plot_maxvals[ch];
                offset += max00 / 20;
                m_plot_y2[ch] = offset;
            }
        }
    }

    int xrange_min = q->xRange().x;
    int xrange_max = q->xRange().y;

    //m_plot_area.setXRange(xrange_min - 1, xrange_max + 1);
    if (M > 0) {
        if (!q->channelFlip()) // ahb, matches above
            q->setYRange(vec2(m_plot_offsets[0] + m_plot_minvals[0] - max00 / 2, m_plot_offsets[M - 1] + m_plot_maxvals[M - 1] + max00 / 2));
        else
            q->setYRange(vec2(m_plot_offsets[M - 1] + m_plot_minvals[M - 1] - max00 / 2, m_plot_offsets[0] + m_plot_maxvals[0] + max00 / 2));
    }

    int NN = xrange_max - xrange_min + 1;
    int msfactor = 1;
    while (NN / msfactor / MULTISCALE_FACTOR > MIN_NUM_PIXELS) {
        msfactor *= MULTISCALE_FACTOR;
    }

    int x1, x2;
    if (msfactor == 1) {
        x1 = xrange_min;
        x2 = xrange_max;
    } else {
        x1 = (xrange_min / msfactor) * msfactor;
        x2 = (xrange_max / msfactor) * msfactor;
    }

    m_data_loader.stopComputation();
    m_data_loader.x1 = x1;
    m_data_loader.x2 = x2;
    m_data_loader.msfactor = msfactor;
    m_data_loader.array_model = m_data;

    //m_data_loader.startComputation();
    /// TODO do we want to do this asynchronously?
    m_data_loader.compute();

    //the label markers
    m_plot_area.clearMarkers();
    int max_range_for_showing_labels = 1e6;

    ///This section is bad.
    float tmprange = (x2 - x1 + 1);
    float tmprange_max = max_range_for_showing_labels;
    float tmprange_min = 10000; //magic, hack, hard-code, bug
    float pp = qMax(0.0F, qMin(1.0F, (1 - (tmprange - tmprange_min) / (tmprange_max - tmprange_min))));
    //pp=pp*pp*pp*pp; //i'm actually proud of this bug
    //int alpha=(int)(pp*255.999);
    //m_plot_area.setMarkerAlpha((int)alpha);
    if ((m_labels) && (pp > 0)) {
        Mda TL = m_labels->getTimepointsLabels(x1, x2);
        int K = TL.N2();
        for (int i = 0; i < K; i++) {
            m_plot_area.addMarker(TL.value(0, i) - 1, TL.value(1, i));
        }
    }
    if ((m_compare_labels) && (pp > 0)) {
        Mda TL = m_compare_labels->getTimepointsLabels(x1, x2);
        int K = TL.N2();
        for (int i = 0; i < K; i++) {
            m_plot_area.addCompareMarker(TL.value(0, i) - 1, TL.value(1, i));
        }
    }

    //split the clips
    if ((x2 - x1 + 1 <= max_range_for_showing_labels / 100) && (q->timepointMapping().totalSize() > 1)) {
        DiskReadMda TM = q->timepointMapping();
        int last_time = 0;
        for (int xx = x1; xx <= x2; xx++) {
            int tt = TM.value(0, xx);
            if ((xx > x1) && (tt != last_time + 1)) {
                m_plot_area.addMarker(xx - 0.5, 0);
            }
            last_time = tt;
        }
    }

    if (m_dim3 > 1) {
        int tmp = (m_max_timepoint + 1) / m_dim3;
        for (int xx = x1; xx <= x2; xx++) {
            if (tmp > 0) {
                if (xx % tmp == 0) {
                    m_plot_area.addMarker(xx - 0.5, 0);
                }
            }
        }
    }

    q->slot_data_loaded();
}

void SSTimeSeriesPlotPrivate::schedule_setup_plot_area()
{
    if (m_setup_plot_area_scheduled)
        return;
    m_setup_plot_area_scheduled = true;
    QTimer::singleShot(0, q, SLOT(slot_setup_plot_area()));
}

void SSTimeSeriesPlotPrivate::draw_control_panel(QPainter* P)
{
    QColor col1(50, 50, 50, 20);
    QColor col2(50, 50, 50, 40);

    m_bar_width = q->width();
    int bar_width = m_bar_width;
    QRect RR1(0, 0, bar_width, m_control_panel_height);
    P->fillRect(RR1, col1);

    Vec2 tmp = q->xRange();
    int t0 = tmp.x, t1 = tmp.y;
    if (t0 < 0) {
        P->fillRect(RR1, col2);
    } else {
        int num_timepoints = m_max_timepoint + 1;
        int xpix0 = bar_width * 1.0 * t0 / num_timepoints;
        int xpix1 = bar_width * 1.0 * t1 / num_timepoints;
        if (xpix1 < xpix0 - 3)
            xpix1 = xpix0 + 3;
        QRect RR2(xpix0, 0, xpix1 - xpix0, m_control_panel_height);
        P->fillRect(RR2, col2);
    }
}

Vec2 SSTimeSeriesPlot::coordToPix(Vec2 coord)
{
    return d->m_plot_area.coordToPix(coord);
}

Vec2 SSTimeSeriesPlot::pixToCoord(Vec2 pix)
{
    return d->m_plot_area.pixToCoord(pix);
}

int SSTimeSeriesPlot::pixToChannel(Vec2 pix)
{
    Vec2 coord = pixToCoord(pix);
    for (int i = 0; i < d->m_plot_offsets.count(); i++) {
        double y1 = d->m_plot_offsets.value(i) + d->m_plot_minvals.value(i);
        double y2 = d->m_plot_offsets.value(i) + d->m_plot_maxvals.value(i);
        if ((y1 <= coord.y) && (coord.y <= y2))
            return i;
    }
    return -1;
}

void SSTimeSeriesPlot::setMargins(int left, int right, int top, int bottom)
{
    d->m_margins[0] = left;
    d->m_margins[1] = right;
    d->m_margins[2] = top;
    d->m_margins[3] = bottom;
}

/*void SSTimeSeriesPlot::setConnectZeros(bool val)
{
	d->m_plot_area.setConnectZeros(val);
}*/

void SSTimeSeriesPlot::slot_replot_needed()
{
    d->schedule_setup_plot_area();
}

void SSTimeSeriesPlot::slot_setup_plot_area()
{
    d->m_setup_plot_area_scheduled = false;
    d->setup_plot_area();
}

void SSTimeSeriesPlot::slot_data_loaded()
{
    //d->m_data_loader.stopComputation(); //because I'm paranoid
    int x1 = d->m_data_loader.x1;
    int x2 = d->m_data_loader.x2;
    int msfactor = d->m_data_loader.msfactor;
    int M = d->m_num_channels;
    Mda data0 = d->m_data_loader.data;
    d->m_plot_area.clearSeries();
    for (int ch = 0; ch < M; ch++) {
        if (msfactor == 1) {
            Mda xvals;
            xvals.allocate(1, x2 - x1 + 1);
            for (int x = x1; x <= x2; x++) {
                xvals.setValue(x, 0, x - x1);
            }
            Mda yvals;
            yvals.allocate(1, x2 - x1 + 1);
            for (int ii = x1; ii <= x2; ii++) {
                double val = data0.value(ch, ii - x1);
                yvals.setValue(val, 0, ii - x1);
            }
            QColor color = d->get_channel_color(ch);
            PlotSeries SS;
            SS.xvals = xvals;
            SS.yvals = yvals;
            SS.color = color;
            SS.offset = d->m_plot_offsets[ch];
            SS.plot_pairs = false;
            QString label0 = d->m_channel_labels.value(ch);
            if (label0.isEmpty())
                label0 = QString("%1").arg(ch + 1);
            SS.name = label0;
            d->m_plot_area.addSeries(SS);
        } else {

            Mda xvals;
            xvals.allocate(1, (x2 - x1) / msfactor * 2);
            for (int x = x1; x < x2; x += msfactor) {
                xvals.setValue(x, 0, (x - x1) / msfactor * 2);
                xvals.setValue(x, 0, (x - x1) / msfactor * 2 + 1);
            }
            Mda yvals;
            yvals.allocate(1, (x2 - x1) / msfactor * 2);
            for (int ii = x1; ii < x2; ii += msfactor) {
                double val1 = data0.value(ch, ((ii - x1) / msfactor) * 2);
                double val2 = data0.value(ch, ((ii - x1) / msfactor) * 2 + 1);
                yvals.setValue(val1, 0, (ii - x1) / msfactor * 2);
                yvals.setValue(val2, 0, (ii - x1) / msfactor * 2 + 1);
            }
            QColor color = d->get_channel_color(ch);
            PlotSeries SS;
            SS.xvals = xvals;
            SS.yvals = yvals;
            SS.color = color;
            SS.offset = d->m_plot_offsets[ch];
            SS.plot_pairs = false;
            SS.name = QString("%1").arg(ch + 1);
            d->m_plot_area.addSeries(SS);
        }
    }
    d->m_image_needs_update = true;
    update();
}

void SSTimeSeriesPlot::setXRange(const Vec2& range)
{
    if ((xRange().x != range.x) || (xRange().y != range.y)) {
        SSAbstractPlot::setXRange(range);

        int xrange_min = xRange().x;
        int xrange_max = xRange().y;
        d->m_plot_area.setXRange(xrange_min - 1, xrange_max + 1);

        d->schedule_setup_plot_area();
    }
}
void SSTimeSeriesPlot::setYRange(const Vec2& range)
{
    if ((yRange().x != range.x) || (yRange().y != range.y)) {
        SSAbstractPlot::setYRange(range);
        d->m_plot_area.setYRange(range.x, range.y);
    }
}

void DataLoaderThread::compute()
{
    data = array_model->loadData(msfactor, x1 / msfactor, x2 / msfactor);
}
