/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
** Created: 5/27/2016
*******************************************************/

#include "spectrogramview.h"
#include <math.h>

#include <QAction>
#include <QIcon>
#include <QImageWriter>
#include <QMessageBox>
#include <QMouseEvent>
#include <QPainter>
#include <QSlider>
#include <mountainprocessrunner.h>
#include <mveegcontext.h>
#include <qlabel.h>
#include <taskprogress.h>
#include <QFileDialog>

class SpectrogramViewCalculator {
public:
    //input
    DiskReadMda32 timeseries;
    int time_resolution = 32;
    QString spectrogram_freq_range = "";

    //output
    DiskReadMda spectrogram;

    void compute();
};

class SpectrogramViewPrivate {
public:
    SpectrogramView* q;
    DiskReadMda m_spectrogram;
    int m_time_resolution = 0; //corresponding to this spectrogram
    double m_window_min = 0, m_window_max = 4;
    double m_brightness_level = 0;
    QString m_spectrogram_freq_range = "";

    SpectrogramViewCalculator m_calculator;

    QColor get_pixel_color(double val);
    Mda get_spectrogram_data(int t1, int t2);

    static void parse_freq_range(int& ret_min, int& ret_max, QString str);
};

SpectrogramView::SpectrogramView(MVAbstractContext* context)
    : MVTimeSeriesViewBase(context)
{
    d = new SpectrogramViewPrivate;
    d->q = this;

    MVEEGContext* c = qobject_cast<MVEEGContext*>(context);
    Q_ASSERT(c);

    {
        QSlider* S = new QSlider;
        S->setOrientation(Qt::Horizontal);
        S->setRange(-100, 100);
        S->setValue(0);
        S->setMaximumWidth(150);
        this->addToolbarControl(new QLabel("brightness:"));
        this->addToolbarControl(S);
        QObject::connect(S, SIGNAL(valueChanged(int)), this, SLOT(slot_brightness_slider_changed(int)));
    }

    {
        QAction* a = new QAction(QIcon(":/image/gear.png"), "Export .csv", this);
        this->addAction(a);
        connect(a, SIGNAL(triggered(bool)), this, SLOT(slot_export_csv()));
    }

    /*
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
    */

    this->recalculateOnOptionChanged("spectrogram_time_resolution");
    this->recalculateOnOptionChanged("spectrogram_freq_range");
    this->recalculateOn(context, SIGNAL(currentTimeseriesChanged()));

    this->recalculate();
}

SpectrogramView::~SpectrogramView()
{
    this->stopCalculation();
    delete d;
}

void SpectrogramView::prepareCalculation()
{
    MVEEGContext* c = qobject_cast<MVEEGContext*>(mvContext());
    Q_ASSERT(c);

    //d->m_layout_needed = true;
    d->m_calculator.timeseries = c->currentTimeseries();
    d->m_calculator.time_resolution = c->option("spectrogram_time_resolution").toInt();
    d->m_calculator.spectrogram_freq_range = c->option("spectrogram_freq_range").toString();

    MVTimeSeriesViewBase::prepareCalculation();
}

void SpectrogramView::runCalculation()
{
    d->m_calculator.compute();

    MVTimeSeriesViewBase::runCalculation();
}

void SpectrogramView::onCalculationFinished()
{
    MVEEGContext* c = qobject_cast<MVEEGContext*>(mvContext());
    Q_ASSERT(c);

    d->m_spectrogram = d->m_calculator.spectrogram;
    d->m_time_resolution = d->m_calculator.time_resolution;
    d->m_spectrogram_freq_range = d->m_calculator.spectrogram_freq_range;

    MVTimeSeriesViewBase::onCalculationFinished();
}

void SpectrogramView::resizeEvent(QResizeEvent* evt)
{
    //d->m_layout_needed = true;
    MVTimeSeriesViewBase::resizeEvent(evt);
}

void SpectrogramView::paintContent(QPainter* painter)
{
    MVEEGContext* c = qobject_cast<MVEEGContext*>(mvContext());
    Q_ASSERT(c);

    int time_resolution = d->m_time_resolution; //corresponding to the output spectrogram
    if (!time_resolution)
        return;

    MVRange timerange = c->currentTimeRange();
    int t1 = ceil(timerange.min / time_resolution);
    int t2 = floor(timerange.max / time_resolution);

    if (t2 <= t1 + 2)
        return;

    Mda data = d->get_spectrogram_data(t1, t2);

    QImage img = QImage(data.N2(), data.N1(), QImage::Format_RGB32);
    for (int n = 0; n < data.N2(); n++) {
        for (int m = 0; m < data.N1(); m++) {
            double val = data.value(m, n);
            QColor col = d->get_pixel_color(val);
            img.setPixel(n, m, col.rgb());
        }
    }

    QRectF geom = this->contentGeometry();
    // TODO: fine adjustment on this based on t1,t2,timerange
    QImage img2 = img.scaled(geom.width(), geom.height(), Qt::IgnoreAspectRatio, Qt::SmoothTransformation);
    painter->drawImage(geom.topLeft(), img2);

    /*
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
    */
}

void SpectrogramView::keyPressEvent(QKeyEvent* evt)
{
    /*
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
    */
}

void SpectrogramView::slot_brightness_slider_changed(int val)
{
    d->m_brightness_level = val;
    this->update();
}

void SpectrogramView::slot_export_csv()
{

    int time_resolution = d->m_time_resolution; //corresponding to the output spectrogram
    if (!time_resolution) {
        QMessageBox::warning(0, "Problem exporting", "Problem exporting. time_resolution = 0");
        return;
    }

    int t1 = 0;
    int t2 = d->m_spectrogram.N2() - 1;

    if (t2 <= t1 + 2) {
        QMessageBox::warning(0, "Problem exporting", "Problem exporting. t2 <= t1+2");
        return;
    }

    Mda data = d->get_spectrogram_data(t1, t2);

    QString txt;
    txt += "timepoint";
    for (int m = 0; m < data.N1(); m++) {
        txt += QString(",Ch%1").arg(m + 1);
    }
    txt += "\n";
    for (int n = 0; n < data.N2(); n++) {
        txt += QString("%1").arg(n * time_resolution);
        for (int m = 0; m < data.N1(); m++) {
            double val = data.value(m, n);
            txt += QString(",%1").arg(val);
        }
        txt += "\n";
    }

    QString fname = QFileDialog::getSaveFileName(this, "Save spectrogram data", "", "*.csv");
    if (fname.isEmpty())
        return;
    if ((!fname.endsWith(".csv")) && (!fname.endsWith(".txt")))
        fname += ".csv";

    TextFile::write(fname, txt);
}

void SpectrogramViewCalculator::compute()
{

    TaskProgress task(TaskProgress::Calculate, "eeg-spectrogram");
    MountainProcessRunner X;
    QString processor_name = "eeg-spectrogram";
    X.setProcessorName(processor_name);

    QMap<QString, QVariant> params;
    params["timeseries"] = timeseries.makePath();
    params["time_resolution"] = time_resolution;
    X.setInputParameters(params);

    QString spectrogram_fname = X.makeOutputFilePath("spectrogram_out");

    X.runProcess();
    spectrogram.setPath(spectrogram_fname);
}

SpectrogramDataFactory::SpectrogramDataFactory(MVMainWindow* mw, QObject* parent)
    : MVAbstractViewFactory(mw, parent)
{
}

QString SpectrogramDataFactory::id() const
{
    return QStringLiteral("open-Spectrogram-data");
}

QString SpectrogramDataFactory::name() const
{
    return tr("Spectrogram Data");
}

QString SpectrogramDataFactory::title() const
{
    return tr("Spectrogram");
}

MVAbstractView* SpectrogramDataFactory::createView(MVAbstractContext* context)
{
    MVEEGContext* c = qobject_cast<MVEEGContext*>(context);
    Q_ASSERT(c);

    SpectrogramView* X = new SpectrogramView(context);
    //QList<int> ks = c->selectedClusters();
    //X->setLabelsToView(ks.toSet());
    return X;
}

QColor SpectrogramViewPrivate::get_pixel_color(double val)
{
    double pct = (val - m_window_min) / (m_window_max - m_window_min);
    double brightness_factor = exp(m_brightness_level / 100 * 5);
    pct *= brightness_factor;
    pct = qMin(1.0, qMax(0.0, pct));
    QColor gray = QColor(pct * 255, pct * 255, pct * 255);
    return gray;
}

Mda SpectrogramViewPrivate::get_spectrogram_data(int t1, int t2)
{
    MVEEGContext* c = qobject_cast<MVEEGContext*>(q->mvContext());
    Q_ASSERT(c);

    int M = m_spectrogram.N1();
    int N2 = t2 - t1 + 1;
    int K = m_spectrogram.N3();

    int ifreq_min, ifreq_max;
    SpectrogramViewPrivate::parse_freq_range(ifreq_min, ifreq_max, m_spectrogram_freq_range);

    Mda X;
    if (!m_spectrogram.readChunk(X, 0, 0, ifreq_min, m_spectrogram.N1(), m_spectrogram.N2(), ifreq_max - ifreq_min + 1)) {
        qWarning() << "Unable to read chunk of spectrogram";
    }

    Mda ret(M, N2);
    for (int n = 0; n < N2; n++) {
        for (int m = 0; m < M; m++) {
            double val = 0;
            for (int ii = ifreq_min; ii <= ifreq_max; ii++) {
                double val2 = X.value(m, t1 + n, ii - ifreq_min);
                if (val2 > val)
                    val = val2;
            }
            ret.setValue(val, m, n);
        }
    }

    return ret;
}

void SpectrogramViewPrivate::parse_freq_range(int& ret_min, int& ret_max, QString str)
{
    QStringList vals = str.split("-");
    if (vals.count() == 1) {
        ret_min = ret_max = (int)str.trimmed().toDouble();
    }
    else if (vals.count() == 2) {
        ret_min = (int)vals[0].trimmed().toDouble();
        ret_max = (int)vals[1].trimmed().toDouble();
    }
    else {
        qWarning() << "Error parsing freq range" << str;
        ret_min = ret_max = 0;
    }
}
