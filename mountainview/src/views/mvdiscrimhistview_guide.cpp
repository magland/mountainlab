#include "mountainprocessrunner.h"
#include "mvdiscrimhistview_guide.h"

#include <QGridLayout>
#include <taskprogress.h>
#include "mlcommon.h"
#include "histogramview.h"
#include <QWheelEvent>
#include "mvmainwindow.h"
#include "actionfactory.h"

struct DiscrimHistogram {
    int k1, k2;
    QVector<double> data1, data2;
};

class MVDiscrimHistViewGuideComputer {
public:
    //input
    QString mlproxy_url;
    DiskReadMda32 timeseries;
    DiskReadMda firings;
    int num_histograms; //old version
    QSet<int> clusters_to_exclude; //old version
    QList<int> cluster_numbers;
    QString method; //old version

    //output
    QList<DiscrimHistogram> histograms;

    void compute();
};

class MVDiscrimHistViewGuidePrivate {
public:
    MVDiscrimHistViewGuide* q;

    int m_num_histograms = 20;

    MVDiscrimHistViewGuideComputer m_computer;
    QList<DiscrimHistogram> m_histograms;

    double m_zoom_factor = 1;

    void set_views();
    QSet<int> get_clusters_to_exclude();
};

MVDiscrimHistViewGuide::MVDiscrimHistViewGuide(MVAbstractContext* context)
    : MVHistogramGrid(context)
{
    d = new MVDiscrimHistViewGuidePrivate;
    d->q = this;

    this->recalculateOn(context, SIGNAL(currentTimeseriesChanged()));
    this->recalculateOn(context, SIGNAL(firingsChanged()), false);
    this->recalculateOn(context, SIGNAL(clusterMergeChanged()), true);
    this->recalculateOn(context, SIGNAL(clusterVisibilityChanged()), true);
    this->recalculateOn(context, SIGNAL(viewMergedChanged()), true);
    this->recalculateOnOptionChanged("discrim_hist_method", true);

    ActionFactory::addToToolbar(ActionFactory::ActionType::ZoomInHorizontal, this, SLOT(slot_zoom_in_horizontal()));
    ActionFactory::addToToolbar(ActionFactory::ActionType::ZoomOutHorizontal, this, SLOT(slot_zoom_out_horizontal()));
    ActionFactory::addToToolbar(ActionFactory::ActionType::PanLeft, this, SLOT(slot_pan_left()));
    ActionFactory::addToToolbar(ActionFactory::ActionType::PanRight, this, SLOT(slot_pan_right()));

    this->recalculate();
}

MVDiscrimHistViewGuide::~MVDiscrimHistViewGuide()
{
    this->stopCalculation();
    delete d;
}

void MVDiscrimHistViewGuide::setNumHistograms(int num)
{
    d->m_num_histograms = num;
}

void MVDiscrimHistViewGuide::prepareCalculation()
{
    MVContext* c = qobject_cast<MVContext*>(mvContext());
    Q_ASSERT(c);

    d->m_computer.mlproxy_url = c->mlProxyUrl();
    d->m_computer.timeseries = c->currentTimeseries();
    d->m_computer.firings = c->firings();
    d->m_computer.num_histograms = d->m_num_histograms;
    d->m_computer.clusters_to_exclude = d->get_clusters_to_exclude();
    d->m_computer.cluster_numbers = c->selectedClusters();
    d->m_computer.method = c->option("discrim_hist_method").toString();
}

void MVDiscrimHistViewGuide::runCalculation()
{
    d->m_computer.compute();
}

double compute_min3(const QList<DiscrimHistogram>& data0)
{
    double ret = 0;
    for (int i = 0; i < data0.count(); i++) {
        QVector<double> tmp = data0[i].data1;
        tmp.append(data0[i].data2);
        for (int j = 0; j < tmp.count(); j++) {
            if (tmp[j] < ret)
                ret = tmp[j];
        }
    }
    return ret;
}

double max3(const QList<DiscrimHistogram>& data0)
{
    double ret = 0;
    for (int i = 0; i < data0.count(); i++) {
        QVector<double> tmp = data0[i].data1;
        tmp.append(data0[i].data2);
        for (int j = 0; j < tmp.count(); j++) {
            if (tmp[j] > ret)
                ret = tmp[j];
        }
    }
    return ret;
}

void MVDiscrimHistViewGuide::onCalculationFinished()
{
    d->m_histograms = d->m_computer.histograms;

    d->set_views();
}

void MVDiscrimHistViewGuideComputer::compute()
{
    TaskProgress task(TaskProgress::Calculate, QString("Discrim Histograms"));

    histograms.clear();

    MountainProcessRunner MPR;
    MPR.setMLProxyUrl(mlproxy_url);

    /*
    MPR.setProcessorName("mv_discrimhist_guide");
    QMap<QString, QVariant> params;
    params["timeseries"] = timeseries.makePath();
    params["firings"] = firings.makePath();
    params["num_histograms"] = num_histograms;
    params["clusters_to_exclude"] = MLUtil::intListToStringList(clusters_to_exclude.toList()).join(",");
    params["method"]=method;
    */

    MPR.setProcessorName("mv_discrimhist_guide2");
    QMap<QString, QVariant> params;
    params["timeseries"] = timeseries.makePath();
    params["firings"] = firings.makePath();
    params["clip_size"] = 50;
    params["add_noise_level"] = 1;
    params["cluster_numbers"] = MLUtil::intListToStringList(cluster_numbers).join(",");
    params["max_comparisons_per_cluster"] = 5;

    MPR.setInputParameters(params);

    QString output_path = MPR.makeOutputFilePath("output");

    MPR.runProcess();

    DiskReadMda output(output_path);
    //output.setRemoteDataType("float32");

    QMap<QString, DiscrimHistogram*> hist_lookup;

    int LL = output.N2();
    for (int i = 0; i < LL; i++) {
        int k1 = output.value(0, i);
        int k2 = output.value(1, i);
        int k0 = output.value(2, i);
        double val = output.value(3, i);
        QString code = QString("%1:%2").arg(k1).arg(k2);
        DiscrimHistogram* HH = hist_lookup.value(code);
        if (!HH) {
            DiscrimHistogram H;
            this->histograms << H;
            DiscrimHistogram* ptr = &this->histograms[this->histograms.count() - 1];
            ptr->k1 = k1;
            ptr->k2 = k2;
            hist_lookup[code] = ptr;
            HH = hist_lookup.value(code);
        }
        {
            if (k0 == k1)
                HH->data1 << val;
            else
                HH->data2 << val;
        }
    }
}

void MVDiscrimHistViewGuide::wheelEvent(QWheelEvent* evt)
{
    Q_UNUSED(evt)
    /*
    if (evt->delta() > 0) {
        slot_zoom_in();
    }
    else if (evt->delta() < 0) {
        slot_zoom_out();
    }
    */
}

void MVDiscrimHistViewGuide::slot_zoom_in_horizontal(double zoom_factor)
{
    QList<HistogramView*> views = this->histogramViews(); //inherited
    for (int i = 0; i < views.count(); i++) {
        views[i]->setXRange(views[i]->xRange() * (1.0 / zoom_factor));
    }
}

void MVDiscrimHistViewGuide::slot_zoom_out_horizontal(double factor)
{
    slot_zoom_in_horizontal(1 / factor);
}

void MVDiscrimHistViewGuide::slot_pan_left(double units)
{
    QList<HistogramView*> views = this->histogramViews(); //inherited
    for (int i = 0; i < views.count(); i++) {
        MVRange range = views[i]->xRange();
        if (range.range()) {
            range = range + units * range.range();
        }
        views[i]->setXRange(range);
    }
}

void MVDiscrimHistViewGuide::slot_pan_right(double units)
{
    slot_pan_left(-units);
}

void MVDiscrimHistViewGuidePrivate::set_views()
{
    MVContext* c = qobject_cast<MVContext*>(q->mvContext());
    Q_ASSERT(c);

    double bin_min = compute_min3(m_histograms);
    double bin_max = max3(m_histograms);
    double max00 = qMax(qAbs(bin_min), qAbs(bin_max));

    int num_bins = 500; //how to choose this?

    QList<HistogramView*> views;
    for (int ii = 0; ii < m_histograms.count(); ii++) {
        int k1 = m_histograms[ii].k1;
        int k2 = m_histograms[ii].k2;
        //if (q->c->clusterIsVisible(k1)) {
        {
            HistogramView* HV = new HistogramView;
            QVector<double> tmp = m_histograms[ii].data1;
            tmp.append(m_histograms[ii].data2);
            HV->setData(tmp);
            HV->setSecondData(m_histograms[ii].data2);
            HV->setColors(c->colors());
            //HV->autoSetBins(50);
            HV->setBinInfo(bin_min, bin_max, num_bins);
            QString title0 = QString("%1/%2").arg(k1).arg(k2);
            HV->setTitle(title0);
            //HV->setDrawVerticalAxisAtZero(true);
            HV->setDrawVerticalAxisAtZero(false);
            HV->setXRange(MVRange(-max00, max00));
            HV->autoCenterXRange();
            HV->setProperty("k1", k1);
            HV->setProperty("k2", k2);
            views << HV;
        }
    }

    q->setHistogramViews(views); //base class
    q->slot_zoom_in_horizontal(2.5); //give it a nice zoom in to start
}

QSet<int> MVDiscrimHistViewGuidePrivate::get_clusters_to_exclude()
{
    MVContext* c = qobject_cast<MVContext*>(q->mvContext());
    Q_ASSERT(c);

    QSet<int> ret;
    QList<int> keys = c->clusterAttributesKeys();
    foreach (int key, keys) {
        if (c->clusterTags(key).contains("rejected"))
            ret.insert(key);
    }
    return ret;
}

MVDiscrimHistGuideFactory::MVDiscrimHistGuideFactory(MVMainWindow* mw, QObject* parent)
    : MVAbstractViewFactory(mw, parent)
{
}

QString MVDiscrimHistGuideFactory::id() const
{
    return QStringLiteral("open-discrim-histograms-guide");
}

QString MVDiscrimHistGuideFactory::name() const
{
    return tr("Discrim Histograms Guide");
}

QString MVDiscrimHistGuideFactory::title() const
{
    return tr("Discrim");
}

MVAbstractView* MVDiscrimHistGuideFactory::createView(MVAbstractContext* context)
{
    MVDiscrimHistViewGuide* X = new MVDiscrimHistViewGuide(context);
    X->setPreferredHistogramWidth(200);
    X->setNumHistograms(80);
    return X;
}

bool MVDiscrimHistGuideFactory::isEnabled(MVAbstractContext* context) const
{
    Q_UNUSED(context)
    return true;
}
