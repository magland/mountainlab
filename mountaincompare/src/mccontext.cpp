#include "mccontext.h"

#include <QMutex>
#include <mountainprocessrunner.h>
#include <taskprogress.h>
#include "mlcommon.h"

/*
class MatchFiringsCalculator {
public:
    //input
    QString mlproxy_url;
    DiskReadMda firings1;
    DiskReadMda firings2;

    //output
    DiskReadMda matched_firings;
    DiskReadMda confusion_matrix;
    DiskReadMda label_map;

    virtual void compute();
};
*/

class MCContextPrivate {
public:
    MCContext* q;
    MVContext* m_mv_context1;
    MVContext* m_mv_context2;
    MVContext* m_mv_context_intersect;

    QMutex m_mutex;
    //bool m_matched_firings_computed = false;
    //MatchFiringsCalculator m_match_firings_calculator;

    /////////////////////////////////////
    DiskReadMda m_firings1;
    DiskReadMda m_firings2;
    DiskReadMda m_matched_firings;
    int m_current_cluster2 = -1;
    QSet<int> m_selected_clusters2;
    DiskReadMda m_confusion_matrix;
    DiskReadMda m_label_map;

    void setup_mv_contexts();
};

MCContext::MCContext()
{
    d = new MCContextPrivate;
    d->q = this;
    d->m_mv_context1 = new MVContext;
    d->m_mv_context2 = new MVContext;

    d->setup_mv_contexts();
}

MCContext::~MCContext()
{
    delete d->m_mv_context1;
    delete d->m_mv_context2;
    delete d;
}

MVContext* MCContext::mvContext1()
{
    return d->m_mv_context1;
}

MVContext* MCContext::mvContext2()
{
    return d->m_mv_context2;
}

void MCContext::setFromMV2FileObject(QJsonObject obj)
{
    d->m_mv_context1->setFromMV2FileObject(obj);
    d->m_mv_context2->setFromMV2FileObject(obj);
    if (obj.contains("firings1")) {
        this->setFirings1(DiskReadMda(obj["firings1"].toObject()));
    }
    if (obj.contains("firings2")) {
        this->setFirings2(DiskReadMda(obj["firings2"].toObject()));
    }
}

QJsonObject MCContext::toMV2FileObject() const
{
    QJsonObject obj = d->m_mv_context1->toMV2FileObject();
    obj["firings1"] = d->m_mv_context1->firings().toPrvObject();
    obj["firings2"] = d->m_mv_context2->firings().toPrvObject();
    return obj;
}

MCEvent MCContext::currentEvent() const
{
    MCEvent ret;
    if (d->m_mv_context1->currentEvent().label >= 0) {
        ret.firings_num = 1;
        ret.time = d->m_mv_context1->currentEvent().time;
        ret.label = d->m_mv_context1->currentEvent().label;
    }
    else {
        ret.firings_num = 2;
        ret.time = d->m_mv_context2->currentEvent().time;
        ret.label = d->m_mv_context2->currentEvent().label;
    }
    return ret;
}

MCCluster MCContext::currentCluster() const
{
    MCCluster ret;
    if (d->m_mv_context2->currentEvent().label >= 0) {
        ret.firings_num = 1;
        ret.num = d->m_mv_context1->currentCluster();
    }
    else {
        ret.firings_num = 2;
        ret.num = d->m_mv_context2->currentCluster();
    }
    return ret;
}

QList<MCCluster> MCContext::selectedClusters() const
{
    QList<MCCluster> ret;
    QList<int> C1 = d->m_mv_context1->selectedClusters();
    QList<int> C2 = d->m_mv_context2->selectedClusters();
    foreach (int c, C1) {
        MCCluster X;
        X.num = c;
        X.firings_num = 1;
        ret << X;
    }
    foreach (int c, C2) {
        MCCluster X;
        X.num = c;
        X.firings_num = 2;
        ret << X;
    }
    return ret;
}

double MCContext::currentTimepoint() const
{
    return d->m_mv_context1->currentTimepoint();
}

MVRange MCContext::currentTimeRange() const
{
    return d->m_mv_context1->currentTimeRange();
}

void MCContext::setCurrentEvent(const MCEvent& evt)
{
    MVEvent evt2;
    evt2.label = evt.label;
    evt2.time = evt.time;
    if (evt.firings_num == 1) {
        d->m_mv_context1->setCurrentEvent(evt2);
        d->m_mv_context2->setCurrentEvent(MVEvent());
    }
    else if (evt.firings_num == 2) {
        d->m_mv_context2->setCurrentEvent(evt2);
        d->m_mv_context1->setCurrentEvent(MVEvent());
    }
    else {
        d->m_mv_context1->setCurrentEvent(MVEvent());
        d->m_mv_context2->setCurrentEvent(MVEvent());
    }
}

void MCContext::setCurrentCluster(const MCCluster& C)
{
    if (C.firings_num == 1) {
        d->m_mv_context1->setCurrentCluster(C.num);
        d->m_mv_context2->setCurrentCluster(-1);
    }
    else if (C.firings_num == 2) {
        d->m_mv_context2->setCurrentCluster(C.num);
        d->m_mv_context1->setCurrentCluster(-1);
    }
    else {
        d->m_mv_context1->setCurrentCluster(-1);
        d->m_mv_context2->setCurrentCluster(-1);
    }
}

bool matches(const QSet<MCCluster>& clusters1, const QSet<MCCluster>& clusters2)
{
    if (clusters1.count() != clusters2.count())
        return false;
    foreach (MCCluster C, clusters1) {
        if (!clusters2.contains(C))
            return false;
    }
    return true;
}

void MCContext::setSelectedClusters(const QList<MCCluster>& Cs)
{
    if (matches(Cs.toSet(), this->selectedClusters().toSet()))
        return;
    QList<int> C1s, C2s;
    foreach (MCCluster C, Cs) {
        if (C.firings_num == 1) {
            C1s << C.num;
        }
        else if (C.firings_num == 2) {
            C2s << C.num;
        }
    }
    d->m_mv_context1->setSelectedClusters(C1s);
    d->m_mv_context2->setSelectedClusters(C2s);
}

void MCContext::setCurrentTimepoint(double tp)
{
    d->m_mv_context1->setCurrentTimepoint(tp);
}

void MCContext::setCurrentTimeRange(const MVRange& range)
{
    d->m_mv_context1->setCurrentTimeRange(range);
}

void MCContext::clickCluster(MCCluster C, Qt::KeyboardModifiers modifiers)
{
    if (C.firings_num == 1) {
        d->m_mv_context1->clickCluster(C.num, modifiers);
        if (!(modifiers & Qt::ControlModifier)) {
            d->m_mv_context2->setSelectedClusters(QList<int>());
        }
    }
    else if (C.firings_num == 2) {
        d->m_mv_context2->clickCluster(C.num, modifiers);
        if (!(modifiers & Qt::ControlModifier)) {
            d->m_mv_context1->setSelectedClusters(QList<int>());
        }
    }
}

QColor MCContext::clusterColor(int k) const
{
    return d->m_mv_context1->clusterColor(k);
}

QColor MCContext::channelColor(int m) const
{
    return d->m_mv_context1->channelColor(m);
}

QColor MCContext::color(QString name, QColor default_color) const
{
    return d->m_mv_context1->color(name, default_color);
}

QMap<QString, QColor> MCContext::colors() const
{
    return d->m_mv_context1->colors();
}

QList<QColor> MCContext::channelColors() const
{
    return d->m_mv_context1->channelColors();
}

QList<QColor> MCContext::clusterColors() const
{
    return d->m_mv_context1->clusterColors();
}

void MCContext::setClusterColors(const QList<QColor>& colors)
{
    d->m_mv_context1->setClusterColors(colors);
    d->m_mv_context2->setClusterColors(colors);
}

void MCContext::setChannelColors(const QList<QColor>& colors)
{
    d->m_mv_context1->setChannelColors(colors);
    d->m_mv_context2->setChannelColors(colors);
}

void MCContext::setColors(const QMap<QString, QColor>& colors)
{
    d->m_mv_context1->setColors(colors);
    d->m_mv_context2->setColors(colors);
}

DiskReadMda32 MCContext::currentTimeseries() const
{
    return d->m_mv_context1->currentTimeseries();
}

DiskReadMda32 MCContext::timeseries(QString name) const
{
    return d->m_mv_context1->timeseries(name);
}

QString MCContext::currentTimeseriesName() const
{
    return d->m_mv_context1->currentTimeseriesName();
}

QStringList MCContext::timeseriesNames() const
{
    return d->m_mv_context1->timeseriesNames();
}

void MCContext::addTimeseries(QString name, DiskReadMda32 timeseries)
{
    d->m_mv_context1->addTimeseries(name, timeseries);
    d->m_mv_context2->addTimeseries(name, timeseries);
}

void MCContext::setCurrentTimeseriesName(QString name)
{
    d->m_mv_context1->setCurrentTimeseriesName(name);
    d->m_mv_context2->setCurrentTimeseriesName(name);
}

DiskReadMda MCContext::firings1()
{
    return d->m_mv_context1->firings();
}

DiskReadMda MCContext::firings2()
{
    return d->m_mv_context2->firings();
}

DiskReadMda MCContext::matchedFirings()
{
    return d->m_matched_firings;
}

void MCContext::setFirings1(const DiskReadMda& F)
{
    d->m_mv_context1->setFirings(F);
}

void MCContext::setFirings2(const DiskReadMda& F)
{
    d->m_mv_context2->setFirings(F);
}

void MCContext::setConfusionMatrix(const DiskReadMda& CM)
{
    d->m_confusion_matrix = CM;
}

void MCContext::setLabelMap(const DiskReadMda& LM)
{
    d->m_label_map = LM;
}

void MCContext::setMatchedFirings(const DiskReadMda& MF)
{
    d->m_matched_firings = MF;
}

void MCContext::setClustersToForceShow(const QList<int>& list)
{
    d->m_mv_context1->setClustersToForceShow(list);
    d->m_mv_context2->setClustersToForceShow(list);
}

/*
void MCContext::computeMatchedFirings()
{
    {
        QMutexLocker locker(&d->m_mutex);
        if (d->m_matched_firings_computed)
            return;
        d->m_matched_firings_computed = true;
    }

    d->m_match_firings_calculator.firings1 = this->firings1();
    d->m_match_firings_calculator.firings2 = this->firings2();
    d->m_match_firings_calculator.compute();
    d->m_confusion_matrix = d->m_match_firings_calculator.confusion_matrix;
    d->m_matched_firings = d->m_match_firings_calculator.matched_firings;
    d->m_label_map = d->m_match_firings_calculator.label_map;
}
*/

double MCContext::sampleRate() const
{
    return d->m_mv_context1->sampleRate();
}

void MCContext::setSampleRate(double sample_rate)
{
    d->m_mv_context1->setSampleRate(sample_rate);
    d->m_mv_context2->setSampleRate(sample_rate);
}

Mda MCContext::confusionMatrix() const
{
    Mda ret;
    if (!d->m_confusion_matrix.readChunk(ret, 0, 0, d->m_confusion_matrix.N1(), d->m_confusion_matrix.N2())) {
        qWarning() << "Unable to read chunk of confusion matrix";
    }
    return ret;
}

QList<int> MCContext::labelMap() const
{
    QList<int> ret;
    for (int i = 0; i < d->m_label_map.totalSize(); i++) {
        ret << d->m_label_map.value(i);
    }
    return ret;
}

void MCContext::slot_context_current_timepoint_changed()
{
    MVContext* context = qobject_cast<MVContext*>(sender());
    if (context == d->m_mv_context1) {
        d->m_mv_context2->setCurrentTimepoint(d->m_mv_context1->currentTimepoint());
    }
    else if (context == d->m_mv_context2) {
        d->m_mv_context1->setCurrentTimepoint(d->m_mv_context2->currentTimepoint());
    }
}

void MCContext::slot_context_current_time_range_changed()
{
    MVContext* context = qobject_cast<MVContext*>(sender());
    if (context == d->m_mv_context1) {
        d->m_mv_context2->setCurrentTimeRange(d->m_mv_context1->currentTimeRange());
    }
    else if (context == d->m_mv_context2) {
        d->m_mv_context1->setCurrentTimeRange(d->m_mv_context2->currentTimeRange());
    }
}

void MCContext::slot_context_current_timeseries_changed()
{
    MVContext* context = qobject_cast<MVContext*>(sender());
    if (context == d->m_mv_context1) {
        d->m_mv_context2->addTimeseries(d->m_mv_context1->currentTimeseriesName(), d->m_mv_context1->currentTimeseries());
        d->m_mv_context2->setCurrentTimeseriesName(d->m_mv_context1->currentTimeseriesName());
    }
    else if (context == d->m_mv_context2) {
        d->m_mv_context1->addTimeseries(d->m_mv_context2->currentTimeseriesName(), d->m_mv_context2->currentTimeseries());
        d->m_mv_context1->setCurrentTimeseriesName(d->m_mv_context2->currentTimeseriesName());
    }
}

/*
void MCContext::computeMergedFirings()
{
    {
        QMutexLocker locker(&d->m_mutex);
        if (d->m_merged_firings_computed)
            return;
    }

    d->m_calculator.mlproxy_url = this->mlProxyUrl();
    d->m_calculator.firings1 = this->firings1();
    d->m_calculator.firings2 = this->firings2();
    d->m_calculator.compute();
    d->m_confusion_matrix = d->m_calculator.confusion_matrix;
    d->m_label_map = d->m_calculator.label_map;
    d->m_firings_merged = d->m_calculator.firings_merged;

    {
        QMutexLocker locker(&d->m_mutex);
        d->m_merged_firings_computed = true;
    }
}

Mda MCContext::confusionMatrix() const
{
    Mda ret;
    d->m_confusion_matrix.readChunk(ret, 0, 0, d->m_confusion_matrix.N1(), d->m_confusion_matrix.N2());
    return ret;
}

QList<int> MCContext::optimalLabelMap() const
{
    QList<int> ret;
    for (int i = 0; i < d->m_label_map.totalSize(); i++) {
        ret << d->m_label_map.value(i);
    }
    return ret;
}

DiskReadMda MCContext::firingsMerged() const
{
    return d->m_firings_merged;
}

DiskReadMda MCContext::firings1()
{
    return d->m_firings1;
}

DiskReadMda MCContext::firings2()
{
    return d->m_firings2;
}

int MCContext::currentCluster2() const
{
    return d->m_current_cluster2;
}

QList<int> MCContext::selectedClusters2() const
{
    QList<int> ret = d->m_selected_clusters2.toList();
    qSort(ret);
    return ret;
}

void MCContext::setFirings1(const DiskReadMda& F)
{
    d->m_firings1 = F;
    emit firings1Changed();
}

void MCContext::setFirings2(const DiskReadMda& F)
{
    d->m_firings2 = F;
    emit firings2Changed();
}

void MCContext::setCurrentCluster2(int k)
{
    if (d->m_current_cluster2 == k)
        return;
    d->m_current_cluster2 = k;
    emit currentCluster2Changed();
}

void MCContext::setSelectedClusters2(const QList<int>& list)
{
    if (list.toSet() == d->m_selected_clusters2)
        return;
    d->m_selected_clusters2 = list.toSet();
    emit selectedClusters2Changed();
}

void MCContext::clickCluster2(int k, Qt::KeyboardModifiers modifiers)
{
    if (modifiers & Qt::ControlModifier) {
        if (k < 0)
            return;
        if (d->m_selected_clusters2.contains(k)) {
            QList<int> tmp = d->m_selected_clusters2.toList();
            tmp.removeAll(k);
            this->setSelectedClusters2(tmp);
        }
        else {
            if (k >= 0) {
                QList<int> tmp = d->m_selected_clusters2.toList();
                tmp << k;
                this->setSelectedClusters2(tmp);
            }
        }
    }
    else {
        //this->setSelectedClusterPairs(QSet<ClusterPair>());
        this->setSelectedClusters2(QList<int>());
        this->setCurrentCluster2(k);
    }
}

void MergeFiringsCalculator::compute()
{
    TaskProgress task(TaskProgress::Calculate, "Compute confusion matrix");
    task.setProgress(0.1);

    MountainProcessRunner MPR;
    MPR.setProcessorName("merge_firings");

    QMap<QString, QVariant> params;
    params["firings1"] = firings1.makePath();
    params["firings2"] = firings2.makePath();
    params["max_matching_offset"] = 30;
    MPR.setInputParameters(params);
    MPR.setMLProxyUrl(mlproxy_url);

    task.log() << "Firings 1/2 dimensions" << firings1.N1() << firings1.N2() << firings2.N1() << firings2.N2();

    QString firings_merged_path = MPR.makeOutputFilePath("firings_merged");
    QString label_map_path = MPR.makeOutputFilePath("label_map");
    QString confusion_matrix_path = MPR.makeOutputFilePath("confusion_matrix");

    MPR.runProcess();

    if (MLUtil::threadInterruptRequested()) {
        task.error(QString("Halted while running process."));
        return;
    }

    firings_merged.setPath(firings_merged_path);
    confusion_matrix.setPath(confusion_matrix_path);
    label_map.setPath(label_map_path);
}

void MCContextPrivate::setup_mv_contexts()
{
    //finish
}
*/

void MCContextPrivate::setup_mv_contexts()
{
    QObject::connect(m_mv_context1, SIGNAL(currentTimeseriesChanged()), q, SIGNAL(currentTimeseriesChanged()));
    QObject::connect(m_mv_context1, SIGNAL(timeseriesNamesChanged()), q, SIGNAL(timeseriesNamesChanged()));
    QObject::connect(m_mv_context1, SIGNAL(currentTimepointChanged()), q, SIGNAL(currentTimepointChanged()));
    QObject::connect(m_mv_context1, SIGNAL(currentTimeRangeChanged()), q, SIGNAL(currentTimeRangeChanged()));
    QObject::connect(m_mv_context1, SIGNAL(clusterColorsChanged(QList<QColor>)), q, SIGNAL(clusterColorsChanged(QList<QColor>)));

    QObject::connect(m_mv_context1, SIGNAL(firingsChanged()), q, SIGNAL(firingsChanged()));
    QObject::connect(m_mv_context2, SIGNAL(firingsChanged()), q, SIGNAL(firingsChanged()));
    QObject::connect(m_mv_context1, SIGNAL(currentEventChanged()), q, SIGNAL(currentEventChanged()));
    QObject::connect(m_mv_context2, SIGNAL(currentEventChanged()), q, SIGNAL(currentEventChanged()));
    QObject::connect(m_mv_context1, SIGNAL(currentClusterChanged()), q, SIGNAL(currentClusterChanged()));
    QObject::connect(m_mv_context2, SIGNAL(currentClusterChanged()), q, SIGNAL(currentClusterChanged()));
    QObject::connect(m_mv_context1, SIGNAL(selectedClustersChanged()), q, SIGNAL(selectedClustersChanged()));
    QObject::connect(m_mv_context2, SIGNAL(selectedClustersChanged()), q, SIGNAL(selectedClustersChanged()));

    QObject::connect(m_mv_context1, SIGNAL(currentTimepointChanged()), q, SLOT(slot_context_current_timepoint_changed()));
    QObject::connect(m_mv_context2, SIGNAL(currentTimepointChanged()), q, SLOT(slot_context_current_timepoint_changed()));
    QObject::connect(m_mv_context1, SIGNAL(currentTimeRangeChanged()), q, SLOT(slot_context_current_time_range_changed()));
    QObject::connect(m_mv_context2, SIGNAL(currentTimeRangeChanged()), q, SLOT(slot_context_current_time_range_changed()));
    QObject::connect(m_mv_context1, SIGNAL(currentTimeseriesChanged()), q, SLOT(slot_context_current_timeseries_changed()), Qt::QueuedConnection);
    QObject::connect(m_mv_context2, SIGNAL(currentTimeseriesChanged()), q, SLOT(slot_context_current_timeseries_changed()), Qt::QueuedConnection);
}

uint qHash(const MCCluster& C)
{
    return qHash(C.toString());
}

/*
void MatchFiringsCalculator::compute()
{
    TaskProgress task(TaskProgress::Calculate, "Compute confusion matrix");
    task.setProgress(0.1);

    MountainProcessRunner MPR;
    MPR.setProcessorName("mountainsort.confusion_matrix");

    QMap<QString, QVariant> params;
    params["firings1"] = firings1.makePath();
    params["firings2"] = firings2.makePath();
    params["max_matching_offset"] = 30;
    MPR.setInputParameters(params);

    task.log() << "Firings 1/2 dimensions" << firings1.N1() << firings1.N2() << firings2.N1() << firings2.N2();

    QString confusion_matrix_path = MPR.makeOutputFilePath("confusion_matrix_out");
    QString matched_firings_path = MPR.makeOutputFilePath("matched_firings_out");
    QString label_map_path = MPR.makeOutputFilePath("label_map_out");

    MPR.runProcess();

    if (MLUtil::threadInterruptRequested()) {
        task.error(QString("Halted while running process."));
        return;
    }

    matched_firings.setPath(matched_firings_path);
    confusion_matrix.setPath(confusion_matrix_path);
    label_map.setPath(label_map_path);
}
*/
