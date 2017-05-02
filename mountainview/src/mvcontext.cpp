/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
** Created: 5/23/2016
*******************************************************/

#include "mvcontext.h"
#include <QAction>
#include <QDebug>
#include <QJsonDocument>
#include <QThread>
#include <cachemanager.h>
#include <taskprogress.h>
#include "mlcommon.h"
#include "mountainprocessrunner.h"
#include "get_sort_indices.h"

struct TimeseriesStruct {
    QString name;
    DiskReadMda32 data;
};

class MVContextPrivate {
public:
    MVContext* q;
    QMap<int, QJsonObject> m_cluster_attributes;
    //needed to do it this way rather than QMap<ClusterPair,QJsonObject>
    //because of a terrible bug that took me hours to conclude that it
    //could not be solve, and is PROBABLY a bug with Qt!!
    QMap<int, QMap<int, QJsonObject> > m_cluster_pair_attributes;
    MVEvent m_current_event;
    int m_current_cluster;
    QList<int> m_selected_clusters;
    double m_current_timepoint;
    MVRange m_current_time_range = MVRange(0, 0); //(0,0) triggers automatic calculation
    QList<QColor> m_cluster_colors;
    QList<QColor> m_channel_colors;
    QMap<QString, TimeseriesStruct> m_timeseries;
    QString m_current_timeseries_name;
    DiskReadMda m_firings;
    DiskReadMda m_firings_subset;
    double m_sample_rate = 0;
    QString m_mlproxy_url; //this is to disappear
    QMap<QString, QColor> m_colors;
    ClusterVisibilityRule m_visibility_rule;
    QJsonObject m_original_object;
    QSet<ClusterPair> m_selected_cluster_pairs;
    ElectrodeGeometry m_electrode_geometry;
    QSet<int> m_clusters_subset;
    bool m_view_merged;
    QSet<int> m_visible_channels;
    QList<double> m_cluster_order_scores;
    QString m_cluster_order_scores_name;
    QSet<int> m_clusters_to_force_show;

    void update_current_and_selected_clusters_according_to_merged();
    void set_default_options();
};

QJsonArray intlist_to_json_array(const QList<int>& X);
QJsonArray doublelist_to_json_array(const QList<double>& X);
QList<int> json_array_to_intlist(const QJsonArray& X);
QList<double> json_array_to_doublelist(const QJsonArray& X);
QJsonArray strlist_to_json_array(const QList<QString>& X);
QList<QString> json_array_to_strlist(const QJsonArray& X);

MVContext::MVContext()
{
    d = new MVContextPrivate;
    d->q = this;
    d->m_current_cluster = 0;
    d->m_current_timepoint = 0;

    d->set_default_options();

    // default colors
    d->m_colors["background"] = QColor(240, 240, 240);
    d->m_colors["frame1"] = QColor(245, 245, 245);
    d->m_colors["info_text"] = QColor(80, 80, 80);
    d->m_colors["view_background"] = QColor(243, 243, 247);
    d->m_colors["view_background_highlighted"] = QColor(210, 230, 250);
    d->m_colors["view_background_selected"] = QColor(220, 240, 250);
    d->m_colors["view_background_hovered"] = QColor(243, 243, 237);
    d->m_colors["view_frame"] = QColor(200, 190, 190);
    d->m_colors["view_frame_selected"] = QColor(100, 80, 80);
    d->m_colors["divider_line"] = QColor(255, 100, 150);
    d->m_colors["calculation-in-progress"] = QColor(130, 130, 140, 50);
    d->m_colors["cluster_label"] = Qt::darkGray;
    d->m_colors["firing_rate_disk"] = Qt::lightGray;
}

MVContext::~MVContext()
{
    delete d;
}

void MVContext::clear()
{
    d->m_cluster_attributes.clear();
    d->m_cluster_pair_attributes.clear();
    d->m_timeseries.clear();
    d->m_firings = DiskReadMda();
    clearOptions();
    d->set_default_options();
}

QJsonObject cluster_attributes_to_object(const QMap<int, QJsonObject>& map)
{
    QJsonObject X;
    QList<int> keys = map.keys();
    foreach (int key, keys) {
        X[QString("%1").arg(key)] = map[key];
    }
    return X;
}

QMap<int, QJsonObject> object_to_cluster_attributes(QJsonObject X)
{
    QMap<int, QJsonObject> ret;
    QStringList keys = X.keys();
    foreach (QString key, keys) {
        ret[key.toInt()] = X[key].toObject();
    }
    return ret;
}

QJsonObject cluster_pair_attributes_to_object(const QMap<int, QMap<int, QJsonObject> >& map)
{
    QJsonObject X;
    QList<int> keys1 = map.keys();
    foreach (int k1, keys1) {
        QList<int> keys2 = map[k1].keys();
        foreach (int k2, keys2) {
            X[QString("%1").arg(ClusterPair(k1, k2).toString())] = map[k1][k2];
        }
    }
    return X;
}

QMap<int, QMap<int, QJsonObject> > object_to_cluster_pair_attributes(QJsonObject X)
{
    QMap<int, QMap<int, QJsonObject> > ret;
    QStringList keys = X.keys();
    foreach (QString key, keys) {
        ClusterPair pair = ClusterPair::fromString(key);
        ret[pair.k1()][pair.k2()] = X[key].toObject();
    }
    return ret;
}

QJsonObject timeseries_map_to_object(const QMap<QString, TimeseriesStruct>& TT)
{
    QJsonObject ret;
    QStringList keys = TT.keys();
    foreach (QString key, keys) {
        QJsonObject obj;
        obj["data"] = TT[key].data.makePath();
        obj["name"] = TT[key].name;
        ret[key] = obj;
    }
    return ret;
}

QJsonObject timeseries_map_to_object_for_mv2(const QMap<QString, TimeseriesStruct>& TT)
{
    QJsonObject ret;
    QStringList keys = TT.keys();
    foreach (QString key, keys) {
        QJsonObject obj;
        obj["data"] = TT[key].data.toPrvObject();
        obj["name"] = TT[key].name;
        ret[key] = obj;
    }
    return ret;
}

QMap<QString, TimeseriesStruct> object_to_timeseries_map(QJsonObject X)
{
    QMap<QString, TimeseriesStruct> ret;
    QStringList keys = X.keys();
    foreach (QString key, keys) {
        QJsonObject obj = X[key].toObject();
        TimeseriesStruct A;
        A.data = DiskReadMda32(obj["data"].toString());
        A.name = obj["name"].toString();
        ret[key] = A;
    }
    return ret;
}

QMap<QString, TimeseriesStruct> object_to_timeseries_map_for_mv2(QJsonObject X)
{
    QMap<QString, TimeseriesStruct> ret;
    QStringList keys = X.keys();
    foreach (QString key, keys) {
        QJsonObject obj = X[key].toObject();
        TimeseriesStruct A;
        A.data = DiskReadMda32(obj["data"].toObject());
        A.name = obj["name"].toString();
        ret[key] = A;
    }
    return ret;
}

QJsonObject MVContext::toMVFileObject() const
{
    // old
    QJsonObject X = d->m_original_object;
    X["cluster_attributes"] = cluster_attributes_to_object(d->m_cluster_attributes);
    X["cluster_pair_attributes"] = cluster_pair_attributes_to_object(d->m_cluster_pair_attributes);
    X["timeseries"] = timeseries_map_to_object(d->m_timeseries);
    X["current_timeseries_name"] = d->m_current_timeseries_name;
    X["firings"] = d->m_firings.makePath();
    X["samplerate"] = d->m_sample_rate;
    //X["options"] = QJsonObject::fromVariantMap(d->m_options);
    X["mlproxy_url"] = d->m_mlproxy_url;
    X["visibility_rule"] = d->m_visibility_rule.toJsonObject();
    X["electrode_geometry"] = d->m_electrode_geometry.toJsonObject();
    X["clusters_subset"] = intlist_to_json_array(d->m_clusters_subset.toList());
    X["view_merged"] = d->m_view_merged;
    X["visible_channels"] = intlist_to_json_array(d->m_visible_channels.toList());
    return X;
}

QJsonObject MVContext::toMV2FileObject() const
{
    QJsonObject X = d->m_original_object;
    X["cluster_attributes"] = cluster_attributes_to_object(d->m_cluster_attributes);
    X["cluster_pair_attributes"] = cluster_pair_attributes_to_object(d->m_cluster_pair_attributes);
    X["timeseries"] = timeseries_map_to_object_for_mv2(d->m_timeseries);
    X["current_timeseries_name"] = d->m_current_timeseries_name;
    X["firings"] = d->m_firings.toPrvObject();
    X["samplerate"] = d->m_sample_rate;
    X["options"] = QJsonObject::fromVariantMap(this->options());
    X["mlproxy_url"] = d->m_mlproxy_url;
    X["visibility_rule"] = d->m_visibility_rule.toJsonObject();
    X["electrode_geometry"] = d->m_electrode_geometry.toJsonObject();
    X["clusters_subset"] = intlist_to_json_array(d->m_clusters_subset.toList());
    X["view_merged"] = d->m_view_merged;
    X["visible_channels"] = intlist_to_json_array(d->m_visible_channels.toList());
    X["cluster_order_scores"] = doublelist_to_json_array(d->m_cluster_order_scores);
    X["cluster_order_scores_name"] = d->m_cluster_order_scores_name;
    return X;
}

void MVContext::setFromMVFileObject(QJsonObject X)
{
    // old
    this->clear();
    d->m_original_object = X; // to preserve unused fields
    d->m_cluster_attributes = object_to_cluster_attributes(X["cluster_attributes"].toObject());
    d->m_cluster_pair_attributes = object_to_cluster_pair_attributes(X["cluster_pair_attributes"].toObject());
    d->m_timeseries = object_to_timeseries_map(X["timeseries"].toObject());
    this->setCurrentTimeseriesName(X["current_timeseries_name"].toString());
    this->setFirings(DiskReadMda(X["firings"].toString()));
    d->m_sample_rate = X["samplerate"].toDouble();
    /*
    if (X.contains("options")) {
        d->m_options = X["options"].toObject().toVariantMap();
    }
    */
    d->m_mlproxy_url = X["mlproxy_url"].toString();
    if (X.contains("visibility_rule")) {
        this->setClusterVisibilityRule(ClusterVisibilityRule::fromJsonObject(X["visibility_rule"].toObject()));
    }
    d->m_electrode_geometry = ElectrodeGeometry::fromJsonObject(X["electrode_geometry"].toObject());
    if (X.contains("clusters_subset")) {
        this->setClustersSubset(json_array_to_intlist(X["clusters_subset"].toArray()).toSet());
    }
    if (X.contains("visible_channels")) {
        this->setVisibleChannels(json_array_to_intlist(X["visible_clusters"].toArray()));
    }
    d->m_view_merged = X["view_merged"].toBool();
    emit this->currentTimeseriesChanged();
    emit this->timeseriesNamesChanged();
    emit this->firingsChanged();
    emit this->clusterMergeChanged();
    emit this->clusterAttributesChanged(0);
    emit this->currentEventChanged();
    emit this->currentClusterChanged();
    emit this->selectedClustersChanged();
    emit this->currentTimepointChanged();
    emit this->currentTimeRangeChanged();
    emit this->optionChanged("");
    emit this->clusterVisibilityChanged();
    emit this->selectedClusterPairsChanged();
    emit this->clusterPairAttributesChanged(ClusterPair(0, 0));
    emit this->viewMergedChanged();
    emit this->visibleChannelsChanged();
}

void MVContext::setFromMV2FileObject(QJsonObject X)
{
    this->clear();
    d->m_original_object = X; // to preserve unused fields
    d->m_cluster_attributes = object_to_cluster_attributes(X["cluster_attributes"].toObject());
    d->m_cluster_pair_attributes = object_to_cluster_pair_attributes(X["cluster_pair_attributes"].toObject());
    d->m_timeseries = object_to_timeseries_map_for_mv2(X["timeseries"].toObject());
    this->setCurrentTimeseriesName(X["current_timeseries_name"].toString());
    this->setFirings(DiskReadMda(X["firings"].toObject()));
    d->m_sample_rate = X["samplerate"].toDouble();
    if (X.contains("options")) {
        this->setOptions(X["options"].toObject().toVariantMap());
    }
    d->m_mlproxy_url = X["mlproxy_url"].toString();
    if (X.contains("visibility_rule")) {
        this->setClusterVisibilityRule(ClusterVisibilityRule::fromJsonObject(X["visibility_rule"].toObject()));
    }
    d->m_electrode_geometry = ElectrodeGeometry::fromJsonObject(X["electrode_geometry"].toObject());
    if (X.contains("clusters_subset")) {
        this->setClustersSubset(json_array_to_intlist(X["clusters_subset"].toArray()).toSet());
    }
    if (X.contains("visible_channels")) {
        this->setVisibleChannels(json_array_to_intlist(X["visible_clusters"].toArray()));
    }
    d->m_view_merged = X["view_merged"].toBool();
    d->m_cluster_order_scores = json_array_to_doublelist(X["cluster_order_scores"].toArray());
    d->m_cluster_order_scores_name = X["cluster_order_scores_name"].toString();

    emit this->currentTimeseriesChanged();
    emit this->timeseriesNamesChanged();
    emit this->firingsChanged();
    emit this->clusterMergeChanged();
    emit this->clusterAttributesChanged(0);
    emit this->currentEventChanged();
    emit this->currentClusterChanged();
    emit this->selectedClustersChanged();
    emit this->currentTimepointChanged();
    emit this->currentTimeRangeChanged();
    emit this->optionChanged("");
    emit this->clusterVisibilityChanged();
    emit this->selectedClusterPairsChanged();
    emit this->clusterPairAttributesChanged(ClusterPair(0, 0));
    emit this->viewMergedChanged();
    emit this->visibleChannelsChanged();
    emit this->clusterOrderChanged();
}

MVEvent MVContext::currentEvent() const
{
    return d->m_current_event;
}

int MVContext::currentCluster() const
{
    return d->m_current_cluster;
}

QList<int> MVContext::selectedClusters() const
{
    return d->m_selected_clusters;
}

double MVContext::currentTimepoint() const
{
    return d->m_current_timepoint;
}

MVRange MVContext::currentTimeRange() const
{
    if ((d->m_current_time_range.min <= 0) && (d->m_current_time_range.max <= 0)) {
        return MVRange(0LL, qMax(0LL, (long long)this->currentTimeseries().N2() - 1));
    }
    return d->m_current_time_range;
}

QList<QColor> MVContext::channelColors() const
{
    return d->m_channel_colors;
}

QList<QColor> MVContext::clusterColors() const
{
    QList<QColor> ret;
    for (int i = 0; i < d->m_cluster_colors.count(); i++) {
        ret << this->clusterColor(i + 1);
    }
    return ret;
}

DiskReadMda32 MVContext::currentTimeseries() const
{
    return timeseries(d->m_current_timeseries_name);
}

DiskReadMda32 MVContext::timeseries(QString name) const
{
    return d->m_timeseries.value(name).data;
}

QString MVContext::currentTimeseriesName() const
{
    return d->m_current_timeseries_name;
}

QColor MVContext::clusterColor(int k) const
{
    if (k < 0)
        return Qt::black;
    if (k == 0)
        return Qt::gray;
    if (d->m_cluster_colors.isEmpty())
        return Qt::black;
    int cluster_color_index_shift = this->option("cluster_color_index_shift", 0).toInt();
    k += cluster_color_index_shift;
    while (k < 1)
        k += d->m_cluster_colors.count();
    return d->m_cluster_colors[(k - 1) % d->m_cluster_colors.count()];
}

QColor MVContext::channelColor(int m) const
{
    if (m < 0)
        return Qt::black;
    if (d->m_channel_colors.isEmpty())
        return Qt::black;
    return d->m_channel_colors[m % d->m_channel_colors.count()];
}

QColor MVContext::color(QString name, QColor default_color) const
{
    return d->m_colors.value(name, default_color);
}

QMap<QString, QColor> MVContext::colors() const
{
    return d->m_colors;
}

QStringList MVContext::timeseriesNames() const
{
    return d->m_timeseries.keys();
}

void MVContext::addTimeseries(QString name, DiskReadMda32 timeseries)
{
    TimeseriesStruct X;
    X.data = timeseries;
    X.name = name;
    d->m_timeseries[name] = X;
    emit this->timeseriesNamesChanged();
    if (name == d->m_current_timeseries_name)
        emit this->currentTimeseriesChanged();
}

DiskReadMda MVContext::firings()
{
    if (d->m_clusters_subset.isEmpty())
        return d->m_firings;
    else {
        return d->m_firings_subset;
    }
}

double MVContext::sampleRate() const
{
    return d->m_sample_rate;
}

void MVContext::setCurrentTimeseriesName(QString name)
{
    if (d->m_current_timeseries_name == name)
        return;
    d->m_current_timeseries_name = name;
    emit this->currentTimeseriesChanged();
}

void MVContext::setFirings(const DiskReadMda& F)
{
    d->m_firings = F;
    emit firingsChanged();
}

void MVContext::setSampleRate(double sample_rate)
{
    /*
    if (!d->m_sample_rate) {
        if (sample_rate<1000) { //assume eeg
            d->m_options["cc_max_dt_msec"]=10000;
            d->m_options["cc_bin_size_msec"]=50;
        }
    }
    */
    d->m_sample_rate = sample_rate;
}

QString MVContext::mlProxyUrl() const
{
    return d->m_mlproxy_url;
}

void MVContext::setMLProxyUrl(QString url)
{
    d->m_mlproxy_url = url;
}

ClusterMerge MVContext::clusterMerge() const
{
    ClusterMerge ret;
    QList<int> keys1 = d->m_cluster_pair_attributes.keys();
    foreach (int k1, keys1) {
        QList<int> keys2 = d->m_cluster_pair_attributes[k1].keys();
        foreach (int k2, keys2) {
            ClusterPair pair(k1, k2);
            if (this->clusterPairTags(pair).contains("merged")) {
                QSet<int> labels;
                labels.insert(pair.k1());
                labels.insert(pair.k2());
                ret.merge(labels);
            }
        }
    }
    return ret;
}

bool MVContext::viewMerged() const
{
    return d->m_view_merged;
}

void MVContext::setViewMerged(bool val)
{
    if (d->m_view_merged == val)
        return;
    d->m_view_merged = val;
    d->update_current_and_selected_clusters_according_to_merged();
    emit this->viewMergedChanged();
}

QJsonObject MVContext::clusterAttributes(int num) const
{
    return d->m_cluster_attributes.value(num);
}

QList<int> MVContext::clusterAttributesKeys() const
{
    return d->m_cluster_attributes.keys();
}

QSet<QString> MVContext::clusterTags(int num) const
{
    return jsonarray2stringset(clusterAttributes(num)["tags"].toArray());
}

QList<QString> MVContext::clusterTagsList(int num) const
{
    QList<QString> ret = clusterTags(num).toList();
    qSort(ret);
    return ret;
}

QSet<QString> MVContext::allClusterTags() const
{
    QSet<QString> ret;
    ret.insert("accepted");
    ret.insert("rejected");
    QList<int> keys = clusterAttributesKeys();
    foreach (int key, keys) {
        QSet<QString> tags0 = clusterTags(key);
        foreach (QString tag, tags0)
            ret.insert(tag);
    }
    return ret;
}

void MVContext::setClusterAttributes(int num, const QJsonObject& obj)
{
    if (d->m_cluster_attributes.value(num) == obj)
        return;
    d->m_cluster_attributes[num] = obj;
    emit this->clusterAttributesChanged(num);
    /// TODO: (HIGH) only emit this if there really was a change
    emit this->clusterVisibilityChanged();
}

void MVContext::setClusterTags(int num, const QSet<QString>& tags)
{
    QJsonObject obj = clusterAttributes(num);
    obj["tags"] = stringset2jsonarray(tags);
    setClusterAttributes(num, obj);
}

QList<int> MVContext::clusterOrder(int max_K) const
{
    QList<int> ret;
    if (d->m_cluster_order_scores.isEmpty()) {
        for (int k = 1; k <= max_K; k++) {
            ret << k;
        }
    }
    QList<int> inds = get_sort_indices(d->m_cluster_order_scores.toVector());

    for (int i = 0; i < inds.count(); i++) {
        ret << inds[i] + 1;
    }
    return ret;
}

QString MVContext::clusterOrderScoresName() const
{
    return d->m_cluster_order_scores_name;
}

QList<double> MVContext::clusterOrderScores() const
{
    return d->m_cluster_order_scores;
}

void MVContext::setClusterOrderScores(QString scores_name, const QList<double>& scores)
{
    if ((scores_name == d->m_cluster_order_scores_name) && (scores == d->m_cluster_order_scores))
        return;
    d->m_cluster_order_scores_name = scores_name;
    d->m_cluster_order_scores = scores;
    emit this->clusterOrderChanged();
}

QJsonObject MVContext::clusterPairAttributes(const ClusterPair& pair) const
{
    return d->m_cluster_pair_attributes.value(pair.k1()).value(pair.k2());
}

QList<ClusterPair> MVContext::clusterPairAttributesKeys() const
{
    QList<ClusterPair> ret;
    QList<int> keys1 = d->m_cluster_pair_attributes.keys();
    foreach (int k1, keys1) {
        QList<int> keys2 = d->m_cluster_pair_attributes[k1].keys();
        foreach (int k2, keys2) {
            ClusterPair pair(k1, k2);
            ret << pair;
        }
    }
    return ret;
}

void MVContext::setClusterPairAttributes(const ClusterPair& pair, const QJsonObject& obj)
{
    if (d->m_cluster_pair_attributes.value(pair.k1()).value(pair.k2()) == obj)
        return;
    ClusterMerge cluster_merge_before = this->clusterMerge();
    d->m_cluster_pair_attributes[pair.k1()][pair.k2()] = obj;
    emit this->clusterPairAttributesChanged(pair);
    if (!(cluster_merge_before == this->clusterMerge())) {
        emit this->clusterMergeChanged();
    }
}

QSet<QString> MVContext::clusterPairTags(const ClusterPair& pair) const
{
    return jsonarray2stringset(clusterPairAttributes(pair)["tags"].toArray());
}

QList<QString> MVContext::clusterPairTagsList(const ClusterPair& pair) const
{
    QList<QString> ret = clusterPairTags(pair).toList();
    qSort(ret);
    return ret;
}

QSet<QString> MVContext::allClusterPairTags() const
{
    QSet<QString> ret;
    ret.insert("merge_candidate");
    ret.insert("merged");
    ret.insert("uncertain");
    QList<ClusterPair> keys = clusterPairAttributesKeys();
    foreach (ClusterPair key, keys) {
        QSet<QString> tags0 = clusterPairTags(key);
        foreach (QString tag, tags0)
            ret.insert(tag);
    }
    return ret;
}

void MVContext::setClusterPairTags(const ClusterPair& pair, const QSet<QString>& tags)
{
    QJsonObject obj = clusterPairAttributes(pair);
    obj["tags"] = stringset2jsonarray(tags);
    setClusterPairAttributes(pair, obj);
}

ClusterVisibilityRule MVContext::clusterVisibilityRule() const
{
    return d->m_visibility_rule;
}

QList<int> MVContext::visibleClusters(int K) const
{
    QList<int> ret;
    for (int k = 1; k <= K; k++) {
        if (this->clusterIsVisible(k)) {
            ret << k;
        }
    }
    return ret;
}

bool MVContext::clusterIsVisible(int k) const
{
    if (this->viewMerged()) {
        if (this->clusterMerge().representativeLabel(k) != k)
            return false;
    }
    return d->m_visibility_rule.isVisible(this, k);
}

void MVContext::setClusterVisibilityRule(const ClusterVisibilityRule& rule)
{
    if (d->m_visibility_rule == rule)
        return;
    d->m_visibility_rule = rule;
    emit this->clusterVisibilityChanged();
}

QList<int> MVContext::clustersToForceShow() const
{
    QList<int> ret = d->m_clusters_to_force_show.toList();
    qSort(ret);
    return ret;
}

void MVContext::setClustersToForceShow(const QList<int>& list)
{
    d->m_clusters_to_force_show = list.toSet();
}

QList<int> MVContext::visibleChannels() const
{
    QList<int> ch = d->m_visible_channels.toList();
    qSort(ch);
    return ch;
}

void MVContext::setVisibleChannels(const QList<int>& channels)
{
    QSet<int> ch = channels.toSet();
    if (ch == d->m_visible_channels)
        return;
    d->m_visible_channels = ch;
    emit visibleChannelsChanged();
}

ElectrodeGeometry MVContext::electrodeGeometry() const
{
    return d->m_electrode_geometry;
}

void MVContext::setElectrodeGeometry(const ElectrodeGeometry& geom)
{
    if (d->m_electrode_geometry == geom)
        return;
    d->m_electrode_geometry = geom;
    emit this->electrodeGeometryChanged();
}

void MVContext::loadClusterMetricsFromFile(QString csv_or_json_file_path)
{
    if (csv_or_json_file_path.endsWith(".csv")) {
        QStringList lines = TextFile::read(csv_or_json_file_path).split("\n", QString::SkipEmptyParts);
        if (lines.isEmpty())
            return;
        QStringList metric_names = lines[0].split(",", QString::SkipEmptyParts);
        for (int i = 1; i < lines.count(); i++) {
            QStringList vals = lines[i].split(",", QString::SkipEmptyParts);
            bool ok;
            int k = vals.value(0).toInt(&ok);
            if (ok) {
                QJsonObject obj = this->clusterAttributes(k);
                QJsonObject metrics = obj["metrics"].toObject();
                for (int j = 1; j < metric_names.count(); j++) {
                    metrics[metric_names[j]] = vals.value(j).toDouble();
                }
                obj["metrics"] = metrics;
                this->setClusterAttributes(k, obj);
            }
        }
    }
    else {
        QString json = TextFile::read(csv_or_json_file_path);
        QJsonObject X = QJsonDocument::fromJson(json.toUtf8()).object();
        QJsonArray clusters = X["clusters"].toArray();
        for (int j = 0; j < clusters.count(); j++) {
            QJsonObject C = clusters[j].toObject();
            int k = C["label"].toDouble();
            QJsonObject C_metrics = C["metrics"].toObject();
            if (k > 0) {
                QJsonObject obj = this->clusterAttributes(k);
                QJsonObject metrics = obj["metrics"].toObject();
                QStringList keys = C_metrics.keys();
                foreach (QString key, keys) {
                    metrics[key] = C_metrics[key];
                }
                obj["metrics"] = metrics;
                this->setClusterAttributes(k, obj);
            }
        }
    }
}

void MVContext::loadClusterPairMetricsFromFile(QString csv_file_path)
{
    QStringList lines = TextFile::read(csv_file_path).split("\n", QString::SkipEmptyParts);
    if (lines.isEmpty())
        return;
    QStringList metric_names = lines[0].split(",", QString::SkipEmptyParts);
    for (int i = 1; i < lines.count(); i++) {
        QStringList vals = lines[i].split(",", QString::SkipEmptyParts);
        bool ok;
        int k1 = vals.value(0).toInt(&ok);
        if (ok) {
            int k2 = vals.value(1).toInt(&ok);
            if (ok) {
                QJsonObject obj = this->clusterPairAttributes(ClusterPair(k1, k2));
                QJsonObject metrics = obj["metrics"].toObject();
                for (int j = 2; j < metric_names.count(); j++) {
                    metrics[metric_names[j]] = vals.value(j).toDouble();
                }
                obj["metrics"] = metrics;
                this->setClusterPairAttributes(ClusterPair(k1, k2), obj);
            }
        }
    }
}

QSet<int> MVContext::clustersSubset() const
{
    return d->m_clusters_subset;
}

class FiringsSubsetCalculator : public QThread {
public:
    //input
    QString mlproxy_url;
    QString firings_path;
    QSet<int> clusters_subset;

    //output
    QString firings_out_path;

    void run()
    {
        TaskProgress task(TaskProgress::Calculate, "Firings subset");
        MountainProcessRunner PR;
        PR.setMLProxyUrl(mlproxy_url);
        PR.setProcessorName("firings_subset");
        QMap<QString, QVariant> params;
        params["firings"] = firings_path;
        QStringList clusters_subset_str;
        QList<int> clusters_subset_list = clusters_subset.toList();
        qSort(clusters_subset_list); //to make it unique
        foreach (int label, clusters_subset_list) {
            clusters_subset_str << QString("%1").arg(label);
        }
        params["labels"] = clusters_subset_str.join(",");
        PR.setInputParameters(params);
        firings_out_path = PR.makeOutputFilePath("firings_out");
        task.log("runProcess");
        PR.runProcess();
        task.log("done with runProcess");
    }
};

void MVContext::setClustersSubset(const QSet<int>& clusters_subset)
{
    if (d->m_clusters_subset == clusters_subset)
        return;
    d->m_clusters_subset = clusters_subset;
    if (!clusters_subset.isEmpty()) {
        FiringsSubsetCalculator* CC = new FiringsSubsetCalculator;
        QObject::connect(CC, SIGNAL(finished()), this, SLOT(slot_firings_subset_calculator_finished()));
        CC->mlproxy_url = d->m_mlproxy_url;
        CC->firings_path = d->m_firings.makePath();
        CC->clusters_subset = d->m_clusters_subset;
        CC->start();
    }
}

void MVContext::setCurrentEvent(const MVEvent& evt)
{
    if (evt == d->m_current_event)
        return;
    d->m_current_event = evt;
    emit currentEventChanged();
    this->setCurrentTimepoint(evt.time);
}

void MVContext::setCurrentCluster(int k)
{
    if (k == d->m_current_cluster)
        return;
    d->m_current_cluster = k;

    QList<int> tmp = selectedClusters();
    if (!tmp.contains(k)) {
        tmp << k;
        this->setSelectedClusters(tmp); //think about this
    }
    d->update_current_and_selected_clusters_according_to_merged();
    emit currentClusterChanged();
}

void MVContext::setSelectedClusters(const QList<int>& ks)
{
    QList<int> ks2 = QList<int>::fromSet(ks.toSet()); //remove duplicates and -1
    ks2.removeAll(-1);
    qSort(ks2);
    if (d->m_selected_clusters == ks2)
        return;
    d->m_selected_clusters = ks2;
    if (!d->m_selected_clusters.contains(d->m_current_cluster)) {
        this->setCurrentCluster(-1);
    }
    d->update_current_and_selected_clusters_according_to_merged();
    emit selectedClustersChanged();
}

void MVContext::setCurrentTimepoint(double tp)
{
    if (tp < 0)
        tp = 0;
    if (tp >= this->currentTimeseries().N2())
        tp = this->currentTimeseries().N2() - 1;
    if (d->m_current_timepoint == tp)
        return;
    d->m_current_timepoint = tp;
    emit currentTimepointChanged();
}

void MVContext::setCurrentTimeRange(const MVRange& range_in)
{
    MVRange range = range_in;
    if (range.max >= this->currentTimeseries().N2()) {
        range = range + (this->currentTimeseries().N2() - 1 - range.max);
    }
    if (range.min < 0) {
        range = range + (0 - range.min);
    }
    if (range.max - range.min < 150) { //don't allow range to be too small
        range.max = range.min + 150;
    }
    if ((range.max >= this->currentTimeseries().N2()) && (range.min == 0)) { //second condition important
        //don't allow it to extend too far
        range.max = this->currentTimeseries().N2() - 1;
    }
    if (d->m_current_time_range == range)
        return;
    d->m_current_time_range = range;
    emit currentTimeRangeChanged();
}

void MVContext::setClusterColors(const QList<QColor>& colors)
{
    d->m_cluster_colors = colors;
    emit clusterColorsChanged(colors);
}

void MVContext::setChannelColors(const QList<QColor>& colors)
{
    d->m_channel_colors = colors;
}

void MVContext::setColors(const QMap<QString, QColor>& colors)
{
    d->m_colors = colors;
}

void MVContext::copySettingsFrom(MVContext* other)
{
    this->setChannelColors(other->channelColors());
    this->setClusterColors(other->clusterColors());
    this->setColors(other->colors());
    this->setMLProxyUrl(other->mlProxyUrl());
    this->setSampleRate(other->sampleRate());
    this->setOptions(other->options());
}

QMap<QString, QJsonObject> MVContext::allPrvObjects()
{
    QMap<QString, QJsonObject> ret;
    ret["firings"] = d->m_firings.toPrvObject();

    QStringList tsnames = this->timeseriesNames();
    foreach (QString tsname, tsnames) {
        ret[tsname] = this->timeseries(tsname).toPrvObject();
    }

    return ret;
}

void MVContext::slot_firings_subset_calculator_finished()
{
    FiringsSubsetCalculator* CC = (FiringsSubsetCalculator*)qobject_cast<QThread*>(sender());
    if (!CC)
        return;
    d->m_firings_subset.setPath(CC->firings_out_path);
    CC->deleteLater();
    emit this->firingsChanged();
}

void MVContext::clickCluster(int k, Qt::KeyboardModifiers modifiers)
{
    if (modifiers & Qt::ControlModifier) {
        if (k < 0)
            return;
        if (d->m_selected_clusters.contains(k)) {
            QList<int> tmp = d->m_selected_clusters;
            tmp.removeAll(k);
            this->setSelectedClusters(tmp);
        }
        else {
            if (k >= 0) {
                QList<int> tmp = d->m_selected_clusters;
                tmp << k;
                this->setSelectedClusters(tmp);
            }
        }
    }
    else {
        this->setSelectedClusterPairs(QSet<ClusterPair>());
        this->setSelectedClusters(QList<int>());
        this->setCurrentCluster(k);
    }
}

QSet<ClusterPair> MVContext::selectedClusterPairs() const
{
    return d->m_selected_cluster_pairs;
}

void MVContext::setSelectedClusterPairs(const QSet<ClusterPair>& pairs)
{
    if (d->m_selected_cluster_pairs == pairs)
        return;
    d->m_selected_cluster_pairs = pairs;
    emit this->selectedClusterPairsChanged();
}

void MVContext::clickClusterPair(const ClusterPair& pair, Qt::KeyboardModifiers modifiers)
{
    if (modifiers & Qt::ControlModifier) {
        if ((pair.k1() <= 0) || (pair.k2() <= 0))
            return;
        QSet<ClusterPair> tmp = d->m_selected_cluster_pairs;
        if (tmp.contains(pair)) {
            tmp.remove(pair);
        }
        else {
            tmp.insert(pair);
        }
        this->setSelectedClusterPairs(tmp);
    }
    else {
        QSet<ClusterPair> tmp;
        tmp.insert(pair);
        this->setSelectedClusterPairs(tmp);
        QList<int> list0;
        list0 << pair.k1() << pair.k2();
        this->setSelectedClusters(list0);
        this->setCurrentCluster(-1);
    }
}

ClusterVisibilityRule::ClusterVisibilityRule()
{
}

ClusterVisibilityRule::ClusterVisibilityRule(const ClusterVisibilityRule& other)
{
    copy_from(other);
}

ClusterVisibilityRule::~ClusterVisibilityRule()
{
}

void ClusterVisibilityRule::operator=(const ClusterVisibilityRule& other)
{
    copy_from(other);
}

bool ClusterVisibilityRule::operator==(const ClusterVisibilityRule& other) const
{
    if (this->view_tags != other.view_tags)
        return false;
    if (this->view_all_tagged != other.view_all_tagged)
        return false;
    if (this->view_all_untagged != other.view_all_untagged)
        return false;
    if (this->hide_rejected != other.hide_rejected)
        return false;
    if (this->use_subset != other.use_subset)
        return false;
    if (this->subset != other.subset)
        return false;
    return true;
}

void ClusterVisibilityRule::copy_from(const ClusterVisibilityRule& other)
{
    this->view_tags = other.view_tags;
    this->view_all_tagged = other.view_all_tagged;
    this->view_all_untagged = other.view_all_untagged;
    this->hide_rejected = other.hide_rejected;

    this->use_subset = other.use_subset;
    this->subset = other.subset;
}

bool ClusterVisibilityRule::isVisible(const MVContext* context, int cluster_num) const
{
    if (this->use_subset) {
        if (!subset.contains(cluster_num))
            return false;
    }

    QSet<QString> tags = context->clusterTags(cluster_num);

    if (hide_rejected) {
        if (tags.contains("rejected"))
            return false;
    }

    if ((view_all_tagged) && (!tags.isEmpty()))
        return true;
    if ((view_all_untagged) && (tags.isEmpty()))
        return true;

    foreach (QString tag, tags) {
        if (view_tags.contains(tag))
            return true;
    }

    return false;
}

QJsonArray intlist_to_json_array(const QList<int>& X)
{
    QJsonArray ret;
    foreach (int x, X) {
        ret.push_back(x);
    }
    return ret;
}

QJsonArray doublelist_to_json_array(const QList<double>& X)
{
    QJsonArray ret;
    foreach (double x, X) {
        ret.push_back(x);
    }
    return ret;
}

QList<int> json_array_to_intlist(const QJsonArray& X)
{
    QList<int> ret;
    for (int i = 0; i < X.count(); i++) {
        ret << X[i].toInt();
    }
    return ret;
}

QList<double> json_array_to_doublelist(const QJsonArray& X)
{
    QList<double> ret;
    for (int i = 0; i < X.count(); i++) {
        ret << X[i].toDouble();
    }
    return ret;
}

QJsonArray strlist_to_json_array(const QList<QString>& X)
{
    QJsonArray ret;
    foreach (QString x, X) {
        ret.push_back(x);
    }
    return ret;
}

QList<QString> json_array_to_strlist(const QJsonArray& X)
{
    QList<QString> ret;
    for (int i = 0; i < X.count(); i++) {
        ret << X[i].toString();
    }
    return ret;
}

QJsonObject ClusterVisibilityRule::toJsonObject() const
{
    QJsonObject obj;
    obj["view_all_tagged"] = this->view_all_tagged;
    obj["view_all_untagged"] = this->view_all_untagged;
    obj["hide_rejected"] = this->hide_rejected;
    obj["view_tags"] = strlist_to_json_array(this->view_tags.toList());
    obj["use_subset"] = this->use_subset;
    obj["subset"] = intlist_to_json_array(this->subset.toList());
    return obj;
}

ClusterVisibilityRule ClusterVisibilityRule::fromJsonObject(const QJsonObject& X)
{
    ClusterVisibilityRule ret;
    ret.view_all_tagged = X["view_all_tagged"].toBool();
    ret.view_all_untagged = X["view_all_untagged"].toBool();
    ret.hide_rejected = X["hide_rejected"].toBool();
    ret.view_tags = QSet<QString>::fromList(json_array_to_strlist(X["view_tags"].toArray()));
    ret.use_subset = X["use_subset"].toBool();
    ret.subset = QSet<int>::fromList(json_array_to_intlist(X["subset"].toArray()));
    return ret;
}

ClusterPair::ClusterPair(int k1, int k2)
{
    set(k1, k2);
}

ClusterPair::ClusterPair(const ClusterPair& other)
{
    set(other.k1(), other.k2());
}

void ClusterPair::set(int k1, int k2)
{
    m_k1 = k1;
    m_k2 = k2;
}

int ClusterPair::k1() const
{
    return m_k1;
}

int ClusterPair::k2() const
{
    return m_k2;
}

void ClusterPair::operator=(const ClusterPair& other)
{
    set(other.k1(), other.k2());
}

bool ClusterPair::operator==(const ClusterPair& other) const
{
    return ((k1() == other.k1()) && (k2() == other.k2()));
}

bool ClusterPair::operator<(const ClusterPair& other) const
{
    if (k1() < other.k1())
        return true;
    if (k2() > other.k2())
        return false;
    return (k2() < other.k2());
}

QString ClusterPair::toString() const
{
    return QString("pair_%1_%2").arg(m_k1).arg(m_k2);
}

ClusterPair ClusterPair::fromString(const QString& str)
{
    QStringList list = str.split("_");
    if (list.count() != 3)
        return ClusterPair(0, 0);
    return ClusterPair(list.value(1).toInt(), list.value(2).toInt());
}

uint qHash(const ClusterPair& pair)
{
    return qHash(pair.toString());
}

QJsonObject ElectrodeGeometry::toJsonObject() const
{
    QJsonArray C;
    for (int j = 0; j < coordinates.count(); j++) {
        QJsonArray X;
        for (int k = 0; k < coordinates[j].count(); k++) {
            X.push_back(coordinates[j][k]);
        }
        C.push_back(X);
    }
    QJsonObject ret;
    ret["coordinates"] = C;
    return ret;
}

bool ElectrodeGeometry::operator==(const ElectrodeGeometry& other)
{
    if (coordinates != other.coordinates)
        return false;
    return true;
}

ElectrodeGeometry ElectrodeGeometry::fromJsonObject(const QJsonObject& obj)
{
    ElectrodeGeometry ret;
    QJsonArray C = obj["coordinates"].toArray();
    for (int i = 0; i < C.count(); i++) {
        QVector<double> tmp;
        QJsonArray X = C[i].toArray();
        for (int j = 0; j < X.count(); j++) {
            tmp << X[j].toDouble();
        }
        ret.coordinates << tmp;
    }
    return ret;
}

ElectrodeGeometry ElectrodeGeometry::loadFromGeomFile(const QString& path)
{
    ElectrodeGeometry ret;
    Mda geom(path);
    if (geom.totalSize() > 10000) {
        qWarning() << "The geometry file is too large for my limited vision." << path << geom.N1() << geom.N2() << geom.N3();
        return ret;
    }
    for (int i = 0; i < geom.N2(); i++) {
        QVector<double> tmp;
        for (int j = 0; j < geom.N1(); j++) {
            tmp << geom.value(j, i);
        }
        ret.coordinates << tmp;
    }
    return ret;
}

void MVContextPrivate::update_current_and_selected_clusters_according_to_merged()
{
    if (m_view_merged) {
        ClusterMerge CM = q->clusterMerge();
        if (m_current_cluster >= 0) {
            m_current_cluster = CM.representativeLabel(m_current_cluster);
        }
        QSet<int> selected_clusters_set = m_selected_clusters.toSet();
        foreach (int k, m_selected_clusters) {
            QList<int> list = CM.getMergeGroup(k);
            foreach (int k2, list) {
                selected_clusters_set.insert(k2);
            }
        }
        m_selected_clusters = selected_clusters_set.toList();
        qSort(m_selected_clusters);
    }
}

void MVContextPrivate::set_default_options()
{
    q->setOption("clip_size", 100);
    q->setOption("cc_max_dt_msec", 100);
    q->setOption("cc_log_time_constant_msec", 1);
    q->setOption("cc_bin_size_msec", 0.5);
    q->setOption("cc_max_est_data_size", "1e4");
    q->setOption("amp_thresh_display", 3);
    q->setOption("discrim_hist_method", "centroid");
}
