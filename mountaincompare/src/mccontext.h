#ifndef MCCONTEXT_H
#define MCCONTEXT_H

#include "mvabstractcontext.h"
#include "mvcontext.h"

///A firing event defined by a time, label, and firings_num
struct MCEvent {
    bool operator==(const MCEvent& other) const
    {
        return ((time == other.time) && (label == other.label) && (firings_num == other.firings_num));
    }

    double time = -1;
    int label = -1;
    int firings_num = 1;
};

struct MCCluster {
    int num = -1;
    int firings_num = 1;
    bool operator==(const MCCluster& other) const
    {
        return ((num == other.num) && (firings_num == other.firings_num));
    }
    QString toString() const
    {
        return QString("%1,%2").arg(num).arg(firings_num);
    }
};

uint qHash(const MCCluster& C);

class MCContextPrivate;
class MCContext : public MVAbstractContext {
    Q_OBJECT
public:
    friend class MCContextPrivate;
    MCContext();
    virtual ~MCContext();

    MVContext* mvContext1();
    MVContext* mvContext2();

    /////////////////////////////////////////////////
    void setFromMV2FileObject(QJsonObject obj) Q_DECL_OVERRIDE;
    QJsonObject toMV2FileObject() const Q_DECL_OVERRIDE;

    /////////////////////////////////////////////////
    MCEvent currentEvent() const;
    MCCluster currentCluster() const;
    QList<MCCluster> selectedClusters() const;
    double currentTimepoint() const;
    MVRange currentTimeRange() const;
    void setCurrentEvent(const MCEvent& evt);
    void setCurrentCluster(const MCCluster& C);
    void setSelectedClusters(const QList<MCCluster>& Cs);
    void setCurrentTimepoint(double tp);
    void setCurrentTimeRange(const MVRange& range);
    void clickCluster(MCCluster C, Qt::KeyboardModifiers modifiers);

    /////////////////////////////////////////////////
    QColor clusterColor(int k) const;
    QColor channelColor(int m) const;
    QColor color(QString name, QColor default_color = Qt::black) const;
    QMap<QString, QColor> colors() const;
    QList<QColor> channelColors() const;
    QList<QColor> clusterColors() const;
    void setClusterColors(const QList<QColor>& colors);
    void setChannelColors(const QList<QColor>& colors);
    void setColors(const QMap<QString, QColor>& colors);

    /////////////////////////////////////////////////
    DiskReadMda32 currentTimeseries() const;
    DiskReadMda32 timeseries(QString name) const;
    QString currentTimeseriesName() const;
    QStringList timeseriesNames() const;
    void addTimeseries(QString name, DiskReadMda32 timeseries);
    void setCurrentTimeseriesName(QString name);

    /////////////////////////////////////////////////
    DiskReadMda firings1();
    DiskReadMda firings2();
    DiskReadMda matchedFirings();
    //DiskReadMda firingsMerged();
    void setFirings1(const DiskReadMda& F);
    void setFirings2(const DiskReadMda& F);
    void setConfusionMatrix(const DiskReadMda& CM);
    void setLabelMap(const DiskReadMda& LM);
    void setMatchedFirings(const DiskReadMda& MF);

    void setClustersToForceShow(const QList<int>& list);

    /////////////////////////////////////////////////
    // these should be set once at beginning
    double sampleRate() const;
    void setSampleRate(double sample_rate);

    /////////////////////////////////////////////////
    Mda confusionMatrix() const;
    QList<int> labelMap() const;

signals:
    void currentTimeseriesChanged();
    void timeseriesNamesChanged();
    void firingsChanged();
    void currentEventChanged();
    void currentClusterChanged();
    void selectedClustersChanged();
    void currentTimepointChanged();
    void currentTimeRangeChanged();
    void clusterColorsChanged(const QList<QColor>&);

private slots:
    void slot_context_current_timepoint_changed();
    void slot_context_current_time_range_changed();
    void slot_context_current_timeseries_changed();

private:
    MCContextPrivate* d;
};

#endif // MCCONTEXT_H
