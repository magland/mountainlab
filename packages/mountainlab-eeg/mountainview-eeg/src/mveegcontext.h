#ifndef MVEEGCONTEXT_H
#define MVEEGCONTEXT_H

#include "mvabstractcontext.h"

#include <diskreadmda.h>

class MVEEGContextPrivate;
class MVEEGContext : public MVAbstractContext {
    Q_OBJECT
public:
    friend class MVEEGContextPrivate;
    MVEEGContext();
    virtual ~MVEEGContext();

    void clear();

    void setSampleRate(double rate);
    double sampleRate() const;

    /////////////////////////////////////////////////
    DiskReadMda currentTimeseries() const;
    DiskReadMda32 timeseries(QString name) const;
    QString currentTimeseriesName() const;
    QStringList timeseriesNames() const;
    void addTimeseries(QString name, DiskReadMda32 timeseries);
    void setCurrentTimeseriesName(QString name);

    /////////////////////////////////////////////////
    double currentTimepoint() const;
    MVRange currentTimeRange() const;
    void setCurrentTimepoint(double tp);
    void setCurrentTimeRange(const MVRange& range);

    void setFromMV2FileObject(QJsonObject obj) Q_DECL_OVERRIDE;
    QJsonObject toMV2FileObject() const Q_DECL_OVERRIDE;

signals:
    void currentTimeseriesChanged();
    void timeseriesNamesChanged();
    void currentTimepointChanged();
    void currentTimeRangeChanged();

private:
    MVEEGContextPrivate* d;
};

#endif // MVEEGCONTEXT_H
