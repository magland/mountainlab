#include "mveegcontext.h"

struct EEGTimeseriesStruct {
    QString name;
    DiskReadMda data;
};

class MVEEGContextPrivate {
public:
    MVEEGContext* q;

    QJsonObject m_original_object; //to preserve unused fields

    QMap<QString, EEGTimeseriesStruct> m_timeseries;
    QString m_current_timeseries_name;
    MVRange m_current_time_range = MVRange(0, 0); //(0,0) triggers automatic calculation
    double m_current_timepoint = 0;

    double m_sample_rate = 0;

    static QMap<QString, EEGTimeseriesStruct> object_to_timeseries_map_for_mv2(QJsonObject X);
    static QJsonObject timeseries_map_to_object_for_mv2(const QMap<QString, EEGTimeseriesStruct>& TT);
};

MVEEGContext::MVEEGContext()
{
    d = new MVEEGContextPrivate;
    d->q = this;

    this->setOption("spectrogram_time_resolution", 32);
    this->setOption("spectrogram_freq_range", "8-12");
    this->setOption("timeseries_freq_range", "*");
}

MVEEGContext::~MVEEGContext()
{
    delete d;
}

void MVEEGContext::clear()
{
    d->m_timeseries.clear();
    clearOptions();
}

void MVEEGContext::setSampleRate(double rate)
{
    d->m_sample_rate = rate;
}

double MVEEGContext::sampleRate() const
{
    return d->m_sample_rate;
}

QMap<QString, EEGTimeseriesStruct> MVEEGContextPrivate::object_to_timeseries_map_for_mv2(QJsonObject X)
{
    QMap<QString, EEGTimeseriesStruct> ret;
    QStringList keys = X.keys();
    foreach (QString key, keys) {
        QJsonObject obj = X[key].toObject();
        EEGTimeseriesStruct A;
        A.data = DiskReadMda(obj["data"].toObject());
        A.name = obj["name"].toString();
        ret[key] = A;
    }
    return ret;
}

QJsonObject MVEEGContextPrivate::timeseries_map_to_object_for_mv2(const QMap<QString, EEGTimeseriesStruct>& TT)
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

void MVEEGContext::setFromMV2FileObject(QJsonObject X)
{
    this->clear();
    d->m_original_object = X; // to preserve unused fields
    d->m_timeseries = MVEEGContextPrivate::object_to_timeseries_map_for_mv2(X["timeseries"].toObject());
    this->setCurrentTimeseriesName(X["current_timeseries_name"].toString());
    //this->setFirings(DiskReadMda(X["firings"].toObject()));
    d->m_sample_rate = X["samplerate"].toDouble();
    if (X.contains("options")) {
        this->setOptions(X["options"].toObject().toVariantMap());
    }
    //d->m_mlproxy_url = X["mlproxy_url"].toString();

    emit this->currentTimeseriesChanged();
    emit this->timeseriesNamesChanged();
    emit this->currentTimepointChanged();
    emit this->currentTimeRangeChanged();
    emit this->optionChanged("");
}

QJsonObject MVEEGContext::toMV2FileObject() const
{
    QJsonObject X = d->m_original_object;
    X["timeseries"] = MVEEGContextPrivate::timeseries_map_to_object_for_mv2(d->m_timeseries);
    X["current_timeseries_name"] = d->m_current_timeseries_name;
    X["samplerate"] = d->m_sample_rate;
    X["options"] = QJsonObject::fromVariantMap(this->options());
    return X;
}

DiskReadMda MVEEGContext::currentTimeseries() const
{
    return timeseries(d->m_current_timeseries_name);
}

DiskReadMda MVEEGContext::timeseries(QString name) const
{
    return d->m_timeseries.value(name).data;
}

QString MVEEGContext::currentTimeseriesName() const
{
    return d->m_current_timeseries_name;
}

QStringList MVEEGContext::timeseriesNames() const
{
    return d->m_timeseries.keys();
}

void MVEEGContext::addTimeseries(QString name, DiskReadMda32 timeseries)
{
    EEGTimeseriesStruct X;
    X.data = timeseries;
    X.name = name;
    d->m_timeseries[name] = X;
    emit this->timeseriesNamesChanged();
    if (name == d->m_current_timeseries_name)
        emit this->currentTimeseriesChanged();
}

void MVEEGContext::setCurrentTimeseriesName(QString name)
{
    if (d->m_current_timeseries_name == name)
        return;
    d->m_current_timeseries_name = name;
    emit this->currentTimeseriesChanged();
}

void MVEEGContext::setCurrentTimepoint(double tp)
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

void MVEEGContext::setCurrentTimeRange(const MVRange& range_in)
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

double MVEEGContext::currentTimepoint() const
{
    return d->m_current_timepoint;
}

MVRange MVEEGContext::currentTimeRange() const
{
    if ((d->m_current_time_range.min <= 0) && (d->m_current_time_range.max <= 0)) {
        return MVRange(0, qMax(0L, this->currentTimeseries().N2() - 1));
    }
    return d->m_current_time_range;
}
