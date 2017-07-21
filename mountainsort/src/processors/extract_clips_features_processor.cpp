/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
** Created: 4/5/2016
*******************************************************/

#include "extract_clips_features_processor.h"
#include "extract_clips.h"
#include "pca.h"

class extract_clips_features_ProcessorPrivate {
public:
    extract_clips_features_Processor* q;
};

extract_clips_features_Processor::extract_clips_features_Processor()
{
    d = new extract_clips_features_ProcessorPrivate;
    d->q = this;

    this->setName("extract_clips_features");
    this->setVersion("0.11");
    this->setInputFileParameters("timeseries", "firings");
    this->setOutputFileParameters("features");
    this->setRequiredParameters("clip_size", "num_features");
    this->setOptionalParameters("subtract_mean");
    this->setOptionalParameters("weighted_pca");
}

extract_clips_features_Processor::~extract_clips_features_Processor()
{
    delete d;
}

bool extract_clips_features_Processor::check(const QMap<QString, QVariant>& params)
{
    if (!this->checkParameters(params))
        return false;
    return true;
}

bool extract_clips_features_Processor::run(const QMap<QString, QVariant>& params)
{
    QString timeseries_path = params["timeseries"].toString();
    QString firings_path = params["firings"].toString();
    QString features_path = params["features"].toString();
    int clip_size = params["clip_size"].toInt();
    int num_features = params["num_features"].toInt();
    int subtract_mean = params.value("subtract_mean", 0).toInt();
    int weighted_pca = params.value("weighted_pca",0).toInt();

    DiskReadMda X(timeseries_path);
    DiskReadMda F(firings_path);
    QVector<double> times;
    //QVector<int> labels;
    for (int i = 0; i < F.N2(); i++) {
        times << F.value(1, i);
        //labels << (int)F.value(2,i);
    }
    Mda clips = extract_clips(X, times, clip_size);
    int M = clips.N1();
    int T = clips.N2();
    int L = clips.N3();
    Mda clips_reshaped(M * T, L);
    int NNN = M * T * L;
    for (int iii = 0; iii < NNN; iii++) {
        clips_reshaped.set(clips.get(iii), iii);
    }
    Mda CC, FF, sigma;
    if (weighted_pca) {
        wpca(CC, FF, sigma, clips_reshaped, num_features, subtract_mean);
    }else {
        pca(CC, FF, sigma, clips_reshaped, num_features, subtract_mean);
    }
    return FF.write32(features_path);
    //return extract_clips_features(timeseries_path,firings_path,features_path,clip_size,num_features);
}
