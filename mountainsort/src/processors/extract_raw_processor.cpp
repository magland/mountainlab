/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
** Created: 4/20/2016
*******************************************************/

#include "extract_raw_processor.h"
#include <QStringList>
#include <QList>
#include <diskreadmda.h>
#include <diskwritemda.h>
#include <QTime>

class extract_raw_ProcessorPrivate {
public:
    extract_raw_Processor* q;
};

extract_raw_Processor::extract_raw_Processor()
{
    d = new extract_raw_ProcessorPrivate;
    d->q = this;

    this->setName("extract_raw");
    this->setVersion("0.12");
    this->setInputFileParameters("timeseries");
    this->setOutputFileParameters("timeseries_out");
    this->setOptionalParameters("t1", "t2", "channels");
}

extract_raw_Processor::~extract_raw_Processor()
{
    delete d;
}

bool extract_raw_Processor::check(const QMap<QString, QVariant>& params)
{
    if (!this->checkParameters(params))
        return false;
    return true;
}

QList<int> str_to_intlist(QString str)
{
    QStringList list = str.split(",");
    QList<int> ret;
    for (int i = 0; i < list.count(); i++) {
        QString val = list[i].trimmed();
        if (!val.isEmpty())
            ret << val.toInt();
    }
    return ret;
}

bool extract_raw_Processor::run(const QMap<QString, QVariant>& params)
{
    QString timeseries_path = params["timeseries"].toString();
    QString timeseries_out_path = params["timeseries_out"].toString();
    long t1 = params.value("t1", "-1").toLongLong();
    long t2 = params.value("t2", "-1").toLongLong();
    QString channels_str = params["channels"].toString();
    QList<int> channels = str_to_intlist(channels_str);

    DiskReadMda X(timeseries_path);
    long M = X.N1();
    long N = X.N2();

    if (t1 < 0) {
        t1 = 1;
        t2 = N;
    }

    if ((t1 < 0) || (t2 < t1) || (t2 > N)) {
        printf("Unexpected input parameters, t1=%ld, t2=%ld, N=%ld\n", t1, t2, N);
        return false;
    }

    if (channels.isEmpty()) {
        for (int m = 1; m <= M; m++)
            channels << m;
    }

    long M2 = channels.count();
    long N2 = t2 - t1 + 1;

    for (int j = 0; j < M2; j++) {
        if ((channels[j] <= 0) || (channels[j] > M)) {
            printf("Unexpected input channel %d (M=%ld)\n", channels[j], M);
            return false;
        }
    }

    DiskWriteMda Y(MDAIO_TYPE_FLOAT32, timeseries_out_path, M2, N2);
    long chunk_size = qMax(1000L, (long)1e6 / M);
    QTime timer;
    timer.start();
    for (long t = 0; t < N2; t += chunk_size) {
        if ((t == 0) || (timer.elapsed() > 5000)) {
            printf("extract raw %ld/%ld (%d%%)\n", t, N2, (int)(t * 100.0 / N2));
            timer.restart();
        }
        long aa = chunk_size;
        if (t + aa > N2)
            aa = N2 - t;
        Mda chunk;
        X.readChunk(chunk, 0, t1 - 1 + t, M, aa);
        Mda chunk2(M2, aa);
        for (int k = 0; k < aa; k++) {
            for (int j = 0; j < channels.count(); j++) {
                chunk2.set(chunk.value(channels[j] - 1, k), j, k);
            }
        }
        Y.writeChunk(chunk2, 0, t);
    }

    return true;
}
