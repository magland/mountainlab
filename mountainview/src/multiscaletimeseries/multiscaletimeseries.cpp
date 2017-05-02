/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
** Created: 5/27/2016
*******************************************************/

#include "multiscaletimeseries.h"
#include "cachemanager.h"
#include "mlcommon.h"

#include <QFile>
#include <diskwritemda.h>
#include <taskprogress.h>
#include <sys/stat.h>
#include <QFileInfo>
#include <QMutex>
#include <diskreadmda32.h>
#include <math.h>
#include "mountainprocessrunner.h"
#include "mlcommon.h"

class MultiScaleTimeSeriesPrivate {
public:
    MultiScaleTimeSeries* q;

    DiskReadMda32 m_data;
    DiskReadMda32 m_multiscale_data;
    bool m_initialized;
    QString m_remote_data_type;
    QString m_mlproxy_url;

    QString get_multiscale_fname();
    bool get_data(Mda32& min, Mda32& max, int t1, int t2, int ds_factor);

    static bool is_power_of_3(int N);
};

MultiScaleTimeSeries::MultiScaleTimeSeries()
{
    d = new MultiScaleTimeSeriesPrivate;
    d->q = this;
    d->m_initialized = false;
    d->m_remote_data_type = "float32";
}

MultiScaleTimeSeries::MultiScaleTimeSeries(const MultiScaleTimeSeries& other)
{
    d = new MultiScaleTimeSeriesPrivate;
    d->q = this;

    d->m_data = other.d->m_data;
    d->m_multiscale_data = other.d->m_multiscale_data;
    d->m_initialized = other.d->m_initialized;
    d->m_mlproxy_url = other.d->m_mlproxy_url;
    d->m_remote_data_type = other.d->m_remote_data_type;
}

MultiScaleTimeSeries::~MultiScaleTimeSeries()
{
    delete d;
}

void MultiScaleTimeSeries::operator=(const MultiScaleTimeSeries& other)
{
    d->m_data = other.d->m_data;
    d->m_multiscale_data = other.d->m_multiscale_data;
    d->m_initialized = other.d->m_initialized;
    d->m_mlproxy_url = other.d->m_mlproxy_url;
    d->m_remote_data_type = other.d->m_remote_data_type;
}

void MultiScaleTimeSeries::setData(const DiskReadMda32& X)
{
    d->m_data = X;
    d->m_multiscale_data = DiskReadMda32();
    d->m_initialized = false;
}

void MultiScaleTimeSeries::setMLProxyUrl(const QString& url)
{
    d->m_mlproxy_url = url;
}

void MultiScaleTimeSeries::initialize()
{
    TaskProgress task("Initializing multiscaletimeseries");
    QString path;
    {
        path = d->m_data.makePath();
        if (path.isEmpty()) {
            qWarning() << "Unable to initialize multiscaletimeseries.... path is empty.";
            task.error() << "Unable to initialize multiscaletimeseries.... path is empty.";
            d->m_initialized = true;
            return;
        }
    }

    MountainProcessRunner MPR;
    QString path_out;
    {
        MPR.setProcessorName("create_multiscale_timeseries");
        QVariantMap params;
        params["timeseries"] = path;
        MPR.setInputParameters(params);
        MPR.setMLProxyUrl(d->m_mlproxy_url);
        path_out = MPR.makeOutputFilePath("timeseries_out");
        //MPR.setDetach(true);
        task.log("Running process");
    }
    MPR.runProcess();
    if (MLUtil::threadInterruptRequested()) {
        return;
    }
    {
        d->m_multiscale_data.setPath(path_out);
        task.log(d->m_data.makePath());
        task.log(d->m_multiscale_data.makePath());
        task.log(QString("%1x%2 -- %3x%4").arg(d->m_multiscale_data.N1()).arg(d->m_multiscale_data.N2()).arg(d->m_data.N1()).arg(d->m_data.N2()));
        d->m_initialized = true;
    }
}

int MultiScaleTimeSeries::N1()
{
    return d->m_data.N1();
}

int MultiScaleTimeSeries::N2()
{
    return d->m_data.N2();
}

bool MultiScaleTimeSeries::getData(Mda32& min, Mda32& max, int t1, int t2, int ds_factor)
{
    return d->get_data(min, max, t1, t2, ds_factor);
}

double MultiScaleTimeSeries::minimum()
{
    int ds_factor = MultiScaleTimeSeries::smallest_power_of_3_larger_than(this->N2() / 3);
    Mda32 min, max;
    this->getData(min, max, 0, 0, ds_factor);
    return min.minimum();
}

double MultiScaleTimeSeries::maximum()
{
    int ds_factor = MultiScaleTimeSeries::smallest_power_of_3_larger_than(this->N2() / 3);
    Mda32 min, max;
    this->getData(min, max, 0, 0, ds_factor);
    return max.maximum();
}

int MultiScaleTimeSeries::smallest_power_of_3_larger_than(int N)
{
    int ret = 1;
    while (ret < N) {
        ret *= 3;
    }
    return ret;
}

bool MultiScaleTimeSeriesPrivate::get_data(Mda32& min, Mda32& max, int t1, int t2, int ds_factor)
{
    int M, N, N2;
    {
        if (!m_initialized) {
            qWarning() << "Cannot get_data. Multiscale timeseries is not initialized!";
            return false;
        }
        M = m_data.N1();
        N2 = m_data.N2();
        N = MultiScaleTimeSeries::smallest_power_of_3_larger_than(N2);

        if ((t2 < 0) || (t1 >= m_data.N2() / ds_factor)) {
            //we are completely out of range, so we return all zeros
            min.allocate(M, (t2 - t1 + 1));
            max.allocate(M, (t2 - t1 + 1));
            return true;
        }
    }

    if ((t1 < 0) || (t2 >= N2 / ds_factor)) {
        //we are somewhat out of range.
        min.allocate(M, t2 - t1 + 1);
        max.allocate(M, t2 - t1 + 1);
        Mda32 min0, max0;
        int s1 = t1, s2 = t2;
        if (s1 < 0)
            s1 = 0;
        if (s2 >= N2 / ds_factor)
            s2 = N2 / ds_factor - 1;
        if (!get_data(min0, max0, s1, s2, ds_factor)) {
            return false;
        }
        if (t1 >= 0) {
            min.setChunk(min0, 0, 0);
            max.setChunk(max0, 0, 0);
        }
        else {
            min.setChunk(min0, 0, -t1);
            max.setChunk(max0, 0, -t1);
        }
        return true;
    }

    if (ds_factor == 1) {
        if (!m_data.readChunk(min, 0, t1, M, t2 - t1 + 1)) {
            qWarning() << "Unable to read chunk of data in multi-scale timeseries (1)";
            return false;
        }
        max = min;
        return true;
    }

    {
        if (!is_power_of_3(ds_factor)) {
            qWarning() << "Invalid ds_factor: " + ds_factor;
            return false;
        }

        //m_multiscale_data.setRemoteDataType(m_remote_data_type);

        int t_offset_min = 0;
        int ds_factor_0 = 3;
        while (ds_factor_0 < ds_factor) {
            t_offset_min += 2 * (N / ds_factor_0);
            ds_factor_0 *= 3;
        }
        int t_offset_max = t_offset_min + N / ds_factor;

        if (!m_multiscale_data.readChunk(min, 0, t1 + t_offset_min, M, t2 - t1 + 1)) {
            qWarning() << "Unable to read chunk of data in multi-scale timeseries (2)";
            return false;
        }
        if (!m_multiscale_data.readChunk(max, 0, t1 + t_offset_max, M, t2 - t1 + 1)) {
            qWarning() << "Unable to read chunk of data in multi-scale timeseries (3)";
            return false;
        }

        if (MLUtil::threadInterruptRequested()) {
            return false;
        }
    }

    return true;
}

bool MultiScaleTimeSeriesPrivate::is_power_of_3(int N)
{
    double val = N;
    while (val > 1) {
        val /= 3;
    }
    return (val == 1);
}
