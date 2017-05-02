/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
** Created: 5/27/2016
*******************************************************/

#ifndef MULTISCALEMDA_H
#define MULTISCALEMDA_H

#include <diskreadmda32.h>

class MultiScaleTimeSeriesPrivate;
class MultiScaleTimeSeries {
public:
    friend class MultiScaleTimeSeriesPrivate;
    MultiScaleTimeSeries();
    MultiScaleTimeSeries(const MultiScaleTimeSeries& other);
    virtual ~MultiScaleTimeSeries();
    void operator=(const MultiScaleTimeSeries& other);
    void setData(const DiskReadMda32& X);
    void setMLProxyUrl(const QString& url);
    void initialize();

    int N1();
    int N2();
    bool getData(Mda32& min, Mda32& max, int t1, int t2, int ds_factor); //returns values at timepoints i1*ds_factor:ds_factor:i2*ds_factor
    double minimum(); //return the global minimum value
    double maximum(); //return the global maximum value

    static int smallest_power_of_3_larger_than(int N);

private:
    MultiScaleTimeSeriesPrivate* d;
};

#endif // MULTISCALEMDA_H
