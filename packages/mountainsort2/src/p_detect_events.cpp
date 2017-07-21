#include "p_detect_events.h"

#include <diskreadmda32.h>
#include <mda.h>
#include "mlcommon.h"

namespace P_detect_events {
    QVector<double> align_events(const QVector<double>& X, const QVector<double>& ptimes, int sign, double detect_interval);
    QVector<double> detect_events(const QVector<double>& X, double detect_threshold, double detect_interval, int sign);
    QVector<double> subsample_events(const QVector<double>& X, double subsample_factor);
}

bool p_detect_events(QString timeseries, QString event_times_out, P_detect_events_opts opts)
{
    DiskReadMda32 X(timeseries);
    bigint M = X.N1();
    bigint N = X.N2();

    printf("Collecting data vector...\n");
    QVector<double> data(N);
    if (opts.central_channel > 0) {
        if (opts.central_channel - 1 >= M) {
            qWarning() << "Central channel is out of range:" << opts.central_channel << M;
            return false;
        }
        for (bigint i = 0; i < N; i++) {
            data[i] = X.value(opts.central_channel - 1, i);
        }
    }
    else {
        for (bigint i = 0; i < N; i++) {
            double best_value = 0;
            bigint best_m = 0;
            for (bigint m = 0; m < M; m++) {
                double val = X.value(m, i);
                if (opts.sign < 0)
                    val = -val;
                if (opts.sign == 0)
                    val = fabs(val);
                if (val > best_value) {
                    best_value = val;
                    best_m = m;
                }
            }
            data[i] = X.value(best_m, i);
        }
    }
    
    //zero mean first 
    double datamean = MLCompute::mean(data);
    for (bigint i = 0; i < N; i++) {
        data[i] = data[i] - datamean;
    }

    QVector<double> event_times;
    if (opts.detect_rms_window > 0) {
        printf("Computing RMS... \n");
        QVector<double> rmsdata(N);  
        rmsdata.fill(0);
        int midwin = opts.detect_rms_window/2;
        /*assert (M == 1) * not for multichannel */
        for (bigint i = 0; i < N; i++) {
            double val = 0;
            for (int wind = 0; wind < opts.detect_rms_window; wind++) {
                val += X.value(1, i + wind) * X.value(1, i + wind);
            }
            rmsdata[i+midwin] = sqrt(val);
        }
        /*double rmsdatamean = MLCompute::mean(rmsdata);
        for (bigint i = 0; i < N; i++) {
            rmsdata[i] = rmsdata[i] - rmsdatamean;
        }*/
        printf("Detecting events...\n");
        
        QVector<double> rms_event_times = P_detect_events::detect_events(rmsdata, opts.detect_threshold, opts.detect_interval, 1); // Power peaks are always +ve
        event_times = P_detect_events::align_events(data, rms_event_times, opts.sign, opts.detect_interval); // Sign now just refers to where we will align the events
        printf("%d events detected.\n", event_times.count());
    } else {
        printf("Detecting events...\n");
        event_times = P_detect_events::detect_events(data, opts.detect_threshold, opts.detect_interval, opts.sign);
        printf("%d events detected.\n", event_times.count());
    }

    printf("%d events detected.\n", event_times.count());

    if ((opts.subsample_factor) && (opts.subsample_factor < 1)) {
        printf("Subsampling by factor %g...\n", opts.subsample_factor);
        event_times = P_detect_events::subsample_events(event_times, opts.subsample_factor);
    }

    printf("Creating result array...\n");
    Mda ret(1, event_times.count());
    for (bigint j = 0; j < event_times.count(); j++) {
        ret.setValue(event_times[j], j);
    }
    printf("Writing result...\n");
    return ret.write64(event_times_out);
}

namespace P_detect_events {
    QVector<double> align_events(const QVector<double>& X, const QVector<double>& ptimes, int sign, double detect_interval)
    {
        int winaround = 10; //hard-coded for now
        int winlen = winaround * 2;
        bigint N = ptimes.size();
        QVector<double> pwin(winlen);
        QVector<double> times;
        double maxv = 0;
        double maxind = 0;
        double maxindprev = 0;
        double val = 0;

        for (bigint n = 0; n < N; n++) {
            bigint Xpos = ptimes[n] - winaround;
            maxv = 0;
            maxind = 0;
            val = 0; 
            // get trace around spike time
            for (int w = 0; w < winlen; w++) {
                val = X[Xpos + w];
                if (sign < 0) // just detect -ve peaks
                    val = -val;
                else if (sign == 0) // detect all peaks
                    val = fabs(val);
                if (val > maxv) { // set peak in trace
                    maxv = val;
                    maxind = Xpos + w;
                }
            }
            if (maxind - maxindprev > detect_interval) {
                times << maxind; 
                maxindprev = maxind; 
            }
        }
        return times;
    }
    
    QVector<double> detect_events(const QVector<double>& X, double detect_threshold, double detect_interval, int sign)
    {
        double mean = MLCompute::mean(X);
        double stdev = MLCompute::stdev(X);
        double threshold2 = detect_threshold * stdev;

        bigint N = X.count();
        QVector<bigint> to_use(N);
        to_use.fill(0);
        bigint last_best_ind = 0;
        double last_best_val = 0;
        for (bigint n = 0; n < N; n++) { // n is current sample 
            double val = X[n]; //- mean); move mean to previous stage
            if (sign < 0) // just detect -ve peaks
                val = -val;
            else if (sign == 0) // detect all peaks
                val = fabs(val);
            if (n - last_best_ind > detect_interval) // look for new peaks
                last_best_val = 0;
            if (val >= threshold2) { // if above detection threshold
                if (last_best_val > 0) { // have we already detected a possible peak
                    if (val > last_best_val) { // is this one bigger than other possible peaks
                        to_use[n] = 1; // set new possible peak
                        to_use[last_best_ind] = 0; // unset previous possible peak
                        last_best_ind = n; // set ind and val for new peak
                        last_best_val = val;
                    }
                }
                else { // this is the first possible peak in this window
                    if (val > 0) {
                        to_use[n] = 1;
                        last_best_ind = n;
                        last_best_val = val;
                    }
                }
            }
        }
        QVector<double> times; // get times from boolean time array
        for (bigint n = 0; n < N; n++) {
            if (to_use[n]) {
                times << n;
            }
        }
        return times;
    }

    double pseudorandomnumber(double i)
    {
        double ret = sin(i + cos(i));
        ret = (ret + 5) - (bigint)(ret + 5);
        return ret;
    }

    QVector<double> subsample_events(const QVector<double>& X, double subsample_factor)
    {
        QVector<double> ret;
        for (bigint i = 0; i < X.count(); i++) {
            double randnum = pseudorandomnumber(i);
            if (randnum <= subsample_factor)
                ret << X[i];
        }
        return ret;
    }
}
