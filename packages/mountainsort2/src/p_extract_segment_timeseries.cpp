/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
** Created: 2/23/2017
*******************************************************/

#include "p_extract_segment_timeseries.h"

#include <diskreadmda32.h>
#include <diskwritemda.h>
#include <mda32.h>

namespace P_extract_segment_timeseries {
QVector<bigint> get_durations(QStringList timeseries_list);
Mda32 extract_channels_from_chunk(const Mda32& chunk, const QList<int>& channels);
}

bool p_extract_segment_timeseries(QString timeseries, QString timeseries_out, bigint t1, bigint t2, const QList<int>& channels)
{
    DiskReadMda32 X(timeseries);
    bigint M = X.N1();
    //bigint N=X.N2();
    bigint N2 = t2 - t1 + 1;

    qDebug().noquote() << QString("Extracting segment timeseries M=%1, N2=%2").arg(M).arg(N2);

    bigint M2 = M;
    if (!channels.isEmpty()) {
        M2 = channels.count();
    }

    //do it this way so we can specify the datatype
    DiskWriteMda Y;
    Y.open(X.mdaioHeader().data_type, timeseries_out, M2, N2);

    //conserve memory as of 3/1/17 -- jfm
    bigint chunk_size = 10000; //conserve memory!
    for (bigint t = 0; t < N2; t += chunk_size) {
        bigint sz = chunk_size;
        if (t + sz > N2)
            sz = N2 - t;
        Mda32 chunk;
        if (!X.readChunk(chunk, 0, t1 + t, M, sz)) {
            qWarning() << "Problem reading chunk.";
            return false;
        }
        if (!channels.isEmpty()) {
            chunk = P_extract_segment_timeseries::extract_channels_from_chunk(chunk, channels);
        }
        if (!Y.writeChunk(chunk, 0, t)) {
            qWarning() << "Problem writing chunk.";
            return false;
        }
    }

    return true;
}

bool p_extract_segment_timeseries_from_concat_list(QStringList timeseries_list, QString timeseries_out, bigint t1, bigint t2, const QList<int>& channels)
{
    QString ts0 = timeseries_list.value(0);
    DiskReadMda32 X(ts0);
    bigint M = X.N1();
    if (!channels.isEmpty()) {
        M = channels.count();
    }
    QVector<bigint> durations = P_extract_segment_timeseries::get_durations(timeseries_list);
    QVector<bigint> start_timepoints, end_timepoints;
    start_timepoints << 0;
    end_timepoints << durations.value(0) - 1;
    for (bigint i = 1; i < durations.count(); i++) {
        start_timepoints << end_timepoints[i - 1] + 1;
        end_timepoints << start_timepoints[i] + durations[i] - 1;
    }

    bigint ii1 = 0;
    while ((ii1 < timeseries_list.count()) && (t1 > end_timepoints[ii1])) {
        ii1++;
    }
    bigint ii2 = ii1;
    while ((ii2 + 1 < timeseries_list.count()) && (start_timepoints[ii2 + 1] <= t2)) {
        ii2++;
    }

    Mda32 out(M, t2 - t1 + 1);
    for (bigint ii = ii1; ii <= ii2; ii++) {
        DiskReadMda32 Y(timeseries_list.value(ii));
        Mda32 chunk;

        //ttA,ttB, the range to read from Y
        //ssA, the position to write it in "out"
        bigint ttA, ttB, ssA;
        if (ii == ii1) {
            bigint offset = t1 - start_timepoints[ii];
            ssA = 0;
            if (ii1 < ii2) {
                ttA = offset;
                ttB = durations[ii] - 1;
            }
            else {
                ttA = offset;
                ttB = ttA + (t2 - t1 + 1) - 1;
            }
        }
        else if (ii == ii2) {
            ttA = 0;
            ttB = t2 - start_timepoints[ii];
            ssA = start_timepoints[ii] - t1;
        }
        else {
            ttA = 0;
            ttB = durations[ii] - 1;
            ssA = start_timepoints[ii] - t1;
        }

        if (!Y.readChunk(chunk, 0, ttA, Y.N1(), ttB - ttA + 1)) {
            qWarning() << "Problem reading chunk in extract_segment_timeseries from concat list";
            return false;
        }
        if (!channels.isEmpty()) {
            chunk = P_extract_segment_timeseries::extract_channels_from_chunk(chunk, channels);
        }
        out.setChunk(chunk, 0, ssA);
    }

    //do it this way so we can specify the datatype
    DiskWriteMda Y;
    if (!Y.open(X.mdaioHeader().data_type, timeseries_out, M, t2 - t1 + 1)) {
        qWarning() << "Error opening file for writing: " << timeseries_out << M << t2 - t1 + 1;
        return false;
    }
    if (!Y.writeChunk(out, 0, 0)) {
        qWarning() << "Problem writing chunk in extract_segment_timeseries";
        return false;
    }
    Y.close();

    return true;
}

namespace P_extract_segment_timeseries {
QVector<bigint> get_durations(QStringList timeseries_list)
{
    QVector<bigint> ret;
    foreach (QString ts, timeseries_list) {
        DiskReadMda32 X(ts);
        ret << X.N2();
    }
    return ret;
}
Mda32 extract_channels_from_chunk(const Mda32& chunk, const QList<int>& channels)
{
    Mda32 ret(channels.count(), chunk.N2());
    for (bigint i = 0; i < chunk.N2(); i++) {
        for (bigint j = 0; j < channels.count(); j++) {
            ret.set(chunk.get(channels[j] - 1, i), j, i);
        }
    }
    return ret;
}
}

bool p_extract_segment_firings(QString firings_path, QString firings_out_path, bigint t1, bigint t2)
{
    Mda firings(firings_path);

    bigint event_times_row = 1;
    if (firings.N1() == 1) { //must be event_times not firings
        event_times_row = 0;
    }

    QVector<bigint> inds;
    for (bigint i = 0; i < firings.N2(); i++) {
        double t0 = firings.value(event_times_row, i);
        if ((t1 <= t0) && (t0 <= t2))
            inds << i;
    }

    Mda ret(firings.N1(), inds.count());
    for (bigint i = 0; i < inds.count(); i++) {
        bigint j = inds[i];
        for (bigint r = 0; r < firings.N1(); r++) {
            ret.setValue(firings.value(r, j), r, i);
        }
        ret.setValue(ret.value(event_times_row, i) - t1, event_times_row, i);
    }
    return ret.write64(firings_out_path);
}
