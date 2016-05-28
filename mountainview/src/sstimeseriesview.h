/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
*******************************************************/

#ifndef SSTIMESERIESVIEW_H
#define SSTIMESERIESVIEW_H

#include <QWidget>
#include "plotarea.h" //for Vec2
#include "sslabelsmodel.h"
#include "ssabstractview.h"

class SSTimeSeriesViewPrivate;
class SSTimeSeriesView : public SSAbstractView {
    Q_OBJECT
public:
    friend class SSTimeSeriesViewPrivate;
    explicit SSTimeSeriesView(QWidget* parent = 0);
    ~SSTimeSeriesView();

    void setData(DiskArrayModel_New* data, bool is_owner = false);
    DiskArrayModel_New* data();
    //void setLabels(Mda *T,Mda *L);
    void setLabels(DiskReadMda* TL, bool is_owner = false);
    void setCompareLabels(DiskReadMda* TL, bool is_owner = false);
    //void setConnectZeros(bool val);
    void setClipMode(bool val);
    bool clipMode();

    void setChannelLabels(const QStringList& labels);
    void setUniformVerticalChannelSpacing(bool val);
    void setTimesLabels(const QList<long>& times, const QList<long>& labels);

    SSLabelsModel* getLabels();

    double currentValue();
    QString viewType();
    void setMarkerLinesVisible(bool val);
    QImage renderImage(int W, int H);

    SSTimeSeriesPlot* plot();

private
slots:
    void slot_request_move_to_timepoint(int t0);
    void slot_context_menu(const QPoint& pos);

signals:
    void requestCenterOnCursor();

protected:
    void keyPressEvent(QKeyEvent*);

protected:
private:
    SSTimeSeriesViewPrivate* d;

signals:

private
slots:
};

#endif // SSTIMESERIESVIEW_H
