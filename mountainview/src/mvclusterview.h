/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
*******************************************************/

#ifndef MVCLUSTERVIEW_H
#define MVCLUSTERVIEW_H

#include <QWidget>
#include "mda.h"
#include "mvutils.h"
#include "affinetransformation.h"

#define MVCV_MODE_HEAT_DENSITY 1
#define MVCV_MODE_LABEL_COLORS 2
#define MVCV_MODE_TIME_COLORS 3
#define MVCV_MODE_AMPLITUDE_COLORS 4

/** \class MVClusterView
 *  \brief A rotatable view of datapoints in 3D space
 *
 *  Several modes are available
 */

class MVClusterViewPrivate;
class MVClusterView : public QWidget {
    Q_OBJECT
public:
    friend class MVClusterViewPrivate;
    MVClusterView(QWidget* parent = 0);
    virtual ~MVClusterView();
    ///Set the 3xL data array. Note if first dimension is >3, it only uses first 3 dims.
    void setData(const Mda& X);
    ///Have we actually set the data yet?
    bool hasData();
    ///Set L times corresponding to datapoints for purpose of currentEvent()
    void setTimes(const QList<double>& times);
    ///Set L labels corresponding to datapoints for purpose of currentEvent()
    void setLabels(const QList<int>& labels);
    ///Set L amplitudes corresponding to datapoints for purpose of currentEvent()
    void setAmplitudes(const QList<double>& amps);
    ///Set the display mode (several options, not yet documented)
    void setMode(int mode);
    ///Select the data point corresponding to the given evt time and label
    void setCurrentEvent(MVEvent evt, bool do_emit = false);
    ///Return the event corresponding to the current data point (see setTimes(), setLabels())
    MVEvent currentEvent();
    ///The index of the current data point, corresponding to the array from setData()
    int currentEventIndex();
    ///The current rotation (affine transformation) for purpose of synchronizing multiple cluster views
    AffineTransformation transformation();
    ////Set the current rotation (affine transformation) for purpose of synchronizing multiple cluster views
    void setTransformation(const AffineTransformation& T);

    QImage renderImage(int W, int H);
signals:
    ///The current data point has changed!
    void currentEventChanged();
    ///The current rotation (affine transformation) has changed - useful for synchronizing multiple cluster views
    void transformationChanged();

private
slots:
    void slot_context_menu(const QPoint& pos);

protected:
    void paintEvent(QPaintEvent* evt);
    void mouseMoveEvent(QMouseEvent* evt);
    void mousePressEvent(QMouseEvent* evt);
    void mouseReleaseEvent(QMouseEvent* evt);
    void wheelEvent(QWheelEvent* evt);

private:
    MVClusterViewPrivate* d;
};

#endif // MVCLUSTERVIEW_H
