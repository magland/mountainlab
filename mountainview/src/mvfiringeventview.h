/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
*******************************************************/

#ifndef MVFIRINGEVENTVIEW_H
#define MVFIRINGEVENTVIEW_H

#include <QWidget>
#include "diskreadmda.h"
#include "mvutils.h"

/** \class MVFiringEventView
 *  \brief Shows the firing rates (in color) of a selected neuron (or cluster) as a function of time. And shows distribution of amplitudes as well.
 */

class MVFiringEventViewPrivate;
class MVFiringEventView : public QWidget {
    Q_OBJECT
public:
    friend class MVFiringEventViewPrivate;
    MVFiringEventView();
    virtual ~MVFiringEventView();
    ///Set the RxL array of firing event info
    void setFirings(const Mda& firings);
    ///The sampling rate, (not sure if used right now)
    void setSampleRate(double ff); //in Hz
    ///The currently selected event by user click (does it work?)
    void setCurrentEvent(MVEvent evt);
    ///Controls display. Not explained in detail right now.
    void setEpochs(const QList<Epoch>& epochs);
    QImage renderImage(int W, int H);
signals:

protected:
    void paintEvent(QPaintEvent* evt);
    void mouseReleaseEvent(QMouseEvent* evt);

private slots:
    void slot_update();
    void slot_context_menu(const QPoint& pos);

private:
    MVFiringEventViewPrivate* d;
};

#endif // MVFiringEventView_H
