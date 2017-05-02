/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
*******************************************************/

#ifndef MVCLIPSVIEW_H
#define MVCLIPSVIEW_H

#include "mvutils.h"

#include <QWidget>
#include <diskreadmda32.h>
#include <mvabstractcontext.h>

/** \class MVClipsView
 *  \brief View a set of clips. Usually each clip contains a single spike.
 */

class MVClipsViewPrivate;
class MVClipsView : public QWidget {
    Q_OBJECT
public:
    friend class MVClipsViewPrivate;
    MVClipsView(MVAbstractContext* context);
    virtual ~MVClipsView();

    void setClips(const DiskReadMda32& clips);
    /// TODO: (MEDIUM) in mvclipsview implement times/labels for purpose of current event and labeling
    //void setTimes(const QVector<double>& times);
    //void setLabels(const QVector<int>& labels);
protected:
    void paintEvent(QPaintEvent* evt);
signals:
private slots:

private:
    MVClipsViewPrivate* d;
};

#endif // MVCLIPSVIEW_H
