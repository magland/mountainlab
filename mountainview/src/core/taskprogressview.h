/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland and Witold Wysota
** Created: 4/30/2016
*******************************************************/

#ifndef TASKPROGRESSVIEW_H
#define TASKPROGRESSVIEW_H

#include <QTreeView>

class TaskProgressViewPrivate;
class TaskProgressView : public QTreeView {
    Q_OBJECT
public:
    friend class TaskProgressViewPrivate;
    TaskProgressView();
    virtual ~TaskProgressView();
public slots:
    void copySelectedToClipboard();
    void showLogMessages(const QModelIndex& index);
protected:
    void updateGeometries();
private:
    TaskProgressViewPrivate* d;
};

#endif // TASKPROGRESSVIEW_H
