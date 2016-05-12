/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
*******************************************************/

#ifndef MVCLUSTERDETAILWIDGET_H
#define MVCLUSTERDETAILWIDGET_H

#include "diskreadmda.h"
#include "mda.h"
#include <QWidget>

/** \class MVClusterDetailWidget
 *  \brief Display a view of each cluster -- mainly the template shapes and some stats
 *
 *  The user may click to change the current cluster, or use the Ctrl/Shift keys to select multiple clusters.
 */

class MVClusterDetailWidgetPrivate;
class MVClusterDetailWidget : public QWidget {
    Q_OBJECT
public:
    friend class MVClusterDetailWidgetPrivate;
    MVClusterDetailWidget(QWidget* parent = 0);
    virtual ~MVClusterDetailWidget();

    void setMscmdServerUrl(const QString& url);

    ///Set the time series, from which the templates and stats will be derived
    void setTimeseries(DiskReadMda& X);

    ///Set the firings info, from which the templates and stats will be derived
    void setFirings(const DiskReadMda& X);

    ///The size of the templates to display
    void setClipSize(int T);

    ///This is important when we have split the clusters into amplitude shells, and we want to group them together
    void setGroupNumbers(const QList<int>& group_numbers);

    ///So we can display the firing rate (events per second)
    void setSampleRate(double freq);

    ///To make the color scheme uniform. TODO: handle this correctly
    void setChannelColors(const QList<QColor>& colors);

    ///To make the color scheme uniform. TODO: handle this correctly
    void setColors(const QMap<QString, QColor>& colors);

    ///The current label number (selected by user)
    int currentK();

    ///The selected label numbers (selected by user)
    QList<int> selectedKs();

    ///Set current label number (affects highlighting)
    void setCurrentK(int k);

    ///Set selected label numbers (affects highlighting)
    void setSelectedKs(const QList<int>& ks);

    ///Create an image of the current view
    QImage renderImage(int W = 0, int H = 0);

protected:
    void paintEvent(QPaintEvent* evt);
    void keyPressEvent(QKeyEvent* evt);
    void mousePressEvent(QMouseEvent* evt);
    void mouseReleaseEvent(QMouseEvent* evt);
    void mouseMoveEvent(QMouseEvent* evt);
    void mouseDoubleClickEvent(QMouseEvent* evt);
    void wheelEvent(QWheelEvent* evt);
signals:
    ///A cluster has been double-clicked (or enter pressed?)
    void signalTemplateActivated();
    ///The current label number has changed
    void signalCurrentKChanged();
    ///The selected label numbers have changed
    void signalSelectedKsChanged();
    ///The user has zoomed in (or out?). TODO: low-priority Is this needed?
    void signalZoomedIn();
private slots:
    void slot_context_menu(const QPoint& pos);
    void slot_calculator_finished();

private:
    MVClusterDetailWidgetPrivate* d;
};

#endif // MVCLUSTERDETAILWIDGET_H
