/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
*******************************************************/

#ifndef MVCLUSTERDETAILWIDGET_H
#define MVCLUSTERDETAILWIDGET_H

#include "diskreadmda.h"
#include "mda.h"
#include <QWidget>
#include "mvabstractview.h"

/** \class ClusterDetailView
 *  \brief Display a view of each cluster -- mainly the template shapes and some stats
 *
 *  The user may click to change the current cluster, or use the Ctrl/Shift keys to select multiple clusters.
 */

class ClusterView;
class ClusterDetailViewPrivate;
class ClusterDetailView : public MVAbstractView {
    Q_OBJECT
public:
    friend class ClusterDetailViewPrivate;
    friend class ClusterView;
    ClusterDetailView(MVAbstractContext* context);
    virtual ~ClusterDetailView();

    void prepareCalculation() Q_DECL_OVERRIDE;
    void runCalculation() Q_DECL_OVERRIDE;
    void onCalculationFinished() Q_DECL_OVERRIDE;

    void zoomAllTheWayOut();

    ///Create an image of the current view
    QImage renderImage(int W = 0, int H = 0);

    QJsonObject exportStaticView() Q_DECL_OVERRIDE;
    void loadStaticView(const QJsonObject& X) Q_DECL_OVERRIDE;

protected:
    void leaveEvent(QEvent*);
    void paintEvent(QPaintEvent* evt);
    void keyPressEvent(QKeyEvent* evt);
    void mousePressEvent(QMouseEvent* evt);
    void mouseReleaseEvent(QMouseEvent* evt);
    void mouseMoveEvent(QMouseEvent* evt);
    void mouseDoubleClickEvent(QMouseEvent* evt);
    void wheelEvent(QWheelEvent* evt);
    void prepareMimeData(QMimeData& mimeData, const QPoint& pos);
signals:
    ///A cluster has been double-clicked (or enter pressed?)
    //void signalTemplateActivated();
private slots:
    //void slot_context_menu(const QPoint& pos);
    //void slot_export_waveforms();
    void slot_export_image();
    void slot_toggle_stdev_shading();
    void slot_zoom_in();
    void slot_zoom_out();
    void slot_vertical_zoom_in();
    void slot_vertical_zoom_out();
    void slot_export_static_view();
    void slot_update_sort_order();
    void slot_view_properties();
    void slot_export_template_waveforms();
    void slot_export_template_waveform_stdevs();

private:
    ClusterDetailViewPrivate* d;
};

#endif // MVCLUSTERDETAILWIDGET_H
