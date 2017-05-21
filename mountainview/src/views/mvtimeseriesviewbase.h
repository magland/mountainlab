/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
** Created: 6/13/2016
*******************************************************/

#ifndef MVTimeSeriesViewBaseBASE_H
#define MVTimeSeriesViewBaseBASE_H

#include <QWidget>
#include <diskreadmda.h>
#include <mvutils.h>
#include "mvabstractview.h"

struct mvtsv_colors {
    QColor marker_color = QColor(200, 0, 0, 120);
    QColor text_color = Qt::black;
    QColor axis_color = Qt::black;
    QColor background_color = Qt::white;
};

struct mvtsv_prefs {
    bool show_current_timepoint = true;
    bool show_time_axis = true;
    bool show_time_axis_ticks = true;
    bool show_time_axis_scale = true;
    bool show_time_axis_scale_only_1ms = false;

    int num_label_levels = 3;
    int label_font_height = 12;
    double mleft = 30, mright = 30;
    double mtop = 40, mbottom = 40;
    bool markers_visible = true;

    mvtsv_colors colors;
};

class MVTimeSeriesViewBasePrivate;
class MVTimeSeriesViewBase : public MVAbstractView {
    Q_OBJECT
public:
    friend class MVTimeSeriesViewBasePrivate;
    MVTimeSeriesViewBase(MVAbstractContext* context);
    virtual ~MVTimeSeriesViewBase();

    virtual void prepareCalculation() Q_DECL_OVERRIDE;
    virtual void runCalculation() Q_DECL_OVERRIDE;
    virtual void onCalculationFinished() Q_DECL_OVERRIDE;

    virtual void paintContent(QPainter* painter) = 0;
    void setNumTimepoints(int N); //called by subclass
    mvtsv_prefs prefs() const;
    void setPrefs(mvtsv_prefs prefs);

    void setLabelsToView(const QSet<int>& labels);
    void setActivated(bool val);
    void setMarkersVisible(bool val);
    void setMargins(double mleft, double mright, double mtop, double mbottom);

    void addSpecialEvents(const QList<MVEvent>& events);

    void setClipSize(int clip_size); //for when viewing clips

    void setTimeRange(MVRange range);
    MVRange timeRange();
    void setSelectedTimeRange(MVRange range);

    double amplitudeFactor() const;

protected:
    void resizeEvent(QResizeEvent* evt);
    void paintEvent(QPaintEvent* evt);
    void mousePressEvent(QMouseEvent* evt);
    void mouseReleaseEvent(QMouseEvent* evt);
    void mouseMoveEvent(QMouseEvent* evt);
    void wheelEvent(QWheelEvent* evt);
    void keyPressEvent(QKeyEvent* evt);

    // I wanted to make this protected, but ran into some trouble
public:
    QRectF contentGeometry();
    double time2xpix(double t) const;
    double xpix2time(double x) const;

signals:
    void clicked();

private slots:
    void slot_scroll_to_current_timepoint();
    void slot_zoom_in();
    void slot_zoom_out();

private:
    MVTimeSeriesViewBasePrivate* d;
};

#endif // MVTimeSeriesViewBaseBASE_H
