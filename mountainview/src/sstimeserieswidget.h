/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
*******************************************************/

#ifndef SSTIMESERIESWIDGET_H
#define SSTIMESERIESWIDGET_H

#include <QMainWindow>
#include "ssabstractview.h"

class SSTimeSeriesWidgetPrivate;
class SSTimeSeriesWidget : public QWidget {
    Q_OBJECT
public:
    friend class SSTimeSeriesWidgetPrivate;
    explicit SSTimeSeriesWidget(QWidget* parent = 0);
    ~SSTimeSeriesWidget();

    void addView(SSAbstractView* V);
    void setClipData(const Mda& X);
    SSAbstractView* view(int index = 0);

private:
    SSTimeSeriesWidgetPrivate* d;

private
slots:
    void slot_current_x_changed();
    void slot_current_channel_changed();
    void slot_x_range_changed();
    void slot_replot_needed();
    void slot_selection_range_changed();
    void slot_view_clicked();

    void slot_center_on_cursor();

    void slot_navigation_instructions();

signals:

public
slots:
};

#endif // SSTIMESERIESWIDGET_H
