/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
*******************************************************/

#ifndef MVOVERVIEW2WIDGET_H
#define MVOVERVIEW2WIDGET_H

#include <QWidget>
#include "mda.h"
#include <QTabWidget>
#include "mvutils.h"

/** \class MVOverview2Widget
 *  \brief The main window (for now) showing an overview of the results of sorting by providing a number of interactive and synchronized views.
 *
 *  Presents user with a rich set of views. Cross-correlograms, raw data, cluster details, rotatable 3D views, firing rate vs. time view, etc.
 */

class MVOverview2WidgetPrivate;
class MVOverview2Widget : public QWidget {
    Q_OBJECT
public:
    friend class MVOverview2WidgetPrivate;
    MVOverview2Widget(QWidget* parent = 0);
    virtual ~MVOverview2Widget();
    ///The path to the timeseries that was sorted. For example, raw, filtered, or pre-processed. Usually all three of these are set, so user can choose between them in dropdown selection box
    Q_INVOKABLE void addTimeseriesPath(const QString& name, const QString& path);
    ///The name of the timeseries being viewed... corresponds to name in addTimeseriesPath()
    Q_INVOKABLE void setCurrentTimeseriesName(const QString& name);
    ///Set the path to the results of sorting.
    Q_INVOKABLE void setFiringsPath(const QString& firings);
    ///The sample rate for the dataset
    Q_INVOKABLE void setSampleRate(float freq);
    ///Open the initial views
    Q_INVOKABLE void setDefaultInitialization();
    ///Corresponds to MVFiringRateView::setEpochs()
    Q_INVOKABLE void setEpochs(const QList<Epoch>& epochs);
    Q_INVOKABLE void setParameterValue(const QString& name, const QVariant& value);
    Q_INVOKABLE QVariant getParameterValue(const QString& name, const QVariant& defaultval);
    Q_INVOKABLE void clickButton(const QString& name);
    Q_INVOKABLE QImage generateImage(const QMap<QString, QVariant>& params);
    Q_INVOKABLE int getMaxLabel();
    Q_INVOKABLE void setMscmdServerUrl(const QString& url);

protected:
    void resizeEvent(QResizeEvent* evt);
    void keyPressEvent(QKeyEvent* evt);

signals:

public slots:

private slots:
    void slot_control_panel_button_clicked(QString str);
    void slot_control_panel_combobox_activated(QString str);
    void slot_auto_correlogram_activated(int k);
    //void slot_templates_clicked();
    void slot_details_current_k_changed();
    void slot_details_selected_ks_changed();
    void slot_details_template_activated();
    void slot_cross_correlogram_current_index_changed();
    void slot_cross_correlogram_selected_indices_changed();
    void slot_clips_view_current_event_changed();
    void slot_clips_widget_current_event_changed();
    void slot_cluster_view_current_event_changed();
    //void slot_cross_correlogram_computer_finished();

private:
    MVOverview2WidgetPrivate* d;
};

/*
class CustomTabWidget : public QTabWidget {
	Q_OBJECT
public:
	MVOverview2Widget *q;
	CustomTabWidget(MVOverview2Widget *q);
protected:
	void mousePressEvent(QMouseEvent *evt);
private slots:
	void slot_tab_close_requested(int num);
	void slot_tab_bar_clicked();
	void slot_tab_bar_double_clicked();
	void slot_switch_to_other_tab_widget();
};
*/

#endif // MVOVERVIEW2WIDGET_H
