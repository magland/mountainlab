/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
** Created: 6/6/2016
*******************************************************/

#ifndef SPIKESPYWIDGET_H
#define SPIKESPYWIDGET_H

#include "diskreadmda.h"
#include "mvabstractcontext.h"

#include <QWidget>
#include <diskreadmda32.h>

struct SpikeSpyViewData {
    DiskReadMda32 timeseries;
    DiskReadMda firings;
};

class SpikeSpyWidgetPrivate;
class SpikeSpyWidget : public QWidget {
    Q_OBJECT
public:
    friend class SpikeSpyWidgetPrivate;
    SpikeSpyWidget(MVAbstractContext* context);
    virtual ~SpikeSpyWidget();
    void addView(const SpikeSpyViewData& data);

private slots:
    void slot_show_tasks();
    void slot_open_mountainview();
    void slot_view_clicked();

private:
    SpikeSpyWidgetPrivate* d;
};

#endif // SPIKESPYWIDGET_H
