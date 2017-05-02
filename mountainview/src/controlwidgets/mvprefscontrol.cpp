#include "mvprefscontrol.h"

#include <QGridLayout>
#include <QLabel>
#include <QLineEdit>
#include <QMap>
#include <QTimer>
#include <mvcontext.h>
#include "mlcommon.h"

class MVPrefsControlPrivate {
public:
    MVPrefsControl* q;
};

MVPrefsControl::MVPrefsControl(MVAbstractContext* context, MVMainWindow* mw)
    : MVAbstractControl(context, mw)
{
    d = new MVPrefsControlPrivate;
    d->q = this;

    QGridLayout* glayout = new QGridLayout;
    int row = 0;
    {
        QWidget* X = this->createStringControl("visible_channels");
        X->setToolTip("Comma-separated list of channels to view");
        connect(context, SIGNAL(visibleChannelsChanged()), this, SLOT(updateControls()));
        glayout->addWidget(new QLabel("Visible channels:"), row, 0);
        glayout->addWidget(X, row, 1);
        row++;
    }
    /*
    {
        QWidget* X = this->createChoicesControl("timeseries_for_spikespray");
        connect(context, SIGNAL(timeseriesNamesChanged()), this, SLOT(updateControls()));
        connect(context, SIGNAL(currentTimeseriesChanged()), this, SLOT(updateControls()));
        glayout->addWidget(new QLabel("For spikespray:"), row, 0);
        glayout->addWidget(X, row, 1);
        row++;
    }
    */
    {
        QWidget* X = this->createIntControl("clip_size");
        X->setToolTip("The clip size used for all the views");
        context->onOptionChanged("clip_size", this, SLOT(updateControls()));
        glayout->addWidget(new QLabel("Clip size:"), row, 0);
        glayout->addWidget(X, row, 1);
        row++;
    }
    {
        QWidget* X = this->createDoubleControl("cc_max_dt_msec");
        X->setToolTip("Maximum dt in milliseconds for cross-correlograms.");
        context->onOptionChanged("cc_max_dt_msec", this, SLOT(updateControls()));
        glayout->addWidget(new QLabel("CC max. dt (msec):"), row, 0);
        glayout->addWidget(X, row, 1);
        row++;
    }
    {
        QWidget* X = this->createDoubleControl("cc_log_time_constant_msec");
        context->onOptionChanged("cc_log_time_constant_msec", this, SLOT(updateControls()));
        glayout->addWidget(new QLabel("CC log time constant (msec):"), row, 0);
        glayout->addWidget(X, row, 1);
        row++;
    }
    {
        QWidget* X = this->createDoubleControl("cc_bin_size_msec");
        context->onOptionChanged("cc_bin_size_msec", this, SLOT(updateControls()));
        glayout->addWidget(new QLabel("CC bin size (msec):"), row, 0);
        glayout->addWidget(X, row, 1);
        row++;
    }
    {
        QWidget* X = this->createDoubleControl("cc_max_est_data_size");
        context->onOptionChanged("cc_max_est_data_size", this, SLOT(updateControls()));
        glayout->addWidget(new QLabel("CC max est data size:"), row, 0);
        glayout->addWidget(X, row, 1);
        row++;
    }
    {
        QWidget* X = this->createDoubleControl("amp_thresh_display");
        context->onOptionChanged("amp_thresh_display", this, SLOT(updateControls()));
        glayout->addWidget(new QLabel("Amp. thresh. (display):"), row, 0);
        glayout->addWidget(X, row, 1);
        row++;
    }
    /*
    {
        QWidget* X = this->createChoicesControl("discrim_hist_method");
        QStringList choices;
        choices << "centroid"
                << "svm";
        this->setChoices("discrim_hist_method", choices);
        context->onOptionChanged("discrim_hist_method", this, SLOT(updateControls()));
        glayout->addWidget(new QLabel("Discrim hist method:"), row, 0);
        glayout->addWidget(X, row, 1);
        row++;
    }
    */
    this->setLayout(glayout);

    updateControls();
}

MVPrefsControl::~MVPrefsControl()
{
    delete d;
}

QString MVPrefsControl::title() const
{
    return "Preferences";
}

void MVPrefsControl::updateContext()
{
    MVContext* c = qobject_cast<MVContext*>(mvContext());
    Q_ASSERT(c);

    //QString ts4ss = this->controlValue("timeseries_for_spikespray").toString();
    //if (ts4ss == "Default timeseries")
    //    ts4ss = "";
    c->setVisibleChannels(MLUtil::stringListToIntList(this->controlValue("visible_channels").toString().split(",", QString::SkipEmptyParts)));
    //c->setOption("timeseries_for_spikespray", ts4ss);
    c->setOption("clip_size", this->controlValue("clip_size").toInt());
    c->setOption("cc_max_dt_msec", this->controlValue("cc_max_dt_msec").toDouble());
    c->setOption("cc_log_time_constant_msec", this->controlValue("cc_log_time_constant_msec").toDouble());
    c->setOption("cc_bin_size_msec", this->controlValue("cc_bin_size_msec").toDouble());
    c->setOption("cc_max_est_data_size", this->controlValue("cc_max_est_data_size").toDouble());
    c->setOption("amp_thresh_display", this->controlValue("amp_thresh_display").toDouble());
    //c->setOption("discrim_hist_method", this->controlValue("discrim_hist_method").toString());
}

void MVPrefsControl::updateControls()
{
    MVContext* c = qobject_cast<MVContext*>(mvContext());
    Q_ASSERT(c);

    QStringList choices = c->timeseriesNames();
    choices.insert(0, "Default timeseries");
    this->setControlValue("visible_channels", MLUtil::intListToStringList(c->visibleChannels()).join(","));
    //this->setChoices("timeseries_for_spikespray", choices);
    this->setControlValue("clip_size", c->option("clip_size").toInt());
    this->setControlValue("cc_max_dt_msec", c->option("cc_max_dt_msec").toDouble());
    this->setControlValue("cc_log_time_constant_msec", c->option("cc_log_time_constant_msec").toDouble());
    this->setControlValue("cc_bin_size_msec", c->option("cc_bin_size_msec").toDouble());
    this->setControlValue("cc_max_est_data_size", c->option("cc_max_est_data_size").toDouble());
    this->setControlValue("amp_thresh_display", c->option("amp_thresh_display").toDouble());
    //this->setControlValue("discrim_hist_method", c->option("discrim_hist_method"));
}
