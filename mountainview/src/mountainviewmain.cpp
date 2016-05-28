#include <QApplication>
#include <QDebug>
#include <qdatetime.h>
#include <QFile>
#include <QFileInfo>
#include <QMessageBox>
#include <QProcess>
#include <QStringList>
#include "textfile.h"
#include "usagetracking.h"
#include "mda.h"
#include <QDesktopServices>
#include <QDesktopWidget>
#include <QDir>
#include <QImageWriter>
#include "commandlineparams.h"
#include "diskarraymodel_new.h"
#include "histogramview.h"
#include "mvlabelcomparewidget.h"
#include "mvoverview2widget.h"
#include "sstimeserieswidget.h"
#include "sstimeseriesview.h"
#include "mvclusterwidget.h"
#include "closemehandler.h"
#include "remotereadmda.h"
#include "taskprogress.h"

//void run_export_instructions(MVOverview2Widget* W, const QStringList& instructions);

/// TODO provide mountainview usage information
/// TODO auto correlograms for selected clusters
/// TODO figure out what to do when #channels and/or #clusters is huge
/// TODO 0.9.1 -- make sure to handle merging with other views, such as clips etc. Make elegant way

int main(int argc, char* argv[])
{
    QApplication a(argc, argv);
    CloseMeHandler::start();

    setbuf(stdout, 0);

    /// Witold I don't want to do this here! It should be in the taskprogress.h. What can I do?
    qRegisterMetaType<TaskInfo>();

    CLParams CLP = commandlineparams(argc, argv);

    QString mv_fname;
    if (CLP.unnamed_parameters.value(0).endsWith(".mv")) {
        mv_fname = CLP.unnamed_parameters.value(0);
    }

    if (CLP.unnamed_parameters.value(0) == "unit_test") {
        QString arg2 = CLP.unnamed_parameters.value(1);
        if (arg2 == "remotereadmda") {
            unit_test_remote_read_mda();
        } else if (arg2 == "remotereadmda2") {
            QString arg3 = CLP.unnamed_parameters.value(2, "http://localhost:8000/firings.mda");
            unit_test_remote_read_mda_2(arg3);
        } else {
            qWarning() << "No such unit test: " + arg2;
        }
        return 0;
    }

    QString mode = CLP.named_parameters.value("mode", "overview2").toString();
    if (mode == "overview2") {
        printf("overview2...\n");
        QString raw_path = CLP.named_parameters["raw"].toString();
        QString pre_path = CLP.named_parameters["pre"].toString();
        QString filt_path = CLP.named_parameters["filt"].toString();
        QString firings_path = CLP.named_parameters["firings"].toString();
        double samplerate = CLP.named_parameters.value("samplerate", 20000).toDouble();
        QString epochs_path = CLP.named_parameters["epochs"].toString();
        QString window_title = CLP.named_parameters["window_title"].toString();
        MVOverview2Widget* W = new MVOverview2Widget;
        W->setMLProxyUrl(CLP.named_parameters.value("mlproxy_url", "").toString());
        {
            W->setWindowTitle(window_title);
            W->show();
            W->move(QApplication::desktop()->screen()->rect().topLeft() + QPoint(200, 200));

            int W0 = 1400, H0 = 1000;
            QRect geom = QApplication::desktop()->geometry();
            if ((geom.width() - 100 < W0) || (geom.height() - 100 < H0)) {
                //W->showMaximized();
                W->resize(geom.width() - 100, geom.height() - 100);
            } else {
                W->resize(W0, H0);
            }

            qApp->processEvents();
        }

        if (!mv_fname.isEmpty()) {
            W->loadMVFile(mv_fname);
        }

        if (!pre_path.isEmpty()) {
            W->addTimeseriesPath("Preprocessed Data", pre_path);
        }
        if (!filt_path.isEmpty()) {
            W->addTimeseriesPath("Filtered Data", filt_path);
        }
        if (!raw_path.isEmpty()) {
            W->addTimeseriesPath("Raw Data", raw_path);
        }

        if (!epochs_path.isEmpty()) {
            QList<Epoch> epochs = read_epochs(epochs_path);
            W->setEpochs(epochs);
        }
        if (window_title.isEmpty())
            window_title = pre_path;
        if (window_title.isEmpty())
            window_title = filt_path;
        if (window_title.isEmpty())
            window_title = raw_path;
        if (!firings_path.isEmpty()) {
            W->setFiringsPath(firings_path);
        }
        if (samplerate) {
            W->setSampleRate(samplerate);
        }

        W->setDefaultInitialization();
    } else if (mode == "spikespy") {
        printf("spikespy...\n");
        QString timeseries_path = CLP.named_parameters["timeseries"].toString();
        QString firings_path = CLP.named_parameters["firings"].toString();
        double samplerate = CLP.named_parameters["samplerate"].toDouble();
        SSTimeSeriesWidget* W = new SSTimeSeriesWidget;
        SSTimeSeriesView* V = new SSTimeSeriesView;
        V->setSampleRate(samplerate);
        DiskArrayModel_New* DAM = new DiskArrayModel_New;
        DAM->setPath(timeseries_path);
        V->setData(DAM, true);
        Mda firings;
        firings.read(firings_path);
        QList<long> times, labels;
        for (int i = 0; i < firings.N2(); i++) {
            times << (long)firings.value(1, i);
            labels << (long)firings.value(2, i);
        }
        V->setTimesLabels(times, labels);
        W->addView(V);
        W->show();
        W->move(QApplication::desktop()->screen()->rect().topLeft() + QPoint(200, 200));
        W->resize(1800, 1200);
    }

    int ret = a.exec();

    printf("Number of files open: %d, number of unfreed mallocs: %d, number of unfreed megabytes: %g\n", jnumfilesopen(), jmalloccount(), (int)jbytesallocated() * 1.0 / 1000000);

    return ret;
}
