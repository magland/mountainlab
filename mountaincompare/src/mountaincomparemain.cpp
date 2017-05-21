#include <QApplication>
#include <QDebug>
#include <qdatetime.h>
#include <QFile>
#include <QFileInfo>
#include <QMessageBox>
#include <QProcess>
#include <QStringList>

#include "usagetracking.h"
#include "mda.h"
#include <QDesktopServices>
#include <QDesktopWidget>
#include <QDir>
#include <QImageWriter>
#include "closemehandler.h"
#include "remotereadmda.h"
#include "taskprogress.h"
#include "clipsviewplugin.h"

#include <QHBoxLayout>
#include <QJsonDocument>
#include <QRunnable>
#include <QSettings>
#include <QThreadPool>
#include <QtConcurrentRun>
#include <mvclipswidget.h>
#include <mvclusterwidget.h>
#include <mvopenviewscontrol.h>
#include <mvtimeseriescontrol.h>
#include <qprocessmanager.h>
#include <tabber.h>

#include "multiscaletimeseries.h"
#include "taskprogressview.h"
#include "mvcontrolpanel2.h"
#include "mvabstractcontrol.h"
#include "initialize_confusion_matrix.h"

#include "mccontext.h"
#include "mcviewfactories.h"
#include "mvmainwindow.h"
#include "clusterdetailplugin.h"

#include "objectregistry.h"
#include "signal.h"

#include <views/confusionmatrixview.h>
#include <QFileDialog>

void set_nice_size(QWidget* W);
QColor brighten(QColor col, int amount);
QList<QColor> generate_colors_ahb();

void sig_handler(int signum)
{
    (void)signum;
    QProcessManager* manager = ObjectRegistry::getObject<QProcessManager>();
    if (manager) {
        manager->closeAll();
    }
    abort();
}

int main(int argc, char* argv[])
{
    QApplication a(argc, argv);
    ObjectRegistry registry;

    //The process manager
    QProcessManager* processManager = new QProcessManager;
    registry.addAutoReleasedObject(processManager);
    signal(SIGINT, sig_handler);
    signal(SIGKILL, sig_handler);
    signal(SIGTERM, sig_handler);

    // make sure task progress monitor is instantiated in the main thread
    TaskManager::TaskProgressMonitor* monitor = TaskManager::TaskProgressMonitor::globalInstance();
    Q_UNUSED(monitor);
    CloseMeHandler::start();

    setbuf(stdout, 0);

    CLParams CLP(argc, argv);

    QStringList keys = CLP.named_parameters.keys();
    foreach (QString key, keys) {
        QVariant val = CLP.named_parameters[key].toString();
        if ((!QFile::exists(val.toString())) && (QFile::exists(val.toString() + ".prv"))) {
            val = val.toString() + ".prv";
            CLP.named_parameters[key] = val;
        }
    }

    QList<QColor> channel_colors;
    QStringList color_strings;
    color_strings
        << "#282828"
        << "#402020"
        << "#204020"
        << "#202070";
    for (int i = 0; i < color_strings.count(); i++)
        channel_colors << QColor(brighten(color_strings[i], 80));

    QList<QColor> label_colors = generate_colors_ahb();

    MCContext* context = new MCContext; //note that the view agent does not get deleted. :(
    context->setChannelColors(channel_colors);
    context->setClusterColors(label_colors);
    MVMainWindow* W = new MVMainWindow(context);

    if (CLP.named_parameters.contains("samplerate")) {
        context->setSampleRate(CLP.named_parameters.value("samplerate", 0).toDouble());
    }
    if (CLP.named_parameters.contains("firings1")) {
        QString firings_path = CLP.named_parameters["firings1"].toString();
        context->setFirings1(DiskReadMda(firings_path));
        W->setWindowTitle(firings_path);
    }
    if (CLP.named_parameters.contains("firings2")) {
        QString firings_path = CLP.named_parameters["firings2"].toString();
        context->setFirings2(DiskReadMda(firings_path));
        W->setWindowTitle(firings_path);
    }
    if (CLP.named_parameters.contains("raw")) {
        QString raw_path = CLP.named_parameters["raw"].toString();
        context->addTimeseries("Raw Data", DiskReadMda32(raw_path));
        context->setCurrentTimeseriesName("Raw Data");
    }
    if (CLP.named_parameters.contains("filt")) {
        QString filt_path = CLP.named_parameters["filt"].toString();
        context->addTimeseries("Filtered Data", DiskReadMda32(filt_path));
        context->setCurrentTimeseriesName("Filtered Data");
    }
    if (CLP.named_parameters.contains("pre")) {
        QString pre_path = CLP.named_parameters["pre"].toString();
        context->addTimeseries("Preprocessed Data", DiskReadMda32(pre_path));
        context->setCurrentTimeseriesName("Preprocessed Data");
    }
    /*
    if (CLP.named_parameters.contains("mlproxy_url")) {
        QString mlproxy_url = CLP.named_parameters.value("mlproxy_url", "").toString();
        context->setMLProxyUrl(mlproxy_url);
    }
    */
    if (CLP.named_parameters.contains("window_title")) {
        QString window_title = CLP.named_parameters["window_title"].toString();
        W->setWindowTitle(window_title);
    }
    if (CLP.named_parameters.contains("geom")) {
        QString geom_path = CLP.named_parameters["geom"].toString();
        ElectrodeGeometry eg = ElectrodeGeometry::loadFromGeomFile(geom_path);
        context->mvContext1()->setElectrodeGeometry(eg);
        context->mvContext2()->setElectrodeGeometry(eg);
    }

    /*
    if (CLP.named_parameters.contains("clusters")) {
        QStringList clusters_subset_str = CLP.named_parameters["clusters"].toString().split(",", QString::SkipEmptyParts);
        QList<int> clusters_subset;
        foreach (QString label, clusters_subset_str) {
            clusters_subset << label.toInt();
        }
        context->setClustersSubset(clusters_subset.toSet());
    }
    */

    //W->loadPlugin(new ClusterDetailPlugin());
    //W->loadPlugin(new ClipsViewPlugin());

    W->registerViewFactory(new ClusterDetail1Factory(W));
    W->registerViewFactory(new ClusterDetail2Factory(W));

    W->registerViewFactory(new MVSpikeSpray1Factory(W));
    W->registerViewFactory(new MVSpikeSpray2Factory(W));

    W->registerViewFactory(new MVTimeSeriesView1Factory(W));
    W->registerViewFactory(new MVTimeSeriesView2Factory(W));
    W->registerViewFactory(new MVTimeSeriesViewIntersectFactory(W));

    W->registerViewFactory(new ConfusionMatrixViewFactory(W));
    //W->registerViewFactory(new CompareClustersFactory(W));

    W->registerViewFactory(new MVPCAFeatures1Factory(W));
    W->registerViewFactory(new MVPCAFeatures2Factory(W));

    W->addControl(new MVOpenViewsControl(context, W), true);
    W->addControl(new MVTimeseriesControl(context->mvContext1(), W), true);

    TaskProgressView* TPV = new TaskProgressView;
    TPV->show();

    Initialize_confusion_matrix* ICM = new Initialize_confusion_matrix; //need pointer to go into lambda expression
    ICM->firings1 = context->firings1().makePath();
    ICM->firings2 = context->firings2().makePath();
    ICM->relabel = CLP.named_parameters.contains("relabel");
    QObject::connect(ICM, &Initialize_confusion_matrix::finishedInGui, [=]() {
        if (ICM->relabel)
            context->setFirings2(DiskReadMda(ICM->firings2_relabeled));
        context->setConfusionMatrix(DiskReadMda(ICM->confusion_matrix));
        context->setMatchedFirings(DiskReadMda(ICM->matched_firings));
        context->setLabelMap(DiskReadMda(ICM->label_map));
        DiskReadMda CM(ICM->confusion_matrix);
        int Kmax=qMax(CM.N1()-1,CM.N2()-1);
        QList<int> clusters_to_force_show;
        for (int kk=1; kk<=Kmax; kk++) {
            clusters_to_force_show << kk;
        }
        context->setClustersToForceShow(clusters_to_force_show);
        W->setCurrentContainerName("north");
        W->openView("open-cluster-details-1");
        W->setCurrentContainerName("south");
        W->openView("open-cluster-details-2");
        delete TPV;
        delete ICM;
        set_nice_size(W);
        W->show();
        W->openView("open-confusion-matrix");
    });
    ICM->start();

    a.processEvents();

    return a.exec();
}

QColor brighten(QColor col, int amount)
{
    const int r = qBound(0, col.red() + amount, 255);
    const int g = qBound(0, col.green() + amount, 255);
    const int b = qBound(0, col.blue() + amount, 255);
    return QColor(r, g, b, col.alpha());
}
/*
 * The following list can be auto-generated using
 * https://github.com/magland/mountainlab_devel/blob/dev-jfm-3/view/ncolorpicker.m
 * Run with no arguments
 */
static float colors_ahb[256][3] = {
    { 0.201, 0.000, 1.000 },
    { 0.467, 1.000, 0.350 },
    { 1.000, 0.000, 0.761 },
    { 0.245, 0.700, 0.656 },
    { 1.000, 0.839, 0.000 },
    { 0.746, 0.350, 1.000 },
    { 0.000, 1.000, 0.059 },
    { 0.700, 0.245, 0.435 },
    { 0.000, 0.555, 1.000 },
    { 0.782, 1.000, 0.350 },
    { 0.894, 0.000, 1.000 },
    { 0.245, 0.700, 0.394 },
    { 1.000, 0.055, 0.000 },
    { 0.350, 0.353, 1.000 },
    { 0.296, 1.000, 0.000 },
    { 0.700, 0.245, 0.641 },
    { 0.000, 1.000, 0.708 },
    { 1.000, 0.747, 0.350 },
    { 0.458, 0.000, 1.000 },
    { 0.262, 0.700, 0.245 },
    { 1.000, 0.000, 0.575 },
    { 0.350, 0.861, 1.000 },
    { 0.855, 1.000, 0.000 },
    { 0.604, 0.245, 0.700 },
    { 0.000, 1.000, 0.207 },
    { 1.000, 0.350, 0.450 },
    { 0.000, 0.225, 1.000 },
    { 0.441, 0.700, 0.245 },
    { 1.000, 0.000, 0.968 },
    { 0.350, 1.000, 0.697 },
    { 1.000, 0.378, 0.000 },
    { 0.374, 0.245, 0.700 },
    { 0.135, 1.000, 0.000 },
    { 1.000, 0.350, 0.811 },
    { 0.000, 1.000, 0.994 },
    { 0.700, 0.670, 0.245 },
    { 0.667, 0.000, 1.000 },
    { 0.350, 1.000, 0.416 },
    { 1.000, 0.000, 0.344 },
    { 0.245, 0.452, 0.700 },
    { 0.590, 1.000, 0.000 },
    { 0.958, 0.350, 1.000 },
    { 0.000, 1.000, 0.384 },
    { 0.700, 0.313, 0.245 },
    { 0.086, 0.000, 1.000 },
    { 0.509, 1.000, 0.350 },
    { 1.000, 0.000, 0.825 },
    { 0.245, 0.700, 0.604 },
    { 1.000, 0.710, 0.000 },
    { 0.692, 0.350, 1.000 },
    { 0.000, 1.000, 0.005 },
    { 0.700, 0.245, 0.477 },
    { 0.000, 0.688, 1.000 },
    { 0.851, 1.000, 0.350 },
    { 0.835, 0.000, 1.000 },
    { 0.245, 0.700, 0.361 },
    { 1.000, 0.000, 0.067 },
    { 0.350, 0.434, 1.000 },
    { 0.371, 1.000, 0.000 },
    { 0.700, 0.245, 0.667 },
    { 0.000, 1.000, 0.606 },
    { 1.000, 0.661, 0.350 },
    { 0.361, 0.000, 1.000 },
    { 0.287, 0.700, 0.245 },
    { 1.000, 0.000, 0.654 },
    { 0.350, 0.944, 1.000 },
    { 0.973, 1.000, 0.000 },
    { 0.573, 0.245, 0.700 },
    { 0.000, 1.000, 0.145 },
    { 1.000, 0.350, 0.522 },
    { 0.000, 0.356, 1.000 },
    { 0.481, 0.700, 0.245 },
    { 0.977, 0.000, 1.000 },
    { 0.350, 1.000, 0.640 },
    { 1.000, 0.247, 0.000 },
    { 0.324, 0.245, 0.700 },
    { 0.196, 1.000, 0.000 },
    { 1.000, 0.350, 0.855 },
    { 0.000, 1.000, 0.875 },
    { 0.700, 0.612, 0.245 },
    { 0.589, 0.000, 1.000 },
    { 0.350, 1.000, 0.380 },
    { 1.000, 0.000, 0.442 },
    { 0.245, 0.513, 0.700 },
    { 0.691, 1.000, 0.000 },
    { 0.922, 0.350, 1.000 },
    { 0.000, 1.000, 0.308 },
    { 0.700, 0.256, 0.245 },
    { 0.000, 0.035, 1.000 },
    { 0.554, 1.000, 0.350 },
    { 1.000, 0.000, 0.884 },
    { 0.245, 0.700, 0.555 },
    { 1.000, 0.578, 0.000 },
    { 0.632, 0.350, 1.000 },
    { 0.050, 1.000, 0.000 },
    { 0.700, 0.245, 0.516 },
    { 0.000, 0.818, 1.000 },
    { 0.925, 1.000, 0.350 },
    { 0.772, 0.000, 1.000 },
    { 0.245, 0.700, 0.332 },
    { 1.000, 0.000, 0.183 },
    { 0.350, 0.517, 1.000 },
    { 0.453, 1.000, 0.000 },
    { 0.700, 0.245, 0.692 },
    { 0.000, 1.000, 0.512 },
    { 1.000, 0.574, 0.350 },
    { 0.256, 0.000, 1.000 },
    { 0.313, 0.700, 0.245 },
    { 1.000, 0.000, 0.727 },
    { 0.350, 1.000, 0.976 },
    { 1.000, 0.903, 0.000 },
    { 0.540, 0.245, 0.700 },
    { 0.000, 1.000, 0.087 },
    { 1.000, 0.350, 0.590 },
    { 0.000, 0.489, 1.000 },
    { 0.524, 0.700, 0.245 },
    { 0.922, 0.000, 1.000 },
    { 0.350, 1.000, 0.587 },
    { 1.000, 0.118, 0.000 },
    { 0.271, 0.245, 0.700 },
    { 0.262, 1.000, 0.000 },
    { 1.000, 0.350, 0.896 },
    { 0.000, 1.000, 0.762 },
    { 0.700, 0.553, 0.245 },
    { 0.503, 0.000, 1.000 },
    { 0.356, 1.000, 0.350 },
    { 1.000, 0.000, 0.533 },
    { 0.245, 0.573, 0.700 },
    { 0.799, 1.000, 0.000 },
    { 0.883, 0.350, 1.000 },
    { 0.000, 1.000, 0.239 },
    { 0.700, 0.245, 0.289 },
    { 0.000, 0.161, 1.000 },
    { 0.604, 1.000, 0.350 },
    { 1.000, 0.000, 0.941 },
    { 0.245, 0.700, 0.510 },
    { 1.000, 0.445, 0.000 },
    { 0.568, 0.350, 1.000 },
    { 0.106, 1.000, 0.000 },
    { 0.700, 0.245, 0.551 },
    { 0.000, 0.945, 1.000 },
    { 1.000, 0.997, 0.350 },
    { 0.704, 0.000, 1.000 },
    { 0.245, 0.700, 0.304 },
    { 1.000, 0.000, 0.292 },
    { 0.350, 0.603, 1.000 },
    { 0.542, 1.000, 0.000 },
    { 0.683, 0.245, 0.700 },
    { 0.000, 1.000, 0.425 },
    { 1.000, 0.489, 0.350 },
    { 0.145, 0.000, 1.000 },
    { 0.341, 0.700, 0.245 },
    { 1.000, 0.000, 0.793 },
    { 0.350, 1.000, 0.900 },
    { 1.000, 0.775, 0.000 },
    { 0.504, 0.245, 0.700 },
    { 0.000, 1.000, 0.032 },
    { 1.000, 0.350, 0.653 },
    { 0.000, 0.622, 1.000 },
    { 0.571, 0.700, 0.245 },
    { 0.865, 0.000, 1.000 },
    { 0.350, 1.000, 0.539 },
    { 1.000, 0.000, 0.006 },
    { 0.245, 0.275, 0.700 },
    { 0.333, 1.000, 0.000 },
    { 1.000, 0.350, 0.934 },
    { 0.000, 1.000, 0.656 },
    { 0.700, 0.493, 0.245 },
    { 0.410, 0.000, 1.000 },
    { 0.392, 1.000, 0.350 },
    { 1.000, 0.000, 0.616 },
    { 0.245, 0.632, 0.700 },
    { 0.914, 1.000, 0.000 },
    { 0.841, 0.350, 1.000 },
    { 0.000, 1.000, 0.175 },
    { 0.700, 0.245, 0.341 },
    { 0.000, 0.290, 1.000 },
    { 0.658, 1.000, 0.350 },
    { 1.000, 0.000, 0.995 },
    { 0.245, 0.700, 0.468 },
    { 1.000, 0.312, 0.000 },
    { 0.499, 0.350, 1.000 },
    { 0.165, 1.000, 0.000 },
    { 0.700, 0.245, 0.584 },
    { 0.000, 1.000, 0.933 },
    { 1.000, 0.916, 0.350 },
    { 0.629, 0.000, 1.000 },
    { 0.245, 0.700, 0.278 },
    { 1.000, 0.000, 0.394 },
    { 0.350, 0.689, 1.000 },
    { 0.639, 1.000, 0.000 },
    { 0.658, 0.245, 0.700 },
    { 0.000, 1.000, 0.346 },
    { 1.000, 0.406, 0.350 },
    { 0.027, 0.000, 1.000 },
    { 0.372, 0.700, 0.245 },
    { 1.000, 0.000, 0.855 },
    { 0.350, 1.000, 0.828 },
    { 1.000, 0.644, 0.000 },
    { 0.464, 0.245, 0.700 },
    { 0.023, 1.000, 0.000 },
    { 1.000, 0.350, 0.710 },
    { 0.000, 0.753, 1.000 },
    { 0.621, 0.700, 0.245 },
    { 0.804, 0.000, 1.000 },
    { 0.350, 1.000, 0.495 },
    { 1.000, 0.000, 0.125 },
    { 0.245, 0.333, 0.700 },
    { 0.411, 1.000, 0.000 },
    { 1.000, 0.350, 0.970 },
    { 0.000, 1.000, 0.558 },
    { 0.700, 0.432, 0.245 },
    { 0.309, 0.000, 1.000 },
    { 0.428, 1.000, 0.350 },
    { 1.000, 0.000, 0.692 },
    { 0.245, 0.689, 0.700 },
    { 1.000, 0.965, 0.000 },
    { 0.796, 0.350, 1.000 },
    { 0.000, 1.000, 0.116 },
    { 0.700, 0.245, 0.390 },
    { 0.000, 0.422, 1.000 },
    { 0.718, 1.000, 0.350 },
    { 0.950, 0.000, 1.000 },
    { 0.245, 0.700, 0.429 },
    { 1.000, 0.182, 0.000 },
    { 0.425, 0.350, 1.000 },
    { 0.228, 1.000, 0.000 },
    { 0.700, 0.245, 0.613 },
    { 0.000, 1.000, 0.817 },
    { 1.000, 0.833, 0.350 },
    { 0.547, 0.000, 1.000 },
    { 0.245, 0.700, 0.253 },
    { 1.000, 0.000, 0.488 },
    { 0.350, 0.776, 1.000 },
    { 0.744, 1.000, 0.000 },
    { 0.632, 0.245, 0.700 },
    { 0.000, 1.000, 0.273 },
    { 1.000, 0.350, 0.374 },
    { 0.000, 0.097, 1.000 },
    { 0.405, 0.700, 0.245 },
    { 1.000, 0.000, 0.913 },
    { 0.350, 1.000, 0.760 },
    { 1.000, 0.511, 0.000 },
    { 0.421, 0.245, 0.700 },
    { 0.078, 1.000, 0.000 },
    { 1.000, 0.350, 0.763 },
    { 0.000, 0.882, 1.000 },
    { 0.674, 0.700, 0.245 },
    { 0.738, 0.000, 1.000 },
    { 0.350, 1.000, 0.454 },
    { 1.000, 0.000, 0.238 },
    { 0.245, 0.392, 0.700 },
    { 0.497, 1.000, 0.000 },
    { 0.994, 0.350, 1.000 },
    { 0.000, 1.000, 0.467 },
    { 0.700, 0.372, 0.245 },
};

QList<QColor> generate_colors_ahb()
{
    QList<QColor> ret;
    for (int i = 0; i < 256; i++) {
        ret << QColor(colors_ahb[i][0] * 255, colors_ahb[i][1] * 255, colors_ahb[i][2] * 255);
    }
    //now we shift/cycle it over by one
    ret.insert(0, ret[ret.count() - 1]);
    ret = ret.mid(0, ret.count() - 1);
    return ret;
}

void set_nice_size(QWidget* W)
{
    int W0 = 1800, H0 = 1200;
    QRect geom = QApplication::desktop()->geometry();
    if (geom.width() - 100 < W0)
        W0 = geom.width() - 100;
    if (geom.height() - 100 < H0)
        H0 = geom.height() - 100;
    W->resize(W0, H0);
}
