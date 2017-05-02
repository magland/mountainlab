#include <QApplication>
#include <QDebug>
#include <qdatetime.h>
#include <QFile>
#include <QFileInfo>
#include <QMessageBox>
#include <QProcess>
#include <QStringList>
#include "cachemanager.h"
#include "resolveprvsdialog.h"
#include "prvgui.h"

#include "usagetracking.h"
#include "mda.h"
#include <QDesktopServices>
#include <QDesktopWidget>
#include <QDir>
#include <QImageWriter>
#include "mvmainwindow.h"
#include "closemehandler.h"
#include "taskprogress.h"

#include <QHBoxLayout>
#include <QJsonDocument>
#include <QRunnable>
#include <QThreadPool>
#include <QtConcurrentRun>
#include <tabber.h>
#include <QFileDialog>
#include <QSettings>
#include <objectregistry.h>
#include <icounter.h>
#include <sslvopenviewscontrol.h>
#include <sslvprefscontrol.h>
#include <sslvtoolscontrol.h>
#include <firingeventplugin.h>
//#include <sslvclustermetricsview.h>
#include <sslvclustermetricsplugin.h>
#include <sslvtimespanplugin.h>
#include <clusterdetailplugin.h>
#include "correlogramplugin.h"

void setup_main_window(MVMainWindow* W);

QColor brighten(QColor col, int amount);
QList<QColor> generate_colors_ahb();
QList<QColor> generate_colors_old(const QColor& bg, const QColor& fg, int noColors);

#include "multiscaletimeseries.h"
//#include "mvfile.h"
#include "taskprogressview.h"
#include "mvcontrolpanel2.h"
#include "mvabstractcontrol.h"
#include "sslvcontext.h"

void set_nice_size(QWidget* W);
bool check_whether_prv_objects_need_to_be_downloaded_or_regenerated(QJsonObject obj);
bool check_whether_prv_objects_need_to_be_downloaded_or_regenerated(QList<PrvRecord> prvs);

void try_to_automatically_download_and_regenerate_prv_objects(QJsonObject obj);
void try_to_automatically_download_and_regenerate_prv_objects(QList<PrvRecord> prvs);

int main(int argc, char* argv[])
{
    QApplication a(argc, argv);
    ObjectRegistry registry;
    CounterManager* counterManager = new CounterManager;
    registry.addAutoReleasedObject(counterManager);

    ObjectRegistry::addAutoReleasedObject(new IIntCounter("allocated_bytes"));
    ObjectRegistry::addAutoReleasedObject(new IIntCounter("freed_bytes"));
    ObjectRegistry::addAutoReleasedObject(new IIntCounter("remote_processing_time"));
    ObjectRegistry::addAutoReleasedObject(new IIntCounter("bytes_downloaded"));
    ObjectRegistry::addAutoReleasedObject(new IIntCounter("bytes_read"));
    ObjectRegistry::addAutoReleasedObject(new IIntCounter("bytes_written"));

    QList<ICounterBase*> counters = ObjectRegistry::getObjects<ICounterBase>();
    counterManager->setCounters(counters);
    counterManager->connect(ObjectRegistry::instance(), &ObjectRegistry::objectAdded, [counterManager](QObject* o) {
      if (ICounterBase *cntr = qobject_cast<ICounterBase*>(o)) {
          counterManager->addCounter(cntr);
      }
    });

    // make sure task progress monitor is instantiated in the main thread
    TaskManager::TaskProgressMonitor* monitor = TaskManager::TaskProgressMonitor::globalInstance();
    registry.addObject(monitor);
    CloseMeHandler::start();

    setbuf(stdout, 0);

    CLParams CLP(argc, argv);
    {
        TaskProgressView TPV;
        TPV.show();

        // add on a .prv if the [file] does not exist, but the [file].prv does
        // TODO this is very sloppy and dangerous... should not be done here. Should rather be done in kron-view
        QStringList keys = CLP.named_parameters.keys();
        foreach (QString key, keys) {
            QVariant val = CLP.named_parameters[key].toString();
            if ((!QFile::exists(val.toString())) && (QFile::exists(val.toString() + ".prv"))) {
                val = val.toString() + ".prv";
                CLP.named_parameters[key] = val;
            }
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

    QString mv2_fname;
    if (CLP.unnamed_parameters.value(0).endsWith(".mv2")) {
        mv2_fname = CLP.unnamed_parameters.value(0);
    }

    QString mode = CLP.named_parameters.value("mode", "default").toString();
    if (mode == "default") {

        SSLVContext* context = new SSLVContext; //note that the view agent does not get deleted. :(
        context->setChannelColors(channel_colors);
        context->setClusterColors(label_colors);
        MVMainWindow* W = new MVMainWindow(context);

        if (mv2_fname.isEmpty()) {
            SSLVContext dc; //dummy context
            if (CLP.named_parameters.contains("samplerate")) {
                dc.setSampleRate(CLP.named_parameters.value("samplerate", 0).toDouble());
            }
            if (CLP.named_parameters.contains("firings")) {
                QString firings_path = CLP.named_parameters["firings"].toString();
                dc.setFirings(DiskReadMda(firings_path));
                W->setWindowTitle(firings_path);
            }
            if (CLP.named_parameters.contains("raw")) {
                QStringList raw_path = MLUtil::toStringList(CLP.named_parameters["raw"]);
                DiskReadMda32 ts(2, raw_path);
                dc.setTimeseries(ts);
            }
            if (CLP.named_parameters.contains("firings")) {
                QString firings_path = CLP.named_parameters["firings"].toString();
                DiskReadMda firings0(firings_path);
                dc.setFirings(firings0);
            }
            if (CLP.named_parameters.contains("window_title")) {
                QString window_title = CLP.named_parameters["window_title"].toString();
                W->setWindowTitle(window_title);
            }
            if (CLP.named_parameters.contains("cluster_metrics")) {
                QString cluster_metrics_path = CLP.named_parameters["cluster_metrics"].toString();
                QJsonObject obj = QJsonDocument::fromJson(TextFile::read(cluster_metrics_path).toUtf8()).object();
                dc.setClusterMetrics(obj);
            }

            QJsonObject mv2 = dc.toMV2FileObject();
            mv2_fname = CacheManager::globalInstance()->makeLocalFile() + ".mv2";
            QString mv2_text = QJsonDocument(mv2).toJson();
            TextFile::write(mv2_fname, mv2_text);
            CacheManager::globalInstance()->setTemporaryFileDuration(mv2_fname, 600);
        }

        if (!mv2_fname.isEmpty()) {
            TaskProgressView TPV;
            TPV.show();
            bool done_checking = false;
            QJsonObject obj;
            while (!done_checking) {
                QString json = TextFile::read(mv2_fname);
                obj = QJsonDocument::fromJson(json.toLatin1()).object();
                if (check_whether_prv_objects_need_to_be_downloaded_or_regenerated(obj)) {
                    ResolvePrvsDialog dlg;
                    if (dlg.exec() == QDialog::Accepted) {
                        if (dlg.choice() == ResolvePrvsDialog::OpenPrvGui) {
                            int exit_code = system(("prv-gui " + mv2_fname).toUtf8().data());
                            Q_UNUSED(exit_code)
                            done_checking = false; //check again
                        }
                        else if (dlg.choice() == ResolvePrvsDialog::AutomaticallyDownloadAndRegenerate) {
                            try_to_automatically_download_and_regenerate_prv_objects(obj);
                            done_checking = false; //check again
                        }
                        else {
                            done_checking = true;
                        }
                    }
                    else {
                        return -1;
                    }
                }
                else
                    done_checking = true;
            }
            context->setFromMV2FileObject(obj);
            context->setMV2FileName(mv2_fname);
        }

        set_nice_size(W);

        W->show();

        setup_main_window(W);
        W->addControl(new OpenViewsControl(context, W), true);
        W->addControl(new SSLVPrefsControl(context, W), true);
        W->addControl(new SSLVToolsControl(context, W), true);

        a.processEvents();

        W->setCurrentContainerName("north");
        W->openView("open-time-span");
        W->setCurrentContainerName("south");
        W->openView("open-cluster-metrics");

        return a.exec();
    }

    int ret = a.exec();

    printf("Number of files open: %d, number of unfreed mallocs: %d, number of unfreed megabytes: %g\n", jnumfilesopen(), jmalloccount(), (int)jbytesallocated() * 1.0 / 1000000);

    return ret;
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

// generate_colors() is adapted from code by...
/*
 * Copyright (c) 2008 Helder Correia
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
*/

QList<QColor> generate_colors_old(const QColor& bg, const QColor& fg, int noColors)
{
    QList<QColor> colors;
    const int HUE_BASE = (bg.hue() == -1) ? 90 : bg.hue();
    int h, s, v;

    for (int i = 0; i < noColors; i++) {
        h = int(HUE_BASE + (360.0 / noColors * i)) % 360;
        s = 240;
        v = int(qMax(bg.value(), fg.value()) * 0.85);

        // take care of corner cases
        const int M = 35;
        if ((h < bg.hue() + M && h > bg.hue() - M)
            || (h < fg.hue() + M && h > fg.hue() - M)) {
            h = ((bg.hue() + fg.hue()) / (i + 1)) % 360;
            s = ((bg.saturation() + fg.saturation() + 2 * i) / 2) % 256;
            v = ((bg.value() + fg.value() + 2 * i) / 2) % 256;
        }

        colors.append(QColor::fromHsv(h, s, v));
    }

    return colors;
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

void setup_main_window(MVMainWindow* W)
{
    W->loadPlugin(new SSLVClusterMetricsPlugin);
    W->loadPlugin(new SSLVTimeSpanPlugin);
    W->loadPlugin(new ClusterDetailPlugin);
    W->loadPlugin(new CorrelogramPlugin);
    W->loadPlugin(new FiringEventPlugin);
    //W->loadPlugin(new TimeseriesPlugin);
    //W->loadPlugin(new SpectrogramPlugin);

    //W->loadPlugin(new ClusterDetailPlugin);
    //W->loadPlugin(new IsolationMatrixPlugin);
    //W->loadPlugin(new ClipsViewPlugin);
    //W->loadPlugin(new ClusterContextMenuPlugin);
    //W->loadPlugin(new CurationProgramPlugin);
    /*
    W->registerViewFactory(new MVTemplatesView3Factory(W));
    W->registerViewFactory(new MVAutoCorrelogramsFactory(W));
    W->registerViewFactory(new MVSelectedAutoCorrelogramsFactory(W));
    W->registerViewFactory(new MVCrossCorrelogramsFactory(W));
    W->registerViewFactory(new MVMatrixOfCrossCorrelogramsFactory(W));
    W->registerViewFactory(new MVSelectedCrossCorrelogramsFactory(W));
    W->registerViewFactory(new MVTimeSeriesDataFactory(W));
    W->registerViewFactory(new MVPCAFeaturesFactory(W));
    W->registerViewFactory(new MVChannelFeaturesFactory(W));
    W->registerViewFactory(new MVSpikeSprayFactory(W));
    W->registerViewFactory(new MVFiringEventsFactory(W));
    W->registerViewFactory(new MVAmplitudeHistograms3Factory(W));
    W->registerViewFactory(new MVDiscrimHistFactory(W));
    W->registerViewFactory(new MVDiscrimHistGuideFactory(W));
    W->registerViewFactory(new MVFireTrackFactory(W));
    */
}

bool check_whether_prv_objects_need_to_be_downloaded_or_regenerated(QJsonObject obj)
{
    QList<PrvRecord> prvs = find_prvs("", obj);
    return check_whether_prv_objects_need_to_be_downloaded_or_regenerated(prvs);
}

QString check_if_on_local_disk(PrvRecord prv)
{
    QString cmd = "prv";
    QStringList args;
    args << "locate";
    args << "--checksum=" + prv.checksum;
    args << "--fcs=" + prv.fcs;
    args << QString("--size=%1").arg(prv.size);
    args << QString("--original_path=%1").arg(prv.original_path);
    args << "--local-only";
    QString output = exec_process_and_return_output(cmd, args);
    return output.split("\n").last();
}

bool check_whether_prv_objects_need_to_be_downloaded_or_regenerated(QList<PrvRecord> prvs)
{
    foreach (PrvRecord prv, prvs) {
        QString path = check_if_on_local_disk(prv);
        if (path.isEmpty()) {
            return true;
        }
    }
    return false;
}

void try_to_automatically_download_and_regenerate_prv_objects(QJsonObject obj)
{
    QList<PrvRecord> prvs = find_prvs("", obj);
    return try_to_automatically_download_and_regenerate_prv_objects(prvs);
}

void system_call_keeping_gui_alive(QString cmd, QString task_label)
{
    if (task_label.isEmpty())
        task_label = "Running " + cmd;
    TaskProgress task(task_label);
    task.log() << cmd;
    QProcess P;
    P.setReadChannelMode(QProcess::MergedChannels);
    P.start(cmd);
    while (!P.waitForFinished(10)) {
        qApp->processEvents();
    }
    task.log() << P.readAll();
}

void try_to_automatically_download_and_regenerate_prv_objects(QList<PrvRecord> prvs)
{
    foreach (PrvRecord prv, prvs) {
        if (check_if_on_local_disk(prv).isEmpty()) {
            QString src_fname = CacheManager::globalInstance()->makeLocalFile() + ".tmp.prv";
            QString dst_fname = src_fname + ".recover";
            TextFile::write(src_fname, QJsonDocument(prv.original_object).toJson());
            QString cmd = QString("prv recover %1 %2").arg(src_fname).arg(dst_fname);
            system_call_keeping_gui_alive(cmd.toLatin1().data(), "Recovering " + prv.label);
            QFile::remove(src_fname);
        }
    }
}
