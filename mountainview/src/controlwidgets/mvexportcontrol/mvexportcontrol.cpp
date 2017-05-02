/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
** Created: 6/27/2016
*******************************************************/

#include "flowlayout.h"
#include "mvexportcontrol.h"
#include "mlcommon.h"

#include <QGridLayout>
#include <QLabel>
#include <QLineEdit>
#include <QMap>
#include <QPushButton>
#include <QTimer>
#include <QFileDialog>
#include <QSettings>
#include <QJsonDocument>
#include <QMessageBox>
#include <QTextBrowser>
#include <QProcess>
#include "taskprogress.h"
#include <QThread>
#include <cachemanager.h>
#include "diskreadmda.h"
#include "mvcontext.h"

class MVExportControlPrivate {
public:
    MVExportControl* q;
};

MVExportControl::MVExportControl(MVAbstractContext* context, MVMainWindow* mw)
    : MVAbstractControl(context, mw)
{
    d = new MVExportControlPrivate;
    d->q = this;

    QFont fnt = this->font();
    fnt.setPixelSize(qMax(14, fnt.pixelSize() - 6));
    this->setFont(fnt);

    FlowLayout* flayout = new FlowLayout;
    this->setLayout(flayout);
    {
        QPushButton* B = new QPushButton("Export .mv2 document");
        connect(B, SIGNAL(clicked(bool)), this, SLOT(slot_export_mv2_document()));
        flayout->addWidget(B);
    }
    {
        QPushButton* B = new QPushButton("Export firings");
        connect(B, SIGNAL(clicked(bool)), this, SLOT(slot_export_firings()));
        flayout->addWidget(B);
    }
    {
        QPushButton* B = new QPushButton("Export curated firings");
        connect(B, SIGNAL(clicked(bool)), this, SLOT(slot_export_curated_firings()));
        flayout->addWidget(B);
    }
    {
        QPushButton* B = new QPushButton("Open PRV manager");
        connect(B, SIGNAL(clicked(bool)), this, SLOT(slot_open_prv_manager()));
        flayout->addWidget(B);
    }

    /*
    {
        QPushButton* B = new QPushButton("Export .mv document");
        connect(B, SIGNAL(clicked(bool)), this, SLOT(slot_export_mv_document()));
        flayout->addWidget(B);
    }
    {
        QPushButton* B = new QPushButton("Export static views (.smv)");
        connect(B, SIGNAL(clicked(bool)), this, SLOT(slot_export_static_views()));
        flayout->addWidget(B);
    }
    {
        QPushButton* B = new QPushButton("Share views on web");
        connect(B, SIGNAL(clicked(bool)), this, SLOT(slot_share_views_on_web()));
        flayout->addWidget(B);
    }
    {
        QPushButton* B = new QPushButton("Export firings.mda");
        connect(B, SIGNAL(clicked(bool)), this, SLOT(slot_export_firings_file()));
        flayout->addWidget(B);
    }

    {
        QPushButton* B = new QPushButton("Export cluster_curation.mda");
        connect(B, SIGNAL(clicked(bool)), this, SLOT(slot_export_cluster_curation_file()));
        flayout->addWidget(B);
    }
    */

    updateControls();
}

MVExportControl::~MVExportControl()
{
    delete d;
}

QString MVExportControl::title() const
{
    return "Export / Upload / Download";
}

void MVExportControl::updateContext()
{
}

void MVExportControl::updateControls()
{
}

class UploadStuffForMV2 : public QThread {
public:
    QStringList all_prv_paths;
    ///TODO: how can we guarantee these guys won't get deleted during this potentially rather lengthy process.
    /// They should really be in memory from the start!
    bool do_ensure_local;
    bool do_ensure_remote;
    bool raw_only;
    QString server;

    void run()
    {
        TaskProgress task(".mv2 -- saving and/or uploading files");
        for (int i = 0; i < all_prv_paths.count(); i++) {
            QString path = all_prv_paths[i];
            task.log() << QString("i = %1 / %2").arg(i + 1).arg(all_prv_paths.count());
            task.setProgress((i * 2 + 0.5) / (all_prv_paths.count() * 2));
            if (do_ensure_local) {
                task.log() << "Ensuring local: " + path;
                if (!ensure_local(path)) {
                    task.error() << "Error ensuring local: " + path;
                    return;
                }
            }
            task.setProgress((i * 2 + 1 + 0.5) / (all_prv_paths.count() * 2));
            if (do_ensure_remote) {
                task.log() << QString("Ensuring remote (%1): %2").arg(server).arg(path);
                if (!ensure_remote(path, server)) {
                    task.error() << QString("Error ensuring remote (%1): %2").arg(server).arg(path);
                    return;
                }
            }
        }
    }

private:
    bool ensure_local(QString prv_path)
    {
        QString cmd = "prv";
        QStringList args;
        args << "ensure-local" << prv_path;
        if (raw_only)
            args << "--raw-only";
        return (QProcess::execute(cmd, args) == 0);
    }

    bool ensure_remote(QString prv_path, QString server)
    {
        QString cmd = "prv";
        QStringList args;
        args << "ensure-remote" << prv_path << "--server=" + server;
        if (raw_only)
            args << "--raw-only";
        return (QProcess::execute(cmd, args) == 0);
    }
};

void MVExportControl::slot_export_mv2_document()
{
    //MVContext* c = qobject_cast<MVContext*>(mvContext());
    //Q_ASSERT(c);

    //QSettings settings("SCDA", "MountainView");
    //QString default_dir = settings.value("default_export_dir", "").toString();
    QString default_dir = QDir::currentPath();
    QString fname = QFileDialog::getSaveFileName(this, "Export mv2 document", default_dir, "*.mv2");
    if (fname.isEmpty())
        return;
    //settings.setValue("default_export_dir", QFileInfo(fname).path());
    if (QFileInfo(fname).suffix() != "mv2")
        fname = fname + ".mv2";
    TaskProgress task("export mv2 document");
    task.log() << "Writing: " + fname;
    QJsonObject obj = mvContext()->toMV2FileObject();
    QString json = QJsonDocument(obj).toJson();
    if (TextFile::write(fname, json)) {
        mvContext()->setMV2FileName(fname);
        task.log() << QString("Wrote %1 kilobytes").arg(json.count() * 1.0 / 1000);
    }
    else {
        task.error("Error writing .mv2 file: " + fname);
    }
}

/*
void MVExportControl::slot_export_mv_document()
{
    //QSettings settings("SCDA", "MountainView");
    //QString default_dir = settings.value("default_export_dir", "").toString();
    QString default_dir = QDir::currentPath();
    QString fname = QFileDialog::getSaveFileName(this, "Export mountainview document", default_dir, "*.mv");
    if (fname.isEmpty())
        return;
    //settings.setValue("default_export_dir", QFileInfo(fname).path());
    if (QFileInfo(fname).suffix() != "mv")
        fname = fname + ".mv";
    QJsonObject obj = this->c->toMVFileObject();
    QString json = QJsonDocument(obj).toJson();
    if (!TextFile::write(fname, json)) {
        TaskProgress task("export mountainview document");
        task.error("Error writing .mv file: " + fname);
    }
    //MVFile ff = mainWindow()->getMVFile();
    //if (!ff.write(fname)) {
    //    TaskProgress task("export mountainview document");
    //    task.error("Error writing .mv file: " + fname);
    //}
}
*/

#include "computationthread.h"
/// TODO fix the DownloadComputer2 and don't require computationthread.h
class DownloadComputer2 : public ComputationThread {
public:
    //inputs
    QString source_path;
    QString dest_path;
    bool use_float64;

    void compute();
};
void DownloadComputer2::compute()
{
    TaskProgress task("Downlading");
    task.setDescription(QString("Downloading %1 to %2").arg(source_path).arg(dest_path));
    DiskReadMda X(source_path);
    Mda Y;
    task.setProgress(0.2);
    task.log(QString("Reading/Downloading %1x%2x%3").arg(X.N1()).arg(X.N2()).arg(X.N3()));
    if (!X.readChunk(Y, 0, 0, 0, X.N1(), X.N2(), X.N3())) {
        if (MLUtil::threadInterruptRequested()) {
            task.error("Halted download: " + source_path);
        }
        else {
            task.error("Failed to readChunk from: " + source_path);
        }
        return;
    }
    task.setProgress(0.8);
    if (use_float64) {
        task.log("Writing 64-bit to " + dest_path);
        Y.write64(dest_path);
    }
    else {
        task.log("Writing 32-bit to " + dest_path);
        Y.write32(dest_path);
    }
}

void export_file(QString source_path, QString dest_path, bool use_float64)
{
    DownloadComputer2* C = new DownloadComputer2;
    C->source_path = source_path;
    C->dest_path = dest_path;
    C->use_float64 = use_float64;
    C->setDeleteOnComplete(true);
    C->startComputation();
}

void MVExportControl::slot_export_firings()
{
    //MVContext* c = qobject_cast<MVContext*>(mvContext());
    //Q_ASSERT(c);

    //QSettings settings("SCDA", "MountainView");
    //QString default_dir = settings.value("default_export_dir", "").toString();
    QString default_dir = QDir::currentPath();
    QString fname = QFileDialog::getSaveFileName(this, "Export original firings", default_dir, "*.mda");
    if (fname.isEmpty())
        return;
    //settings.setValue("default_export_dir", QFileInfo(fname).path());
    if (QFileInfo(fname).suffix() != "mda")
        fname = fname + ".mda";

    QJsonObject obj = mvContext()->toMV2FileObject();
    DiskReadMda firings(obj["firings"].toObject());
    //DiskReadMda firings = c->firings();
    export_file(firings.makePath(), fname, true);
}

QString get_local_path_of_firings_file_or_current_path(const DiskReadMda& X)
{
    Q_UNUSED(X)
    /*QString path = X.makePath();
    if (!path.isEmpty()) {
        if (!path.startsWith("http:")) {
            if (QFile::exists(path)) {
                return QFileInfo(path).path();
            }
        }
    }*/
    return QDir::currentPath();
}

void MVExportControl::slot_export_curated_firings()
{
    MVContext* c = qobject_cast<MVContext*>(mvContext());
    Q_ASSERT(c);

    QString default_dir = get_local_path_of_firings_file_or_current_path(c->firings());
    QString fname = QFileDialog::getSaveFileName(this, "Export curated firings file", default_dir + "/firings.curated.mda", "*.mda");
    if (fname.isEmpty())
        return;
    if (QFileInfo(fname).suffix() != "mda")
        fname = fname + ".mda";

    DiskReadMda F = c->firings();

    QVector<int> labels;
    for (long i = 0; i < F.N2(); i++) {
        int label = F.value(2, i);
        labels << label;
    }

    int K = MLCompute::max(labels);

    QSet<int> accepted_clusters;
    QMap<int, int> merge_map;
    for (int k = 1; k <= K; k++) {
        if (!c->clusterTags(k).contains("rejected")) {
            accepted_clusters.insert(k);
            merge_map[k] = c->clusterMerge().representativeLabel(k);
        }
    }

    QVector<long> inds_to_use;
    for (long i = 0; i < F.N2(); i++) {
        int label = F.value(2, i);
        if (accepted_clusters.contains(label)) {
            inds_to_use << i;
        }
    }

    Mda F2(F.N1(), inds_to_use.count());
    for (long j = 0; j < inds_to_use.count(); j++) {
        long i = inds_to_use[j];
        for (int a = 0; a < F.N1(); a++) {
            F2.set(F.value(a, i), a, j);
        }
        int label = F.value(2, i);
        label = merge_map[label];
        F2.setValue(label, 2, j);
    }

    if (!F2.write64(fname)) {
        QMessageBox::warning(0, "Problem exporting firings curated file", "Unable to write file: " + fname);
    }
}

void MVExportControl::slot_open_prv_manager()
{
    MVContext* c = qobject_cast<MVContext*>(mvContext());
    Q_ASSERT(c);

    QString fname = c->mv2FileName();
    if (fname.isEmpty()) {
        fname = CacheManager::globalInstance()->makeLocalFile() + ".prv";
        QJsonObject obj = c->toMV2FileObject();
        QString json = QJsonDocument(obj).toJson();
        if (!TextFile::write(fname, json)) {
            TaskProgress errtask;
            errtask.error("Error writing .mv2 file: " + fname);
            return;
        }
        CacheManager::globalInstance()->setTemporaryFileDuration(fname, 600);
    }
    int exit_code = system(("prv-gui " + fname).toUtf8().data());
    Q_UNUSED(exit_code)
}

/*
void MVExportControl::slot_export_cluster_curation_file()
{
    //first row is the cluster number
    //second row is 0 if not accepted, 1 if accepted
    //third row is the merge label (i.e., smallest cluster number in the merge group)
    QSet<int> clusters_set = c->clustersSubset();
    if (clusters_set.isEmpty()) {
        QList<int> keys = c->clusterAttributesKeys();
        foreach (int key, keys) {
            clusters_set.insert(key);
        }
        QList<ClusterPair> pairs = c->clusterPairAttributesKeys();
        foreach (ClusterPair pair, pairs) {
            clusters_set.insert(pair.k1());
            clusters_set.insert(pair.k2());
        }
    }
    QList<int> clusters = clusters_set.toList();
    qSort(clusters);
    int num = clusters.count();
    Mda cluster_curation(3, num);
    for (int i = 0; i < num; i++) {
        int accepted = 0;
        if (c->clusterTags(clusters[i]).contains("accepted"))
            accepted = 1;
        int merge_label = c->clusterMerge().representativeLabel(clusters[i]);
        cluster_curation.setValue(clusters[i], 0, i);
        cluster_curation.setValue(accepted, 1, i);
        cluster_curation.setValue(merge_label, 2, i);
    }

    //QSettings settings("SCDA", "MountainView");
    //QString default_dir = settings.value("default_export_dir", "").toString();
    QString default_dir = QDir::currentPath();
    QString fname = QFileDialog::getSaveFileName(this, "Export cluster curation array", default_dir, "*.mda");
    if (fname.isEmpty())
        return;
    //settings.setValue("default_export_dir", QFileInfo(fname).path());
    if (QFileInfo(fname).suffix() != "mda")
        fname = fname + ".mda";

    if (!cluster_curation.write32(fname)) {
        QMessageBox::warning(0, "Problem exporting cluster curation array", "Unable to write file: " + fname);
    }
}
*/

/*
void MVExportControl::slot_export_static_views()
{
    //QSettings settings("SCDA", "MountainView");
    //QString default_dir = settings.value("default_export_dir", "").toString();
    QString default_dir = QDir::currentPath();
    QString fname = QFileDialog::getSaveFileName(this, "Export static views", default_dir, "*.smv");
    if (fname.isEmpty())
        return;
    //settings.setValue("default_export_dir", QFileInfo(fname).path());
    if (QFileInfo(fname).suffix() != "smv")
        fname = fname + ".smv";
    QJsonObject obj = this->mainWindow()->exportStaticViews();
    if (!TextFile::write(fname, QJsonDocument(obj).toJson())) {
        qWarning() << "Unable to write file: " + fname;
    }
}

void MVExportControl::slot_share_views_on_web()
{
    QTextBrowser* tb = new QTextBrowser;
    tb->setOpenExternalLinks(true);
    tb->setAttribute(Qt::WA_DeleteOnClose);
    QString html;
    html += "<p>";
    html += "<p>First export a .smv file. Then upload it on this website:</p>";
    QString url0 = "http://datalaunch.org/upload";
    html += "<p><a href=\"" + url0 + "\">" + url0 + "</a></p>";
    html += "<p>Note that for now, only templates and cross-correlograms views may be shared.</p>";
    tb->show();
    tb->resize(600, 400);
    tb->setHtml(html);
}
*/
