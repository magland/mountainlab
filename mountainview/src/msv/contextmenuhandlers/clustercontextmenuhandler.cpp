#include "clustercontextmenuhandler.h"
#include "mvmainwindow.h"
#include <QAction>
#include <QJsonDocument>
#include <QMenu>
#include <QSet>
#include <QSignalMapper>
#include "cachemanager.h"
#include <QApplication>
#include <QProcess>
#include <mvcontext.h>

MVClusterContextMenuHandler::MVClusterContextMenuHandler(MVMainWindow* mw, QObject* parent)
    : QObject(parent)
    , MVAbstractContextMenuHandler(mw)
{
}

bool MVClusterContextMenuHandler::canHandle(const QMimeData& md) const
{
    return md.hasFormat("application/x-msv-clusters");
}

QList<QAction*> MVClusterContextMenuHandler::actions(const QMimeData& md)
{
    QSet<int> clusters;
    QDataStream ds(md.data("application/x-msv-clusters"));
    ds >> clusters;
    QList<QAction*> actions;

    MVMainWindow* mw = this->mainWindow();
    MVAbstractContext* context = mw->mvContext();

    MVContext* c = qobject_cast<MVContext*>(context);
    Q_ASSERT(c);

    /// Witold, is there a reason that the clusterList is being added as data to some of the actions below? If not, let's remove this.
    QVariantList clusterList;
    foreach (int cc, clusters)
        clusterList << cc;
    int first_cluster = clusters.values().first();

    //TAGS
    {
        actions << addTagMenu(clusters);
        actions << removeTagMenu(clusters);
        QAction* action = new QAction("Clear tags", 0);
        connect(action, &QAction::triggered, [clusters, c]() {
            foreach (int cluster_number, clusters) {
                c->setClusterTags(cluster_number, QSet<QString>());
            }
        });
        actions << action;
    }

    //Separator
    {
        QAction* action = new QAction(0);
        action->setSeparator(true);
        actions << action;
    }

    //MERGE
    {
        {
            QAction* A = new QAction(0);
            A->setEnabled(clusters.count() >= 2);
            A->setText("Merge selected");
            QObject::connect(A, &QAction::triggered, [clusters, c]() {
                QList<int> clusters_list=clusters.toList();
                for (int i1=0; i1<clusters_list.count(); i1++) {
                    for (int i2=0; i2<clusters_list.count(); i2++) {
                        if (i1!=i2) {
                            ClusterPair pair(clusters_list[i1],clusters_list[i2]);
                            QSet<QString> tags=c->clusterPairTags(pair);
                            tags.insert("merged");
                            c->setClusterPairTags(pair,tags);
                            {
                                //if we are merging then the user most likely wants to view as merged
                                c->setViewMerged(true);
                            }
                        }
                    }
                }
            });
            actions << A;
        }
        {
            QAction* A = new QAction(0);
            A->setEnabled(can_unmerge_selected_clusters(c, clusters));
            A->setText("Unmerge selected");
            QObject::connect(A, &QAction::triggered, [clusters, c]() {
                QList<ClusterPair> pairs=c->clusterPairAttributesKeys();
                for (int i=0; i<pairs.count(); i++) {
                    QSet<QString> tags=c->clusterPairTags(pairs[i]);
                    if (tags.contains("merged")) {
                        if ((clusters.contains(pairs[i].k1()))||(clusters.contains(pairs[i].k2()))) {
                            tags.remove("merged");
                        }
                    }
                    c->setClusterPairTags(pairs[i],tags);
                }
            });
            actions << A;
        }
    }

    //Separator
    {
        QAction* action = new QAction(0);
        action->setSeparator(true);
        actions << action;
    }

    //Views
    foreach (MVAbstractViewFactory* factory, this->mainWindow()->viewFactories()) {
        actions.append(factory->actions(md));
    }

    //Separator
    {
        QAction* action = new QAction(0);
        action->setSeparator(true);
        actions << action;
    }

    //CROSS-CORRELOGRAMS
    {
        {
            QAction* action = new QAction("Auto-correlograms", 0);
            action->setToolTip("Open auto-correlograms for these clusters");
            action->setEnabled(clusters.count() >= 1);
            action->setData(first_cluster);
            connect(action, &QAction::triggered, [mw]() {
                mw->openView("open-selected-auto-correlograms");
            });
            actions << action;
        }
        {
            QAction* action = new QAction("Cross-correlograms", 0);
            action->setToolTip("Open cross-correlograms associated with this cluster");
            action->setEnabled(clusters.count() == 1);
            action->setData(first_cluster);
            connect(action, &QAction::triggered, [mw]() {
                mw->openView("open-cross-correlograms");
            });
            actions << action;
        }
        {
            QAction* action = new QAction("Matrix of cross-correlograms", 0);
            action->setEnabled(clusters.count() >= 2);
            action->setData(clusterList);
            connect(action, &QAction::triggered, [mw]() {
                mw->openView("open-matrix-of-cross-correlograms");
            });
            actions << action;
        }
    }

    //Separator
    {
        QAction* action = new QAction(0);
        action->setSeparator(true);
        actions << action;
    }

    //DISCRIMINATION HISTOGRAMS
    {
        QAction* action = new QAction("Discrim Hist", 0);
        action->setToolTip("Open discrimination histograms for these clusters");
        action->setEnabled(clusters.count() >= 2);

        connect(action, &QAction::triggered, [mw]() {
            mw->openView("open-discrim-histograms");
        });
        actions << action;
    }

    //SPIKE SPRAY
    {
        QAction* action = new QAction("Spike spray", 0);
        action->setToolTip("Open spike spray for these clusters");
        action->setEnabled(clusters.count() >= 1);
        action->setData(clusterList);
        connect(action, &QAction::triggered, [mw]() {
            mw->openView("open-spike-spray");
        });
        actions << action;
    }

    //CLIPS
    {
        QAction* action = new QAction("Clips", 0);
        action->setToolTip("View all spike clips for these clusters");
        action->setEnabled(clusters.count() >= 1);
        action->setData(clusterList);
        connect(action, &QAction::triggered, [mw]() {
            mw->openView("open-clips");
        });
        actions << action;
    }

    //Separator
    {
        QAction* action = new QAction(0);
        action->setSeparator(true);
        actions << action;
    }

    //EXTRACT SELECTED CLUSTERS
    {
        QAction* action = new QAction("Extract selected clusters in new window", 0);
        action->setToolTip("Extract selected clusters in new window");
        action->setEnabled(clusters.count() >= 1);
        connect(action, SIGNAL(triggered(bool)), this, SLOT(slot_extract_selected_clusters()));
        //connect(action, &QAction::triggered, [mw]() {
        //    mw->extractSelectedClusters();
        //});
        actions << action;
    }

    return actions;
}

void MVClusterContextMenuHandler::slot_extract_selected_clusters()
{
    MVContext* c = qobject_cast<MVContext*>(mainWindow()->mvContext());
    Q_ASSERT(c);

    QString tmp_fname = CacheManager::globalInstance()->makeLocalFile() + ".mv2";
    //QJsonObject obj = this->mainWindow()->mvContext()->toMVFileObject();
    QJsonObject obj = this->mainWindow()->mvContext()->toMV2FileObject();
    QString json = QJsonDocument(obj).toJson();
    TextFile::write(tmp_fname, json);
    CacheManager::globalInstance()->setTemporaryFileDuration(tmp_fname, 600);
    QString exe = qApp->applicationFilePath();
    QStringList args;
    QList<int> clusters = c->selectedClusters();
    QStringList clusters_str;
    foreach (int cluster, clusters) {
        clusters_str << QString("%1").arg(cluster);
    }

    args << tmp_fname << "--clusters=" + clusters_str.join(",");
    qDebug().noquote() << "EXECUTING: " + exe + " " + args.join(" ");
    QProcess::startDetached(exe, args);
}

bool MVClusterContextMenuHandler::can_unmerge_selected_clusters(MVContext* context, const QSet<int>& clusters)
{
    ClusterMerge CM = context->clusterMerge();
    foreach (int k, clusters) {
        if (CM.representativeLabel(k) != k)
            return true;
        if (CM.getMergeGroup(k).count() > 1)
            return true;
    }
    return false;
}

QAction* MVClusterContextMenuHandler::addTagMenu(const QSet<int>& clusters) const
{
    MVContext* c = qobject_cast<MVContext*>(mainWindow()->mvContext());
    Q_ASSERT(c);

    QMenu* M = new QMenu;
    M->setTitle("Add tag");
    QSignalMapper* mapper = new QSignalMapper(M);
    foreach (const QString& tag, validTags()) {
        QAction* a = new QAction(tag, M);
        a->setData(tag);
        M->addAction(a);
        mapper->setMapping(a, tag);
        connect(a, SIGNAL(triggered()), mapper, SLOT(map()));
    }
    /// Witold, I am a bit nervous about passing the context pointer into the lambda function below
    connect(mapper, static_cast<void (QSignalMapper::*)(const QString&)>(&QSignalMapper::mapped),
        [clusters, c](const QString& tag) {
        foreach(int cluster, clusters) {
            QSet<QString> clTags = c->clusterTags(cluster);
            clTags.insert(tag);
            c->setClusterTags(cluster, clTags);
        }
        });

    return M->menuAction();
}

QAction* MVClusterContextMenuHandler::removeTagMenu(const QSet<int>& clusters) const
{
    MVContext* c = qobject_cast<MVContext*>(mainWindow()->mvContext());
    Q_ASSERT(c);

    QSet<QString> tags_set;
    foreach (int cluster_number, clusters) {
        QJsonObject attributes = c->clusterAttributes(cluster_number);
        QJsonArray tags = attributes["tags"].toArray();
        for (int i = 0; i < tags.count(); i++) {
            tags_set.insert(tags[i].toString());
        }
    }
    QStringList tags_list = tags_set.toList();
    qSort(tags_list);
    QMenu* M = new QMenu;
    M->setTitle("Remove tag");
    QSignalMapper* mapper = new QSignalMapper(M);
    foreach (const QString& tag, tags_list) {
        QAction* a = new QAction(tag, M);
        a->setData(tag);
        M->addAction(a);
        mapper->setMapping(a, tag);
        connect(a, SIGNAL(triggered()), mapper, SLOT(map()));
    }
    connect(mapper, static_cast<void (QSignalMapper::*)(const QString&)>(&QSignalMapper::mapped),
        [clusters, c](const QString& tag) {
        foreach(int cluster, clusters) {
            QSet<QString> clTags = c->clusterTags(cluster);
            clTags.remove(tag);
            c->setClusterTags(cluster, clTags);
        }
        });
    M->setEnabled(!tags_list.isEmpty());
    return M->menuAction();
}

QStringList MVClusterContextMenuHandler::validTags() const
{
    MVContext* c = qobject_cast<MVContext*>(mainWindow()->mvContext());
    Q_ASSERT(c);

    /// TODO (LOW) these go in a configuration file
    QSet<QString> set = c->allClusterTags();
    set << "accepted"
        << "rejected"
        << "noise"
        << "mua"
        << "artifact";
    QStringList ret = set.toList();
    qSort(ret);
    return ret;
}
