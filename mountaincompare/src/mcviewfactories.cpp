#include "mccontext.h"
#include "mcviewfactories.h"

#include "clusterdetailview.h"
#include <views/compareclusterview.h>
#include <mvtimeseriesview2.h>
#include <QMessageBox>
#include <mvspikesprayview.h>
#include <mvclusterwidget.h>

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
ClusterDetail1Factory::ClusterDetail1Factory(MVMainWindow* mw, QObject* parent)
    : MVAbstractViewFactory(mw, parent)
{
}

QString ClusterDetail1Factory::id() const
{
    return QStringLiteral("open-cluster-details-1");
}

QString ClusterDetail1Factory::name() const
{
    return tr("Cluster Details 1");
}

QString ClusterDetail1Factory::title() const
{
    return tr("Details 1");
}

MVAbstractViewFactory::PreferredOpenLocation ClusterDetail1Factory::preferredOpenLocation() const
{
    return PreferredOpenLocation::North;
}

MVAbstractView* ClusterDetail1Factory::createView(MVAbstractContext* context)
{
    MCContext* mc_context = qobject_cast<MCContext*>(context);
    if (!mc_context)
        return 0;

    ClusterDetailView* X = new ClusterDetailView(mc_context->mvContext1());
    return X;
}

ClusterDetail2Factory::ClusterDetail2Factory(MVMainWindow* mw, QObject* parent)
    : MVAbstractViewFactory(mw, parent)
{
}

QString ClusterDetail2Factory::id() const
{
    return QStringLiteral("open-cluster-details-2");
}

QString ClusterDetail2Factory::name() const
{
    return tr("Cluster Details 2");
}

QString ClusterDetail2Factory::title() const
{
    return tr("Details 2");
}

MVAbstractViewFactory::PreferredOpenLocation ClusterDetail2Factory::preferredOpenLocation() const
{
    return PreferredOpenLocation::South;
}

MVAbstractView* ClusterDetail2Factory::createView(MVAbstractContext* context)
{
    MCContext* mc_context = qobject_cast<MCContext*>(context);
    if (!mc_context)
        return 0;

    ClusterDetailView* X = new ClusterDetailView(mc_context->mvContext2());
    return X;
}

////////////////////////////////////////////////////////////////////////
MVTimeSeriesView1Factory::MVTimeSeriesView1Factory(MVMainWindow* mw, QObject* parent)
    : MVAbstractViewFactory(mw, parent)
{
}

QString MVTimeSeriesView1Factory::id() const
{
    return QStringLiteral("open-timeseries-data-1");
}

QString MVTimeSeriesView1Factory::name() const
{
    return tr("Timeseries Data 1");
}

QString MVTimeSeriesView1Factory::title() const
{
    return tr("Timeseries 1");
}

MVAbstractViewFactory::PreferredOpenLocation MVTimeSeriesView1Factory::preferredOpenLocation() const
{
    return PreferredOpenLocation::North;
}

MVAbstractView* MVTimeSeriesView1Factory::createView(MVAbstractContext* context)
{
    MCContext* c = qobject_cast<MCContext*>(context);
    Q_ASSERT(c);

    MVTimeSeriesView2* X = new MVTimeSeriesView2(c->mvContext1());
    QList<int> ks = c->mvContext1()->selectedClusters();
    X->setLabelsToView(ks.toSet());
    return X;
}

MVTimeSeriesView2Factory::MVTimeSeriesView2Factory(MVMainWindow* mw, QObject* parent)
    : MVAbstractViewFactory(mw, parent)
{
}

QString MVTimeSeriesView2Factory::id() const
{
    return QStringLiteral("open-timeseries-data-2");
}

QString MVTimeSeriesView2Factory::name() const
{
    return tr("Timeseries Data 2");
}

QString MVTimeSeriesView2Factory::title() const
{
    return tr("Timeseries 2");
}

MVAbstractViewFactory::PreferredOpenLocation MVTimeSeriesView2Factory::preferredOpenLocation() const
{
    return PreferredOpenLocation::South;
}

MVAbstractView* MVTimeSeriesView2Factory::createView(MVAbstractContext* context)
{
    MCContext* c = qobject_cast<MCContext*>(context);
    Q_ASSERT(c);

    MVTimeSeriesView2* X = new MVTimeSeriesView2(c->mvContext2());
    QList<int> ks = c->mvContext2()->selectedClusters();
    X->setLabelsToView(ks.toSet());
    return X;
}

MVTimeSeriesViewIntersectFactory::MVTimeSeriesViewIntersectFactory(MVMainWindow* mw, QObject* parent)
    : MVAbstractViewFactory(mw, parent)
{
}

QString MVTimeSeriesViewIntersectFactory::id() const
{
    return QStringLiteral("open-timeseries-data-intersect");
}

QString MVTimeSeriesViewIntersectFactory::name() const
{
    return tr("Timeseries Data Intersect");
}

QString MVTimeSeriesViewIntersectFactory::title() const
{
    return tr("Timeseries Intersect");
}

MVAbstractViewFactory::PreferredOpenLocation MVTimeSeriesViewIntersectFactory::preferredOpenLocation() const
{
    return PreferredOpenLocation::North;
}

MVAbstractView* MVTimeSeriesViewIntersectFactory::createView(MVAbstractContext* context)
{
    MCContext* c = qobject_cast<MCContext*>(context);
    Q_ASSERT(c);

    MVTimeSeriesView2* X = new MVTimeSeriesView2(c->mvContext1());
    //QList<int> ks = c->mvContext1()->selectedClusters();
    QSet<int> ks1 = c->mvContext1()->selectedClusters().toSet();
    QSet<int> ks2 = c->mvContext2()->selectedClusters().toSet();
    if (ks1.isEmpty())
        ks1.insert(0);
    if (ks2.isEmpty())
        ks2.insert(0);
    //X->setLabelsToView(ks.toSet());
    DiskReadMda matched_firings = c->matchedFirings();
    QList<MVEvent> events; //change MVEvent to MVSpecialEvent so we can put on a special label
    for (bigint i = 0; i < matched_firings.N2(); i++) {
        int k1 = matched_firings.value(2, i);
        int k2 = matched_firings.value(3, i);
        if ((ks1.contains(k1)) && (ks2.contains(k2))) {
            double t0 = matched_firings.value(1, i);
            events << MVEvent(t0, 0);
        }
    }
    X->addSpecialEvents(events);
    return X;
}

////////////////////////////////////////////////////////////////////////
MVSpikeSpray1Factory::MVSpikeSpray1Factory(MVMainWindow* mw, QObject* parent)
    : MVAbstractViewFactory(mw, parent)
{
}

QString MVSpikeSpray1Factory::id() const
{
    return QStringLiteral("open-spike-spray-1");
}

QString MVSpikeSpray1Factory::name() const
{
    return tr("Spike Spray 1");
}

QString MVSpikeSpray1Factory::title() const
{
    return tr("Spike Spray 1");
}

MVAbstractViewFactory::PreferredOpenLocation MVSpikeSpray1Factory::preferredOpenLocation() const
{
    return PreferredOpenLocation::North;
}

MVAbstractView* MVSpikeSpray1Factory::createView(MVAbstractContext* context)
{
    MCContext* c = qobject_cast<MCContext*>(context);
    Q_ASSERT(c);

    MVContext* cc = c->mvContext1();

    QList<int> ks = cc->selectedClusters();
    if (ks.isEmpty())
        ks = cc->clusterVisibilityRule().subset.toList();
    qSort(ks);
    if (ks.isEmpty()) {
        QMessageBox::warning(0, "Unable to open spike spray", "You must select at least one cluster.");
        return Q_NULLPTR;
    }
    MVSpikeSprayView* X = new MVSpikeSprayView(cc);
    X->setLabelsToUse(ks.toSet());
    return X;
}

MVSpikeSpray2Factory::MVSpikeSpray2Factory(MVMainWindow* mw, QObject* parent)
    : MVAbstractViewFactory(mw, parent)
{
}

QString MVSpikeSpray2Factory::id() const
{
    return QStringLiteral("open-spike-spray-2");
}

QString MVSpikeSpray2Factory::name() const
{
    return tr("Spike Spray 2");
}

QString MVSpikeSpray2Factory::title() const
{
    return tr("Spike Spray 2");
}

MVAbstractViewFactory::PreferredOpenLocation MVSpikeSpray2Factory::preferredOpenLocation() const
{
    return PreferredOpenLocation::South;
}

MVAbstractView* MVSpikeSpray2Factory::createView(MVAbstractContext* context)
{
    MCContext* c = qobject_cast<MCContext*>(context);
    Q_ASSERT(c);

    MVContext* cc = c->mvContext2();

    QList<int> ks = cc->selectedClusters();
    if (ks.isEmpty())
        ks = cc->clusterVisibilityRule().subset.toList();
    qSort(ks);
    if (ks.isEmpty()) {
        QMessageBox::warning(0, "Unable to open spike spray", "You must select at least one cluster.");
        return Q_NULLPTR;
    }
    MVSpikeSprayView* X = new MVSpikeSprayView(cc);
    X->setLabelsToUse(ks.toSet());
    return X;
}

////////////////////////////////////////////////////////////////////////

MVPCAFeatures1Factory::MVPCAFeatures1Factory(MVMainWindow* mw, QObject* parent)
    : MVAbstractViewFactory(mw, parent)
{
}

QString MVPCAFeatures1Factory::id() const
{
    return QStringLiteral("open-pca-features-1");
}

QString MVPCAFeatures1Factory::name() const
{
    return tr("PCA Features 1");
}

QString MVPCAFeatures1Factory::title() const
{
    return tr("PCA features 1");
}

MVAbstractView* MVPCAFeatures1Factory::createView(MVAbstractContext* context)
{
    MCContext* c = qobject_cast<MCContext*>(context);
    Q_ASSERT(c);

    QList<MCCluster> clusters = c->selectedClusters();
    QList<int> ks;
    foreach (MCCluster C, clusters) {
        if (C.firings_num == 1)
            ks << C.num;
    }
    qSort(ks);
    if (ks.count() == 0) {
        QMessageBox::information(0, "Unable to open clusters", "You must select at least one cluster.");
        return Q_NULLPTR;
    }
    MVClusterWidget* X = new MVClusterWidget(c->mvContext1());
    X->setLabelsToUse(ks);
    X->setFeatureMode("pca");
    return X;
}

bool MVPCAFeatures1Factory::isEnabled(MVAbstractContext* context) const
{
    MCContext* c = qobject_cast<MCContext*>(context);
    Q_ASSERT(c);

    QList<MCCluster> clusters = c->selectedClusters();
    QList<int> ks;
    foreach (MCCluster C, clusters) {
        if (C.firings_num == 1)
            ks << C.num;
    }

    return (!ks.isEmpty());
}

////////////////////////////////////////////////////////////////////////

MVPCAFeatures2Factory::MVPCAFeatures2Factory(MVMainWindow* mw, QObject* parent)
    : MVAbstractViewFactory(mw, parent)
{
}

QString MVPCAFeatures2Factory::id() const
{
    return QStringLiteral("open-pca-features-2");
}

QString MVPCAFeatures2Factory::name() const
{
    return tr("PCA Features 2");
}

QString MVPCAFeatures2Factory::title() const
{
    return tr("PCA features 2");
}

MVAbstractView* MVPCAFeatures2Factory::createView(MVAbstractContext* context)
{
    MCContext* c = qobject_cast<MCContext*>(context);
    Q_ASSERT(c);

    QList<MCCluster> clusters = c->selectedClusters();
    QList<int> ks;
    foreach (MCCluster C, clusters) {
        if (C.firings_num == 2)
            ks << C.num;
    }
    qSort(ks);
    if (ks.count() == 0) {
        QMessageBox::information(0, "Unable to open clusters", "You must select at least one cluster.");
        return Q_NULLPTR;
    }
    MVClusterWidget* X = new MVClusterWidget(c->mvContext2());
    X->setLabelsToUse(ks);
    X->setFeatureMode("pca");

    return X;
}

bool MVPCAFeatures2Factory::isEnabled(MVAbstractContext* context) const
{
    MCContext* c = qobject_cast<MCContext*>(context);
    Q_ASSERT(c);

    QList<MCCluster> clusters = c->selectedClusters();
    QList<int> ks;
    foreach (MCCluster C, clusters) {
        if (C.firings_num == 2)
            ks << C.num;
    }

    return (!ks.isEmpty());
}

////////////////////////////////////////////////////////////////////////
CompareClustersFactory::CompareClustersFactory(MVMainWindow* mw, QObject* parent)
    : MVAbstractViewFactory(mw, parent)
{
}

QString CompareClustersFactory::id() const
{
    return QStringLiteral("open-compare-clusters");
}

QString CompareClustersFactory::name() const
{
    return tr("Compare Clusters");
}

QString CompareClustersFactory::title() const
{
    return tr("Compare Clusters");
}

MVAbstractViewFactory::PreferredOpenLocation CompareClustersFactory::preferredOpenLocation() const
{
    return PreferredOpenLocation::Floating;
}

MVAbstractView* CompareClustersFactory::createView(MVAbstractContext* context)
{
    MCContext* mc_context = qobject_cast<MCContext*>(context);
    if (!mc_context)
        return 0;

    CompareClusterView* X = new CompareClusterView(mc_context);

    return X;
}
