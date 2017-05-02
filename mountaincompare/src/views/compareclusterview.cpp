/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
** Created: 8/18/2016
*******************************************************/

#include "compareclusterview.h"
#include "mvclusterview.h"

#include <QHBoxLayout>
#include <mountainprocessrunner.h>
#include <taskprogress.h>

class CompareClusterViewCalculator {
public:
    //input
    DiskReadMda32 timeseries;
    DiskReadMda firings_merged;
    QList<int> ks_1;
    int clip_size;

    //output
    Mda features;
    QVector<int> labels;

    void compute();
};

class CompareClusterViewPrivate {
public:
    CompareClusterView* q;
    MVClusterView* m_view1;
    MVClusterView* m_view2;
    QList<MVClusterView*> m_all_views;

    CompareClusterViewCalculator m_calculator;
};

CompareClusterView::CompareClusterView(MCContext* context)
    : MVAbstractView(context)
{
    d = new CompareClusterViewPrivate;
    d->q = this;

    d->m_view1 = new MVClusterView(context);
    d->m_view1->setMode(MVCV_MODE_HEAT_DENSITY);
    d->m_view2 = new MVClusterView(context);
    d->m_view2->setMode(MVCV_MODE_LABEL_COLORS);
    d->m_all_views << d->m_view1 << d->m_view2;

    QHBoxLayout* hlayout = new QHBoxLayout;
    this->setLayout(hlayout);

    hlayout->addWidget(d->m_view1);
    hlayout->addWidget(d->m_view2);

    recalculate();
}

CompareClusterView::~CompareClusterView()
{
    this->stopCalculation();
    delete d;
}

void CompareClusterView::prepareCalculation()
{
    d->m_calculator.timeseries = mcContext()->currentTimeseries();
    //d->m_calculator.firings_merged = mcContext()->firingsMerged();
    //d->m_calculator.ks_1 = mcContext()->selectedClusters();
    d->m_calculator.clip_size = mcContext()->option("clip_size").toInt();
}

void CompareClusterView::runCalculation()
{
    d->m_calculator.compute();
}

void CompareClusterView::onCalculationFinished()
{
    Mda X = d->m_calculator.features;
    double max_abs_val = qMax(qAbs(X.minimum()), qAbs(X.maximum()));
    AffineTransformation T;
    T.setIdentity();
    if (max_abs_val) {
        double factor = 1 / 1.2;
        T.scale(1.0 / max_abs_val * factor, 1.0 / max_abs_val * factor, 1.0 / max_abs_val * factor);
    }
    foreach (MVClusterView* V, d->m_all_views) {
        V->setLabels(d->m_calculator.labels);
        V->setData(X);
        V->setTransformation(T);
    }
}

MCContext* CompareClusterView::mcContext()
{
    return qobject_cast<MCContext*>(mvContext());
}

DiskReadMda extract_firings_subset(DiskReadMda firings, QList<int> ks, int labels_row_number)
{
    TaskProgress task("extract_firings_subset");
    task.log() << "Dimensions of firings:" << firings.N1() << firings.N2();
    task.log() << "ks:" << ks;
    task.log() << "labels_row_number:" << labels_row_number;
    MountainProcessRunner MPR;
    MPR.setProcessorName("mv_subfirings");
    QMap<QString, QVariant> params;
    params["firings"] = firings.makePath();
    params["labels"] = MLUtil::intListToStringList(ks).join(",");
    ;
    MPR.setInputParameters(params);
    QString firings_out_path = MPR.makeOutputFilePath("firings_out");
    MPR.runProcess();
    DiskReadMda ret(firings_out_path);

    task.log() << "Dimensions of firings_out:" << ret.N1() << ret.N2();
    return ret;
}

DiskReadMda compute_firings_features(DiskReadMda32 timeseries, DiskReadMda firings, int npca, int clip_size)
{
    TaskProgress task("compute_firings_features");
    task.log() << "timeseries dimensions:" << timeseries.N1() << timeseries.N2();
    task.log() << "firings dimensions:" << firings.N1() << firings.N2();
    task.log() << "npca/clip_size:" << npca << clip_size;
    MountainProcessRunner MPR;
    MPR.setProcessorName("extract_clips_features");
    QMap<QString, QVariant> params;
    params["timeseries"] = timeseries.makePath();
    params["firings"] = firings.makePath();
    params["clip_size"] = clip_size;
    params["num_features"] = npca;
    MPR.setInputParameters(params);
    QString features_path = MPR.makeOutputFilePath("features");
    MPR.runProcess();
    DiskReadMda ret(features_path);
    task.log() << "Dimensions of features (output):" << ret.N1() << ret.N2();
    return ret;
}

QVector<int> extract_labels(DiskReadMda firings, int labels_row_num)
{
    QVector<int> ret(firings.N2());
    for (int i = 0; i < firings.N2(); i++) {
        ret[i] = firings.value(labels_row_num, i);
    }
    return ret;
}

void CompareClusterViewCalculator::compute()
{
    DiskReadMda firings_merged_subset = extract_firings_subset(firings_merged, ks_1, 2);

    DiskReadMda features0;
    features0 = compute_firings_features(timeseries, firings_merged_subset, 3, clip_size);
    labels = extract_labels(firings_merged_subset, 3);
    if (!features0.readChunk(features, 0, 0, features0.N1(), features0.N2())) {
        qWarning() << "Unable to read chunk in CompareClusterViewCalculator";
    }
}
