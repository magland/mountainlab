/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
** Created: 7/29/2016
*******************************************************/

#include "isolationmatrixview.h"

#include <QHBoxLayout>
#include <QLabel>
#include <QRadioButton>
#include <mountainprocessrunner.h>
#include <taskprogress.h>
#include "actionfactory.h"
#include "matrixview.h"
#include "get_sort_indices.h"

struct CMVControlBar {
    QWidget* widget;
    QMap<QString, QRadioButton*> permutation_buttons;
    CMVControlBar(IsolationMatrixView* q)
    {
        widget = new QWidget;
        widget->setFixedHeight(50);
        widget->setFont(QFont("Arial", 12));
        QHBoxLayout* hlayout = new QHBoxLayout;
        widget->setLayout(hlayout);

        hlayout->addWidget(new QLabel("Permutation:"));
        {
            QRadioButton* B = new QRadioButton("None");
            B->setProperty("name", "none");
            hlayout->addWidget(B);
            permutation_buttons[B->property("name").toString()] = B;
        }
        {
            QRadioButton* B = new QRadioButton("Row");
            B->setProperty("name", "row");
            hlayout->addWidget(B);
            permutation_buttons[B->property("name").toString()] = B;
        }
        {
            QRadioButton* B = new QRadioButton("Column");
            B->setProperty("name", "column");
            hlayout->addWidget(B);
            permutation_buttons[B->property("name").toString()] = B;
        }
        {
            QRadioButton* B = new QRadioButton("Both1");
            B->setProperty("name", "both1");
            hlayout->addWidget(B);
            permutation_buttons[B->property("name").toString()] = B;
        }
        {
            QRadioButton* B = new QRadioButton("Both2");
            B->setProperty("name", "both2");
            hlayout->addWidget(B);
            permutation_buttons[B->property("name").toString()] = B;
        }
        foreach (QRadioButton* B, permutation_buttons) {
            QObject::connect(B, SIGNAL(clicked(bool)), q, SLOT(slot_permutation_mode_button_clicked()));
        }

        hlayout->addStretch();

        permutation_buttons["column"]->setChecked(true);
    }
    IsolationMatrixView::PermutationMode permutationMode()
    {
        QString str;
        foreach (QRadioButton* B, permutation_buttons) {
            if (B->isChecked())
                str = B->property("name").toString();
        }
        if (str == "none")
            return IsolationMatrixView::NoPermutation;
        if (str == "row")
            return IsolationMatrixView::RowPermutation;
        if (str == "column")
            return IsolationMatrixView::ColumnPermutation;
        if (str == "both1")
            return IsolationMatrixView::BothRowBasedPermutation;
        if (str == "both2")
            return IsolationMatrixView::BothColumnBasedPermutation;
        return IsolationMatrixView::NoPermutation;
    }
};

class IsolationMatrixViewCalculator {
public:
    //input
    QString mlproxy_url;
    DiskReadMda32 timeseries;
    DiskReadMda firings;
    int clip_size = 50;
    double add_noise_level = 0.25;
    QList<int> cluster_numbers;

    //output
    DiskReadMda confusion_matrix;

    virtual void compute();
};

class IsolationMatrixViewPrivate {
public:
    IsolationMatrixView* q;

    IsolationMatrixViewCalculator m_calculator;

    Mda m_confusion_matrix;
    QList<int> m_optimal_label_map;
    MatrixView* m_matrix_view; //raw numbers
    MatrixView* m_matrix_view_rn; //row normalized
    MatrixView* m_matrix_view_cn; //column normalized
    QList<MatrixView*> m_all_matrix_views;
    CMVControlBar* m_control_bar;

    Mda row_normalize(const Mda& A);
    Mda column_normalize(const Mda& A);
    void update_permutations();
    void set_current_clusters(int k1, int k2);
};

IsolationMatrixView::IsolationMatrixView(MVContext* mvcontext)
    : MVAbstractView(mvcontext)
{
    d = new IsolationMatrixViewPrivate;
    d->q = this;

    d->m_matrix_view = new MatrixView;
    d->m_matrix_view->setMode(MatrixView::CountsMode);
    d->m_matrix_view->setTitle("Isolation matrix");
    d->m_matrix_view_rn = new MatrixView;
    d->m_matrix_view_rn->setMode(MatrixView::PercentMode);
    d->m_matrix_view_rn->setTitle("Row-normalized");
    d->m_matrix_view_cn = new MatrixView;
    d->m_matrix_view_cn->setMode(MatrixView::PercentMode);
    d->m_matrix_view_cn->setTitle("Column-normalized");

    d->m_all_matrix_views << d->m_matrix_view << d->m_matrix_view_cn << d->m_matrix_view_rn;

    d->m_control_bar = new CMVControlBar(this);

    QHBoxLayout* hlayout = new QHBoxLayout;
    hlayout->addWidget(d->m_matrix_view);
    hlayout->addWidget(d->m_matrix_view_rn);
    hlayout->addWidget(d->m_matrix_view_cn);

    QVBoxLayout* vlayout = new QVBoxLayout;
    this->setLayout(vlayout);
    vlayout->addLayout(hlayout);
    vlayout->addWidget(d->m_control_bar->widget);

    if (!mvContext()) {
        qCritical() << "mvContext is null" << __FILE__ << __LINE__;
        return;
    }

    foreach (MatrixView* MV, d->m_all_matrix_views) {
        QObject::connect(MV, SIGNAL(currentElementChanged()), this, SLOT(slot_matrix_view_current_element_changed()));
    }

    this->recalculateOn(mvContext(), SIGNAL(firingsChanged()), false);

    //Important to do a queued connection here! because we are changing two things at the same time
    QObject::connect(mvContext(), SIGNAL(currentClusterChanged()), this, SLOT(slot_update_current_elements_based_on_context()), Qt::QueuedConnection);

    this->recalculate();
}

IsolationMatrixView::~IsolationMatrixView()
{
    this->stopCalculation();
    delete d;
}

void IsolationMatrixView::prepareCalculation()
{
    if (!mvContext())
        return;

    MVContext* c = qobject_cast<MVContext*>(mvContext());
    Q_ASSERT(c);

    d->m_calculator.clip_size = 50;
    d->m_calculator.add_noise_level = 0.25;
    d->m_calculator.timeseries = c->currentTimeseries();
    d->m_calculator.firings = c->firings();
    d->m_calculator.cluster_numbers = c->selectedClusters();
}

void IsolationMatrixView::runCalculation()
{
    d->m_calculator.compute();
}

void IsolationMatrixView::onCalculationFinished()
{
    DiskReadMda confusion_matrix = d->m_calculator.confusion_matrix;
    int A1 = confusion_matrix.N1();
    int A2 = confusion_matrix.N2();
    if (!confusion_matrix.readChunk(d->m_confusion_matrix, 0, 0, A1, A2)) {
        qWarning() << "Unable to read chunk for confusion matrix in isolation matrix view";
        return;
    }

    d->m_optimal_label_map.clear();
    for (int k = 1; k <= 1000; k++) {
        d->m_optimal_label_map << k;
    }

    d->m_matrix_view->setMatrix(d->m_confusion_matrix);
    d->m_matrix_view->setValueRange(0, d->m_confusion_matrix.maximum());
    d->m_matrix_view_rn->setMatrix(d->row_normalize(d->m_confusion_matrix));
    d->m_matrix_view_cn->setMatrix(d->column_normalize(d->m_confusion_matrix));

    QStringList row_labels, col_labels;
    for (int m = 0; m < A1 - 1; m++) {
        int cnum = m + 1;
        if (m < d->m_calculator.cluster_numbers.count())
            cnum = d->m_calculator.cluster_numbers[m];
        row_labels << QString("%1").arg(cnum);
    }
    for (int n = 0; n < A2 - 1; n++) {
        int cnum = n + 1;
        if (n < d->m_calculator.cluster_numbers.count())
            cnum = d->m_calculator.cluster_numbers[n];
        col_labels << QString("%1").arg(cnum);
    }

    foreach (MatrixView* MV, d->m_all_matrix_views) {
        MV->setLabels(row_labels, col_labels);
    }

    d->update_permutations();
    slot_update_current_elements_based_on_context();
}

void IsolationMatrixView::keyPressEvent(QKeyEvent* evt)
{
    Q_UNUSED(evt)
    /*
    if (evt->key() == Qt::Key_Up) {
        slot_vertical_zoom_in();
    }
    if (evt->key() == Qt::Key_Down) {
        slot_vertical_zoom_out();
    }
    */
}

void IsolationMatrixView::prepareMimeData(QMimeData& mimeData, const QPoint& pos)
{
    /*
    QByteArray ba;
    QDataStream ds(&ba, QIODevice::WriteOnly);
    ds << mvContext()->selectedClusters();
    mimeData.setData("application/x-msv-clusters", ba); // selected cluster data
    */

    MVAbstractView::prepareMimeData(mimeData, pos); // call base class implementation
}

void IsolationMatrixView::slot_permutation_mode_button_clicked()
{
    d->update_permutations();
}

void IsolationMatrixView::slot_matrix_view_current_element_changed()
{
    MatrixView* MV = qobject_cast<MatrixView*>(sender());
    if (!MV)
        return;

    MVContext* c = qobject_cast<MVContext*>(mvContext());
    Q_ASSERT(c);

    QPoint a = MV->currentElement();
    if ((a.x() >= 0) && (a.y() >= 0)) {
        c->clickCluster(a.x() + 1, Qt::NoModifier);
    }
    else {
        c->setCurrentCluster(-1);
    }
}

void IsolationMatrixView::slot_update_current_elements_based_on_context()
{
    MVContext* c = qobject_cast<MVContext*>(mvContext());
    Q_ASSERT(c);

    int k1 = c->currentCluster();
    if (k1 <= 0) {
        foreach (MatrixView* MV, d->m_all_matrix_views) {
            MV->setCurrentElement(QPoint(-1, -1));
        }
    }
    else {
        foreach (MatrixView* MV, d->m_all_matrix_views) {
            MV->setCurrentElement(QPoint(k1 - 1, k1 - 1));
        }
    }
}

Mda IsolationMatrixViewPrivate::row_normalize(const Mda& A)
{
    Mda B = A;
    int M = B.N1();
    int N = B.N2();
    for (int m = 0; m < M; m++) {
        double sum = 0;
        for (int n = 0; n < N; n++) {
            sum += B.value(m, n);
        }
        if (sum) {
            for (int n = 0; n < N; n++) {
                B.setValue(B.value(m, n) / sum, m, n);
            }
        }
    }
    return B;
}

Mda IsolationMatrixViewPrivate::column_normalize(const Mda& A)
{
    Mda B = A;
    int M = B.N1();
    int N = B.N2();
    for (int n = 0; n < N; n++) {
        double sum = 0;
        for (int m = 0; m < M; m++) {
            sum += B.value(m, n);
        }
        if (sum) {
            for (int m = 0; m < M; m++) {
                B.setValue(B.value(m, n) / sum, m, n);
            }
        }
    }
    return B;
}

void IsolationMatrixViewPrivate::update_permutations()
{
    int M = m_confusion_matrix.N1();
    int N = m_confusion_matrix.N2();

    QVector<int> perm_rows;
    QVector<int> perm_cols;

    IsolationMatrixView::PermutationMode permutation_mode = m_control_bar->permutationMode();

    if (permutation_mode == IsolationMatrixView::RowPermutation) {
        perm_rows.fill(-1, M);
        for (int i = 0; i < M - 1; i++) {
            int val = m_optimal_label_map.value(i);
            if ((val - 1 >= 0) && (val - 1 < M - 1)) {
                perm_rows[i] = val - 1;
            }
            else {
                perm_rows[i] = -1; //will be filled in later
            }
        }
        for (int i = 0; i < M - 1; i++) {
            if (perm_rows[i] == -1) {
                for (int j = 0; j < M - 1; j++) {
                    if (!perm_rows.contains(j)) {
                        perm_rows[i] = j;
                        break;
                    }
                }
            }
        }
        perm_rows[M - 1] = M - 1; //unclassified row
    }

    if (permutation_mode == IsolationMatrixView::ColumnPermutation) {
        perm_cols.fill(-1, N);
        for (int i = 0; i < N - 1; i++) {
            int val = m_optimal_label_map.indexOf(i + 1);
            if ((val >= 0) && (val < N - 1)) {
                perm_cols[i] = val;
            }
            else {
                perm_cols[i] = -1; //will be filled in later
            }
        }
        for (int i = 0; i < N - 1; i++) {
            if (perm_cols[i] == -1) {
                for (int j = 0; j < N - 1; j++) {
                    if (!perm_cols.contains(j)) {
                        perm_cols[i] = j;
                        break;
                    }
                }
            }
        }
        perm_cols[N - 1] = N - 1; //unclassified column
    }

    if (permutation_mode == IsolationMatrixView::BothRowBasedPermutation) {
        Mda A = row_normalize(m_confusion_matrix);

        QVector<double> diag_entries(M - 1);
        for (int m = 0; m < M - 1; m++) {
            if (m_optimal_label_map.value(m) - 1 >= 0) {
                diag_entries[m] = A.value(m, m_optimal_label_map[m] - 1);
            }
            else {
                diag_entries[m] = 0;
            }
        }
        QList<int> sort_inds = get_sort_indices(diag_entries);
        perm_rows.fill(-1, M);
        perm_cols.fill(-1, N);

        for (int ii = 0; ii < sort_inds.count(); ii++) {
            int m = sort_inds[sort_inds.count() - 1 - ii];
            perm_rows[m] = ii;
            int n = m_optimal_label_map.value(m) - 1;
            if (n >= 0) {
                perm_cols[n] = ii;
            }
        }
        for (int i = 0; i < N - 1; i++) {
            if ((perm_cols[i] == -1) || (perm_cols[i] >= N)) {
                for (int j = 0; j < N - 1; j++) {
                    if (!perm_cols.contains(j)) {
                        perm_cols[i] = j;
                        break;
                    }
                }
            }
        }
        perm_rows[M - 1] = M - 1; //unclassified row
        perm_cols[N - 1] = N - 1; //unclassified column
    }

    if (permutation_mode == IsolationMatrixView::BothColumnBasedPermutation) {
        Mda A = column_normalize(m_confusion_matrix);

        QVector<double> diag_entries(N - 1);
        for (int n = 0; n < N - 1; n++) {
            if (m_optimal_label_map.indexOf(n + 1) >= 0) { //something maps to it
                diag_entries[n] = A.value(m_optimal_label_map.indexOf(n + 1), n);
            }
            else {
                diag_entries[n] = 0;
            }
        }
        QList<int> sort_inds = get_sort_indices(diag_entries);
        perm_rows.fill(-1, M);
        perm_cols.fill(-1, N);

        for (int ii = 0; ii < sort_inds.count(); ii++) {
            int n = sort_inds[sort_inds.count() - 1 - ii];
            perm_cols[n] = ii;
            int m = m_optimal_label_map.indexOf(n + 1);
            if (m >= 0) {
                perm_rows[m] = ii;
            }
        }
        for (int i = 0; i < M - 1; i++) {
            if ((perm_rows[i] == -1) || (perm_rows[i] >= M)) {
                for (int j = 0; j < M - 1; j++) {
                    if (!perm_rows.contains(j)) {
                        perm_rows[i] = j;
                        break;
                    }
                }
            }
        }
        perm_rows[M - 1] = M - 1; //unclassified row
        perm_cols[N - 1] = N - 1; //unclassified column
    }

    if (perm_rows.isEmpty()) {
        for (int i = 0; i < M; i++)
            perm_rows << i;
    }
    if (perm_cols.isEmpty()) {
        for (int i = 0; i < N; i++)
            perm_cols << i;
    }

    //put the zero rows/columns at the end
    QList<int> row_sums, col_sums;
    for (int i = 0; i < M; i++)
        row_sums << 0;
    for (int i = 0; i < N; i++)
        col_sums << 0;
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            row_sums[i] += m_confusion_matrix.value(i, j);
            col_sums[j] += m_confusion_matrix.value(i, j);
        }
    }
    for (int i = 0; i < M; i++) {
        if (row_sums[i] == 0) {
            for (int k = 0; k < M; k++) {
                if (perm_rows[k] > perm_rows[i])
                    perm_rows[k]--;
            }
            perm_rows[i] = -1;
        }
    }
    for (int i = 0; i < N; i++) {
        if (col_sums[i] == 0) {
            for (int k = 0; k < N; k++) {
                if (perm_cols[k] > perm_cols[i])
                    perm_cols[k]--;
            }
            perm_cols[i] = -1;
        }
    }

    foreach (MatrixView* MV, m_all_matrix_views) {
        MV->setIndexPermutations(perm_rows, perm_cols);
    }
}

void IsolationMatrixViewPrivate::set_current_clusters(int k1, int k2)
{
    if ((k1 < 0) || (k2 < 0)) {
        foreach (MatrixView* V, m_all_matrix_views) {
            V->setCurrentElement(QPoint(-1, -1));
        }
    }
    else {
        foreach (MatrixView* V, m_all_matrix_views) {
            V->setCurrentElement(QPoint(k1 - 1, k2 - 1));
        }
    }
}

void IsolationMatrixViewCalculator::compute()
{
    MountainProcessRunner MPR;
    MPR.setProcessorName("");
    MPR.setProcessorName("noise_nearest");
    QMap<QString, QVariant> params;
    params["timeseries"] = timeseries.makePath();
    params["firings"] = firings.makePath();
    params["clip_size"] = clip_size;
    params["add_noise_level"] = add_noise_level;
    params["cluster_numbers"] = MLUtil::intListToStringList(cluster_numbers).join(",");
    MPR.setInputParameters(params);
    MPR.setMLProxyUrl(mlproxy_url);

    QString confusion_matrix_fname = MPR.makeOutputFilePath("confusion_matrix");

    MPR.runProcess();
    confusion_matrix.setPath(confusion_matrix_fname);
}
