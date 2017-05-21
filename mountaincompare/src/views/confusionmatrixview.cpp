/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
** Created: 7/29/2016
*******************************************************/

#include "confusionmatrixview.h"

#include <QFileDialog>
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
    CMVControlBar(ConfusionMatrixView* q)
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
    ConfusionMatrixView::PermutationMode permutationMode()
    {
        QString str;
        foreach (QRadioButton* B, permutation_buttons) {
            if (B->isChecked())
                str = B->property("name").toString();
        }
        if (str == "none")
            return ConfusionMatrixView::NoPermutation;
        if (str == "row")
            return ConfusionMatrixView::RowPermutation;
        if (str == "column")
            return ConfusionMatrixView::ColumnPermutation;
        if (str == "both1")
            return ConfusionMatrixView::BothRowBasedPermutation;
        if (str == "both2")
            return ConfusionMatrixView::BothColumnBasedPermutation;
        return ConfusionMatrixView::NoPermutation;
    }
};

class ConfusionMatrixViewPrivate {
public:
    ConfusionMatrixView* q;

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

ConfusionMatrixView::ConfusionMatrixView(MVAbstractContext* mvcontext)
    : MVAbstractView(mvcontext)
{
    d = new ConfusionMatrixViewPrivate;
    d->q = this;

    d->m_matrix_view = new MatrixView;
    d->m_matrix_view->setMode(MatrixView::CountsMode);
    d->m_matrix_view->setTitle("Confusion matrix");
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

    if (!mcContext()) {
        qCritical() << "mcContext is null" << __FILE__ << __LINE__;
        return;
    }

    foreach (MatrixView* MV, d->m_all_matrix_views) {
        QObject::connect(MV, SIGNAL(currentElementChanged()), this, SLOT(slot_matrix_view_current_element_changed()));
    }

    this->recalculateOn(mcContext(), SIGNAL(firingsChanged()), false);
    this->recalculateOn(mcContext(), SIGNAL(firings2Changed()), false);

    {
        QAction* A = new QAction(QString("Export .csv"), this);
        //A->setProperty("action_type", "toolbar");
        QObject::connect(A, SIGNAL(triggered(bool)), this, SLOT(slot_export_csv()));
        this->addAction(A);
    }
    {
        QAction* A = new QAction(QString("Export .mda"), this);
        //A->setProperty("action_type", "toolbar");
        QObject::connect(A, SIGNAL(triggered(bool)), this, SLOT(slot_export_mda()));
        this->addAction(A);
    }

    //Important to do a queued connection here! because we are changing two things at the same time
    QObject::connect(mcContext(), SIGNAL(currentClusterChanged()), this, SLOT(slot_update_current_elements_based_on_context()), Qt::QueuedConnection);
    //QObject::connect(mcContext(), SIGNAL(currentCluster2Changed()), this, SLOT(slot_update_current_elements_based_on_context()), Qt::QueuedConnection);

    QObject::connect(mcContext(), SIGNAL(selectedClustersChanged()), this, SLOT(slot_update_selected()));

    this->recalculate();
}

ConfusionMatrixView::~ConfusionMatrixView()
{
    this->stopCalculation();
    delete d;
}

void ConfusionMatrixView::prepareCalculation()
{
    if (!mcContext())
        return;
}

void ConfusionMatrixView::runCalculation()
{
    //mcContext()->computeMatchedFirings();
}

void ConfusionMatrixView::onCalculationFinished()
{
    DiskReadMda confusion_matrix = mcContext()->confusionMatrix();
    int A1 = confusion_matrix.N1();
    int A2 = confusion_matrix.N2();
    if (!confusion_matrix.readChunk(d->m_confusion_matrix, 0, 0, A1, A2)) {
        qWarning() << "Unable to read chunk of confusion matrix in ConfusionMatrixView";
        return;
    }
    d->m_optimal_label_map = mcContext()->labelMap();

    d->m_matrix_view->setMatrix(d->m_confusion_matrix);
    d->m_matrix_view->setValueRange(0, d->m_confusion_matrix.maximum());
    d->m_matrix_view_rn->setMatrix(d->row_normalize(d->m_confusion_matrix));
    d->m_matrix_view_cn->setMatrix(d->column_normalize(d->m_confusion_matrix));

    QStringList row_labels, col_labels;
    for (int m = 0; m < A1 - 1; m++) {
        row_labels << QString("%1").arg(m + 1);
    }
    for (int n = 0; n < A2 - 1; n++) {
        col_labels << QString("%1").arg(n + 1);
    }

    foreach (MatrixView* MV, d->m_all_matrix_views) {
        MV->setLabels(row_labels, col_labels);
    }

    d->update_permutations();
    slot_update_current_elements_based_on_context();
    slot_update_selected();
}

MCContext* ConfusionMatrixView::mcContext()
{
    return qobject_cast<MCContext*>(mvContext());
}

void ConfusionMatrixView::keyPressEvent(QKeyEvent* evt)
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

void ConfusionMatrixView::prepareMimeData(QMimeData& mimeData, const QPoint& pos)
{
    /*
    QByteArray ba;
    QDataStream ds(&ba, QIODevice::WriteOnly);
    ds << mvContext()->selectedClusters();
    mimeData.setData("application/x-msv-clusters", ba); // selected cluster data
    */

    MVAbstractView::prepareMimeData(mimeData, pos); // call base class implementation
}

void ConfusionMatrixView::slot_permutation_mode_button_clicked()
{
    d->update_permutations();
}

void ConfusionMatrixView::slot_matrix_view_current_element_changed()
{
    MatrixView* MV = qobject_cast<MatrixView*>(sender());
    if (!MV)
        return;
    QPoint a = MV->currentElement();
    if ((a.x() >= 0) && (a.y() >= 0)) {
        MCCluster C1, C2;
        C1.firings_num = 1;
        C1.num = a.x() + 1;
        C2.firings_num = 2;
        C2.num = a.y() + 1;
        QList<MCCluster> list;
        list << C1 << C2;
        mcContext()->setSelectedClusters(list);
        //mcContext()->clickCluster(a.x() + 1, Qt::NoModifier);
        //mcContext()->clickCluster2(a.y() + 1, Qt::NoModifier);
    }
    else {
        //mcContext()->setCurrentCluster(-1);
        //mcContext()->setCurrentCluster2(-1);
    }
}

void ConfusionMatrixView::slot_update_current_elements_based_on_context()
{
    /*
    int k1 = mcContext()->currentCluster();
    int k2 = mcContext()->currentCluster2();
    if (k1 <= 0) {
        foreach (MatrixView* MV, d->m_all_matrix_views) {
            MV->setCurrentElement(QPoint(-1, -1));
        }
    }
    else {
        if (k2 <= 0) {
            k2 = 1;
        }
        foreach (MatrixView* MV, d->m_all_matrix_views) {
            MV->setCurrentElement(QPoint(k1 - 1, k2 - 1));
        }
    }
    */
}

void ConfusionMatrixView::slot_export_csv()
{
    Mda CM = d->m_matrix_view->matrix();
    QStringList RL = d->m_matrix_view->rowLabels();
    QStringList CL = d->m_matrix_view->columnLabels();
    QVector<int> RP = d->m_matrix_view->rowIndexPermutationInv();
    QVector<int> CP = d->m_matrix_view->columnIndexPermutationInv();

    QString txt;
    for (int r = 0; r < CM.N1(); r++) {
        int r2 = RP.value(r);
        if ((r2 >= 0) && (r2 < CM.N1() - 1)) {
            txt += RL.value(r2); //row label
            for (int c = 0; c < CM.N2(); c++) {
                int c2 = CP.value(c);
                if ((c2 >= 0) && (c2 < CM.N2() - 1))
                    txt += QString(",%1").arg(CM.value(r2, c2)); //entry
            }
            txt += QString(",%1").arg(CM.value(r2, CM.N2() - 1)); //unclassified column
            txt += "\n";
        }
    }

    txt += "0"; //unclassified row
    for (int c = 0; c < CM.N2(); c++) {
        int c2 = CP.value(c);
        if ((c2 >= 0) && (c2 < CM.N2() - 1)) {
            txt += QString(",%1").arg(CM.value(CM.N1() - 1, c2));
        }
    }
    txt += ",0\n"; //nothing is doubly-unclassified

    //column labels
    txt += "0";
    for (int c = 0; c < CM.N2() - 1; c++) {
        int c2 = CP.value(c);
        if ((c2 >= 0) && (c2 < CM.N2() - 1)) {
            txt += QString(",%1").arg(CL.value(c2)); //column label
        }
    }
    txt += ",0\n";

    QString fname = QFileDialog::getSaveFileName(this, "Save confusion matrix data", "", "*.csv");
    if (fname.isEmpty())
        return;
    TextFile::write(fname, txt);
}

void ConfusionMatrixView::slot_export_mda()
{
    Mda CM = d->m_matrix_view->matrix();

    QString fname = QFileDialog::getSaveFileName(this, "Save confusion matrix data", "", "*.mda");
    if (fname.isEmpty())
        return;
    CM.write64(fname);
}

void ConfusionMatrixView::slot_update_selected()
{
    QList<MCCluster> cluster_list = mcContext()->selectedClusters();
    QSet<QPoint> pts;
    foreach (MCCluster C, cluster_list) {
        if (C.firings_num == 1) {
            for (int i = 0; i < d->m_confusion_matrix.N2(); i++) {
                QPoint pt(C.num - 1, i);
                pts.insert(pt);
            }
        }
        else if (C.firings_num == 2) {
            for (int i = 0; i < d->m_confusion_matrix.N1(); i++) {
                QPoint pt(i, C.num - 1);
                pts.insert(pt);
            }
        }
    }

    foreach (MatrixView* V, d->m_all_matrix_views) {
        V->setSelectedElements(pts.toList());
    }
}

ConfusionMatrixViewFactory::ConfusionMatrixViewFactory(MVMainWindow* mw, QObject* parent)
    : MVAbstractViewFactory(mw, parent)
{
}

QString ConfusionMatrixViewFactory::id() const
{
    return QStringLiteral("open-confusion-matrix");
}

QString ConfusionMatrixViewFactory::name() const
{
    return tr("Confusion Matrix");
}

QString ConfusionMatrixViewFactory::title() const
{
    return tr("Confusion Matrix");
}

MVAbstractViewFactory::PreferredOpenLocation ConfusionMatrixViewFactory::preferredOpenLocation() const
{
    return PreferredOpenLocation::Floating;
}

MVAbstractView* ConfusionMatrixViewFactory::createView(MVAbstractContext* context)
{
    ConfusionMatrixView* X = new ConfusionMatrixView(context);

    if (!X->mcContext()) {
        qCritical() << "mcContext is null" << __FILE__ << __LINE__;
        delete X;
        return 0;
    }

    return X;
}

Mda ConfusionMatrixViewPrivate::row_normalize(const Mda& A)
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

Mda ConfusionMatrixViewPrivate::column_normalize(const Mda& A)
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

void ConfusionMatrixViewPrivate::update_permutations()
{
    int M = m_confusion_matrix.N1();
    int N = m_confusion_matrix.N2();

    QVector<int> perm_rows;
    QVector<int> perm_cols;

    ConfusionMatrixView::PermutationMode permutation_mode = m_control_bar->permutationMode();

    if (permutation_mode == ConfusionMatrixView::RowPermutation) {
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

    if (permutation_mode == ConfusionMatrixView::ColumnPermutation) {
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

    if (permutation_mode == ConfusionMatrixView::BothRowBasedPermutation) {
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
            if ((n >= 0) && (n < N - 1)) {
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

    if (permutation_mode == ConfusionMatrixView::BothColumnBasedPermutation) {
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
            if ((m >= 0) && (m < M - 1)) {
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

void ConfusionMatrixViewPrivate::set_current_clusters(int k1, int k2)
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
