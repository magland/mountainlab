/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
** Created: 7/5/2016
*******************************************************/

#include <QLabel>
#include <QVBoxLayout>
#include <QPushButton>
#include <QMessageBox>
#include <flowlayout.h>
#include <QLineEdit>
#include <mountainprocessrunner.h>
#include <taskprogress.h>
#include <QThread>
#include "guidev1.h"

#include "mv_compute_templates.h"

class GuideV1Private {
public:
    GuideV1* q;

    MVContext* m_context;
    MVMainWindow* m_main_window;

    QWizardPage* make_page_1();
    QWizardPage* make_page_2();
    QWizardPage* make_page_3();
    QWizardPage* make_page_4();

    QAbstractButton* make_instructions_button(QString text, QString instructions);
    QAbstractButton* make_open_view_button(QString text, QString view_id, QString container_name);

    bool m_cluster_channel_matrix_computed = false;
    DiskReadMda m_cluster_channel_matrix;
    QLineEdit* m_channel_edit; //this goes in separate class for page 4
    QAbstractButton* m_prev_channel_button;
    QAbstractButton* m_next_channel_button;

    void show_instructions(QString title, QString instructions);
    void compute_cluster_channel_matrix();
    void update_buttons();
};

GuideV1::GuideV1(MVContext* mvcontext, MVMainWindow* mw)
{
    d = new GuideV1Private;
    d->q = this;
    d->m_context = mvcontext;
    d->m_main_window = mw;

    this->addPage(d->make_page_1());
    this->addPage(d->make_page_2());
    this->addPage(d->make_page_3());
    this->addPage(d->make_page_4());

    this->resize(800, 600);

    this->setWindowTitle("Guide Version 1");
}

GuideV1::~GuideV1()
{
    delete d;
}

void GuideV1::slot_button_clicked()
{
    QString action = sender()->property("action").toString();
    if (action == "open_view") {
        d->m_main_window->setCurrentContainerName(sender()->property("container-name").toString());
        d->m_main_window->openView(sender()->property("view-id").toString());
    }
    else if (action == "show_instructions") {
        d->show_instructions(sender()->property("title").toString(), sender()->property("instructions").toString());
    }
}

QSet<int> find_clusters_to_use(int ch, const DiskReadMda& cluster_channel_matrix)
{
    int K = cluster_channel_matrix.N1();
    int M = cluster_channel_matrix.N2();
    QVector<double> maxvals;
    for (int k = 0; k < K; k++) {
        double maxval = 0;
        for (int m = 0; m < M; m++) {
            maxval = qMax(maxval, cluster_channel_matrix.value(k, m));
        }
        maxvals << maxval;
    }

    QSet<int> ret;
    for (int k = 0; k < K; k++) {
        if (cluster_channel_matrix.value(k, ch) > maxvals[k] * 0.8) {
            ret.insert(k + 1);
        }
    }
    return ret;
}

void GuideV1::slot_next_channel(int offset)
{
    if (!d->m_cluster_channel_matrix_computed) {
        d->m_channel_edit->setText("Calculating....");
        d->compute_cluster_channel_matrix();
        return;
    }
    int ch = d->m_channel_edit->text().toInt();
    ch += offset;
    d->m_channel_edit->setText(QString("%1").arg(ch));
    ClusterVisibilityRule rule = d->m_context->clusterVisibilityRule();
    rule.use_subset = true;
    rule.subset = find_clusters_to_use(ch - 1, d->m_cluster_channel_matrix);
    d->m_context->setClusterVisibilityRule(rule);
    d->m_context->setSelectedClusters(QList<int>());

    d->update_buttons();
}

void GuideV1::slot_previous_channel()
{
    slot_next_channel(-1);
}

QWizardPage* GuideV1Private::make_page_1()
{
    QWizardPage* page = new QWizardPage;
    page->setTitle("Cluster detail view");

    QLabel* label = new QLabel(TextFile::read(":/guides/guidev1/page_1.txt"));
    label->setWordWrap(true);

    QVBoxLayout* layout = new QVBoxLayout;
    layout->addWidget(label);
    page->setLayout(layout);

    FlowLayout* flayout = new FlowLayout;
    flayout->addWidget(make_open_view_button("Cluster details", "open-cluster-details", "north"));
    layout->addLayout(flayout);

    return page;
}

QWizardPage* GuideV1Private::make_page_2()
{
    QWizardPage* page = new QWizardPage;
    page->setTitle("Event Filter");

    QLabel* label = new QLabel(TextFile::read(":/guides/guidev1/page_2.txt"));
    label->setWordWrap(true);

    QVBoxLayout* layout = new QVBoxLayout;
    layout->addWidget(label);
    page->setLayout(layout);

    //FlowLayout* flayout = new FlowLayout;
    //flayout->addWidget(make_open_view_button("Cluster details", "open-cluster-details", "north"));
    //layout->addLayout(flayout);

    return page;
}

QWizardPage* GuideV1Private::make_page_3()
{
    QWizardPage* page = new QWizardPage;
    page->setTitle("Auto-Correlograms");

    QLabel* label = new QLabel(TextFile::read(":/guides/guidev1/page_3.txt"));
    label->setWordWrap(true);

    QVBoxLayout* layout = new QVBoxLayout;
    layout->addWidget(label);
    page->setLayout(layout);

    FlowLayout* flayout = new FlowLayout;
    flayout->addWidget(make_open_view_button("Auto-correlograms", "open-auto-correlograms", "south"));
    layout->addLayout(flayout);

    return page;
}

QWizardPage* GuideV1Private::make_page_4()
{
    /// TODO, page_4 should be it's own class
    QWizardPage* page = new QWizardPage;
    page->setTitle("Channel Groups");

    QLabel* label = new QLabel(TextFile::read(":/guides/guidev1/page_4.txt"));
    label->setWordWrap(true);

    QHBoxLayout* hlayout = new QHBoxLayout;
    {
        hlayout->addWidget(new QLabel("Channel"));
        m_channel_edit = new QLineEdit;
        m_channel_edit->setText("");
        hlayout->addWidget(m_channel_edit);
        {
            QPushButton* B = new QPushButton("Prev channel");
            QObject::connect(B, SIGNAL(clicked(bool)), q, SLOT(slot_previous_channel()));
            hlayout->addWidget(B);
            m_prev_channel_button = B;
        }
        {
            QPushButton* B = new QPushButton("Next channel");
            QObject::connect(B, SIGNAL(clicked(bool)), q, SLOT(slot_next_channel()));
            hlayout->addWidget(B);
            m_next_channel_button = B;
        }
        hlayout->addStretch();
    }
    this->update_buttons();

    FlowLayout* flayout = new FlowLayout;
    flayout->addWidget(make_open_view_button("Cluster details", "open-cluster-details", "north"));
    flayout->addWidget(make_open_view_button("Discrim histograms", "open-discrim-histograms", "south"));
    flayout->addWidget(make_open_view_button("Auto-correlograms", "open-auto-correlograms", "south"));
    flayout->addWidget(make_open_view_button("Cross-correlograms", "open-matrix-of-cross-correlograms", "south"));
    flayout->addWidget(make_open_view_button("PCA features", "open-pca-features", "south"));
    flayout->addWidget(make_open_view_button("Channel features", "open-channel-features", "south"));
    flayout->addWidget(make_open_view_button("Spike spray", "open-spike-spray", "south"));
    flayout->addWidget(make_open_view_button("Firing events", "open-firing-events", "south"));
    flayout->addWidget(make_open_view_button("Amplitude histograms", "open-amplitude-histograms", "south"));

    QVBoxLayout* layout = new QVBoxLayout;
    layout->addWidget(label);
    layout->addLayout(hlayout);
    layout->addLayout(flayout);
    page->setLayout(layout);

    return page;
}

QAbstractButton* GuideV1Private::make_instructions_button(QString text, QString instructions)
{
    QPushButton* B = new QPushButton(text);
    B->setProperty("action", "show_instructions");
    B->setProperty("instructions", instructions);
    QObject::connect(B, SIGNAL(clicked(bool)), q, SLOT(slot_button_clicked()));
    return B;
}

QAbstractButton* GuideV1Private::make_open_view_button(QString text, QString view_id, QString container_name)
{
    QPushButton* B = new QPushButton(text);
    B->setProperty("action", "open_view");
    B->setProperty("view-id", view_id);
    B->setProperty("container-name", container_name);
    QObject::connect(B, SIGNAL(clicked(bool)), q, SLOT(slot_button_clicked()));
    return B;
}

void GuideV1Private::show_instructions(QString title, QString instructions)
{
    QMessageBox::information(q, title, instructions);
}

class compute_cluster_channel_matrix_thread : public QThread {
public:
    //input
    DiskReadMda32 timeseries;
    DiskReadMda firings;
    QString mlproxy_url;
    int clip_size;

    //output
    DiskReadMda templates_out;
    DiskReadMda cluster_channel_matrix;

    void run();
};

void GuideV1Private::compute_cluster_channel_matrix()
{
    compute_cluster_channel_matrix_thread* T = new compute_cluster_channel_matrix_thread;
    QObject::connect(T, SIGNAL(finished()), q, SLOT(slot_cluster_channel_matrix_computed()));
    QObject::connect(T, SIGNAL(finished()), T, SLOT(deleteLater()));
    T->firings = m_context->firings();
    T->timeseries = m_context->currentTimeseries();
    T->clip_size = m_context->option("clip_size").toInt();
    T->mlproxy_url = m_context->mlProxyUrl();

    T->start();
}

void GuideV1Private::update_buttons()
{
    int ch = m_channel_edit->text().toInt();
    m_prev_channel_button->setEnabled(ch > 1);
    m_next_channel_button->setEnabled((ch == 0) || (ch + 1 <= m_cluster_channel_matrix.N2()));
}

void GuideV1::slot_cluster_channel_matrix_computed()
{
    compute_cluster_channel_matrix_thread* T = (compute_cluster_channel_matrix_thread*)qobject_cast<QThread*>(sender());
    if (!T)
        return;

    d->m_cluster_channel_matrix_computed = true;
    d->m_cluster_channel_matrix = T->cluster_channel_matrix;
    this->slot_next_channel();
}

void compute_cluster_channel_matrix_thread::run()
{
    TaskProgress task(TaskProgress::Calculate, "mp_compute_templates_stdevs");
    task.log("mlproxy_url: " + mlproxy_url);
    MountainProcessRunner X;
    QString processor_name = "mv_compute_templates";
    X.setProcessorName(processor_name);

    QMap<QString, QVariant> params;
    params["timeseries"] = timeseries.makePath();
    params["firings"] = firings.makePath();
    params["clip_size"] = clip_size;
    X.setInputParameters(params);
    X.setMLProxyUrl(mlproxy_url);

    QString templates_fname = X.makeOutputFilePath("templates");
    QString stdevs_fname = X.makeOutputFilePath("stdevs");

    X.runProcess();
    templates_out.setPath(templates_fname);
    //templates_out.setRemoteDataType("float32_q8");

    int M = templates_out.N1();
    int T = templates_out.N2();
    int K = templates_out.N3();
    Mda XX(K, M);
    for (int k = 0; k < K; k++) {
        Mda W;
        templates_out.readChunk(W, 0, 0, k, M, T, 1);
        for (int m = 0; m < M; m++) {
            double sumsqr = 0;
            for (int t = 0; t < T; t++) {
                double tmp = W.value(m, t);
                sumsqr += tmp * tmp;
            }
            XX.setValue(sqrt(sumsqr), k, m);
        }
    }
    cluster_channel_matrix = XX;
}
