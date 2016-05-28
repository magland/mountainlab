#include "mvoverview2widget.h"
#include "diskreadmda.h"
#include "sstimeseriesview.h"
#include "sstimeserieswidget.h"
#include "mvcrosscorrelogramswidget2.h"
#include "mvoverview2widgetcontrolpanel.h"
#ifdef USE_LAPACK
#include "get_pca_features.h"
#else
#include "get_principal_components.h"
#endif
#include "get_sort_indices.h"
#include "mvclusterdetailwidget.h"
#include "mvclipsview.h"
#include "mvclusterwidget.h"
#include "mvfiringeventview.h"
#include "extract_clips.h"
#include "tabber.h"
#include "computationthread.h"
#include "mountainprocessrunner.h"
#include "mvclipswidget.h"
#include "taskprogressview.h"
#include "mvcontrolpanel.h"
#include "taskprogress.h"
#include "clustermerge.h"
#include "mvviewagent.h"
#include "mvstatusbar.h"

#include <QHBoxLayout>
#include <QMessageBox>
#include <QSplitter>
#include <QTime>
#include <QTimer>
#include <math.h>
#include <QProgressDialog>
#include "msmisc.h"
#include "mvutils.h"
#include <QColor>
#include <QStringList>
#include <QSet>
#include <QKeyEvent>
#include <QFileDialog>
#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonArray>
#include <QInputDialog>
#include <QAbstractButton>
#include "textfile.h"

/// TODO put styles in central place?
#define MV_STATUS_BAR_HEIGHT 30

class MVOverview2WidgetPrivate {
public:
    MVOverview2Widget* q;
    QMap<QString, QString> m_timeseries_paths;
    QString m_current_timeseries_name;
    DiskReadMda m_timeseries;
    DiskReadMda m_firings_original;
    Mda m_firings_split;
    DiskReadMda m_firings;
    QList<Epoch> m_epochs;
    QList<int> m_original_cluster_numbers;
    QList<int> m_original_cluster_offsets;
    int m_current_k;
    QSet<int> m_selected_ks;
    float m_samplerate;
    MVEvent m_current_event;
    //QString m_mscmdserver_url;
    QString m_mlproxy_url;
    QString m_mv_fname;

    MVControlPanel* m_control_panel_new;
    TaskProgressView* m_task_progress_view;

    QSplitter* m_splitter1, *m_splitter2, *m_left_splitter;
    TabberTabWidget* m_tabs1, *m_tabs2;
    Tabber* m_tabber;
    QProgressDialog* m_progress_dialog;

    Mda m_cross_correlograms_data;
    //Mda m_templates_data;

    MVViewAgent m_view_agent;

    QList<QColor> m_channel_colors;
    QMap<QString, QColor> m_colors;

    void create_cross_correlograms_data();
    //void create_templates_data();

    void update_sizes();
    //void update_templates();
    void update_all_widgets();
    void update_cluster_details();
    void update_clips();
    void update_cluster_views();
    void update_firing_event_views();
    void do_shell_split_and_event_filter();
    void add_tab(QWidget* W, QString label);

    MVCrossCorrelogramsWidget2* open_auto_correlograms();
    MVCrossCorrelogramsWidget2* open_cross_correlograms(int k);
    MVCrossCorrelogramsWidget2* open_matrix_of_cross_correlograms();
    //void open_templates();
    MVClusterDetailWidget* open_cluster_details();
    void open_timeseries();
    void open_clips();
    void open_clusters();
    void open_firing_events();
    void find_nearby_events();

    void annotate_selected();
    void merge_selected();
    void unmerge_selected();

    void update_cross_correlograms();
    void update_timeseries_views();
    void move_to_timepoint(double tp);
    void update_widget(QWidget* W);

    void set_cross_correlograms_current_index(int index);
    void set_cross_correlograms_selected_indices(const QList<int>& indices);
    void set_templates_current_number(int kk);
    void set_templates_selected_numbers(const QList<int>& kks);

    void set_times_labels_for_timeseries_widget(SSTimeSeriesWidget* WW);

    QList<QWidget*> get_all_widgets();
    TabberTabWidget* tab_widget_of(QWidget* W);

    void remove_widgets_of_type(QString widget_type);

    Mda compute_centroid(Mda& clips);
    Mda compute_geometric_median(Mda& clips, int num_iterations);
    void compute_geometric_median(int M, int N, double* output, double* input, int num_iterations);

    void set_progress(QString title, QString text, float frac);
    void set_current_event(MVEvent evt);

    long cc_max_dt_timepoints();

    void export_mountainview_document();
    void export_original_firings();
    void export_filtered_firings();
    void export_file(QString source_path, QString dest_path, bool use_float64);

    QString make_absolute_path(QString path); //use basepath of m_mv_fname if path is relative
    QString current_timeseries_path();

    QVariant get_cluster_attribute(int k, QString attr);
    void set_cluster_attribute(int k, QString attr, QVariant val);
    void set_button_enabled(QString name, bool val);

    //void start_cross_correlograms_computer();
};

QColor brighten(QColor col, int amount)
{
    int r = col.red() + amount;
    int g = col.green() + amount;
    int b = col.blue() + amount;
    if (r > 255)
        r = 255;
    if (r < 0)
        r = 0;
    if (g > 255)
        g = 255;
    if (g < 0)
        g = 0;
    if (b > 255)
        b = 255;
    if (b < 0)
        b = 0;
    return QColor(r, g, b, col.alpha());
}

MVOverview2Widget::MVOverview2Widget(QWidget* parent)
    : QWidget(parent)
{
    d = new MVOverview2WidgetPrivate;
    d->q = this;

    d->m_current_k = 0;
    d->m_samplerate = 20000;

    d->m_progress_dialog = 0;
    d->m_current_event.time = -1;
    d->m_current_event.label = -1;

    d->m_control_panel_new = new MVControlPanel;
    connect(d->m_control_panel_new, SIGNAL(userAction(QString)), this, SLOT(slot_control_panel_user_action(QString)));

    QSplitter* splitter1 = new QSplitter;
    splitter1->setOrientation(Qt::Horizontal);
    d->m_splitter1 = splitter1;

    QSplitter* splitter2 = new QSplitter;
    splitter2->setOrientation(Qt::Vertical);
    d->m_splitter2 = splitter2;

    QScrollArea* CP = new QScrollArea;
    CP->setWidget(d->m_control_panel_new);
    CP->setWidgetResizable(true);

    d->m_task_progress_view = new TaskProgressView;

    QSplitter* left_splitter = new QSplitter(Qt::Vertical);
    left_splitter->addWidget(CP);
    left_splitter->addWidget(d->m_task_progress_view);
    d->m_left_splitter = left_splitter;

    splitter1->addWidget(left_splitter);
    splitter1->addWidget(splitter2);

    d->m_tabber = new Tabber;
    d->m_tabs2 = d->m_tabber->createTabWidget("south");
    d->m_tabs1 = d->m_tabber->createTabWidget("north");

    splitter2->addWidget(d->m_tabs1);
    splitter2->addWidget(d->m_tabs2);

    MVStatusBar* status_bar = new MVStatusBar();
    status_bar->setFixedHeight(MV_STATUS_BAR_HEIGHT);

    QVBoxLayout* vlayout = new QVBoxLayout;
    vlayout->setSpacing(0);
    vlayout->setMargin(0);
    vlayout->addWidget(splitter1);
    vlayout->addWidget(status_bar);
    this->setLayout(vlayout);

    QStringList color_strings;
    color_strings
        << "#282828"
        << "#402020"
        << "#204020"
        << "#202070";
    for (int i = 0; i < color_strings.count(); i++)
        d->m_channel_colors << QColor(brighten(color_strings[i], 80));

    d->m_colors["background"] = QColor(240, 240, 240);
    d->m_colors["frame1"] = QColor(245, 245, 245);
    d->m_colors["info_text"] = QColor(80, 80, 80);
    d->m_colors["view_background"] = QColor(245, 245, 245);
    d->m_colors["view_background_highlighted"] = QColor(210, 230, 250);
    d->m_colors["view_background_selected"] = QColor(220, 240, 250);
    d->m_colors["view_background_hovered"] = QColor(240, 245, 240);

    d->m_colors["view_frame_selected"] = QColor(50, 20, 20);
    d->m_colors["divider_line"] = QColor(255, 100, 150);

    slot_update_buttons();

    //connect(&d->m_cross_correlogram_computer,SIGNAL(computationFinished()),this,SLOT(slot_cross_correlogram_computer_finished()));
}

MVOverview2Widget::~MVOverview2Widget()
{
    delete d;
}

void MVOverview2Widget::addTimeseriesPath(const QString& name, const QString& path)
{
    d->m_timeseries_paths[name] = path;
    QStringList choices = d->m_timeseries_paths.keys();
    qSort(choices);
    d->m_control_panel_new->setTimeseriesChoices(choices);
    if (d->m_timeseries_paths.count() == 1) {
        this->setCurrentTimeseriesName(name);
    }
}

void MVOverview2Widget::setCurrentTimeseriesName(const QString& name)
{
    {
        MVViewOptions opts = d->m_control_panel_new->viewOptions();
        opts.timeseries = name;
        d->m_control_panel_new->setViewOptions(opts);
    }

    d->update_timeseries_views();
    d->update_cluster_details();
    d->update_clips();
    d->update_cluster_views();
    d->update_firing_event_views();
}

void MVOverview2Widget::setFiringsPath(const QString& firings)
{
    d->m_firings_original.setPath(d->make_absolute_path(firings));
    d->do_shell_split_and_event_filter();
    d->update_cross_correlograms();
    d->update_cluster_details();
    d->update_timeseries_views();
    //d->start_cross_correlograms_computer();
    slot_update_buttons();
}

void MVOverview2Widget::setSampleRate(float freq)
{
    d->m_samplerate = freq;
    //d->start_cross_correlograms_computer();
}

void MVOverview2Widget::setDefaultInitialization()
{
    //d->open_templates();
    d->open_cluster_details();
    d->m_tabber->switchCurrentContainer();
    d->open_auto_correlograms();
}

void MVOverview2Widget::setEpochs(const QList<Epoch>& epochs)
{
    d->m_epochs = epochs;
}

/*
QImage MVOverview2Widget::generateImage(const QMap<QString, QVariant>& params)
{
    QString type0 = params.value("type").toString();
    if (type0 == "templates") {
        MVClusterDetailWidget* X = d->open_cluster_details();
        return X->renderImage();
    }
    else if (type0 == "auto_correlograms") {
        MVCrossCorrelogramsWidget2* X = d->open_auto_correlograms();
        return X->renderImage();
    }
    else if (type0 == "cross_correlograms") {
        int k = params.value("k", 0).toInt();
        MVCrossCorrelogramsWidget2* X = d->open_cross_correlograms(k);
        return X->renderImage();
    }
    else {
        qWarning() << "Unknown type in generateImage: " << type0;
        return QImage();
    }
}
*/

int MVOverview2Widget::getMaxLabel()
{
    int ret = 0;
    for (long i = 0; i < d->m_firings.N2(); i++) {
        int label = (int)d->m_firings.value(i);
        if (label > ret)
            ret = label;
    }
    return ret;
}

void MVOverview2Widget::setMLProxyUrl(const QString& url)
{
    d->m_mlproxy_url = url;
}

void MVOverview2Widget::setClusterMerge(ClusterMerge CM)
{
    d->m_view_agent.setClusterMerge(CM);
}

/*
void MVOverview2Widget::setMscmdServerUrl(const QString& url)
{
    d->m_mscmdserver_url = url;
}
*/

void MVOverview2Widget::loadMVFile(const QString& mv_fname)
{
    d->m_mv_fname = mv_fname;

    TaskProgress task("loading .mv file: " + mv_fname);
    QString json = read_text_file(mv_fname);
    QJsonParseError err;
    QJsonObject obj = QJsonDocument::fromJson(json.toLatin1(), &err).object();
    if (err.error != QJsonParseError::NoError) {
        task.error("Error parsing .mv file: " + mv_fname);
        task.log("JSON parse error: " + err.errorString());
        return;
    }

    //important to do this first
    if (obj.contains("mlproxy_url")) {
        QString url = obj["mlproxy_url"].toString();
        task.log(QString("mlproxy_url = %1").arg(url));
        this->setMLProxyUrl(url);
    }

    QString mv_version = obj["mv_version"].toString();
    task.log("MV version: " + mv_version);

    if (obj.contains("firings")) {
        QString path = obj["firings"].toString();
        task.log("Setting firings path: " + path);
        this->setFiringsPath(path);
    }
    if (obj.contains("timeseries")) {
        QJsonArray ts = obj["timeseries"].toArray();
        for (int i = 0; i < ts.count(); i++) {
            QJsonObject tsobj = ts[i].toObject();
            QString name = tsobj["name"].toString();
            QString path = tsobj["path"].toString();
            task.log("Adding timeseries " + name + ": " + path);
            this->addTimeseriesPath(name, path);
        }
    }

    if (obj.contains("samplerate")) {
        double rate = obj["samplerate"].toDouble();
        task.log(QString("samplerate = %1").arg(rate));
        this->setSampleRate(rate);
    }

    if (obj.contains("view_options")) {
        QJsonObject obj0 = obj["view_options"].toObject();
        task.log("VIEW OPTIONS:");
        task.log(QJsonDocument(obj0).toJson(QJsonDocument::Compact));
        MVViewOptions opts = MVViewOptions::fromJsonObject(obj0);
        d->m_control_panel_new->setViewOptions(opts);
    }

    if (obj.contains("event_filter")) {
        QJsonObject obj0 = obj["event_filter"].toObject();
        task.log("EVENT FILTER:");
        task.log(QJsonDocument(obj0).toJson(QJsonDocument::Compact));
        MVEventFilter filter = MVEventFilter::fromJsonObject(obj0);
        d->m_control_panel_new->setEventFilter(filter);
    }

    if (obj.contains("annotations")) {
        if (obj["annotations"].toObject().contains("cluster_attributes")) {
            QJsonObject obj2 = obj["annotations"].toObject()["cluster_attributes"].toObject();
            QMap<int, QJsonObject> CA;
            QStringList keys = obj2.keys();
            foreach(QString key, keys)
            {
                bool ok;
                int num = key.toInt(&ok);
                if (ok) {
                    CA[num] = obj2[key].toObject();
                }
            }
            d->m_view_agent.setClusterAttributes(CA);
        } else
            d->m_view_agent.setClusterAttributes(QMap<int, QJsonObject>());
        if (obj["annotations"].toObject().contains("cluster_merge")) {
            QJsonArray CM = obj["annotations"].toObject()["cluster_merge"].toArray();
            QString json = QJsonDocument(CM).toJson(QJsonDocument::Compact);
            d->m_view_agent.setClusterMerge(ClusterMerge::fromJson(json));
        } else {
            d->m_view_agent.setClusterMerge(ClusterMerge());
        }
    }

    d->do_shell_split_and_event_filter();
}

void MVOverview2Widget::saveMVFile(const QString& mv_fname)
{
    TaskProgress task("saving .mv file: " + mv_fname);

    QJsonObject obj;

    obj["mv_version"] = 0.1;

    obj["firings"] = d->m_firings_original.path();
    QJsonArray ts;
    {
        QStringList keys = d->m_timeseries_paths.keys();
        foreach(QString key, keys)
        {
            QJsonObject tsobj;
            tsobj["name"] = key;
            tsobj["path"] = d->m_timeseries_paths[key];
            ts << tsobj;
        }
        obj["timeseries"] = ts;
    }

    obj["samplerate"] = d->m_samplerate;

    obj["view_options"] = d->m_control_panel_new->viewOptions().toJsonObject();
    obj["event_filter"] = d->m_control_panel_new->eventFilter().toJsonObject();

    obj["mlproxy_url"] = d->m_mlproxy_url;

    QJsonObject cluster_attributes;
    QMap<int, QJsonObject> CA = d->m_view_agent.clusterAttributes();
    {
        QList<int> keys = CA.keys();
        foreach(int key, keys)
        {
            if (!CA[key].isEmpty())
                cluster_attributes[QString("%1").arg(key)] = CA[key];
        }
    }

    QJsonObject annotations;
    annotations["cluster_attributes"] = cluster_attributes;
    annotations["cluster_merge"] = QJsonDocument::fromJson(d->m_view_agent.clusterMerge().toJson().toLatin1()).array();
    obj["annotations"] = annotations;

    if (!write_text_file(mv_fname, QJsonDocument(obj).toJson())) {
        task.error("Error writing .mv file: " + mv_fname);
    }
}

void MVOverview2Widget::resizeEvent(QResizeEvent* evt)
{
    Q_UNUSED(evt)
    d->update_sizes();
}

void MVOverview2Widget::keyPressEvent(QKeyEvent* evt)
{
    if ((evt->key() == Qt::Key_W) && (evt->modifiers() & Qt::ControlModifier)) {
        this->close();
    } else if (evt->key() == Qt::Key_A) {
        d->annotate_selected();
    } else if (evt->key() == Qt::Key_M) {
        d->merge_selected();
    } else if (evt->key() == Qt::Key_U) {
        d->unmerge_selected();
    } else
        evt->ignore();
}

/*
void MVOverview2Widget::slot_control_panel_button_clicked(QString str)
{
    if (str == "update_cross_correlograms") {
        //d->start_cross_correlograms_computer();
        d->update_cross_correlograms();
    }
    else if (str == "update_templates") {
        //	d->update_templates();
        d->update_cluster_details();
        d->update_clips();
    }
    else if (str == "update_cluster_details") {
        d->update_cluster_details();
    }
    else if ((str == "update_shell_split") || (str == "use_shell_split")) {
        d->do_shell_split_and_event_filter();
        d->remove_widgets_of_type("cross_correlograms");
        d->remove_widgets_of_type("matrix_of_cross_correlograms");
        d->remove_widgets_of_type("clips");
        d->remove_widgets_of_type("find_nearby_events");
        d->remove_widgets_of_type("clusters");
        d->remove_widgets_of_type("firing_events");
        d->update_cluster_details();
        d->update_cross_correlograms();
    }
    else if ((str == "update_event_filter") || (str == "use_event_filter")) {
        d->do_shell_split_and_event_filter();
        //d->do_event_filter();
        //d->start_cross_correlograms_computer();
        d->update_all_widgets();
    }
    else if (str == "open_auto_correlograms") {
        d->open_auto_correlograms();
    }
    else if (str == "open_matrix_of_cross_correlograms") {
        d->open_matrix_of_cross_correlograms();
    }
    //else if (str=="open_templates") {
    //    d->open_templates();
    //}
    else if (str == "open_cluster_details") {
        d->open_cluster_details();
    }
    else if (str == "open_timeseries") {
        d->open_timeseries();
    }
    else if (str == "open_clips") {
        d->open_clips();
    }
    else if (str == "open_clusters") {
        d->open_clusters();
    }
    else if (str == "open_firing_events") {
        d->open_firing_events();
    }
    else if (str == "find_nearby_events") {
        d->find_nearby_events();
    }
    else if (str == "template_method") {
        d->update_cluster_details();
    }
}
*/

void MVOverview2Widget::slot_control_panel_user_action(QString str)
{
    if ((str == "apply_shell_splitting") || (str == "apply_filter")) {
        d->do_shell_split_and_event_filter();
        d->update_all_widgets();
    } else if (str == "update_all_open_views") {
        d->update_all_widgets();
    } else if (str == "open-cluster-details") {
        d->open_cluster_details();
    } else if (str == "open-auto-correlograms") {
        d->open_auto_correlograms();
    } else if (str == "open-matrix-of-cross-correlograms") {
        d->open_matrix_of_cross_correlograms();
    } else if (str == "open-timeseries-data") {
        d->open_timeseries();
    } else if (str == "open-clips") {
        d->open_clips();
    } else if (str == "open-clusters") {
        d->open_clusters();
    } else if (str == "open-firing-events") {
        d->open_firing_events();
    } else if (str == "find-nearby-events") {
        d->find_nearby_events();
    } else if (str == "annotate_selected") {
        d->annotate_selected();
    } else if (str == "merge_selected") {
        d->merge_selected();
    } else if (str == "unmerge_selected") {
        d->unmerge_selected();
    } else if (str == "export_mountainview_document") {
        d->export_mountainview_document();
    } else if (str == "export_original_firings") {
        d->export_original_firings();
    } else if (str == "export_filtered_firings") {
        d->export_filtered_firings();
    } else {
        TaskProgress task(str);
        task.error("user action not yet implemented.");
    }
}

void MVOverview2Widget::slot_auto_correlogram_activated(int index)
{
    TabberTabWidget* TW = d->tab_widget_of((QWidget*)sender());
    d->m_tabber->setCurrentContainer(TW);
    d->m_tabber->switchCurrentContainer();
    d->open_cross_correlograms(index + 1);
}

void MVOverview2Widget::slot_details_current_k_changed()
{
    MVClusterDetailWidget* X = (MVClusterDetailWidget*)sender();
    int index = X->currentK() - 1;
    d->m_current_k = X->currentK();
    d->set_cross_correlograms_current_index(index);
}

void MVOverview2Widget::slot_details_selected_ks_changed()
{
    MVClusterDetailWidget* X = (MVClusterDetailWidget*)sender();
    QList<int> ks = X->selectedKs();
    QList<int> indices;
    foreach(int k, ks)
    indices << k - 1;
    d->m_selected_ks = ks.toSet();
    d->set_cross_correlograms_selected_indices(indices);
    slot_update_buttons();
}

void MVOverview2Widget::slot_details_template_activated()
{
    MVClusterDetailWidget* X = (MVClusterDetailWidget*)sender();
    int k = X->currentK();
    if (k < 0)
        return;
    TabberTabWidget* TW = d->tab_widget_of((QWidget*)sender());
    d->m_tabber->setCurrentContainer(TW);
    d->m_tabber->switchCurrentContainer();
    d->open_clips();
}

void MVOverview2Widget::slot_cross_correlogram_current_index_changed()
{
    MVCrossCorrelogramsWidget2* X = (MVCrossCorrelogramsWidget2*)sender();
    d->m_current_k = X->currentLabel1();
    d->set_cross_correlograms_current_index(X->currentIndex());
    d->set_templates_current_number(X->currentLabel2());
}

void MVOverview2Widget::slot_cross_correlogram_selected_indices_changed()
{
    MVCrossCorrelogramsWidget2* X = (MVCrossCorrelogramsWidget2*)sender();
    d->m_selected_ks = X->selectedLabels1().toSet();
    d->set_cross_correlograms_selected_indices(X->selectedIndices());
    d->set_templates_selected_numbers(X->selectedLabels2());
}

void MVOverview2Widget::slot_clips_view_current_event_changed()
{
    MVClipsView* V = (MVClipsView*)sender();
    MVEvent evt = V->currentEvent();
    d->set_current_event(evt);
}

void MVOverview2Widget::slot_clips_widget_current_event_changed()
{
    MVClipsWidget* W = (MVClipsWidget*)sender();
    MVEvent evt = W->currentEvent();
    d->set_current_event(evt);
}

void MVOverview2Widget::slot_cluster_view_current_event_changed()
{
    MVClusterWidget* W = (MVClusterWidget*)sender();
    MVEvent evt = W->currentEvent();
    d->set_current_event(evt);
}

void MVOverview2Widget::slot_update_buttons()
{
    bool has_peaks = (d->m_firings.value(0, 3) != 0); //for now we just test the very first one (might be problematic)
    bool something_selected = (!d->m_selected_ks.isEmpty());

    d->set_button_enabled("open-cluster-details", true);
    d->set_button_enabled("open-auto-correlograms", true);
    d->set_button_enabled("open-matrix-of-cross-correlograms", something_selected);
    d->set_button_enabled("open-timeseries-data", !d->m_timeseries.path().startsWith("http"));
    d->set_button_enabled("open-clips", something_selected);
    d->set_button_enabled("open-clusters", something_selected);
    d->set_button_enabled("open-firing-events", (something_selected) && (has_peaks));
    d->set_button_enabled("find-nearby-events", d->m_selected_ks.count() >= 2);

    d->set_button_enabled("annotate_selected", something_selected);
    d->set_button_enabled("merge_selected", d->m_selected_ks.count() >= 2);
    d->set_button_enabled("unmerge_selected", something_selected);
    d->set_button_enabled("export_mountainview_document", true);
    d->set_button_enabled("export_original_firings", true);
    d->set_button_enabled("export_filtered_firings", true);
}

void MVOverview2WidgetPrivate::update_sizes()
{
    float W0 = q->width();
    float H0 = q->height() - MV_STATUS_BAR_HEIGHT;

    int W1 = W0 / 3;
    if (W1 < 150)
        W1 = 150;
    if (W1 > 400)
        W1 = 400;
    int W2 = W0 - W1;

    int H1 = H0 / 2;
    int H2 = H0 / 2;
    //int H3=H0-H1-H2;

    {
        QList<int> sizes;
        sizes << W1 << W2;
        m_splitter1->setSizes(sizes);
    }
    {
        QList<int> sizes;
        sizes << H1 << H2;
        m_splitter2->setSizes(sizes);
    }
    {
        int tv_height;
        if (H0 > 1200) {
            tv_height = 300;
        }
        if (H0 > 900) {
            tv_height = 200;
        } else {
            tv_height = 100;
        }
        int cp_height = H0 - tv_height;
        QList<int> sizes;
        sizes << cp_height << tv_height;
        m_left_splitter->setSizes(sizes);
    }
}

void MVOverview2WidgetPrivate::update_all_widgets()
{
    QList<QWidget*> list = get_all_widgets();
    foreach(QWidget * W, list)
    {
        update_widget(W);
    }
}

void MVOverview2WidgetPrivate::update_cluster_details()
{
    QList<QWidget*> list = get_all_widgets();
    foreach(QWidget * W, list)
    {
        if (W->property("widget_type") == "cluster_details") {
            update_widget(W);
        }
    }
}

void MVOverview2WidgetPrivate::update_clips()
{
    QList<QWidget*> list = get_all_widgets();
    foreach(QWidget * W, list)
    {
        if (W->property("widget_type") == "clips") {
            update_widget(W);
        }
        if (W->property("widget_type") == "find_nearby_events") {
            update_widget(W);
        }
    }
}

void MVOverview2WidgetPrivate::update_cluster_views()
{
    QList<QWidget*> list = get_all_widgets();
    foreach(QWidget * W, list)
    {
        if (W->property("widget_type") == "clusters") {
            update_widget(W);
        }
    }
}

void MVOverview2WidgetPrivate::update_firing_event_views()
{
    QList<QWidget*> list = get_all_widgets();
    foreach(QWidget * W, list)
    {
        if (W->property("widget_type") == "firing_events") {
            update_widget(W);
        }
    }
}

double get_max(QList<double>& list)
{
    double ret = list.value(0);
    for (int i = 0; i < list.count(); i++) {
        if (list[i] > ret)
            ret = list[i];
    }
    return ret;
}

double get_min(QList<double>& list)
{
    double ret = list.value(0);
    for (int i = 0; i < list.count(); i++) {
        if (list[i] < ret)
            ret = list[i];
    }
    return ret;
}

void MVOverview2WidgetPrivate::do_shell_split_and_event_filter()
{
    TaskProgress task("do_shell_split_and_event_filter()");

    MountainProcessRunner MT;
    QString processor_name = "mv_firings_filter";

    MT.setProcessorName(processor_name);

    QMap<QString, QVariant> params;

    MVEventFilter evt_filter = m_control_panel_new->eventFilter();

    params["use_shell_split"] = evt_filter.use_shell_split;
    params["shell_width"] = evt_filter.shell_increment;
    params["min_per_shell"] = evt_filter.min_per_shell;
    params["use_event_filter"] = evt_filter.use_event_filter;
    //if (evt_filter.min_detectability_score) {
    params["min_detectability_score"] = evt_filter.min_detectability_score;
    //}
    //if (evt_filter.max_outlier_score) {
    params["max_outlier_score"] = evt_filter.max_outlier_score;
    //}
    params["min_amplitude"] = 0;

    params["firings"] = m_firings_original.makePath();

    QStringList debug_keys = params.keys();
    foreach(QString key, debug_keys)
    {
        task.log(QString("%1 = %2").arg(key).arg(params[key].toString()));
    }

    MT.setInputParameters(params);
    //MT.setMscmdServerUrl(m_mscmdserver_url);
    MT.setMLProxyUrl(m_mlproxy_url);

    QString firings_out = MT.makeOutputFilePath("firings_out");
    QString original_cluster_numbers_out = MT.makeOutputFilePath("original_cluster_numbers");
    /// TODO this should be called in a separate thread MT.runProcess
    MT.runProcess(0);
    m_firings.setPath(firings_out);
    task.log("m_firings path = " + firings_out);
    task.log(QString("m_firings.N1()=%1 m_firings.N2()=%2").arg(m_firings.N1()).arg(m_firings.N2()));
    m_original_cluster_numbers.clear();
    m_original_cluster_offsets.clear();
    DiskReadMda AA(original_cluster_numbers_out);
    int offset = 0;
    for (int i = 0; i < AA.totalSize(); i++) {
        if (AA.value(i) != AA.value(i - 1)) {
            offset = 0;
        }
        offset++;
        m_original_cluster_numbers << AA.value(i);
        m_original_cluster_offsets << offset;
    }

    this->set_templates_current_number(-1);
    this->set_templates_selected_numbers(QList<int>());
}

void MVOverview2WidgetPrivate::add_tab(QWidget* W, QString label)
{
    W->setFocusPolicy(Qt::StrongFocus);
    //current_tab_widget()->addTab(W,label);
    //current_tab_widget()->setCurrentIndex(current_tab_widget()->count()-1);
    m_tabber->addWidget(m_tabber->currentContainerName(), label, W);
    W->setProperty("tab_label", label); //won't be needed in future, once Tabber is fully implemented
}

MVCrossCorrelogramsWidget2* MVOverview2WidgetPrivate::open_auto_correlograms()
{
    MVCrossCorrelogramsWidget2* X = new MVCrossCorrelogramsWidget2;
    X->setViewAgent(&m_view_agent);
    X->setProperty("widget_type", "auto_correlograms");
    add_tab(X, "Auto-Correlograms");
    QObject::connect(X, SIGNAL(indexActivated(int)), q, SLOT(slot_auto_correlogram_activated(int)));
    QObject::connect(X, SIGNAL(currentIndexChanged()), q, SLOT(slot_cross_correlogram_current_index_changed()));
    QObject::connect(X, SIGNAL(selectedIndicesChanged()), q, SLOT(slot_cross_correlogram_selected_indices_changed()));
    update_widget(X);
    return X;
}

MVCrossCorrelogramsWidget2* MVOverview2WidgetPrivate::open_cross_correlograms(int k)
{
    MVCrossCorrelogramsWidget2* X = new MVCrossCorrelogramsWidget2;
    X->setViewAgent(&m_view_agent);
    X->setProperty("widget_type", "cross_correlograms");
    X->setProperty("kk", k);
    add_tab(X, QString("CC for %1(%2)").arg(m_original_cluster_numbers.value(k)).arg(m_original_cluster_offsets.value(k)));
    QObject::connect(X, SIGNAL(currentIndexChanged()), q, SLOT(slot_cross_correlogram_current_index_changed()));
    QObject::connect(X, SIGNAL(selectedIndicesChanged()), q, SLOT(slot_cross_correlogram_selected_indices_changed()));
    update_widget(X);
    return X;
}

QStringList int_list_to_string_list(const QList<int>& list)
{
    QStringList list2;
    for (int i = 0; i < list.count(); i++)
        list2 << QString("%1").arg(list[i]);
    return list2;
}

QList<int> string_list_to_int_list(const QList<QString>& list)
{
    QList<int> list2;
    for (int i = 0; i < list.count(); i++)
        list2 << list[i].toInt();
    return list2;
}

MVCrossCorrelogramsWidget2* MVOverview2WidgetPrivate::open_matrix_of_cross_correlograms()
{
    MVCrossCorrelogramsWidget2* X = new MVCrossCorrelogramsWidget2;
    X->setViewAgent(&m_view_agent);
    X->setProperty("widget_type", "matrix_of_cross_correlograms");
    QList<int> ks = m_selected_ks.toList();
    qSort(ks);
    if (ks.isEmpty())
        return X;
    X->setProperty("ks", int_list_to_string_list(ks));
    add_tab(X, QString("CC Matrix"));
    update_widget(X);
    return X;
}

//void MVOverview2WidgetPrivate::open_templates()
//{
//	SSTimeSeriesView *X=new SSTimeSeriesView;
//	X->initialize();
//	X->setProperty("widget_type","templates");
//	add_tab(X,QString("Templates"));
//    QObject::connect(X,SIGNAL(currentXChanged()),q,SLOT(slot_templates_clicked()));
//	update_widget(X);
//}

MVClusterDetailWidget* MVOverview2WidgetPrivate::open_cluster_details()
{
    MVClusterDetailWidget* X = new MVClusterDetailWidget;
    //X->setMscmdServerUrl(m_mscmdserver_url);
    X->setViewAgent(&m_view_agent);
    X->setMLProxyUrl(m_mlproxy_url);
    X->setChannelColors(m_channel_colors);
    DiskReadMda TT(current_timeseries_path());
    X->setTimeseries(TT);
    //X->setFirings(DiskReadMda(m_firings)); //done in update_widget
    X->setSampleRate(m_samplerate);
    QObject::connect(X, SIGNAL(signalTemplateActivated()), q, SLOT(slot_details_template_activated()));
    QObject::connect(X, SIGNAL(signalCurrentKChanged()), q, SLOT(slot_details_current_k_changed()));
    QObject::connect(X, SIGNAL(signalSelectedKsChanged()), q, SLOT(slot_details_selected_ks_changed()));
    X->setProperty("widget_type", "cluster_details");
    add_tab(X, QString("Details"));
    update_widget(X);
    return X;
}

void MVOverview2WidgetPrivate::open_timeseries()
{
    SSTimeSeriesWidget* X = new SSTimeSeriesWidget;
    SSTimeSeriesView* V = new SSTimeSeriesView;
    V->initialize();
    V->setSampleRate(m_samplerate);
    X->addView(V);
    X->setProperty("widget_type", "timeseries");
    add_tab(X, QString("Timeseries"));
    update_widget(X);
}

void MVOverview2WidgetPrivate::open_clips()
{
    QList<int> ks = m_selected_ks.toList();
    qSort(ks);
    if (ks.count() == 0) {
        QMessageBox::information(q, "Unable to open clips", "You must select at least one cluster.");
        return;
    }

    MVClipsWidget* X = new MVClipsWidget;
    //X->setMscmdServerUrl(m_mscmdserver_url);
    X->setMLProxyUrl(m_mlproxy_url);
    X->setProperty("widget_type", "clips");
    X->setProperty("ks", int_list_to_string_list(ks));
    q->connect(X, SIGNAL(currentEventChanged()), q, SLOT(slot_clips_widget_current_event_changed()));
    QString tab_title = "Clips";
    if (ks.count() == 1) {
        int kk = ks[0];
        tab_title = QString("Clips %1(%2)").arg(m_original_cluster_numbers.value(kk)).arg(m_original_cluster_offsets.value(kk));
    }

    add_tab(X, tab_title);
    update_widget(X);

    /*
    MVClipsView* X = MVClipsView::newInstance();
    X->setProperty("widget_type", "clips");
    X->setProperty("ks", int_list_to_string_list(ks));
    q->connect(X, SIGNAL(currentEventChanged()), q, SLOT(slot_clips_view_current_event_changed()));
    QString tab_title = "Clips";
    if (ks.count() == 1) {
        int kk = ks[0];
        tab_title = QString("Clips %1(%2)").arg(m_original_cluster_numbers.value(kk)).arg(m_original_cluster_offsets.value(kk));
    }
    add_tab(X, tab_title);
    update_widget(X);
    X->setXRange(vec2(0, 5000));
    */
}

void MVOverview2WidgetPrivate::open_clusters()
{
    QList<int> ks = m_selected_ks.toList();
    qSort(ks);
    if (ks.count() == 0) {
        QMessageBox::information(q, "Unable to open clusters", "You must select at least one cluster.");
        return;
    }
    MVClusterWidget* X = new MVClusterWidget;
    //X->setMscmdServerUrl(m_mscmdserver_url);
    X->setMLProxyUrl(m_mlproxy_url);
    X->setProperty("widget_type", "clusters");
    X->setProperty("ks", int_list_to_string_list(ks));
    q->connect(X, SIGNAL(currentEventChanged()), q, SLOT(slot_cluster_view_current_event_changed()));
    add_tab(X, QString("Clusters"));
    update_widget(X);
}

void MVOverview2WidgetPrivate::open_firing_events()
{
    QList<int> ks = m_selected_ks.toList();
    qSort(ks);
    if (ks.count() == 0) {
        QMessageBox::information(q, "Unable to open firing events", "You must select at least one cluster.");
        return;
    }
    MVFiringEventView* X = new MVFiringEventView;
    X->setProperty("widget_type", "firing_events");
    X->setProperty("ks", int_list_to_string_list(ks));
    //q->connect(X,SIGNAL(currentEventChanged()),q,SLOT(slot_firing_events_view_current_event_changed()));
    add_tab(X, QString("Firing Events"));
    update_widget(X);
}

void MVOverview2WidgetPrivate::find_nearby_events()
{
    QList<int> ks = m_selected_ks.toList();
    qSort(ks);
    if (ks.count() < 2) {
        QMessageBox::information(q, "Problem finding nearby events", "You must select at least two clusters.");
        return;
    }

    MVClipsView* X = MVClipsView::newInstance();
    X->setProperty("widget_type", "find_nearby_events");
    X->setProperty("ks", int_list_to_string_list(ks));
    q->connect(X, SIGNAL(currentEventChanged()), q, SLOT(slot_clips_view_current_event_changed()));
    QString tab_title = "Nearby Clips";
    if (ks.count() == 2) {
        int kk1 = ks[0];
        int kk2 = ks[1];
        tab_title = QString("Clips %1(%2) : Clips %3(%4)")
                        .arg(m_original_cluster_numbers.value(kk1))
                        .arg(m_original_cluster_offsets.value(kk1))
                        .arg(m_original_cluster_numbers.value(kk2))
                        .arg(m_original_cluster_offsets.value(kk2));
    }
    add_tab(X, tab_title);
    update_widget(X);
    X->setXRange(vec2(0, 5000));
}

void MVOverview2WidgetPrivate::annotate_selected()
{
    if (m_selected_ks.isEmpty())
        return;
    QList<int> ks = m_selected_ks.toList();
    QString common_assessment;
    for (int i = 0; i < ks.count(); i++) {
        QString aa = get_cluster_attribute(ks[i], "assessment").toString();
        if (i == 0)
            common_assessment = aa;
        else {
            if (common_assessment != aa)
                common_assessment = "";
        }
    }
    bool ok;
    QString new_assessment = QInputDialog::getText(0, "Set cluster assessment", "Cluster assessment:", QLineEdit::Normal, common_assessment, &ok);
    if (ok) {
        for (int i = 0; i < ks.count(); i++) {
            set_cluster_attribute(ks[i], "assessment", new_assessment);
        }
    }
}

void MVOverview2WidgetPrivate::merge_selected()
{
    ClusterMerge CM = m_view_agent.clusterMerge();
    CM.merge(m_selected_ks);
    m_view_agent.setClusterMerge(CM);
}

void MVOverview2WidgetPrivate::unmerge_selected()
{
    ClusterMerge CM = m_view_agent.clusterMerge();
    CM.unmerge(m_selected_ks);
    m_view_agent.setClusterMerge(CM);
}

void MVOverview2WidgetPrivate::update_cross_correlograms()
{
    QList<QWidget*> widgets = get_all_widgets();
    foreach(QWidget * W, widgets)
    {
        QString widget_type = W->property("widget_type").toString();
        if ((widget_type == "auto_correlograms") || (widget_type == "cross_correlograms")) {
            update_widget(W);
        }
    }
}

void MVOverview2WidgetPrivate::update_timeseries_views()
{
    QList<QWidget*> widgets = get_all_widgets();
    foreach(QWidget * W, widgets)
    {
        QString widget_type = W->property("widget_type").toString();
        if (widget_type == "timeseries") {
            update_widget(W);
        }
    }
}

void MVOverview2WidgetPrivate::move_to_timepoint(double tp)
{
    QList<QWidget*> widgets = get_all_widgets();
    foreach(QWidget * W, widgets)
    {
        QString widget_type = W->property("widget_type").toString();
        if (widget_type == "timeseries") {
            SSTimeSeriesWidget* V = (SSTimeSeriesWidget*)W;
            V->view(0)->setCurrentTimepoint(tp);
        }
    }
}

void subtract_features_mean(Mda& F)
{
    if (F.N2() == 0)
        return;
    double mean[F.N1()];
    for (int i = 0; i < F.N1(); i++)
        mean[i] = 0;
    for (int i = 0; i < F.N2(); i++) {
        for (int j = 0; j < F.N1(); j++) {
            mean[j] += F.value(j, i);
        }
    }
    for (int i = 0; i < F.N1(); i++)
        mean[i] /= F.N2();
    for (int i = 0; i < F.N2(); i++) {
        for (int j = 0; j < F.N1(); j++) {
            F.setValue(F.value(j, i) - mean[j], j, i);
        }
    }
}

void normalize_features(Mda& F, bool individual_channels)
{
    if (individual_channels) {
        for (int m = 0; m < F.N1(); m++) {
            double sumsqr = 0;
            double sum = 0;
            int NN = F.N2();
            for (int i = 0; i < NN; i++) {
                sumsqr += F.value(m, i) * F.value(m, i);
                sum += F.value(m, i);
            }
            double norm = 1;
            if (NN >= 2)
                norm = sqrt((sumsqr - sum * sum / NN) / (NN - 1));
            for (int i = 0; i < NN; i++)
                F.setValue(F.value(m, i) / norm, m, i);
        }
    } else {
        double sumsqr = 0;
        double sum = 0;
        int NN = F.totalSize();
        for (int i = 0; i < NN; i++) {
            sumsqr += F.get(i) * F.get(i);
            sum += F.get(i);
        }
        double norm = 1;
        if (NN >= 2)
            norm = sqrt((sumsqr - sum * sum / NN) / (NN - 1));
        for (int i = 0; i < NN; i++)
            F.set(F.get(i) / norm, i);
    }
}

void MVOverview2WidgetPrivate::update_widget(QWidget* W)
{
    QString widget_type = W->property("widget_type").toString();
    if (widget_type == "auto_correlograms") {
        TaskProgress task("update auto correlograms");
        MVCrossCorrelogramsWidget2* WW = (MVCrossCorrelogramsWidget2*)W;
        WW->setSampleRate(m_samplerate);
        WW->setMaxDtTimepoints(cc_max_dt_timepoints());
        WW->setColors(m_colors);
        WW->setFirings(m_firings);
        task.log(QString("m_samplerate = %1").arg(m_samplerate));
        //WW->setCrossCorrelogramsData(DiskReadMda(m_cross_correlograms_data));
        QStringList text_labels;
        QList<int> labels1, labels2;
        for (int i = 1; i < m_original_cluster_numbers.count(); i++) {
            labels1 << i;
            labels2 << i;
            if ((i == 1) || (m_original_cluster_numbers[i] != m_original_cluster_numbers[i - 1])) {
                text_labels << QString("Auto %1").arg(m_original_cluster_numbers[i]);
            } else
                text_labels << "";
        }
        //WW->setTextLabels(labels);
        WW->setLabelPairs(labels1, labels2, text_labels);
        //WW->updateWidget();
    } else if (widget_type == "cross_correlograms") {
        MVCrossCorrelogramsWidget2* WW = (MVCrossCorrelogramsWidget2*)W;
        int k = W->property("kk").toInt();
        WW->setColors(m_colors);
        WW->setSampleRate(m_samplerate);
        WW->setMaxDtTimepoints(cc_max_dt_timepoints());
        WW->setFirings(DiskReadMda(m_firings));
        //WW->setCrossCorrelogramsData(DiskReadMda(m_cross_correlograms_data));
        //WW->setBaseLabel(k);
        QStringList text_labels;
        QList<int> labels1, labels2;
        for (int i = 1; i < m_original_cluster_numbers.count(); i++) {
            labels1 << k;
            labels2 << i;
            if ((i == 1) || (m_original_cluster_numbers[i] != m_original_cluster_numbers[i - 1])) {
                text_labels << QString("Cross %1(%2)/%3").arg(m_original_cluster_numbers[k + 1]).arg(m_original_cluster_offsets.value(k + 1)).arg(m_original_cluster_numbers[i]);
            } else
                text_labels << "";
        }
        //WW->setTextLabels(labels);
        WW->setLabelPairs(labels1, labels2, text_labels);
        //WW->updateWidget();
    } else if (widget_type == "matrix_of_cross_correlograms") {
        MVCrossCorrelogramsWidget2* WW = (MVCrossCorrelogramsWidget2*)W;
        QList<int> ks = string_list_to_int_list(W->property("ks").toStringList());
        WW->setColors(m_colors);
        WW->setSampleRate(m_samplerate);
        WW->setMaxDtTimepoints(cc_max_dt_timepoints());
        WW->setFirings(DiskReadMda(m_firings));
        //WW->setCrossCorrelogramsData(DiskReadMda(m_cross_correlograms_data));
        //WW->setLabelNumbers(ks);
        QStringList text_labels;
        QList<int> labels1, labels2;
        //text_labels << "";
        for (int a1 = 0; a1 < ks.count(); a1++) {
            QString str1 = QString("%1(%2)").arg(m_original_cluster_numbers[ks[a1]]).arg(m_original_cluster_offsets[ks[a1]]);
            for (int a2 = 0; a2 < ks.count(); a2++) {
                QString str2 = QString("%1(%2)").arg(m_original_cluster_numbers[ks[a2]]).arg(m_original_cluster_offsets[ks[a2]]);
                text_labels << QString("%1/%2").arg(str1).arg(str2);
                labels1 << ks.value(a1);
                labels2 << ks.value(a2);
            }
        }
        //WW->setTextLabels(labels);
        WW->setLabelPairs(labels1, labels2, text_labels);
        //WW->updateWidget();
    }
    /*else if (widget_type=="templates") {
        printf("Setting templates data...\n");
        SSTimeSeriesView *WW=(SSTimeSeriesView *)W;
        Mda TD=m_templates_data;
        DiskArrayModel_New *MM=new DiskArrayModel_New;
        MM->setFromMda(TD);
        int KK=TD.N3();
        QList<long> times,labels;
        int last_k=-1;
        for (int kk=1; kk<=KK; kk++) {
            int k=m_original_cluster_numbers.value(kk);
            if (k!=last_k) {
                times << (long)(TD.N2()*(kk-1+0.5));
                labels << k;
            }
            last_k=k;
        }
        WW->setData(MM,true);
        WW->setTimesLabels(times,labels);
        WW->setMarkerLinesVisible(false);
        printf(".\n");
    }*/
    else if (widget_type == "cluster_details") {
        MVClusterDetailWidget* WW = (MVClusterDetailWidget*)W;
        //int clip_size = m_control_panel->getParameterValue("clip_size").toInt();
        TaskProgress task("Update cluster details");
        int clip_size = m_control_panel_new->viewOptions().clip_size;
        WW->setColors(m_colors);
        DiskReadMda TT(current_timeseries_path());
        WW->setTimeseries(TT);
        WW->setClipSize(clip_size);
        WW->setFirings(DiskReadMda(m_firings));
        WW->setGroupNumbers(m_original_cluster_numbers);
        WW->zoomAllTheWayOut();
        task.log(QString("clip_size=%1, m_firings.N1()=%2, m_firings.N2()=%3").arg(clip_size).arg(m_firings.N1()).arg(m_firings.N2()));
    } else if (widget_type == "clips") {
        MVClipsWidget* WW = (MVClipsWidget*)W;
        int clip_size = m_control_panel_new->viewOptions().clip_size;
        QList<int> ks = string_list_to_int_list(WW->property("ks").toStringList());
        DiskReadMda TT(current_timeseries_path());
        WW->setTimeseries(TT);
        WW->setClipSize(clip_size);
        WW->setFirings(m_firings);
        WW->setLabelsToUse(ks);
    } else if (widget_type == "find_nearby_events") {
        printf("Extracting clips...\n");
        MVClipsView* WW = (MVClipsView*)W;
        QList<int> ks = string_list_to_int_list(WW->property("ks").toStringList());

        QList<int> labels;
        QList<double> times;

        for (int n = 0; n < m_firings.N2(); n++) {
            times << m_firings.value(1, n);
            labels << (int)m_firings.value(2, n);
        }

        int clip_size = m_control_panel_new->viewOptions().clip_size;

        QList<long> inds;
        if (widget_type == "clips") {
            for (int ik = 0; ik < ks.count(); ik++) {
                int kk = ks[ik];
                for (int n = 0; n < labels.count(); n++) {
                    if (labels[n] == kk) {
                        inds << n;
                    }
                }
            }
        } else if (widget_type == "find_nearby_events") {
            long max_time = (long)(compute_max(times) + 1);
            Mda occupied(ks.count(), max_time + 1);
            for (int ik = 0; ik < ks.count(); ik++) {
                int kk = ks[ik];
                for (int n = 0; n < labels.count(); n++) {
                    if (labels[n] == kk) {
                        occupied.setValue(1, ik, (long)times[n]);
                    }
                }
            }
            int radius = clip_size / 2 - 1;
            int kk1 = ks.value(0);
            for (long n = 0; n < labels.count(); n++) {
                if (labels[n] == kk1) {
                    bool okay = true;
                    for (int ik = 1; (ik < ks.count()) && (okay); ik++) {
                        bool found = false;
                        for (int dt = -radius; dt <= radius; dt++) {
                            if (occupied.value(ik, (long)times[n] + dt)) {
                                found = true;
                            }
                        }
                        if (!found)
                            okay = false;
                    }
                    if (okay)
                        inds << n;
                }
            }
        }

        QList<double> times_kk;
        QList<int> labels_kk;
        for (long i = 0; i < inds.count(); i++) {
            long n = inds[i];
            times_kk << times[n];
            labels_kk << labels[n];
        }

        DiskReadMda TT(current_timeseries_path());
        Mda clips = extract_clips(TT, times_kk, clip_size);
        printf("Setting data...\n");
        //DiskArrayModel_New *DAM=new DiskArrayModel_New;
        //DAM->setFromMda(clips);
        //WW->setData(DAM,true);
        WW->setClips(clips);
        WW->setTimes(times_kk);
        WW->setLabels(labels_kk);
        printf(".\n");
    } else if (widget_type == "clusters") {
        MVClusterWidget* WW = (MVClusterWidget*)W;
        int clip_size = m_control_panel_new->viewOptions().clip_size;
        QList<int> ks = string_list_to_int_list(WW->property("ks").toStringList());
        DiskReadMda TT(current_timeseries_path());
        WW->setTimeseries(TT);
        WW->setClipSize(clip_size);
        WW->setFirings(m_firings);
        WW->setLabelsToUse(ks);
    } else if (widget_type == "firing_events") {
        MVFiringEventView* WW = (MVFiringEventView*)W;
        QList<int> ks = string_list_to_int_list(WW->property("ks").toStringList());
        QSet<int> ks_set = ks.toSet();

        QList<double> times, amplitudes;
        QList<int> labels;
        for (int i = 0; i < m_firings.N2(); i++) {
            int label = (int)m_firings.value(2, i);
            if (ks_set.contains(label)) {
                times << m_firings.value(1, i);
                amplitudes << m_firings.value(3, i);
                labels << label;
            }
        }
        Mda firings2;
        firings2.allocate(4, times.count());
        for (int i = 0; i < times.count(); i++) {
            firings2.setValue(times[i], 1, i);
            firings2.setValue(labels[i], 2, i);
            firings2.setValue(amplitudes[i], 3, i);
        }
        WW->setFirings(firings2);
        WW->setSampleRate(m_samplerate);
        WW->setEpochs(m_epochs);
    } else if (widget_type == "timeseries") {
        SSTimeSeriesWidget* WW = (SSTimeSeriesWidget*)W;
        DiskArrayModel_New* X = new DiskArrayModel_New;
        X->setPath(current_timeseries_path());
        ((SSTimeSeriesView*)(WW->view()))->setData(X, true);
        set_times_labels_for_timeseries_widget(WW);
    }
}

void MVOverview2WidgetPrivate::set_cross_correlograms_current_index(int index)
{
    QList<QWidget*> widgets = get_all_widgets();
    foreach(QWidget * W, widgets)
    {
        QString widget_type = W->property("widget_type").toString();
        if ((widget_type == "auto_correlograms") || (widget_type == "cross_correlograms")) {
            MVCrossCorrelogramsWidget2* WW = (MVCrossCorrelogramsWidget2*)W;
            WW->setCurrentIndex(index);
        }
    }
}

void MVOverview2WidgetPrivate::set_cross_correlograms_selected_indices(const QList<int>& indices)
{
    QList<QWidget*> widgets = get_all_widgets();
    foreach(QWidget * W, widgets)
    {
        QString widget_type = W->property("widget_type").toString();
        if ((widget_type == "auto_correlograms") || (widget_type == "cross_correlograms")) {
            MVCrossCorrelogramsWidget2* WW = (MVCrossCorrelogramsWidget2*)W;
            WW->setSelectedIndices(indices);
        }
    }
}

void MVOverview2WidgetPrivate::set_templates_current_number(int kk)
{
    QList<QWidget*> widgets = get_all_widgets();
    foreach(QWidget * W, widgets)
    {
        QString widget_type = W->property("widget_type").toString();
        if (widget_type == "cluster_details") {
            MVClusterDetailWidget* WW = (MVClusterDetailWidget*)W;
            WW->setCurrentK(kk);
        }
    }
}

void MVOverview2WidgetPrivate::set_templates_selected_numbers(const QList<int>& kks)
{
    QList<QWidget*> widgets = get_all_widgets();
    foreach(QWidget * W, widgets)
    {
        QString widget_type = W->property("widget_type").toString();
        if (widget_type == "cluster_details") {
            MVClusterDetailWidget* WW = (MVClusterDetailWidget*)W;
            WW->setSelectedKs(kks);
        }
    }
}

void MVOverview2WidgetPrivate::set_times_labels_for_timeseries_widget(SSTimeSeriesWidget* WW)
{
    QList<long> times, labels;
    for (int n = 0; n < m_firings_original.N2(); n++) {
        long label0 = (long)m_firings_original.value(2, n);
        if ((m_selected_ks.isEmpty()) || (m_selected_ks.contains(label0))) {
            times << (long)m_firings_original.value(1, n);
            labels << label0;
        }
    }
    SSTimeSeriesView* V = (SSTimeSeriesView*)WW->view();
    V->setTimesLabels(times, labels);

    /*
    QList<QWidget*> widgets = get_all_widgets();
    foreach (QWidget* W, widgets) {
        QString widget_type = W->property("widget_type").toString();
        if (widget_type == "timeseries") {
            SSTimeSeriesWidget* WW = (SSTimeSeriesWidget*)W;
            SSTimeSeriesView* V = (SSTimeSeriesView*)WW->view();
            V->setTimesLabels(times, labels);
        }
    }
    */
}

QList<QWidget*> MVOverview2WidgetPrivate::get_all_widgets()
{
    return m_tabber->allWidgets();
    /*
    QList<QWidget *> ret;
    for (int i=0; i<m_tabs1->count(); i++) {
        ret << m_tabs1->widget(i);
    }
    for (int i=0; i<m_tabs2->count(); i++) {
        ret << m_tabs2->widget(i);
    }
    return ret;
    */
}

TabberTabWidget* MVOverview2WidgetPrivate::tab_widget_of(QWidget* W)
{
    for (int i = 0; i < m_tabs1->count(); i++) {
        if (m_tabs1->widget(i) == W)
            return m_tabs1;
    }
    for (int i = 0; i < m_tabs2->count(); i++) {
        if (m_tabs2->widget(i) == W)
            return m_tabs2;
    }
    return m_tabs1;
}

void MVOverview2WidgetPrivate::remove_widgets_of_type(QString widget_type)
{
    QList<QWidget*> list = get_all_widgets();
    foreach(QWidget * W, list)
    {
        if (W->property("widget_type") == widget_type) {
            delete W;
        }
    }
}

Mda MVOverview2WidgetPrivate::compute_centroid(Mda& clips)
{
    int M = clips.N1();
    int T = clips.N2();
    int NC = clips.N3();
    Mda ret;
    ret.allocate(M, T);
    double* retptr = ret.dataPtr();
    double* clipsptr = clips.dataPtr();
    for (int i = 0; i < NC; i++) {
        int iii = i * M * T;
        int jjj = 0;
        for (int t = 0; t < T; t++) {
            for (int m = 0; m < M; m++) {
                retptr[jjj] += clipsptr[iii];
                iii++;
                jjj++;
            }
        }
    }
    if (NC) {
        for (int t = 0; t < T; t++) {
            for (int m = 0; m < M; m++) {
                double val = ret.value(m, t) / NC;
                ret.setValue(val, m, t);
            }
        }
    }
    return ret;
}

Mda MVOverview2WidgetPrivate::compute_geometric_median(Mda& clips, int num_iterations)
{
    int M = clips.N1();
    int T = clips.N2();
    int L = clips.N3();
    double* clipsptr = clips.dataPtr();
    Mda ret;
    ret.allocate(M, T);
    if (L == 0)
        return ret;

    int num_features = 6;
    Mda features;
    features.allocate(num_features, L);
    double* featuresptr = features.dataPtr();
#ifdef USE_LAPACK
    get_pca_features(M * T, L, num_features, featuresptr, clipsptr);
#else
    get_pca_features_2(M * T, L, num_features, featuresptr, clipsptr);
#endif

    Mda geomedian;
    geomedian.allocate(num_features, 1);
    double* geomedianptr = geomedian.dataPtr();
    compute_geometric_median(num_features, L, geomedianptr, featuresptr, num_iterations);

    QList<double> dists;
    {
        int kkk = 0;
        for (int j = 0; j < L; j++) {
            double sumsqr = 0;
            for (int a = 0; a < num_features; a++) {
                double diff = featuresptr[kkk] - geomedianptr[a];
                sumsqr += diff * diff;
                kkk++;
            }
            dists << sqrt(sumsqr);
        }
    }

    QList<double> dists_sorted = dists;
    qSort(dists_sorted);
    double cutoff = dists_sorted[(int)(L * 0.3)];
    QList<int> inds_to_use;
    for (int j = 0; j < L; j++) {
        if (dists[j] <= cutoff)
            inds_to_use << j;
    }

    Mda clips2;
    clips2.allocate(M, T, inds_to_use.count());
    for (int j = 0; j < inds_to_use.count(); j++) {
        int jj = inds_to_use[j];
        for (int t = 0; t < T; t++) {
            for (int m = 0; m < M; m++) {
                clips2.setValue(clips.value(m, t, jj), m, t, j);
            }
        }
    }

    ret = compute_centroid(clips2);
    return ret;
}

void MVOverview2WidgetPrivate::compute_geometric_median(int M, int N, double* output, double* input, int num_iterations)
{
    double* weights = (double*)malloc(sizeof(double) * N);
    double* dists = (double*)malloc(sizeof(double) * N);
    for (int j = 0; j < N; j++)
        weights[j] = 1;
    for (int it = 1; it <= num_iterations; it++) {
        float sumweights = 0;
        for (int j = 0; j < N; j++)
            sumweights += weights[j];
        if (sumweights)
            for (int j = 0; j < N; j++)
                weights[j] /= sumweights;
        for (int i = 0; i < M; i++)
            output[i] = 0;
        {
            //compute output
            int kkk = 0;
            for (int j = 0; j < N; j++) {
                int i = 0;
                for (int m = 0; m < M; m++) {
                    output[i] += weights[j] * input[kkk];
                    i++;
                    kkk++;
                }
            }
        }
        {
            //compute dists
            int kkk = 0;
            for (int j = 0; j < N; j++) {
                int i = 0;
                double sumsqr = 0;
                for (int m = 0; m < M; m++) {
                    double diff = output[i] - input[kkk];
                    i++;
                    kkk++;
                    sumsqr += diff * diff;
                }
                dists[j] = sqrt(sumsqr);
            }
        }
        {
            //compute weights
            for (int j = 0; j < N; j++) {
                if (dists[j])
                    weights[j] = 1 / dists[j];
                else
                    weights[j] = 0;
            }
        }
    }
    free(dists);
    free(weights);
}

void MVOverview2WidgetPrivate::set_progress(QString title, QString text, float frac)
{
    if (!m_progress_dialog) {
        m_progress_dialog = new QProgressDialog;
        m_progress_dialog->setCancelButton(0);
    }
    static QTime* timer = 0;
    if (!timer) {
        timer = new QTime;
        timer->start();
        m_progress_dialog->show();
        m_progress_dialog->repaint();
        qApp->processEvents();
    }
    if (timer->elapsed() > 500) {
        timer->restart();
        if (!m_progress_dialog->isVisible()) {
            m_progress_dialog->show();
        }
        m_progress_dialog->setLabelText(text);
        m_progress_dialog->setWindowTitle(title);
        m_progress_dialog->setValue((int)(frac * 100));
        m_progress_dialog->repaint();
        qApp->processEvents();
    }
    if (frac >= 1) {
        delete m_progress_dialog;
        m_progress_dialog = 0;
    }
}

void MVOverview2WidgetPrivate::set_current_event(MVEvent evt)
{
    if ((m_current_event.time == evt.time) && (m_current_event.label == evt.label)) {
        return;
    }
    m_current_event = evt;
    QList<QWidget*> widgets = get_all_widgets();
    foreach(QWidget * W, widgets)
    {
        QString widget_type = W->property("widget_type").toString();
        if (widget_type == "clips") {
            MVClipsWidget* WW = (MVClipsWidget*)W;
            WW->setCurrentEvent(evt);
        } else if (widget_type == "find_nearby_events") {
            MVClipsView* WW = (MVClipsView*)W;
            WW->setCurrentEvent(evt);
        } else if (widget_type == "clusters") {
            MVClusterWidget* WW = (MVClusterWidget*)W;
            WW->setCurrentEvent(evt);
        } else if (widget_type == "firing_events") {
            MVFiringEventView* WW = (MVFiringEventView*)W;
            WW->setCurrentEvent(evt);
        } else if (widget_type == "cluster_details") {
            MVClusterDetailWidget* WW = (MVClusterDetailWidget*)W;
            if (evt.label > 0) {
                WW->setCurrentK(evt.label);
            }
        } else if (widget_type == "timeseries") {
            SSTimeSeriesWidget* WW = (SSTimeSeriesWidget*)W;
            SSTimeSeriesView* VV = (SSTimeSeriesView*)WW->view(0);
            if (evt.time >= 0) {
                VV->setCurrentX(evt.time);
            }
        }
    }
}

long MVOverview2WidgetPrivate::cc_max_dt_timepoints()
{
    //return (int)(m_control_panel->getParameterValue("max_dt").toInt() * m_samplerate / 1000);
    return m_control_panel_new->viewOptions().cc_max_dt_msec * m_samplerate / 1000;
}

class DownloadComputer : public ComputationThread {
public:
    //inputs
    QString source_path;
    QString dest_path;
    bool use_float64;

    void compute();
};
void DownloadComputer::compute()
{
    TaskProgress task("Downlading");
    task.setDescription(QString("Downloading %1 to %2").arg(source_path).arg(dest_path));
    DiskReadMda X(source_path);
    X.setHaltAgent(this);
    Mda Y;
    task.setProgress(0.2);
    task.log(QString("Reading/Downloading %1x%2x%3").arg(X.N1()).arg(X.N2()).arg(X.N3()));
    if (!X.readChunk(Y, 0, 0, 0, X.N1(), X.N2(), X.N3())) {
        if (this->stopRequested()) {
            task.error("Halted download: " + source_path);
        } else {
            task.error("Failed to readChunk from: " + source_path);
        }
        return;
    }
    task.setProgress(0.8);
    if (use_float64) {
        task.log("Writing 64-bit to " + dest_path);
        Y.write64(dest_path);
    } else {
        task.log("Writing 32-bit to " + dest_path);
        Y.write32(dest_path);
    }
}

void MVOverview2WidgetPrivate::export_mountainview_document()
{
    QString default_dir = "";
    QString fname = QFileDialog::getSaveFileName(q, "Export mountainview document", default_dir, "*.mv");
    if (QFileInfo(fname).suffix() != "mv")
        fname = fname + ".mv";
    q->saveMVFile(fname);
}

void MVOverview2WidgetPrivate::export_original_firings()
{
    QString default_dir = "";
    QString fname = QFileDialog::getSaveFileName(q, "Export original firings", default_dir, "*.mda");
    if (QFileInfo(fname).suffix() != "mda")
        fname = fname + ".mda";
    if (!fname.isEmpty()) {
        export_file(m_firings_original.path(), fname, true);
    }
}

void MVOverview2WidgetPrivate::export_filtered_firings()
{
    QString default_dir = "";
    QString fname = QFileDialog::getSaveFileName(q, "Export filtered firings", default_dir, "*.mda");
    if (QFileInfo(fname).suffix() != "mda")
        fname = fname + ".mda";
    if (!fname.isEmpty()) {
        export_file(m_firings.path(), fname, true);
    }
}

void MVOverview2WidgetPrivate::export_file(QString source_path, QString dest_path, bool use_float64)
{
    DownloadComputer* C = new DownloadComputer;
    C->source_path = source_path;
    C->dest_path = dest_path;
    C->use_float64 = use_float64;
    C->setDeleteOnComplete(true);
    C->startComputation();
}

QString MVOverview2WidgetPrivate::make_absolute_path(QString path)
{
    if (path.startsWith("http"))
        return path;
    if (QFileInfo(path).isAbsolute())
        return path;
    if (!m_mv_fname.isEmpty()) {
        return QFileInfo(m_mv_fname).path() + "/" + path;
    }
    return path;
}

QString MVOverview2WidgetPrivate::current_timeseries_path()
{
    return make_absolute_path(m_timeseries_paths.value(m_control_panel_new->viewOptions().timeseries));
}

QVariant MVOverview2WidgetPrivate::get_cluster_attribute(int k, QString attr)
{
    return m_view_agent.clusterAttributes().value(k).value(attr).toVariant();
}

void MVOverview2WidgetPrivate::set_cluster_attribute(int k, QString attr, QVariant val)
{
    QMap<int, QJsonObject> CA = m_view_agent.clusterAttributes();
    CA[k][attr] = QJsonValue::fromVariant(val);
    m_view_agent.setClusterAttributes(CA);
}

void MVOverview2WidgetPrivate::set_button_enabled(QString name, bool val)
{
    QAbstractButton* B = m_control_panel_new->findButton(name);
    if (B)
        B->setEnabled(val);
}
