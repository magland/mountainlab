/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
** Created: 2/22/2017
*******************************************************/

#include "mlcommon.h"

#include <QJsonArray>
#include <QJsonDocument>
#include <QCoreApplication>
#include "mountainsort2_main.h"
#include "p_extract_neighborhood_timeseries.h"
#include "p_extract_segment_timeseries.h"
#include "p_bandpass_filter.h"
#include "p_detect_events.h"
#include "p_extract_clips.h"
#include "p_sort_clips.h"
#include "p_consolidate_clusters.h"
#include "p_create_firings.h"
#include "p_combine_firings.h"
#include "p_fit_stage.h"
#include "p_whiten.h"
#include "p_apply_timestamp_offset.h"
#include "p_link_segments.h"
#include "p_cluster_metrics.h"
#include "p_split_firings.h"
#include "p_concat_timeseries.h"
#include "p_concat_firings.h"
#include "p_compute_templates.h"
#include "p_load_test.h"
#include "p_compute_amplitudes.h"
#include "p_extract_time_interval.h"
#include "p_isolation_metrics.h"
#include "p_generate_background_dataset.h"

#include "omp.h"
#include "p_confusion_matrix.h"

QJsonObject get_spec()
{
    QJsonArray processors;

    {
        ProcessorSpec X("mountainsort.extract_neighborhood_timeseries", "0.1");
        X.addInputs("timeseries");
        X.addOutputs("timeseries_out");
        X.addRequiredParameters("channels");
        processors.push_back(X.get_spec());
    }
    {
        ProcessorSpec X("mountainsort.extract_segment_timeseries", "0.1");
        X.addInputs("timeseries");
        X.addOutputs("timeseries_out");
        X.addRequiredParameters("t1", "t2");
        X.addOptionalParameter("channels", "Comma-separated list of channels to extract", "");
        processors.push_back(X.get_spec());
    }
    {
        ProcessorSpec X("mountainsort.extract_segment_firings", "0.11");
        X.addInputs("firings");
        X.addOutputs("firings_out");
        X.addRequiredParameters("t1", "t2");
        processors.push_back(X.get_spec());
    }
    {
        ProcessorSpec X("mountainsort.bandpass_filter", "0.18");
        X.addInputs("timeseries");
        X.addOutputs("timeseries_out");
        X.addRequiredParameters("samplerate", "freq_min", "freq_max");
        X.addOptionalParameter("freq_wid", "", 1000);
        X.addOptionalParameter("quantization_unit", "", 0);
        X.addOptionalParameter("testcode", "", "");
        processors.push_back(X.get_spec());
    }
    {
        ProcessorSpec X("mountainsort.whiten", "0.1");
        X.addInputs("timeseries");
        X.addOutputs("timeseries_out");
        //X.addRequiredParameters();
        processors.push_back(X.get_spec());
    }
    {
        ProcessorSpec X("mountainsort.compute_whitening_matrix", "0.11");
        X.addInputs("timeseries_list");
        X.addOutputs("whitening_matrix_out");
        X.addOptionalParameter("channels");
        //X.addRequiredParameters();
        processors.push_back(X.get_spec());
    }
    {
        ProcessorSpec X("mountainsort.whiten_clips", "0.1");
        X.addInputs("clips", "whitening_matrix");
        X.addOutputs("clips_out");
        //X.addRequiredParameters();
        processors.push_back(X.get_spec());
    }
    {
        ProcessorSpec X("mountainsort.apply_whitening_matrix", "0.1");
        X.addInputs("timeseries", "whitening_matrix");
        X.addOutputs("timeseries_out");
        //X.addRequiredParameters();
        processors.push_back(X.get_spec());
    }
    {
        ProcessorSpec X("mountainsort.detect_events", "0.13");
        X.addInputs("timeseries");
        X.addOutputs("event_times_out");
        X.addRequiredParameters("central_channel", "detect_threshold", "detect_interval", "sign");
        X.addOptionalParameter("subsample_factor", "", 1);
        processors.push_back(X.get_spec());
    }
    {
        ProcessorSpec X("mountainsort.extract_clips", "0.11");
        X.addInputs("timeseries", "event_times");
        X.addOutputs("clips_out");
        X.addRequiredParameters("clip_size");
        X.addOptionalParameter("channels");
        processors.push_back(X.get_spec());
    }
    {
        ProcessorSpec X("mountainsort.compute_templates", "0.11");
        X.addInputs("timeseries", "firings");
        X.addOutputs("templates_out");
        X.addRequiredParameters("clip_size");
        X.addOptionalParameter("clusters", "Comma-separated list of clusters to inclue", "");
        processors.push_back(X.get_spec());
    }
    {
        ProcessorSpec X("mountainsort.sort_clips", "0.11a");
        X.addInputs("clips");
        X.addOutputs("labels_out");
        //X.addRequiredParameters();
        X.addOptionalParameter("isocut_threshold", "", 1);
        processors.push_back(X.get_spec());
    }
    {
        ProcessorSpec X("mountainsort.reorder_labels", "0.1");
        X.addInputs("templates", "firings");
        X.addOutputs("firings_out");
        //X.addRequiredParameters();
        processors.push_back(X.get_spec());
    }
    {
        ProcessorSpec X("mountainsort.consolidate_clusters", "0.1");
        X.addInputs("clips", "labels");
        X.addOptionalInputs("amplitudes");
        X.addOutputs("labels_out");
        X.addRequiredParameters("central_channel");
        processors.push_back(X.get_spec());
    }
    {
        ProcessorSpec X("mountainsort.create_firings", "0.1");
        X.addInputs("event_times", "labels");
        X.addOptionalInputs("amplitudes");
        X.addOutputs("firings_out");
        X.addRequiredParameters("central_channel");
        processors.push_back(X.get_spec());
    }
    {
        ProcessorSpec X("mountainsort.combine_firings", "0.1");
        X.addInputs("firings_list");
        X.addOutputs("firings_out");
        //X.addRequiredParameters();
        X.addOptionalParameters("increment_labels", "", "true");
        processors.push_back(X.get_spec());
    }
    {
        ProcessorSpec X("mountainsort.fit_stage", "0.1");
        X.addInputs("timeseries", "firings");
        X.addOutputs("firings_out");
        //X.addRequiredParameters();
        processors.push_back(X.get_spec());
    }
    {
        ProcessorSpec X("mountainsort.apply_timestamp_offset", "0.1");
        X.addInputs("firings");
        X.addOutputs("firings_out");
        X.addRequiredParameters("timestamp_offset");
        processors.push_back(X.get_spec());
    }
    {
        ProcessorSpec X("mountainsort.link_segments", "0.1");
        X.addInputs("firings", "firings_prev", "Kmax_prev");
        X.addOutputs("firings_out", "Kmax_out");
        X.addOutputs("firings_subset_out", "Kmax_out");
        X.addRequiredParameters("t1", "t2", "t1_prev", "t2_prev");
        processors.push_back(X.get_spec());
    }
    {
        ProcessorSpec X("mountainsort.cluster_metrics", "0.11");
        X.addInputs("timeseries", "firings");
        X.addOutputs("cluster_metrics_out");
        X.addRequiredParameters("samplerate");
        processors.push_back(X.get_spec());
    }
    {
        ProcessorSpec X("mountainsort.isolation_metrics", "0.13");
        X.addInputs("timeseries", "firings");
        X.addOutputs("metrics_out");
        X.addOptionalOutputs("pair_metrics_out");
        processors.push_back(X.get_spec());
    }
    {
        ProcessorSpec X("mountainsort.combine_cluster_metrics", "0.1");
        X.addInputs("metrics_list");
        X.addOutputs("metrics_out");
        processors.push_back(X.get_spec());
    }
    {
        ProcessorSpec X("mountainsort.split_firings", "0.14");
        X.addInputs("timeseries_list", "firings");
        X.addOutputs("firings_out_list");
        //X.addRequiredParameters();
        processors.push_back(X.get_spec());
    }
    {
        ProcessorSpec X("mountainsort.concat_timeseries", "0.11");
        X.addInputs("timeseries_list");
        X.addOutputs("timeseries_out");
        //X.addRequiredParameters();
        processors.push_back(X.get_spec());
    }
    {
        ProcessorSpec X("mountainsort.concat_firings", "0.13");
        X.addInputs("firings_list");
        X.addOptionalInputs("timeseries_list");
        X.addOutputs("firings_out");
        X.addOptionalOutputs("timeseries_out");
        //X.addRequiredParameters();
        processors.push_back(X.get_spec());
    }
    {
        ProcessorSpec X("mountainsort.concat_event_times", "0.1");
        X.addInputs("event_times_list");
        X.addOutputs("event_times_out");
        //X.addRequiredParameters();
        processors.push_back(X.get_spec());
    }
    {
        ProcessorSpec X("mountainsort.load_test", "0.1");
        X.addOutputs("stats_out");
        X.addOptionalParameters("num_read_bytes", "num_write_bytes", "num_cpu_ops");
        processors.push_back(X.get_spec());
    }
    {
        ProcessorSpec X("mountainsort.misc_test", "0.1");
        X.addInput("dir");
        X.addOutputs("info_out");
        processors.push_back(X.get_spec());
    }
    {
        ProcessorSpec X("mountainsort.compute_amplitudes", "0.1");
        X.addInputs("timeseries", "event_times");
        X.addOutputs("amplitudes_out");
        X.addRequiredParameters("central_channel");
        processors.push_back(X.get_spec());
    }
    {
        ProcessorSpec X("mountainsort.extract_time_interval", "0.1");
        X.addInputs("timeseries_list", "firings");
        X.addOutputs("timeseries_out", "firings_out");
        X.addRequiredParameters("t1", "t2");
        processors.push_back(X.get_spec());
    }
    {
        ProcessorSpec X("mountainsort.confusion_matrix", "0.13");
        X.addInputs("firings1", "firings2");
        X.addOutputs("confusion_matrix_out");
        X.addOptionalOutputs("matched_firings_out", "label_map_out", "firings2_relabeled_out", "firings2_relabel_map_out");
        X.addOptionalParameter("max_matching_offset", "", 30);
        X.addOptionalParameter("relabel_firings2", "", "false");
        processors.push_back(X.get_spec());
    }
    {
        ProcessorSpec X("mountainsort.generate_background_dataset", "0.1");
        X.addInputs("timeseries", "event_times");
        X.addOutputs("timeseries_out");
        processors.push_back(X.get_spec());
    }

    QJsonObject ret;
    ret["processors"] = processors;
    return ret;
}

int main(int argc, char* argv[])
{
    QCoreApplication app(argc, argv);

    CLParams CLP(argc, argv);

    QString arg1 = CLP.unnamed_parameters.value(0);

    if (arg1 == "spec") {
        QJsonObject spec = get_spec();
        QString json = QJsonDocument(spec).toJson(QJsonDocument::Indented);
        printf("%s\n", json.toUtf8().data());
        return 0;
    }

    bool ret = false;

    if (CLP.named_parameters.contains("_request_num_threads")) {
        int num_threads = CLP.named_parameters.value("_request_num_threads", 0).toInt();
        if (num_threads) {
            qDebug().noquote() << "Setting num threads:" << num_threads;
            omp_set_num_threads(num_threads);
        }
    }

    if (arg1 == "mountainsort.extract_neighborhood_timeseries") {
        QString timeseries = CLP.named_parameters["timeseries"].toString();
        QString timeseries_out = CLP.named_parameters["timeseries_out"].toString();
        QStringList channels_str = CLP.named_parameters["channels"].toString().split(",", QString::SkipEmptyParts);
        QList<int> channels = MLUtil::stringListToIntList(channels_str);
        ret = p_extract_neighborhood_timeseries(timeseries, timeseries_out, channels);
    }
    else if (arg1 == "mountainsort.extract_segment_timeseries") {
        QStringList timeseries_list = MLUtil::toStringList(CLP.named_parameters["timeseries"]);
        QString timeseries_out = CLP.named_parameters["timeseries_out"].toString();
        bigint t1 = CLP.named_parameters["t1"].toDouble(); //to double to handle scientific notation
        bigint t2 = CLP.named_parameters["t2"].toDouble(); //to double to handle scientific notation
        QStringList channels_str = CLP.named_parameters["channels"].toString().split(",", QString::SkipEmptyParts);
        QList<int> channels = MLUtil::stringListToIntList(channels_str);
        if (timeseries_list.count() <= 1) {
            ret = p_extract_segment_timeseries(timeseries_list.value(0), timeseries_out, t1, t2, channels);
        }
        else {
            ret = p_extract_segment_timeseries_from_concat_list(timeseries_list, timeseries_out, t1, t2, channels);
        }
    }
    else if (arg1 == "mountainsort.extract_segment_firings") {
        QString firings = CLP.named_parameters["firings"].toString();
        QString firings_out = CLP.named_parameters["firings_out"].toString();
        bigint t1 = CLP.named_parameters["t1"].toDouble(); //to double to handle scientific notation
        bigint t2 = CLP.named_parameters["t2"].toDouble(); //to double to handle scientific notation
        ret = p_extract_segment_firings(firings, firings_out, t1, t2);
    }
    else if (arg1 == "mountainsort.bandpass_filter") {
        QString timeseries = CLP.named_parameters["timeseries"].toString();
        QString timeseries_out = CLP.named_parameters["timeseries_out"].toString();
        Bandpass_filter_opts opts;
        opts.samplerate = CLP.named_parameters["samplerate"].toDouble();
        opts.freq_min = CLP.named_parameters["freq_min"].toDouble();
        opts.freq_max = CLP.named_parameters["freq_max"].toDouble();
        opts.freq_wid = CLP.named_parameters.value("freq_wid", 1000).toDouble();
        opts.quantization_unit = CLP.named_parameters.value("quantization_unit").toDouble();
        opts.testcode = CLP.named_parameters.value("testcode", "").toString();
        ret = p_bandpass_filter(timeseries, timeseries_out, opts);
    }
    else if (arg1 == "mountainsort.whiten") {
        QString timeseries = CLP.named_parameters["timeseries"].toString();
        QString timeseries_out = CLP.named_parameters["timeseries_out"].toString();
        Whiten_opts opts;
        ret = p_whiten(timeseries, timeseries_out, opts);
    }
    else if (arg1 == "mountainsort.compute_whitening_matrix") {
        QStringList timeseries_list = MLUtil::toStringList(CLP.named_parameters["timeseries_list"]);
        QString whitening_matrix_out = CLP.named_parameters["whitening_matrix_out"].toString();
        QStringList channels_str = CLP.named_parameters["channels"].toString().split(",", QString::SkipEmptyParts);
        QList<int> channels = MLUtil::stringListToIntList(channels_str);
        Whiten_opts opts;
        ret = p_compute_whitening_matrix(timeseries_list, channels, whitening_matrix_out, opts);
    }
    else if (arg1 == "mountainsort.whiten_clips") {
        QString clips = CLP.named_parameters["clips"].toString();
        QString whitening_matrix = CLP.named_parameters["whitening_matrix"].toString();
        QString clips_out = CLP.named_parameters["clips_out"].toString();
        Whiten_opts opts;
        ret = p_whiten_clips(clips, whitening_matrix, clips_out, opts);
    }
    else if (arg1 == "mountainsort.apply_whitening_matrix") {
        QString timeseries = CLP.named_parameters["timeseries"].toString();
        QString whitening_matrix = CLP.named_parameters["whitening_matrix"].toString();
        QString timeseries_out = CLP.named_parameters["timeseries_out"].toString();
        Whiten_opts opts;
        ret = p_apply_whitening_matrix(timeseries, whitening_matrix, timeseries_out, opts);
    }
    else if (arg1 == "mountainsort.detect_events") {
        QString timeseries = CLP.named_parameters["timeseries"].toString();
        QString event_times_out = CLP.named_parameters["event_times_out"].toString();
        P_detect_events_opts opts;
        opts.central_channel = CLP.named_parameters["central_channel"].toInt();
        opts.detect_threshold = CLP.named_parameters["detect_threshold"].toDouble();
        opts.detect_interval = CLP.named_parameters["detect_interval"].toDouble();
        opts.sign = CLP.named_parameters["sign"].toInt();
        opts.subsample_factor = CLP.named_parameters["subsample_factor"].toDouble();
        ret = p_detect_events(timeseries, event_times_out, opts);
    }
    else if (arg1 == "mountainsort.extract_clips") {
        QStringList timeseries_list = MLUtil::toStringList(CLP.named_parameters["timeseries"]);
        QString event_times = CLP.named_parameters["event_times"].toString();
        QString clips_out = CLP.named_parameters["clips_out"].toString();
        QStringList channels_str = CLP.named_parameters["channels"].toString().split(",", QString::SkipEmptyParts);
        QList<int> channels = MLUtil::stringListToIntList(channels_str);
        ret = p_extract_clips(timeseries_list, event_times, channels, clips_out, CLP.named_parameters);
    }
    else if (arg1 == "mountainsort.compute_templates") {
        QStringList timeseries_list = MLUtil::toStringList(CLP.named_parameters["timeseries"]);
        QString firings = CLP.named_parameters["firings"].toString();
        QString templates_out = CLP.named_parameters["templates_out"].toString();
        int clip_size = CLP.named_parameters["clip_size"].toInt();
        QList<int> clusters = MLUtil::stringListToIntList(CLP.named_parameters["clusters"].toString().split(",", QString::SkipEmptyParts));
        ret = p_compute_templates(timeseries_list, firings, templates_out, clip_size, clusters);
    }
    else if (arg1 == "mountainsort.sort_clips") {
        QString clips = CLP.named_parameters["clips"].toString();
        QString labels_out = CLP.named_parameters["labels_out"].toString();
        Sort_clips_opts opts;
        ret = p_sort_clips(clips, labels_out, opts);
    }
    else if (arg1 == "mountainsort.reorder_labels") {
        QString templates = CLP.named_parameters["templates"].toString();
        QString firings = CLP.named_parameters["firings"].toString();
        QString firings_out = CLP.named_parameters["firings_out"].toString();
        ret = p_reorder_labels(templates, firings, firings_out);
    }
    else if (arg1 == "mountainsort.consolidate_clusters") {
        QString clips = CLP.named_parameters["clips"].toString();
        QString labels = CLP.named_parameters["labels"].toString();
        QString labels_out = CLP.named_parameters["labels_out"].toString();
        Consolidate_clusters_opts opts;
        opts.central_channel = CLP.named_parameters["central_channel"].toInt();
        ret = p_consolidate_clusters(clips, labels, labels_out, opts);
    }
    else if (arg1 == "mountainsort.create_firings") {
        QString event_times = CLP.named_parameters["event_times"].toString();
        QString labels = CLP.named_parameters["labels"].toString();
        QString amplitudes = CLP.named_parameters["amplitudes"].toString();
        QString firings_out = CLP.named_parameters["firings_out"].toString();
        int central_channel = CLP.named_parameters["central_channel"].toInt();
        ret = p_create_firings(event_times, labels, amplitudes, firings_out, central_channel);
    }
    else if (arg1 == "mountainsort.combine_firings") {
        QStringList firings_list = MLUtil::toStringList(CLP.named_parameters["firings_list"]);
        QString firings_out = CLP.named_parameters["firings_out"].toString();
        QString tmp = CLP.named_parameters.value("increment_labels", "").toString();
        bool increment_labels = (tmp == "true");
        ret = p_combine_firings(firings_list, firings_out, increment_labels);
    }
    else if (arg1 == "mountainsort.fit_stage") {
        QString timeseries = CLP.named_parameters["timeseries"].toString();
        QString firings = CLP.named_parameters["firings"].toString();
        QString firings_out = CLP.named_parameters["firings_out"].toString();
        Fit_stage_opts opts;
        ret = p_fit_stage(timeseries, firings, firings_out, opts);
    }
    else if (arg1 == "mountainsort.apply_timestamp_offset") {
        QString firings = CLP.named_parameters["firings"].toString();
        QString firings_out = CLP.named_parameters["firings_out"].toString();
        double timestamp_offset = CLP.named_parameters["timestamp_offset"].toDouble();
        ret = p_apply_timestamp_offset(firings, firings_out, timestamp_offset);
    }
    else if (arg1 == "mountainsort.link_segments") {
        QString firings = CLP.named_parameters["firings"].toString();
        QString firings_prev = CLP.named_parameters["firings_prev"].toString();
        QString Kmax_prev = CLP.named_parameters["Kmax_prev"].toString();

        QString firings_out = CLP.named_parameters["firings_out"].toString();
        QString firings_subset_out = CLP.named_parameters["firings_subset_out"].toString();
        QString Kmax_out = CLP.named_parameters["Kmax_out"].toString();

        double t1 = CLP.named_parameters["t1"].toDouble();
        double t2 = CLP.named_parameters["t2"].toDouble();
        double t1_prev = CLP.named_parameters["t1_prev"].toDouble();
        double t2_prev = CLP.named_parameters["t2_prev"].toDouble();
        ret = p_link_segments(firings, firings_prev, Kmax_prev, firings_out, firings_subset_out, Kmax_out, t1, t2, t1_prev, t2_prev);
    }
    else if (arg1 == "mountainsort.cluster_metrics") {
        QString timeseries = CLP.named_parameters["timeseries"].toString();
        QString firings = CLP.named_parameters["firings"].toString();
        QString cluster_metrics_out = CLP.named_parameters["cluster_metrics_out"].toString();
        Cluster_metrics_opts opts;
        opts.samplerate = CLP.named_parameters["samplerate"].toDouble();
        ret = p_cluster_metrics(timeseries, firings, cluster_metrics_out, opts);
    }
    else if (arg1 == "mountainsort.isolation_metrics") {
        QStringList timeseries_list = MLUtil::toStringList(CLP.named_parameters["timeseries"]);
        QString firings = CLP.named_parameters["firings"].toString();
        QString metrics_out = CLP.named_parameters["metrics_out"].toString();
        QString pair_metrics_out = CLP.named_parameters["pair_metrics_out"].toString();
        P_isolation_metrics_opts opts;
        ret = p_isolation_metrics(timeseries_list, firings, metrics_out, pair_metrics_out, opts);
    }
    else if (arg1 == "mountainsort.combine_cluster_metrics") {
        QStringList metrics_list = MLUtil::toStringList(CLP.named_parameters["metrics_list"]);
        QString metrics_out = CLP.named_parameters["metrics_out"].toString();
        ret = p_combine_cluster_metrics(metrics_list, metrics_out);
    }
    else if (arg1 == "mountainsort.split_firings") {
        QStringList timeseries_list = MLUtil::toStringList(CLP.named_parameters["timeseries_list"]);
        QString firings = CLP.named_parameters["firings"].toString();
        QStringList firings_out_list = MLUtil::toStringList(CLP.named_parameters["firings_out_list"]);
        ret = p_split_firings(timeseries_list, firings, firings_out_list);
    }
    else if (arg1 == "mountainsort.concat_timeseries") {
        QStringList timeseries_list = MLUtil::toStringList(CLP.named_parameters["timeseries_list"]);
        QString timeseries_out = CLP.named_parameters["timeseries_out"].toString();
        ret = p_concat_timeseries(timeseries_list, timeseries_out);
    }
    else if (arg1 == "mountainsort.concat_firings") {
        QStringList timeseries_list = MLUtil::toStringList(CLP.named_parameters["timeseries_list"]);
        QStringList firings_list = MLUtil::toStringList(CLP.named_parameters["firings_list"]);
        QString timeseries_out = CLP.named_parameters["timeseries_out"].toString();
        QString firings_out = CLP.named_parameters["firings_out"].toString();
        ret = p_concat_firings(timeseries_list, firings_list, timeseries_out, firings_out);
    }
    else if (arg1 == "mountainsort.concat_event_times") {
        QStringList event_times_list = MLUtil::toStringList(CLP.named_parameters["event_times_list"]);
        QString event_times_out = CLP.named_parameters["event_times_out"].toString();
        ret = p_concat_event_times(event_times_list, event_times_out);
    }
    else if (arg1 == "mountainsort.load_test") {
        QString stats_out = CLP.named_parameters["stats_out"].toString();
        P_load_test_opts opts;
        opts.num_cpu_ops = CLP.named_parameters["num_cpu_ops"].toDouble();
        opts.num_read_bytes = CLP.named_parameters["num_read_bytes"].toDouble();
        opts.num_write_bytes = CLP.named_parameters["num_write_bytes"].toDouble();
        ret = p_load_test(stats_out, opts);
    }
    else if (arg1 == "mountainsort.misc_test") {
        QString dir = CLP.named_parameters["dir"].toString();
        QString info_out = CLP.named_parameters["info_out"].toString();
        ret = p_misc_test(dir, info_out);
    }
    else if (arg1 == "mountainsort.compute_amplitudes") {
        P_compute_amplitudes_opts opts;
        QString timeseries = CLP.named_parameters["timeseries"].toString();
        QString event_times = CLP.named_parameters["event_times"].toString();
        opts.central_channel = CLP.named_parameters["central_channel"].toInt();
        QString amplitudes_out = CLP.named_parameters["amplitudes_out"].toString();
        ret = p_compute_amplitudes(timeseries, event_times, amplitudes_out, opts);
    }
    else if (arg1 == "mountainsort.extract_time_interval") {
        P_compute_amplitudes_opts opts;
        QStringList timeseries_list = MLUtil::toStringList(CLP.named_parameters["timeseries_list"]);
        QString firings = CLP.named_parameters["firings"].toString();
        QString timeseries_out = CLP.named_parameters["timeseries_out"].toString();
        QString firings_out = CLP.named_parameters["firings_out"].toString();
        bigint t1 = CLP.named_parameters["t1"].toDouble();
        bigint t2 = CLP.named_parameters["t2"].toDouble();
        ret = p_extract_time_interval(timeseries_list, firings, timeseries_out, firings_out, t1, t2);
    }
    else if (arg1 == "mountainsort.confusion_matrix") {
        P_confusion_matrix_opts opts;
        QString firings1 = CLP.named_parameters["firings1"].toString();
        QString firings2 = CLP.named_parameters["firings2"].toString();
        QString confusion_matrix_out = CLP.named_parameters["confusion_matrix_out"].toString();
        QString matched_firings_out = CLP.named_parameters.value("matched_firings_out").toString();
        QString label_map_out = CLP.named_parameters.value("label_map_out").toString();
        QString firings2_relabeled_out = CLP.named_parameters.value("firings2_relabeled_out").toString();
        QString firings2_relabel_map_out = CLP.named_parameters.value("firings2_relabel_map_out").toString();
        if (CLP.named_parameters.contains("max_matching_offset")) {
            opts.max_matching_offset = CLP.named_parameters.value("max_matching_offset").toInt();
        }
        opts.relabel_firings2 = (CLP.named_parameters.value("relabel_firings2", "false").toString() == "true");
        ret = p_confusion_matrix(firings1, firings2, confusion_matrix_out, matched_firings_out, label_map_out, firings2_relabeled_out, firings2_relabel_map_out, opts);
    }
    else if (arg1 == "mountainsort.generate_background_dataset") {
        QString timeseries = CLP.named_parameters["timeseries"].toString();
        QString event_times = CLP.named_parameters["event_times"].toString();
        QString timeseries_out = CLP.named_parameters["timeseries_out"].toString();
        P_generate_background_dataset_opts opts;
        ret = p_generate_background_dataset(timeseries, event_times, timeseries_out, opts);
    }
    else {
        qWarning() << "Unexpected processor name: " + arg1;
        return -1;
    }

    if (!ret)
        return -1;

    return 0;
}

QJsonObject ProcessorSpecFile::get_spec()
{
    QJsonObject ret;
    ret["name"] = name;
    ret["description"] = description;
    ret["optional"] = optional;
    return ret;
}

QJsonObject ProcessorSpecParam::get_spec()
{
    QJsonObject ret;
    ret["name"] = name;
    ret["description"] = description;
    ret["optional"] = optional;
    if (default_value.isValid())
        ret["default_value"] = QJsonValue::fromVariant(default_value);
    return ret;
}

void ProcessorSpec::addInputs(QString name1, QString name2, QString name3, QString name4, QString name5)
{
    addInput(name1);
    if (!name2.isEmpty())
        addInput(name2);
    if (!name3.isEmpty())
        addInput(name3);
    if (!name4.isEmpty())
        addInput(name4);
    if (!name5.isEmpty())
        addInput(name5);
}

void ProcessorSpec::addOptionalInputs(QString name1, QString name2, QString name3, QString name4, QString name5)
{
    addInput(name1, "", true);
    if (!name2.isEmpty())
        addInput(name2, "", true);
    if (!name3.isEmpty())
        addInput(name3, "", true);
    if (!name4.isEmpty())
        addInput(name4, "", true);
    if (!name5.isEmpty())
        addInput(name5, "", true);
}

void ProcessorSpec::addOutputs(QString name1, QString name2, QString name3, QString name4, QString name5)
{
    addOutput(name1);
    if (!name2.isEmpty())
        addOutput(name2);
    if (!name3.isEmpty())
        addOutput(name3);
    if (!name4.isEmpty())
        addOutput(name4);
    if (!name5.isEmpty())
        addOutput(name5);
}

void ProcessorSpec::addOptionalOutputs(QString name1, QString name2, QString name3, QString name4, QString name5)
{
    addOutput(name1, "", true);
    if (!name2.isEmpty())
        addOutput(name2, "", true);
    if (!name3.isEmpty())
        addOutput(name3, "", true);
    if (!name4.isEmpty())
        addOutput(name4, "", true);
    if (!name5.isEmpty())
        addOutput(name5, "", true);
}

void ProcessorSpec::addRequiredParameters(QString name1, QString name2, QString name3, QString name4, QString name5)
{
    addRequiredParameter(name1);
    if (!name2.isEmpty())
        addRequiredParameter(name2);
    if (!name3.isEmpty())
        addRequiredParameter(name3);
    if (!name4.isEmpty())
        addRequiredParameter(name4);
    if (!name5.isEmpty())
        addRequiredParameter(name5);
}

void ProcessorSpec::addOptionalParameters(QString name1, QString name2, QString name3, QString name4, QString name5)
{
    addOptionalParameter(name1);
    if (!name2.isEmpty())
        addOptionalParameter(name2);
    if (!name3.isEmpty())
        addOptionalParameter(name3);
    if (!name4.isEmpty())
        addOptionalParameter(name4);
    if (!name5.isEmpty())
        addOptionalParameter(name5);
}

void ProcessorSpec::addInput(QString name, QString description, bool optional)
{
    ProcessorSpecFile X;
    X.name = name;
    X.description = description;
    X.optional = optional;
    inputs.append(X);
}

void ProcessorSpec::addOutput(QString name, QString description, bool optional)
{
    ProcessorSpecFile X;
    X.name = name;
    X.description = description;
    X.optional = optional;
    outputs.append(X);
}

void ProcessorSpec::addRequiredParameter(QString name, QString description)
{
    ProcessorSpecParam X;
    X.name = name;
    X.description = description;
    X.optional = false;
    parameters.append(X);
}

void ProcessorSpec::addOptionalParameter(QString name, QString description, QVariant default_value)
{
    ProcessorSpecParam X;
    X.name = name;
    X.description = description;
    X.optional = true;
    X.default_value = default_value;
    parameters.append(X);
}

QJsonObject ProcessorSpec::get_spec()
{
    QJsonObject ret;
    ret["name"] = processor_name;
    ret["version"] = version;
    ret["description"] = description;
    QJsonArray inputs0;
    for (int i = 0; i < inputs.count(); i++) {
        inputs0.push_back(inputs[i].get_spec());
    }
    ret["inputs"] = inputs0;
    QJsonArray outputs0;
    for (int i = 0; i < outputs.count(); i++) {
        outputs0.push_back(outputs[i].get_spec());
    }
    ret["outputs"] = outputs0;
    QJsonArray parameters0;
    for (int i = 0; i < parameters.count(); i++) {
        parameters0.push_back(parameters[i].get_spec());
    }
    ret["parameters"] = parameters0;
    ret["exe_command"] = qApp->applicationFilePath() + " " + processor_name + " $(arguments)";
    return ret;
}
