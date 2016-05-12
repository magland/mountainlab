QT += core
QT -= gui
#CONFIG += release
#CONFIG -= debug
CONFIG -= app_bundle #Please apple, don't make a bundle today :)

QT += qml

DESTDIR = ../bin
OBJECTS_DIR = ../build
MOC_DIR=../build
TARGET = mountainsort
TEMPLATE = app

QMAKE_POST_LINK += cp $$PWD/../bin/mountainsort $$PWD/../../mountainprocess/processors/mountainsort.mp

INCLUDEPATH += utils core processors mda unit_tests 3rdparty isosplit ../../common/cachemanager

HEADERS += \
    utils/get_command_line_params.h \
    core/msprocessmanager.h \
    core/mountainsort_version.h \
    utils/textfile.h \
    core/msprocessor.h \
    processors/example_processor.h \
    processors/bandpass_filter_processor.h \
    processors/bandpass_filter0.h \
    unit_tests/unit_tests.h \
    msprefs.h \
    processors/detect_processor.h \
    processors/detect.h \
    processors/whiten_processor.h \
    processors/whiten.h \
    utils/get_sort_indices.h \
    utils/eigenvalue_decomposition.h \
    utils/matrix_mda.h \
    processors/branch_cluster_v2_processor.h \
    processors/branch_cluster_v2.h \
    isosplit/isosplit2.h \
    isosplit/isocut.h \
    isosplit/jisotonic.h \
    processors/extract_clips.h \
    utils/msmisc.h \
    processors/remove_duplicate_clusters_processor.h \
    processors/remove_duplicate_clusters.h \
    processors/compute_outlier_scores_processor.h \
    processors/compute_outlier_scores.h \
    processors/copy_processor.h \
    processors/mda2txt_processor.h \
    processors/mask_out_artifacts_processor.h \
    processors/mask_out_artifacts.h \
    utils/get_pca_features.h \
    processors/adjust_times_processor.h \
    processors/adjust_times.h \
    processors/fit_stage_processor.h \
    processors/fit_stage.h \
    utils/compute_templates_0.h \
    processors/compute_templates.h \
    processors/compute_templates_processor.h \
    processors/mv_firings_filter_processor.h \
    processors/mv_firings_filter.h \
    processors/extract_clips_processor.h \
    processors/mv_subfirings_processor.h \
    ../../common/cachemanager/cachemanager.h \
    processors/extract_clips_features_processor.h \
    processors/compute_detectability_scores_processor.h \
    processors/compute_detectability_scores.h \
    processors/merge_labels_processor.h \
    processors/merge_labels.h \
    processors/filter_events.h \
    processors/filter_events_processor.h \
    processors/confusion_matrix.h \
    core/msscriptcontroller.h \
    processors/extract_raw_processor.h \
    processors/merge_across_channels.h \
    processors/merge_across_channels_processor.h \
    processors/geom2adj_processor.h

SOURCES += utils/get_command_line_params.cpp \
    core/msprocessmanager.cpp \
    core/mountainsort_version.cpp \
    utils/textfile.cpp \
    core/msprocessor.cpp \
    processors/example_processor.cpp \
    processors/bandpass_filter_processor.cpp \
    processors/bandpass_filter0.cpp \
    unit_tests/unit_tests.cpp \
    processors/detect_processor.cpp \
    processors/detect.cpp \
    processors/whiten_processor.cpp \
    processors/whiten.cpp \
    utils/get_sort_indices.cpp \
    utils/eigenvalue_decomposition.cpp \
    utils/matrix_mda.cpp \
    processors/branch_cluster_v2_processor.cpp \
    processors/branch_cluster_v2.cpp \
    isosplit/isosplit2.cpp \
    isosplit/isocut.cpp \
    isosplit/jisotonic.cpp \
    processors/extract_clips.cpp \
    utils/msmisc.cpp \
    processors/remove_duplicate_clusters_processor.cpp \
    processors/remove_duplicate_clusters.cpp \
    processors/compute_outlier_scores_processor.cpp \
    processors/compute_outlier_scores.cpp \
    processors/copy_processor.cpp \
    processors/mda2txt_processor.cpp \
    processors/mask_out_artifacts_processor.cpp \
    processors/mask_out_artifacts.cpp \
    utils/get_pca_features.cpp \
    processors/adjust_times_processor.cpp \
    processors/adjust_times.cpp \
    processors/fit_stage_processor.cpp \
    processors/fit_stage.cpp \
    utils/compute_templates_0.cpp \
    processors/compute_templates.cpp \
    processors/compute_templates_processor.cpp \
    processors/mv_firings_filter_processor.cpp \
    processors/mv_firings_filter.cpp \
    processors/extract_clips_processor.cpp \
    processors/mv_subfirings_processor.cpp \
    ../../common/cachemanager/cachemanager.cpp \
    processors/extract_clips_features_processor.cpp \
    processors/compute_detectability_scores_processor.cpp \
    processors/compute_detectability_scores.cpp \
    processors/merge_labels_processor.cpp \
    processors/merge_labels.cpp \
    processors/filter_events.cpp \
    processors/filter_events_processor.cpp \
    processors/confusion_matrix.cpp \
    processors/confusion_matrix_processor.cpp \
    core/msscriptcontroller.cpp \
    processors/extract_raw_processor.cpp \
    processors/fit_stage_new.cpp \
    processors/merge_across_channels.cpp \
    processors/merge_across_channels_processor.cpp \
    processors/geom2adj_processor.cpp

DEFINES += USE_REMOTE_MDA
INCLUDEPATH += ../../common/mda
DEPENDPATH += ../../common/mda
VPATH += ../../common/mda
HEADERS += remotereadmda.h diskreadmda.h diskwritemda.h usagetracking.h mda.h mdaio.h
SOURCES += remotereadmda.cpp diskreadmda.cpp diskwritemda.cpp usagetracking.cpp mda.cpp mdaio.cpp

HEADERS += utils/get_principal_components.h
SOURCES += utils/get_principal_components.cpp

INCLUDEPATH += ../../common
DEPENDPATH += ../../common
VPATH += ../../common
HEADERS += mlutils.h
SOURCES += mlutils.cpp

DISTFILES += \
    ../version.txt


#LAPACK
#On Ubuntu: sudo apt-get install liblapacke-dev
#On CentOS: sudo yum install lapack-devel.i686
INCLUDEPATH += /usr/include/lapacke #this was needed on CentOS
DEFINES += USE_LAPACK
LIBS += -llapack -llapacke

#FFTW
LIBS += -fopenmp -lfftw3 -lfftw3_threads

#OPENMP
!macx {
  QMAKE_LFLAGS += -fopenmp
  QMAKE_CXXFLAGS += -fopenmp
}
#-std=c++11   # AHB removed since not in GNU gcc 4.6.3

#QJSON
HEADERS += 3rdparty/qjson.h
SOURCES += 3rdparty/qjson.cpp
INCLUDEPATH += 3rdparty/qjson
DEPENDPATH += 3rdparty/qjson
VPATH += 3rdparty/qjson
HEADERS += serializer.h serializerrunnable.h parser.h parserrunnable.h json_scanner.h json_parser.hh
SOURCES += serializer.cpp serializerrunnable.cpp parser.cpp parserrunnable.cpp json_scanner.cpp json_parser.cc

#tests
test {
    QT += testlib
    CONFIG += testcase
    TARGET = mountainsort_test
    DEPENDPATH += unit_tests
    SOURCES += unit_tests/testMda.cpp	\
	unit_tests/testMain.cpp	\
	unit_tests/testMdaIO.cpp
    HEADERS += unit_tests/testMda.h \
	unit_tests/testMdaIO.h
} else {
    SOURCES += mountainsortmain.cpp
}
