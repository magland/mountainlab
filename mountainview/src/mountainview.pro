QT += core gui script network

DEFINES += USE_NETWORK
CONFIG += c++11

CONFIG -= app_bundle #Please apple, don't make a bundle today

#we are no longer supporting qt4 for mountainview
#greaterThan(QT_MAJOR_VERSION, 4): QT += widgets
QT += widgets

DESTDIR = ../bin
OBJECTS_DIR = ../build
MOC_DIR=../build
TARGET = mountainview
TEMPLATE = app

HEADERS += \
    histogramview.h \
    mvstatisticswidget.h \
    mvlabelcomparewidget.h \
    mvcdfview.h \
    static_neuron_colors.h \
    mvutils.h \
    mvoverview2widget.h \
    correlationmatrixview.h \
    mvoverview2widgetcontrolpanel.h \
    mvclusterdetailwidget.h \
    mvclipsview.h \
    mvclusterview.h \
    mvclusterwidget.h \
    mvfiringeventview.h \
    diskarraymodel_new.h \
    sstimeserieswidget.h \
    sstimeseriesview.h \
    sstimeseriesplot.h \
    ssabstractview.h \
    ssabstractplot.h \
    plotarea.h \
    sslabelsmodel.h \
    sslabelsmodel1.h \
    sscommon.h \
    tabber.h \
    imagesavedialog.h \
    closemehandler.h \
    computationthread.h \
    set_progress.h \
    mvcrosscorrelogramswidget2.h \
    mountainprocessrunner.h \
    mvclipswidget.h \
    taskprogressview.h \
    mvcontrolpanel.h \
    flowlayout.h \
    clustermerge.h \
    mvviewagent.h \
    mvstatusbar.h
SOURCES += mountainviewmain.cpp \
    histogramview.cpp \
    mvstatisticswidget.cpp \
    mvlabelcomparewidget.cpp \
    mvcdfview.cpp \
    static_neuron_colors.cpp \
    mvutils.cpp \
    mvoverview2widget.cpp \
    correlationmatrixview.cpp \
    mvoverview2widgetcontrolpanel.cpp \
    mvclusterdetailwidget.cpp \
    mvclipsview.cpp \
    mvclusterview.cpp \
    mvclusterwidget.cpp \
    mvfiringeventview.cpp \
    diskarraymodel_new.cpp \
    sstimeserieswidget.cpp \
    sstimeseriesview.cpp \
    sstimeseriesplot.cpp \
    ssabstractview.cpp \
    ssabstractplot.cpp \
    plotarea.cpp \
    sslabelsmodel1.cpp \
    sscommon.cpp \
    tabber.cpp \
    imagesavedialog.cpp \
    closemehandler.cpp \
    computationthread.cpp \
    set_progress.cpp \
    mvcrosscorrelogramswidget2.cpp \
    mountainprocessrunner.cpp \
    mvclipswidget.cpp \
    taskprogressview.cpp \
    mvcontrolpanel.cpp \
    flowlayout.cpp \
    clustermerge.cpp \
    mvviewagent.cpp \
    mvstatusbar.cpp

INCLUDEPATH += ../../common/commandlineparams
VPATH += ../../common/commandlineparams
HEADERS += commandlineparams.h
SOURCES += commandlineparams.cpp

INCLUDEPATH += ../../common/utils
VPATH += ../../common/utils
HEADERS += textfile.h haltagent.h taskprogress.h toolbuttonmenu.h
SOURCES += textfile.cpp haltagent.cpp taskprogress.cpp toolbuttonmenu.cpp

INCLUDEPATH += ../../mountainsort/src/utils
DEPENDPATH += ../../mountainsort/src/utils
VPATH += ../../mountainsort/src/utils
HEADERS += get_sort_indices.h msmisc.h
SOURCES += get_sort_indices.cpp msmisc.cpp
HEADERS += get_pca_features.h get_principal_components.h eigenvalue_decomposition.h
SOURCES += get_pca_features.cpp get_principal_components.cpp eigenvalue_decomposition.cpp
HEADERS += affinetransformation.h
SOURCES += affinetransformation.cpp
HEADERS += compute_templates_0.h
SOURCES += compute_templates_0.cpp

DEFINES += USE_REMOTE_MDA
INCLUDEPATH += ../../common/mda
DEPENDPATH += ../../common/mda
VPATH += ../../common/mda
HEADERS += remotereadmda.h diskreadmda.h diskwritemda.h usagetracking.h mda.h mdaio.h
SOURCES += remotereadmda.cpp diskreadmda.cpp diskwritemda.cpp usagetracking.cpp mda.cpp mdaio.cpp

INCLUDEPATH += ../../mountainsort/src/processors
DEPENDPATH += ../../mountainsort/src/processors
VPATH += ../../mountainsort/src/processors
HEADERS += extract_clips.h
SOURCES += extract_clips.cpp

INCLUDEPATH += ../../common/cachemanager
DEPENDPATH += ../../common/cachemanager
VPATH += ../../common/cachemanager
HEADERS += cachemanager.h
SOURCES += cachemanager.cpp

INCLUDEPATH += ../../common
DEPENDPATH += ../../common
VPATH += ../../common
HEADERS += mlutils.h
SOURCES += mlutils.cpp

RESOURCES += mountainview.qrc

#LAPACK
#On Ubuntu: sudo apt-get install liblapacke-dev
#On CentOS: sudo yum install lapack-devel.i686
#INCLUDEPATH += /usr/include/lapacke #this was needed on CentOS
#DEFINES += USE_LAPACK
#LIBS += -llapack -llapacke

#OPENMP
!macx {
  QMAKE_LFLAGS += -fopenmp
  QMAKE_CXXFLAGS += -fopenmp
}
#-std=c++11   # AHB removed since not in GNU gcc 4.6.3

DISTFILES += \
    ../../server/mountainviewserver/mountainviewserver.py \
    ../../server/mscmdserver/mscmdserver.py \
    ../../server/mdaserver/mdaserver.py \
    ../../server/mountainviewserver/mountainviewserver.example.cfg \
    ../../server/mdaserver/mdaserver.example.cfg \
    ../../server/mscmdserver/mscmdserver.example.cfg
