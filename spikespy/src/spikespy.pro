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
TARGET = spikespy
TEMPLATE = app

HEADERS += \
    sscontroller.h \
    sslabelview.h \
    sstimeseriesplot_prev.h \
    closemehandler.h \
    sstimeserieswidget_prev.h \
    sscommon.h \
    ssabstractview_prev.h \
    ssabstractplot_prev.h \
    diskarraymodel_prev.h \
    sslabelsmodel1_prev.h \
    memorymda_prev.h \
    sstimeseriesview_prev.h \
    mdaobject.h \
    sslabelsmodel_prev.h \
    diskreadmdaold.h \
    sslabelplot_prev.h \
    diskwritemdaold.h \
    extractclipsdialog.h \
    plotarea_prev.h
SOURCES += ssmain.cpp \
    sscontroller.cpp \
    sslabelview.cpp \
    sstimeseriesplot_prev.cpp \
    closemehandler.cpp \
    sstimeserieswidget_prev.cpp \
    sscommon.cpp \
    ssabstractview_prev.cpp \
    ssabstractplot_prev.cpp \
    diskarraymodel_prev.cpp \
    sslabelsmodel1_prev.cpp \
    memorymda_prev.cpp \
    sstimeseriesview_prev.cpp \
    mdaobject.cpp \
    diskreadmdaold.cpp \
    sslabelplot_prev.cpp \
    diskwritemdaold.cpp \
    extractclipsdialog.cpp \
    plotarea_prev.cpp

INCLUDEPATH += ../../common
DEPENDPATH += ../../common
VPATH += ../../common
HEADERS += mlutils.h
SOURCES += mlutils.cpp

INCLUDEPATH += ../../common/commandlineparams
VPATH += ../../common/commandlineparams
HEADERS += commandlineparams.h
SOURCES += commandlineparams.cpp

INCLUDEPATH += ../../common/utils
VPATH += ../../common/utils
HEADERS += textfile.h taskprogress.h
SOURCES += textfile.cpp taskprogress.cpp

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

INCLUDEPATH += ../../common/cachemanager
DEPENDPATH += ../../common/cachemanager
VPATH += ../../common/cachemanager
HEADERS += cachemanager.h
SOURCES += cachemanager.cpp

INCLUDEPATH += ../../mountainsort/src/processors
DEPENDPATH += ../../mountainsort/src/processors
VPATH += ../../mountainsort/src/processors
HEADERS += extract_clips.h
SOURCES += extract_clips.cpp

INCLUDEPATH += ../../mountainsort/src/core
DEPENDPATH += ../../mountainsort/src/core
VPATH += ../../mountainsort/src/core
HEADERS += mscachemanager.h
SOURCES += mscachemanager.cpp

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

FORMS += \
    extractclipsdialog.ui

