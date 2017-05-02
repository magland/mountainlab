QT += core gui network

QT += qml

CONFIG += c++11
QMAKE_CXXFLAGS += -Wno-reorder #qaccordion

CONFIG -= app_bundle #Please apple, don't make a bundle today

include(../../../../mvcommon/mvcommon.pri)
include(../../../../mlcommon/mlcommon.pri)
include(../../../../mlcommon/mda.pri)
include(../../../../mlcommon/taskprogress.pri)

QT += widgets
QT+=concurrent

DESTDIR = ../../bin
OBJECTS_DIR = ../build
MOC_DIR=../build
TARGET = mountainview-eeg
TEMPLATE = app

INCLUDEPATH += plugins views
VPATH += plugins views
HEADERS += clustermetricsplugin.h clustermetricsview.h
SOURCES += clustermetricsplugin.cpp clustermetricsview.cpp

HEADERS += timeseriesplugin.h timeseriesview.h
SOURCES += timeseriesplugin.cpp timeseriesview.cpp
INCLUDEPATH += views/timeseriesview
VPATH += views/timeseriesview
HEADERS += mvtimeseriesrendermanager.h mvtimeseriesviewbase.h
SOURCES += mvtimeseriesrendermanager.cpp mvtimeseriesviewbase.cpp

HEADERS += spectrogramplugin.h spectrogramview.h
SOURCES += spectrogramplugin.cpp spectrogramview.cpp

INCLUDEPATH += msv/contextmenuhandlers
VPATH += msv/contextmenuhandlers

SOURCES += mveegmain.cpp

HEADERS += mveegcontext.h
SOURCES += mveegcontext.cpp

INCLUDEPATH += controlwidgets
VPATH += controlwidgets
HEADERS += openviewscontrol.h eegprefscontrol.h
SOURCES += openviewscontrol.cpp eegprefscontrol.cpp

INCLUDEPATH += ../../../../mountainview/src/multiscaletimeseries
VPATH += ../../../../mountainview/src/multiscaletimeseries
HEADERS += multiscaletimeseries.cpp
SOURCES += multiscaletimeseries.cpp

INCLUDEPATH += ../../../../prv-gui/src
HEADERS += ../../../../prv-gui/src/prvgui.h
SOURCES += ../../../../prv-gui/src/prvgui.cpp

#INCLUDEPATH += ../../../../mountainview/src/controlwidgets
#VPATH += ../../../../mountainview/src/controlwidgets
#HEADERS += resolveprvsdialog.h
#SOURCES += resolveprvsdialog.cpp
#FORMS += resolveprvsdialog.ui

#TODO: Do we need openmp for mountainview?
#OPENMP
!macx {
  QMAKE_LFLAGS += -fopenmp
  QMAKE_CXXFLAGS += -fopenmp
}
#-std=c++11   # AHB removed since not in GNU gcc 4.6.3

