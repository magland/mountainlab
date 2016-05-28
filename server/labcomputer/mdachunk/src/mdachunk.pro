#-------------------------------------------------
#
# Project created by QtCreator 2016-03-30T06:59:21
#
#-------------------------------------------------

QT       += core network

QT       -= gui

DESTDIR = ../bin
OBJECTS_DIR = ../build
MOC_DIR=../build
TARGET = mdachunk
CONFIG   += console
CONFIG   -= app_bundle

TEMPLATE = app

SOURCES += main.cpp

INCLUDEPATH += ../../../../common/commandlineparams
VPATH += ../../../../common/commandlineparams
HEADERS += commandlineparams.h
SOURCES += commandlineparams.cpp

INCLUDEPATH += ../../../../common/utils
VPATH += ../../../../common/utils
HEADERS += textfile.h taskprogress.h
SOURCES += textfile.cpp taskprogress.cpp

DEFINES += USE_REMOTE_MDA
INCLUDEPATH += ../../../../common/mda
DEPENDPATH += ../../../../common/mda
VPATH += ../../../../common/mda
HEADERS += remotereadmda.h diskreadmda.h diskwritemda.h usagetracking.h mda.h mdaio.h
SOURCES += remotereadmda.cpp diskreadmda.cpp diskwritemda.cpp usagetracking.cpp mda.cpp mdaio.cpp

INCLUDEPATH += ../../../../common/cachemanager
DEPENDPATH += ../../../../common/cachemanager
VPATH += ../../../../common/cachemanager
HEADERS += cachemanager.h
SOURCES += cachemanager.cpp

INCLUDEPATH += ../../../../common
DEPENDPATH += ../../../../common
VPATH += ../../../../common
HEADERS += mlutils.h
SOURCES += mlutils.cpp

INCLUDEPATH += ../../../../mountainsort/src/utils
DEPENDPATH += ../../../../mountainsort/src/utils
VPATH += ../../../../mountainsort/src/utils
HEADERS += msmisc.h
SOURCES += msmisc.cpp

INCLUDEPATH += ../../../../mountainsort/src/core
DEPENDPATH += ../../../../mountainsort/src/core
VPATH += ../../../../mountainsort/src/core
HEADERS += mscachemanager.h
SOURCES += mscachemanager.cpp


