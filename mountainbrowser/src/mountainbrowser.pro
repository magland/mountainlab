QT += core gui
QT += webkitwidgets
QT -= gui
CONFIG -= app_bundle #Please apple, don't make a bundle today :)

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets #We do want to support Qt5, but there is no reason not to use Qt4

DEFINES += USE_NETWORK
CONFIG += c++11

DESTDIR = ../bin
OBJECTS_DIR = ../build
MOC_DIR=../build
TARGET = mountainbrowser
TEMPLATE = app

HEADERS += mountainbrowsermain.h \
    mbcontroller.h
SOURCES += mountainbrowsermain.cpp \
    mbcontroller.cpp

INCLUDEPATH += ../../mountainsort/src/utils
DEPENDPATH += ../../mountainsort/src/utils
VPATH += ../../mountainsort/src/utils
HEADERS += get_command_line_params.h textfile.h msmisc.h
SOURCES += get_command_line_params.cpp textfile.cpp msmisc.cpp

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

#QJSON
INCLUDEPATH += ../../mountainsort/src/3rdparty
DEPENDPATH += ../../mountainsort/src/3rdparty
VPATH += ../../mountainsort/src/3rdparty
INCLUDEPATH += ../../mountainsort/src/3rdparty/qjson
DEPENDPATH += ../../mountainsort/src/3rdparty/qjson
VPATH += ../../mountainsort/src/3rdparty/qjson
HEADERS += qjson.h
SOURCES += qjson.cpp
HEADERS += serializer.h serializerrunnable.h parser.h parserrunnable.h json_scanner.h json_parser.hh
SOURCES += serializer.cpp serializerrunnable.cpp parser.cpp parserrunnable.cpp json_scanner.cpp json_parser.cc

DEFINES += USE_REMOTE_MDA
INCLUDEPATH += ../../common/mda
DEPENDPATH += ../../common/mda
VPATH += ../../common/mda
HEADERS += remotereadmda.h diskreadmda.h diskwritemda.h usagetracking.h mda.h mdaio.h
SOURCES += remotereadmda.cpp diskreadmda.cpp diskwritemda.cpp usagetracking.cpp mda.cpp mdaio.cpp


DISTFILES += \
    experiments.json \
    html/style.css \
    franklab.json \
    ../../server/mountainbrowserserver/mountainbrowserserver.py \
    ../../server/mountainbrowserserver/mountainbrowserserver.example.cfg

RESOURCES += \
    mountainbrowser.qrc


