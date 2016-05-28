QT += core
QT -= gui
CONFIG -= app_bundle #Please apple, don't make a bundle today :)

QT += qml

DESTDIR = ../bin
OBJECTS_DIR = ../build
MOC_DIR=../build
TARGET = mountainprocess
TEMPLATE = app

INCLUDEPATH += utils core mda unit_tests 3rdparty

INCLUDEPATH += ../../common/mda
VPATH += ../../common/mda
HEADERS += mda.h mdaio.h usagetracking.h diskreadmda.h diskwritemda.h \
    mpdaemon.h \
    mpdaemoninterface.h
SOURCES += mda.cpp mdaio.cpp usagetracking.cpp diskreadmda.cpp diskwritemda.cpp \
    mpdaemon.cpp \
    mpdaemoninterface.cpp

INCLUDEPATH += ../../common
DEPENDPATH += ../../common
VPATH += ../../common
HEADERS += mlutils.h
SOURCES += mlutils.cpp

#DEFINES += USE_REMOTE_MDA
#HEADERS += remotereadmda.h
#SOURCES += remotereadmda.cpp

INCLUDEPATH += ../../common/cachemanager
VPATH += ../../common/cachemanager
HEADERS += cachemanager.h
SOURCES += cachemanager.cpp

INCLUDEPATH += ../../common/commandlineparams
VPATH += ../../common/commandlineparams
HEADERS += commandlineparams.h
SOURCES += commandlineparams.cpp

INCLUDEPATH += ../../common/utils
VPATH += ../../common/utils
HEADERS += textfile.h haltagent.h taskprogress.h
SOURCES += textfile.cpp haltagent.cpp taskprogress.cpp

HEADERS += \
    processmanager.h \
    scriptcontroller.h \
    unit_tests/unit_tests.h

SOURCES += \
    processmanager.cpp \
    scriptcontroller.cpp \
    unit_tests/unit_tests.cpp

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
    SOURCES += mountainprocessmain.cpp
}
