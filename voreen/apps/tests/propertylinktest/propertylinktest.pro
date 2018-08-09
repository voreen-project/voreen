TARGET = propertylinktest
TEMPLATE = app
LANGUAGE = C++

CONFIG -= qt
CONFIG += console

# Include local configuration
include(../../config.tex)

# Include common configuration
include(../../commonconf.pri)

include(../voreenapp.pri)

# HEADERS += ../../ext/tgt/qt/qtcanvas.h

SOURCES += propertylinktest.cpp 

unix {
  DEFINES += LINUX
}
