TARGET = opencltest
TEMPLATE = app
LANGUAGE = C++

CONFIG -= qt
CONFIG += console

# Include local configuration
include(../../config.txt)

# Include common configuration
include(../../commonconf.pri)

include(../voreenapp.pri)

SOURCES += opencltest.cpp

unix {
  DEFINES += LINUX

  LIBS += -lOpenCL
}
