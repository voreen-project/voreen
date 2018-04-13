####################################################
# Project file for .xml to .dat flow data converter
####################################################
TARGET      = xml2dat
TEMPLATE    = app
LANGUAGE    = C++

CONFIG += console
CONFIG -= qt

include(../../../../config.txt)

include($${VRN_HOME}/commonconf.pri)
include($${VRN_HOME}/apps/voreenapp.pri)

contains(DEFINES, VRN_PRECOMPILE_HEADER) {
  qt {
    PRECOMPILED_HEADER = $${VRN_HOME}/apps/pch_qtapp.h
  }
  else {
    PRECOMPILED_HEADER = $${VRN_HOME}/pch.h
  }
  CONFIG += precompile_header
}

SOURCES += xml2dat.cpp

