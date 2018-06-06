####################################################
# Project file for the Voreen voltool app
####################################################
TARGET = voltool
TEMPLATE = app
LANGUAGE = C++

CONFIG += console
CONFIG -= qt

# check qmake version
QMAKE_VERS = $$[QMAKE_VERSION]
QMAKE_VERSION_CHECK = $$find(QMAKE_VERS, "^[234]\\.")
isEmpty(QMAKE_VERSION_CHECK) {
   error("Your qmake version '$$QMAKE_VERS' is too old, qmake from Qt 4 is required!")
}

# include config
!exists(../../config.txt) {
  error("config.txt not found! copy config-default.txt to config.txt and edit!")
}
include(../../config.txt)

# Include common configuration
include(../../commonconf.pri)

# Set output directory of the executable
VRN_APP_DIRECTORY = "$${VRN_HOME}/bin"

# Include generic app configuration
include(../voreenapp.pri)

win32 {
 contains(DEFINES, VRN_WITH_DEVIL) {
     LIBS += "$${DEVIL_DIR}/lib/ILU.lib"
  }
}

unix {
  contains(DEFINES, VRN_WITH_DEVIL) {
    LIBS += -lILU
  }
}

SOURCES	+= voltool.cpp \
           commands_grad.cpp \
           commands_convert.cpp \
           commands_create.cpp \
           commands_modify.cpp \
           commands_registration.cpp

HEADERS +=  commands_grad.h \
            commands_convert.h \
            commands_create.h \
            commands_modify.h \
            commands_registration.h

exists(voltool-internal.pri) : include(voltool-internal.pri)
