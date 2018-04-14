####################################################
# Project file for the Voreen voltool app
####################################################
TARGET = voltool
TEMPLATE	= app
LANGUAGE	= C++

# Include local configuration
#include(../../config.txt)

# Include common configuration
#include(../../commonconf.pri)

#include(../voreenapp.pri)

CONFIG += console

INCLUDEPATH += ../../ext
INCLUDEPATH += ../../include

contains(DEFINES, VRN_WITH_DEVIL) {
  LIBS += -lILU
}

DEFINES += VRN_NO_OPENGL
DEFINES += VRN_WITH_PVM

SOURCES	+=  voltool.cpp \
            ../../src/core/calc/histogram.cpp \
            ../../src/core/io/datvolumereader.cpp \
            ../../src/core/io/datvolumewriter.cpp \
            ../../src/core/io/ddsbase.cpp \
            ../../src/core/io/interfilevolumereader.cpp \
            ../../src/core/io/nrrdvolumereader.cpp \
            ../../src/core/io/nrrdvolumewriter.cpp \
            ../../src/core/io/quadhidacvolumereader.cpp \
            ../../src/core/io/pvmvolumereader.cpp \
            ../../src/core/io/rawvolumereader.cpp \
            ../../src/core/io/siemensreader.cpp \
            ../../src/core/io/textfilereader.cpp \
            ../../src/core/io/tuvvolumereader.cpp \
            ../../src/core/io/vevovolumereader.cpp \
            ../../src/core/io/volumeserializer.cpp \
            ../../src/core/io/volumeserializerpopulator.cpp \
            ../../src/core/io/volumereader.cpp \
            ../../src/core/io/volumewriter.cpp \
            ../../src/core/datastructures/volume/modality.cpp \
            ../../src/core/datastructures/volume/volumecontainer.cpp \
            ../../src/core/datastructures/volume/volume.cpp \
            ../../src/core/datastructures/volume/volumemetadata.cpp \
            ../../ext/tgt/logmanager.cpp

HEADERS +=  ../../include/voreen/core/calc/histogram.h \
            ../../include/voreen/core/io/datvolumereader.h \
            ../../include/voreen/core/io/datvolumewriter.h \
            ../../include/voreen/core/io/interfilevolumereader.h \
            ../../include/voreen/core/io/nrrdvolumereader.h \
            ../../include/voreen/core/io/nrrdvolumewriter.h \
            ../../include/voreen/core/io/quadhidacvolumereader.h \
            ../../include/voreen/core/io/pvmvolumereader.h \
            ../../include/voreen/core/io/rawvolumereader.h \
            ../../include/voreen/core/io/siemensreader.h \
            ../../include/voreen/core/io/textfilereader.h \
            ../../include/voreen/core/io/tuvvolumereader.h \
            ../../include/voreen/core/io/vevovolumereader.h \
            ../../include/voreen/core/io/volumeserializer.h \
            ../../include/voreen/core/io/volumeserializerpopulator.h \
            ../../include/voreen/core/io/volumereader.h \
            ../../include/voreen/core/io/volumewriter.h \
            ../../include/voreen/core/datastructures/volume/modality.h \
            ../../include/voreen/core/datastructures/volume/volumeatomic.h \
            ../../include/voreen/core/datastructures/volume/volumecontainer.h \
            ../../include/voreen/core/datastructures/volume/volume.h \
            ../../include/voreen/core/datastructures/volume/volumemetadata.h \
            ../../ext/tgt/logmanager.h

# add files, which are not available in the snapshot release
# these files will be added to the snapshot, when they are cleaned up
!contains(DEFINES, VRN_SNAPSHOT) {
SOURCES	+= vqhelper.cpp \
	       aopreprocess.cpp

HEADERS +=  vqhelper.h \
	        aopreprocess.h 
}
