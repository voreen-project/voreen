
# add files, which are not available in the snapshot release
# these files will be added to the snapshot, when they are cleaned up
# this file is not included in the snapshot release

SOURCES	+= commands_motion.cpp \
           commands_streaming.cpp \
           commands_seginfo.cpp

HEADERS +=  commands_motion.h \
            commands_seginfo.h \
            commands_streaming.h
