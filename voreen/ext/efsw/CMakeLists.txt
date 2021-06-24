################################################################################
# Project file for the efsw library
################################################################################
PROJECT(efsw)
CMAKE_MINIMUM_REQUIRED(VERSION 3.5.1 FATAL_ERROR)
INCLUDE(../../cmake/commonconf.cmake)

MESSAGE(STATUS "Configuring efsw Library")

# sources
SET(EFSW_SOURCES
    src/efsw/Debug.cpp
    src/efsw/DirectorySnapshot.cpp
    src/efsw/DirectorySnapshotDiff.cpp
    src/efsw/DirWatcherGeneric.cpp
    src/efsw/FileInfo.cpp
    src/efsw/FileSystem.cpp
    src/efsw/FileWatcher.cpp
    src/efsw/FileWatcherCWrapper.cpp
    src/efsw/FileWatcherFSEvents.cpp
    src/efsw/FileWatcherGeneric.cpp
    src/efsw/FileWatcherImpl.cpp
    src/efsw/FileWatcherInotify.cpp
    src/efsw/FileWatcherKqueue.cpp
    src/efsw/FileWatcherWin32.cpp
    src/efsw/Log.cpp
    src/efsw/Mutex.cpp
    src/efsw/String.cpp
    src/efsw/System.cpp
    src/efsw/Thread.cpp
    src/efsw/Watcher.cpp
    src/efsw/WatcherFSEvents.cpp
    src/efsw/WatcherGeneric.cpp
    src/efsw/WatcherInotify.cpp
    src/efsw/WatcherKqueue.cpp
    src/efsw/WatcherWin32.cpp
)
    
if (WIN32)
    list (APPEND EFSW_SOURCES
          src/efsw/platform/win/FileSystemImpl.cpp
          src/efsw/platform/win/MutexImpl.cpp
          src/efsw/platform/win/SystemImpl.cpp
          src/efsw/platform/win/ThreadImpl.cpp)
else ()
    list (APPEND EFSW_SOURCES
          src/efsw/platform/posix/FileSystemImpl.cpp
          src/efsw/platform/posix/MutexImpl.cpp
          src/efsw/platform/posix/SystemImpl.cpp
          src/efsw/platform/posix/ThreadImpl.cpp)
endif()          

if (APPLE)
    list (REMOVE_ITEM EFSW_SOURCES
          "src/efsw/WatcherInotify.cpp"
          "src/efsw/WatcherWin32.cpp"
          "src/efsw/FileWatcherInotify.cpp"
          "src/efsw/FileWatcherWin32.cpp")

    exec_program(uname ARGS -v  OUTPUT_VARIABLE OSX_VERSION)
       string(REGEX MATCH "[0-9]+" OSX_VERSION ${OSX_VERSION})
       if (NOT OSX_VERSION GREATER 9)
        add_definitions(-DEFSW_FSEVENTS_NOT_SUPPORTED)
       endif()  
elseif (WIN32)
    list (REMOVE_ITEM EFSW_SOURCES 
          "src/efsw/WatcherKqueue.cpp"
          "src/efsw/WatcherFSEvents.cpp"
          "src/efsw/WatcherInotify.cpp"
          "src/efsw/FileWatcherKqueue.cpp"
          "src/efsw/FileWatcherInotify.cpp"
          "src/efsw/FileWatcherFSEvents.cpp")
elseif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    list (REMOVE_ITEM EFSW_SOURCES 
          "src/efsw/WatcherKqueue.cpp"
          "src/efsw/WatcherFSEvents.cpp"
          "src/efsw/WatcherWin32.cpp"
          "src/efsw/FileWatcherKqueue.cpp"
          "src/efsw/FileWatcherWin32.cpp"
          "src/efsw/FileWatcherFSEvents.cpp")

    if (NOT EXISTS "/usr/include/sys/inotify.h" AND NOT EXISTS "/usr/local/include/sys/inotify.h")
        add_definitions(-DEFSW_INOTIFY_NOSYS)
    endif()
elseif (${CMAKE_SYSTEM_NAME} MATCHES "FreeBSD")
    list (REMOVE_ITEM EFSW_SOURCES
          "src/efsw/WatcherInotify.cpp"
          "src/efsw/WatcherWin32.cpp"
          "src/efsw/WatcherFSEvents.cpp"
          "src/efsw/FileWatcherInotify.cpp"
          "src/efsw/FileWatcherWin32.cpp"
          "src/efsw/FileWatcherFSEvents.cpp")
endif()


################################################################################
# define library target
################################################################################
ADD_LIBRARY(efsw ${EFSW_SOURCES})
if (MSVC)
    add_definitions(-D_SCL_SECURE_NO_WARNINGS -DEFSW_DYNAMIC -DEFSW_EXPORTS)
else ()
    add_definitions(-Wall -Wno-long-long)
endif()

if (${CMAKE_BUILD_TYPE} MATCHES "Debug")
    add_definitions(-DDEBUG)
elseif (${CMAKE_BUILD_TYPE} MATCHES "Release")
    add_definitions(-DNDEBUG)
endif()

if (APPLE)
    set(MAC_LIBS "-framework CoreFoundation" "-framework CoreServices")
    target_link_libraries(efsw ${MAC_LIBS})
elseif (NOT (${CMAKE_SYSTEM_NAME} MATCHES "Haiku") AND NOT WIN32)
    target_link_libraries(efsw pthread)
endif()

include_directories(include src)

# group source files by subdir
IF(VRN_GROUP_SOURCE_FILES)
    LIST(APPEND RemovePathPrefixes "ext/efsw")
    DEFINE_SOURCE_GROUPS_FROM_SUBDIR(EFSW_SOURCES ${VRN_HOME} RemovePathPrefixes)
ENDIF()

################################################################################
# deployment
################################################################################
IF(VRN_ADD_INSTALL_TARGET)
    INSTALL(TARGETS efsw 
        RUNTIME DESTINATION .
        LIBRARY DESTINATION .
    )
ENDIF()

