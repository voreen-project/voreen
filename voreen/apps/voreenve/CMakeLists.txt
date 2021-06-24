################################################################################
# Project file for the VoreenVE Application
################################################################################
PROJECT(voreenve)
CMAKE_MINIMUM_REQUIRED(VERSION 3.5.1 FATAL_ERROR)
INCLUDE(../../cmake/commonconf.cmake)

MESSAGE(STATUS "Configuring VoreenVE Application")

# sources
SET(VRN_VE_SOURCES
    main.cpp
    voreenvemainwindow.cpp
    voreenvesplashscreen.cpp
    voreenveapplication.cpp    
    voreenvestartupwizard.cpp
    voreenveaboutbox.cpp
)

# headers
SET(VRN_VE_HEADERS
    voreenvemainwindow.h
    voreenvesplashscreen.h
    voreenveapplication.h
    voreenvestartupwizard.h
    voreenveaboutbox.h
)

# headers that do not need to be run through moc
SET(VRN_VE_HEADERS_NONMOC
)

# add module sources and headers
LIST(APPEND VRN_VE_SOURCES ${VRN_MODULE_VE_SOURCES})
LIST(APPEND VRN_VE_HEADERS ${VRN_MODULE_VE_HEADERS})
LIST(APPEND VRN_VE_HEADERS_NONMOC ${VRN_MODULE_VE_HEADERS_NONMOC})

# find Qt library
SET(VRN_QT_COMPONENTS Gui Widgets Svg)

IF(QT_USE_QTNETWORK)
    LIST(APPEND VRN_QT_COMPONENTS Network)
ENDIF()
FIND_PACKAGE(Qt5Core ${VRN_REQUIRED_QT_VERSION} REQUIRED)
FIND_PACKAGE(Qt5 COMPONENTS ${VRN_QT_COMPONENTS} REQUIRED)
# create list of libraries
set(QT_LIBRARIES Qt5::Core Qt5::Gui Qt5::Widgets Qt5::Svg ${QT5SVG_LIBRARIES})
set(QT_INCLUDE_DIRECTORIES ${Qt5Core_INCLUDES} ${Qt5Gui_INCLUDES} ${Qt5Widgets_INCLUDES}  ${Qt5Svg_INCLUDES})

IF(QT_USE_QTNETWORK)
    LIST(APPEND QT_LIBRARIES Qt5::Network)
    LIST(APPEND QT_INCLUDE_DIRECTORIES ${Qt5Network_INCLUDES})
ENDIF()

################################################################################
# define executable target
################################################################################

# Qt moc
QT5_WRAP_CPP(VRN_VE_MOC_SOURCES ${VRN_VE_HEADERS})
LIST(APPEND VRN_VE_SOURCES ${VRN_VE_MOC_SOURCES})

# Qt resources
QT5_ADD_RESOURCES(VRN_VE_RESOURCES_RCC "${VRN_HOME}/apps/voreenve/voreenve.qrc" ${VRN_MODULE_VE_RESOURCES})

# add msvc pch
IF (VRN_PRECOMPILED_HEADER AND MSVC)
    ADD_MSVC_PRECOMPILED_HEADER("pch_ve.h" "pch_ve.cpp" VRN_VE_SOURCES)
ENDIF()

# define target
IF(WIN32 AND VRN_DEPLOYMENT)
    MESSAGE(STATUS "* Linking against Windows subsystem (deployment build)")
    
    # WIN32 symbol causes the executable to link against the Windows subsystem (thus not showing console window)
    ADD_EXECUTABLE(voreenve WIN32
        ${VRN_VE_SOURCES} ${VRN_VE_HEADERS} ${VRN_VE_HEADERS_NONMOC} 
        ${VRN_VE_RESOURCES_RCC} 
        "${VRN_HOME}/resource/voreenve/icons/winicon.rc") 

    # explicitly define application entry point to prevent linker error when using Windows subsystem 
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /ENTRY:mainCRTStartup") 
ELSE()
    ADD_EXECUTABLE(voreenve
        ${VRN_VE_SOURCES} ${VRN_VE_HEADERS} ${VRN_VE_HEADERS_NONMOC} 
        ${VRN_VE_RESOURCES_RCC} 
        "${VRN_HOME}/resource/voreenve/icons/winicon.rc") 
ENDIF()
ADD_DEFINITIONS(${VRN_DEFINITIONS} ${VRN_MODULE_DEFINITIONS})
INCLUDE_DIRECTORIES(${VRN_INCLUDE_DIRECTORIES})
TARGET_LINK_LIBRARIES(voreenve tgt voreen_core voreen_qt ${VRN_EXTERNAL_LIBRARIES})

# group source files by subdir
IF(VRN_GROUP_SOURCE_FILES)
    LIST(APPEND RemovePathPrefixes "")
    DEFINE_SOURCE_GROUPS_FROM_SUBDIR(VRN_VE_SOURCES ${VRN_HOME} RemovePathPrefixes)
    DEFINE_SOURCE_GROUPS_FROM_SUBDIR(VRN_VE_HEADERS ${VRN_HOME} RemovePathPrefixes)
    DEFINE_SOURCE_GROUPS_FROM_SUBDIR(VRN_VE_HEADERS_NONMOC ${VRN_HOME} RemovePathPrefixes)
ENDIF()

# add Qt dependency
INCLUDE_DIRECTORIES(${QT_INCLUDE_DIRECTORIES})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}) #< location of generated ui headers in out-of-source builds
ADD_DEFINITIONS(${QT_DEFINITIONS})
TARGET_LINK_LIBRARIES(voreenve ${QT_LIBRARIES})

# add gcc pch
IF (VRN_PRECOMPILED_HEADER AND CMAKE_COMPILER_IS_GNUCXX)
    ADD_GCC_PRECOMPILED_HEADER(voreenve pch_ve.h)
ENDIF()


################################################################################
# deployment
################################################################################
IF(VRN_ADD_INSTALL_TARGET)
    INSTALL(TARGETS voreenve
        RUNTIME DESTINATION .
    )
ENDIF()    
