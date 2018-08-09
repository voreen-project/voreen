################################################################################
# Project file for the tgt library
################################################################################
PROJECT(tgt)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.0 FATAL_ERROR)
INCLUDE(../../cmake/commonconf.cmake)

MESSAGE(STATUS "Configuring tgt Library")

# configuration
IF(TGT_WITH_WMI)
    ADD_DEFINITIONS("-DTGT_WITH_WMI")
    ADD_DEFINITIONS("-D_WIN32_DCOM") #< for CoInitializeSecurity() to be defined
    LIST(APPEND VRN_EXTERNAL_LIBRARIES "WbemUuid")
ENDIF(TGT_WITH_WMI)    

# headers
FILE(GLOB TGT_HEADERS *.h event/*.h navigation/*.h)

# sources
SET(TGT_SOURCES
    assert.cpp
    bspline.cpp
    camera.cpp
    catmullromspline.cpp
    cpucapabilities.cpp
    curve.cpp
    exception.cpp
    filesystem.cpp
    font.cpp
    fontmanager.cpp
    framebufferobject.cpp
    frustum.cpp
    glcanvas.cpp
    glcontextbase.cpp
    glcontextmanager.cpp
    gpucapabilities.cpp
    init.cpp
    light.cpp
    logmanager.cpp
    matrixstack.cpp
    naturalcubicspline.cpp
    painter.cpp
    shadermanager.cpp
    spline.cpp
    stopwatch.cpp
    tesselator.cpp
    texture.cpp
    texture2darray.cpp
    texturemanager.cpp
    texturereader.cpp
    texturereaderdevil.cpp
    texturereadertga.cpp
    textureunit.cpp
    tgt_gl.cpp
    timer.cpp
#event
    event/eventhandler.cpp
    event/eventlistener.cpp
    event/keyevent.cpp
    event/mouseevent.cpp
    event/timeevent.cpp
    event/touchevent.cpp
    event/touchpoint.cpp
#navigation 
    navigation/navigation.cpp
    navigation/trackball.cpp
#immediate mode
    immediatemode/immediatemode.cpp
)

LIST(APPEND TGT_HEADERS
    immediatemode/immediatemode.h
    )

IF(NOT VRN_OPENGL_COMPATIBILITY_PROFILE)
    LIST(APPEND TGT_SOURCES
    gl_core.cpp
    immediatemode/gl.cpp
    )
    LIST(APPEND TGT_HEADERS
    immediatemode/gl.h
    )
ENDIF()

################################################################################
# define library target
################################################################################
ADD_LIBRARY(tgt ${TGT_SOURCES} ${TGT_HEADERS})
ADD_DEFINITIONS(${VRN_DEFINITIONS} ${VRN_MODULE_DEFINITIONS})
IF(MSVC)
    ADD_DEFINITIONS("-DTGT_BUILD_DLL")
ENDIF()
INCLUDE_DIRECTORIES(${VRN_INCLUDE_DIRECTORIES})
TARGET_LINK_LIBRARIES(tgt ${VRN_EXTERNAL_LIBRARIES})

# group source files by subdir
IF(VRN_GROUP_SOURCE_FILES)
    LIST(APPEND RemovePathPrefixes "ext/tgt")
    DEFINE_SOURCE_GROUPS_FROM_SUBDIR(TGT_SOURCES ${VRN_HOME} RemovePathPrefixes)
    DEFINE_SOURCE_GROUPS_FROM_SUBDIR(TGT_HEADERS ${VRN_HOME} RemovePathPrefixes)
ENDIF()

################################################################################
# deployment
################################################################################
IF(VRN_ADD_INSTALL_TARGET)
    INSTALL(TARGETS tgt
        RUNTIME DESTINATION .
        LIBRARY DESTINATION .
    )
ENDIF()

