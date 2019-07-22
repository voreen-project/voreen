# Try to find freeGLUT library and include path. Once done this will define:
# GLUT_FOUND
# GLUT_DEFINITIONS
# GLUT_INCLUDE_DIR
# GLUT_LIBRARIES (containing both debug and release libraries on win32)
# win32: GLUT_DEBUG_LIBRARY, GLUT_RELEASE_LIBRARY, GLUT_DEBUG_DLLS, GLUT_RELEASE_DLLS

IF(WIN32)
    SET(GLUT_DIR "${VRN_HOME}/ext/freeglut" CACHE PATH "If freeglut is not found, set this path")
   
    SET(GLUT_DEFINITIONS "-DGLUT_NO_LIB_PRAGMA")
       
    SET(GLUT_INCLUDE_DIR "${GLUT_DIR}")

    # set debug and release library
    IF(VRN_MSVC2015)
        SET(GLUT_DEBUG_LIBRARY      "${GLUT_DIR}/lib/msvc2015/freeglutd.lib")
        SET(GLUT_RELEASE_LIBRARY    "${GLUT_DIR}/lib/msvc2015/freeglut.lib")
        
        SET(GLUT_DEBUG_DLL          "${GLUT_DIR}/lib/msvc2015/freeglutd.dll")
        SET(GLUT_RELEASE_DLL        "${GLUT_DIR}/lib/msvc2015/freeglut.dll")
    ELSEIF(VRN_MSVC2017)
        SET(GLUT_DEBUG_LIBRARY      "${GLUT_DIR}/lib/msvc2017/freeglutd.lib")
        SET(GLUT_RELEASE_LIBRARY    "${GLUT_DIR}/lib/msvc2017/freeglut.lib")
        
        SET(GLUT_DEBUG_DLL          "${GLUT_DIR}/lib/msvc2017/freeglutd.dll")
        SET(GLUT_RELEASE_DLL        "${GLUT_DIR}/lib/msvc2017/freeglut.dll")
    ENDIF()

    IF (GLUT_DEBUG_LIBRARY AND GLUT_RELEASE_LIBRARY)
        SET(GLUT_LIBRARIES debug ${GLUT_DEBUG_LIBRARY} optimized ${GLUT_RELEASE_LIBRARY})
    ENDIF()
    
    IF(GLUT_INCLUDE_DIR AND GLUT_LIBRARIES)
        SET(GLUT_FOUND TRUE)
    ELSE()
        SET(GLUT_FOUND FALSE)
    ENDIF()
    
ELSE(WIN32)
    FIND_PACKAGE(GLUT REQUIRED)
    # we don't need the Xmu and Xi libraries, which have been added to ${GLUT_LIBRARIES} (see module)
    SET(GLUT_LIBRARIES ${GLUT_glut_LIBRARY})
ENDIF(WIN32)

UNSET(GlutVRN_DIR)
MARK_AS_ADVANCED(GLUT_DIR GLUT_INCLUDE_DIR GLUT_LIBRARIES GLUT_DEBUG_LIBRARY GLUT_RELEASE_LIBRARY GLUT_DEBUG_DLL GLUT_RELEASE_DLL)
