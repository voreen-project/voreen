# Try to find GLEW library and include path. Once done this will define:
# GLEW_FOUND
# GLEW_DEFINITIONS
# GLEW_INCLUDE_DIR
# GLEW_LIBRARY (containing both debug and release libraries on win32)
# win32: GLEW_LIBRARY_DEBUG, GLEW_LIBRARY_RELEASE, GLEW_DLL_DEBUG, GLEW_DLL_RELEASE, GLEW_LICENSE_FILE

IF (WIN32)
    SET(GLEW_DIR "${VRN_HOME}/ext/glew" CACHE PATH "If glew is not found, set this path")
   
    SET(GLEW_INCLUDE_DIR "${GLEW_DIR}/include")

    # set debug and release library
    IF(VRN_MSVC)
		SET(GLEW_LIBRARY_DEBUG      "${GLEW_DIR}/lib/glew32.lib")
		SET(GLEW_DLL_DEBUG          "${GLEW_DIR}/lib/glew32.dll")
		SET(GLEW_LIBRARY_RELEASE    "${GLEW_DIR}/lib/glew32.lib")
		SET(GLEW_DLL_RELEASE        "${GLEW_DIR}/lib/glew32.dll")
    ELSEIF(VRN_MINGW)
        SET(GLEW_DEFINITIONS "-DGLEW_STATIC")
        SET(GLEW_LIBRARY_DEBUG      "${GLEW_DIR}/lib/glew32s.lib")
        SET(GLEW_LIBRARY_RELEASE    "${GLEW_DIR}/lib/glew32s.lib")
    ENDIF()

    IF (GLEW_LIBRARY_DEBUG AND GLEW_LIBRARY_RELEASE)
        SET(GLEW_LIBRARY debug ${GLEW_LIBRARY_DEBUG} optimized ${GLEW_LIBRARY_RELEASE})
    ENDIF(GLEW_LIBRARY_DEBUG AND GLEW_LIBRARY_RELEASE)
    
    SET(GLEW_LICENSE_FILE "${GLEW_DIR}/license.txt")
    
ELSE (WIN32)
    FIND_PATH( 
        GLEW_INCLUDE_DIR 
        NAMES GL/glew.h
        PATHS ${GLEW_DIR}/include /usr/include /usr/local/include /sw/include /opt/local/include
        DOC "The directory where GL/glew.h resides"
    )
    FIND_LIBRARY(
        GLEW_LIBRARY
        NAMES GLEW glew
        PATHS ${GLEW_DIR}/lib /usr/lib64 /usr/lib /usr/local/lib64 /usr/local/lib /sw/lib /opt/local/lib /usr/lib/x86_64-linux-gnu
        DOC "The GLEW library"
    )
ENDIF (WIN32)

IF(GLEW_INCLUDE_DIR AND GLEW_LIBRARY)
    SET(GLEW_FOUND TRUE)
ELSE(GLEW_INCLUDE_DIR AND GLEW_LIBRARY)
    SET(GLEW_FOUND FALSE)
ENDIF(GLEW_INCLUDE_DIR AND GLEW_LIBRARY)

UNSET(GlewVRN_DIR)
MARK_AS_ADVANCED(GLEW_DIR GlewVRN_DIR GLEW_INCLUDE_DIR GLEW_LIBRARY GLEW_LIBRARY_DEBUG GLEW_LIBRARY_RELEASE GLEW_DLL_DEBUG GLEW_DLL_RELEASE)
