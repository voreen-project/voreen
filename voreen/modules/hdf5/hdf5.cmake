
################################################################################
# External dependency: HDF5 library
################################################################################

IF(WIN32)
    SET(HDF5_INCLUDE_DIR ${MOD_DIR}/ext/hdf5/include)
    SET(MOD_INCLUDE_DIRECTORIES ${HDF5_INCLUDE_DIR})
    
    # Define used libs:
    LIST(APPEND HDF5_LIB_NAMES 
        "libhdf5" "libhdf5_cpp" "libhdf5_hl" "libhdf5_hl_cpp" "libhdf5_tools" "libszip" "libzlib"
    )
    
    # set debug and release libraries
    IF(VRN_MSVC2012)            
		SET(HDF5_LIB_DIR "${MOD_DIR}/ext/hdf5/lib/msvc2012")
	ELSEIF(VRN_MSVC2013)            
		SET(HDF5_LIB_DIR "${MOD_DIR}/ext/hdf5/lib/msvc2013")
	ELSEIF(VRN_MSVC2015)            
		SET(HDF5_LIB_DIR "${MOD_DIR}/ext/hdf5/lib/msvc2015")
	ELSEIF(VRN_MSVC2017)            
		SET(HDF5_LIB_DIR "${MOD_DIR}/ext/hdf5/lib/msvc2017")
	ENDIF()
    
    # set libraries
    FOREACH(elem ${HDF5_LIB_NAMES})
        LIST(APPEND MOD_DEBUG_LIBRARIES    "${HDF5_LIB_DIR}/${elem}_D.lib")
        LIST(APPEND MOD_RELEASE_LIBRARIES  "${HDF5_LIB_DIR}/${elem}.lib")
    ENDFOREACH()
    
    SET(MOD_INSTALL_FILES
        ${MOD_DIR}/ext/hdf5/COPYING
    )

ELSEIF(UNIX)

    SET(VRN_USE_HDF5_VERSION "1.8" CACHE STRING "HDF5 version")
    SET_PROPERTY(CACHE VRN_USE_HDF5_VERSION PROPERTY STRINGS "1.8" "1.10")

    SET(HDF5_FOUND FALSE)
    IF(${VRN_USE_HDF5_VERSION} MATCHES "1.8")
    
        MESSAGE(STATUS "Trying to find HDF5 libraries (in between versions 1.8.13 and 1.10.0) with C++ support...")

        # First: Try to find the needed headers and libraries in some hard coded locations corresponding to
        # a definite 1.8 library. This is useful if HDF5 version 1.10 is the default on a system, but
        # version 1.8 can be installed in a non-standard location. This is similar to the procedure for ffmpeg
        # packages (at least on Arch Linux).
        FIND_PATH(HDF5_INCLUDE_DIR_MANUAL H5Cpp.h
            PATHS
            /usr/include/hdf5_18/
            /usr/local/include/hdf5_18/
            /opt/local/include/hdf5_18/
            NO_DEFAULT_PATH
        )

        FIND_LIBRARY(HDF5_LIB_C_MANUAL
            NAMES hdf5
            PATHS
            /usr/lib/hdf5_18/
            /usr/local/lib/hdf5_18/
            /opt/local/lib/hdf5_18/
            NO_DEFAULT_PATH
        )

        FIND_LIBRARY(HDF5_LIB_CPP_MANUAL
            NAMES hdf5_cpp
            PATHS
            /usr/lib/hdf5_18/
            /usr/local/lib/hdf5_18/
            /opt/local/lib/hdf5_18/
            NO_DEFAULT_PATH
        )

        IF(HDF5_INCLUDE_DIR_MANUAL AND HDF5_LIB_C_MANUAL AND HDF5_LIB_CPP_MANUAL)
            SET(HDF5_LIBRARIES
                ${HDF5_LIB_C_MANUAL}
                ${HDF5_LIB_CPP_MANUAL}
            )
            SET(HDF5_INCLUDE_DIRS
                ${HDF5_INCLUDE_DIR_MANUAL}
            )

            SET(HDF5_VERSION 1.8)
            SET(HDF5_FOUND TRUE)
        ENDIF()
    ENDIF()
    
    IF(NOT HDF5_FOUND)
        # Only if we did not find HDF5 at the hardcoded locations we resort to the find_package script.
        find_package(HDF5 1.8.13 COMPONENTS C CXX REQUIRED)
    ENDIF()
    
    IF(HDF5_FOUND)
        MESSAGE(STATUS "Found HDF5 version ${HDF5_VERSION}")
        IF(${VRN_USE_HDF5_VERSION} MATCHES "1.8" AND NOT HDF5_VERSION VERSION_LESS 1.10.0)
            MESSAGE(FATAL_ERROR "Could not find HDF5 1.8, but more recent version instead.")
        ELSEIF(${VRN_USE_HDF5_VERSION} MATCHES "1.10" AND HDF5_VERSION VERSION_LESS 1.10.0)
            MESSAGE(FATAL_ERROR "Could not find HDF5 1.10, but less recent version instead.")
        ENDIF()

        SET(MOD_INCLUDE_DIRECTORIES ${HDF5_INCLUDE_DIRS})
        MESSAGE(STATUS "Include Directories: " ${HDF5_INCLUDE_DIRS})
        SET(MOD_LIBRARIES ${HDF5_LIBRARIES})
        MESSAGE(STATUS "Libraries: " ${HDF5_LIBRARIES})
    ELSE()
        MESSAGE(FATAL_ERROR "HDF5 library could not be found")
    ENDIF()
ENDIF()


################################################################################
# Core module resources
################################################################################
SET(MOD_CORE_MODULECLASS HDF5Module)

SET(MOD_CORE_SOURCES
    ${MOD_DIR}/io/hdf5volumereader.cpp
    ${MOD_DIR}/io/hdf5volumewriter.cpp
    ${MOD_DIR}/io/volumediskhdf5.cpp
    ${MOD_DIR}/io/hdf5filevolume.cpp
    ${MOD_DIR}/utils/hdf5utils.cpp
)

SET(MOD_CORE_HEADERS
    ${MOD_DIR}/io/hdf5volumereader.h
    ${MOD_DIR}/io/hdf5volumewriter.h
    ${MOD_DIR}/io/volumediskhdf5.h
    ${MOD_DIR}/io/hdf5filevolume.h
    ${MOD_DIR}/utils/hdf5utils.h
)
