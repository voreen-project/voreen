if(VRN_MSVC2012)
    MESSAGE(FATAL_ERROR "Big Data Image Processing Module does NOT compile with MSVC 2012 - use any compiler supporting the c++11 standard completely")
ENDIF()

IF(NOT VRN_MODULE_HDF5)
    MESSAGE(FATAL_ERROR "Big Data Image Processing Module requires HDF5 Module")
ENDIF()

################################################################################
# External dependency: lz4 library
################################################################################
MESSAGE(STATUS "Trying to find lz4 libraries")
FIND_PACKAGE(LZ4)
IF(LZ4_FOUND)
    MESSAGE(STATUS "  - Found lz4 library")

    MESSAGE(STATUS "Include Directories: " ${LZ4_INCLUDE_DIR})
    MESSAGE(STATUS "Libraries: " ${LZ4_LIBRARIES})
    LIST(APPEND MOD_INCLUDE_DIRECTORIES ${LZ4_INCLUDE_DIR})
    LIST(APPEND MOD_LIBRARIES ${LZ4_LIBRARIES})
ELSE()
    MESSAGE(FATAL_ERROR "Could not find lz4 Library.")
ENDIF()

################################################################################
# Core module resources
################################################################################
SET(MOD_CORE_MODULECLASS BigDataImageProcessingModule)

SET(MOD_CORE_SOURCES
    ${MOD_DIR}/datastructures/lz4slicevolume.cpp
    ${MOD_DIR}/io/lz4slicevolumefilereader.cpp
    ${MOD_DIR}/io/volumedisklz4.cpp
    ${MOD_DIR}/processors/connectedcomponentanalysis.cpp
    ${MOD_DIR}/processors/segmentationquantification.cpp
    ${MOD_DIR}/processors/volumebricksource.cpp
    ${MOD_DIR}/processors/volumebricksave.cpp
    ${MOD_DIR}/processors/binarymedian.cpp
    ${MOD_DIR}/processors/volumeresampletransformation.cpp

    # Volumefiltering
    ${MOD_DIR}/volumefiltering/slicereader.cpp
    ${MOD_DIR}/volumefiltering/parallelvolumefilter.cpp
    ${MOD_DIR}/volumefiltering/gaussianfilter.cpp
    ${MOD_DIR}/volumefiltering/medianfilter.cpp
    ${MOD_DIR}/volumefiltering/binarymedianfilter.cpp
)
IF(VRN_MODULE_PLOTTING)
    LIST(APPEND MOD_CORE_SOURCES
        ${MOD_DIR}/processors/segmentationslicedensity.cpp
    )
ENDIF()
SET(MOD_CORE_HEADERS
    ${MOD_DIR}/datastructures/lz4slicevolume.h
    ${MOD_DIR}/io/lz4slicevolumefilereader.h
    ${MOD_DIR}/io/volumedisklz4.h
    ${MOD_DIR}/processors/connectedcomponentanalysis.h
    ${MOD_DIR}/processors/segmentationquantification.h
    ${MOD_DIR}/processors/volumebricksource.h
    ${MOD_DIR}/processors/volumebricksave.h
    ${MOD_DIR}/processors/binarymedian.cpp
    ${MOD_DIR}/processors/volumeresampletransformation.h

    # Volumefiltering
    ${MOD_DIR}/volumefiltering/slicereader.h
    ${MOD_DIR}/volumefiltering/volumefilter.h
    ${MOD_DIR}/volumefiltering/parallelvolumefilter.h
    ${MOD_DIR}/volumefiltering/gaussianfilter.h
    ${MOD_DIR}/volumefiltering/medianfilter.h
    ${MOD_DIR}/volumefiltering/binarymedianfilter.h
)
IF(VRN_MODULE_PLOTTING)
    LIST(APPEND MOD_CORE_HEADERS
        ${MOD_DIR}/processors/segmentationslicedensity.h
    )
ENDIF()
