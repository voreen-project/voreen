if(VRN_MSVC2012)
    MESSAGE(FATAL_ERROR "Big Data Image Processing Module does NOT compile with MSVC 2012 - use any compiler supporting the c++11 standard completely")
ENDIF()

IF(NOT VRN_MODULE_HDF5)
    MESSAGE(FATAL_ERROR "Big Data Image Processing Module requires HDF5 Module")
ENDIF()

################################################################################
# Core module resources
################################################################################
SET(MOD_CORE_MODULECLASS BigDataImageProcessingModule)

SET(MOD_CORE_SOURCES
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
