if(VRN_MSVC2012 OR VRN_MSVC2013)
    MESSAGE(FATAL_ERROR "Big Data Image Processing Extra Module does NOT compile with MSVC 2012 and MSVC 2013 - use MSVC 2015 or higher")
ENDIF()

IF(NOT VRN_MODULE_BIGDATAIMAGEPROCESSING)
    MESSAGE(FATAL_ERROR "Big Data Image Processing Module requires Big Data Image Processing Module")
ENDIF()

IF(NOT VRN_MODULE_PLOTTING)
    MESSAGE(FATAL_ERROR "Big Data Image Processing Extra Module requires Plotting Module")
ENDIF()

################################################################################
# Core module resources
################################################################################
SET(MOD_CORE_MODULECLASS BigDataImageProcessingExtraModule)

SET(MOD_CORE_SOURCES
    ${MOD_DIR}/processors/fatcellquantification.cpp
    ${MOD_DIR}/processors/segmentationquantification.cpp
    ${MOD_DIR}/processors/segmentationslicedensity.cpp
    ${MOD_DIR}/processors/volumebricksave.cpp
    ${MOD_DIR}/processors/volumebricksource.cpp
)
ENDIF()
SET(MOD_CORE_HEADERS
    ${MOD_DIR}/processors/fatcellquantification.h
    ${MOD_DIR}/processors/segmentationquantification.h
    ${MOD_DIR}/processors/segmentationslicedensity.h
    ${MOD_DIR}/processors/volumebricksave.h
    ${MOD_DIR}/processors/volumebricksource.h
)
ENDIF()
