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
    ${MOD_DIR}/processors/binarymedian.cpp
    ${MOD_DIR}/processors/fatcellquantification.cpp
    ${MOD_DIR}/processors/largetestdatagenerator.cpp
    ${MOD_DIR}/processors/segmentationslicedensity.cpp
    ${MOD_DIR}/processors/segmentationquantification.cpp
    ${MOD_DIR}/processors/volumebricksave.cpp
    ${MOD_DIR}/processors/volumebricksource.cpp
    ${MOD_DIR}/processors/volumecomparison.cpp
    ${MOD_DIR}/processors/volumedetectionlinker.cpp
)
SET(MOD_CORE_HEADERS
    ${MOD_DIR}/processors/binarymedian.h
    ${MOD_DIR}/processors/fatcellquantification.h
    ${MOD_DIR}/processors/largetestdatagenerator.h
    ${MOD_DIR}/processors/segmentationslicedensity.cpp
    ${MOD_DIR}/processors/segmentationslicedensity.h
    ${MOD_DIR}/processors/segmentationquantification.cpp
    ${MOD_DIR}/processors/volumebricksave.h
    ${MOD_DIR}/processors/volumebricksource.h
    ${MOD_DIR}/processors/volumecomparison.h
    ${MOD_DIR}/processors/volumedetectionlinker.h
)
