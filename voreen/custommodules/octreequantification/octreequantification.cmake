################################################################################
# Core module resources 
################################################################################

IF(NOT VRN_MODULE_ULTRAMICROSCOPYDEPLOYMENT)
    MESSAGE(FATAL_ERROR "Octree Quantification Module requires UltramicroscopyDeployment Module")
ENDIF()
IF(NOT VRN_MODULE_PLOTTING)
    MESSAGE(FATAL_ERROR "Octree Quantification Module requires Plotting Module")
ENDIF()
IF(NOT VRN_MODULE_RANDOMWALKER)
    MESSAGE(FATAL_ERROR "Octree Quantification Module requires Randomwalker Module")
ENDIF()
IF(NOT VRN_MODULE_STAGING)
    MESSAGE(FATAL_ERROR "Octree Quantification Module requires Staging Module")
ENDIF()


SET(MOD_CORE_MODULECLASS OctreeQuantificationModule)

SET(MOD_CORE_SOURCES
    ${MOD_DIR}/processors/boundingboxratiolinker.cpp
    ${MOD_DIR}/processors/channelshiftratiolinker.cpp
    ${MOD_DIR}/processors/octreesegmentationquantification.cpp
    ${MOD_DIR}/processors/slicenumberratiolinker.cpp
    ${MOD_DIR}/processors/thresholdtotransfunclinker.cpp
    ${MOD_DIR}/processors/transfuncchannellinker.cpp
    ${MOD_DIR}/util/octreequantificationnodequeue.cpp
    ${MOD_DIR}/util/octreequantificationthread.cpp
    ${MOD_DIR}/util/octreequantificationresults.cpp
)

SET(MOD_CORE_HEADERS
    ${MOD_DIR}/processors/boundingboxratiolinker.h
    ${MOD_DIR}/processors/channelshiftratiolinker.h
    ${MOD_DIR}/processors/octreesegmentationquantification.h
    ${MOD_DIR}/processors/slicenumberratiolinker.h
    ${MOD_DIR}/processors/thresholdtotransfunclinker.h
    ${MOD_DIR}/processors/transfuncchannellinker.h
    ${MOD_DIR}/util/octreequantificationnodequeue.h
    ${MOD_DIR}/util/octreequantificationthread.h
    ${MOD_DIR}/util/octreequantificationresults.h
)
   
