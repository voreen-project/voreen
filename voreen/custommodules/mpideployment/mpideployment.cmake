################################################################################
# MPIDeployment module resources
################################################################################

IF(NOT VRN_MODULE_BASE)
    MESSAGE(FATAL_ERROR "MPIDeployment Module requires Base Module")
ENDIF()
#IF(NOT VRN_MODULE_DEPRECATED)
#    MESSAGE (FATAL_ERROR "MPIDeployment Module requires Deprecated Module")
#ENDIF()
IF(NOT VRN_MODULE_DEVIL)
    MESSAGE (FATAL_ERROR "MPIDeployment Module requires Devil Module")
ENDIF()
IF(NOT VRN_MODULE_FFMPEG)
    MESSAGE (FATAL_ERROR "MPIDeployment Module requires Ffmpeg Module")
ENDIF()
#IF(NOT VRN_MODULE_GDCM)
#    MESSAGE(FATAL_ERROR "MPIDeployment Module requires GDCM Module")
#ENDIF()
IF(NOT VRN_MODULE_HDF5)
    MESSAGE(FATAL_ERROR "MPIDeployment Module requires HDF5 Module")
ENDIF()
IF(NOT VRN_MODULE_OPENCL)
    MESSAGE(FATAL_ERROR "MPIDeployment Module requires OpenCL Module")
ENDIF()
IF(NOT VRN_MODULE_OPENMP)
    MESSAGE(FATAL_ERROR "MPIDeployment Module requires OpenMP Module")
ENDIF()
IF(NOT VRN_MODULE_PLOTTING)
    MESSAGE(FATAL_ERROR "MPIDeployment Module requires Plotting Module")
ENDIF()
IF(NOT VRN_MODULE_RANDOMWALKER)
    MESSAGE(FATAL_ERROR "MPIDeployment Module requires RandomWalker Module")
ENDIF()
IF(NOT VRN_MODULE_STAGING)
    MESSAGE(FATAL_ERROR "MPIDeployment Module requires Staging Module")
ENDIF()
IF(NOT VRN_MODULE_TIFF)
    MESSAGE(FATAL_ERROR "MPIDeployment Module requires TIFF Module")
ENDIF()

SET(MOD_CORE_MODULECLASS MPIDeploymentModule)

# deployment
SET(MOD_INSTALL_DIRECTORIES
#    ${MOD_DIR}/workspaces
)

SET(MOD_CORE_SOURCES
    #Processors
    ${MOD_DIR}/processors/transfuncalphachannelanimation.cpp
)

SET(MOD_CORE_HEADERS
    #Processors
    ${MOD_DIR}/processors/transfuncalphachannelanimation.h
)
