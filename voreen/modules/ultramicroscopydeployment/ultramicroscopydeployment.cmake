################################################################################
# UltramicroscopyDeployment module resources
################################################################################

IF(NOT VRN_MODULE_BASE)
    MESSAGE(FATAL_ERROR "UltramicroscopyDeployment Module requires Base Module")
ENDIF()
IF(NOT VRN_MODULE_DEVIL)
    MESSAGE (FATAL_ERROR "UltramicroscopyDeployment Module requires Devil Module")
ENDIF()
IF(NOT VRN_MODULE_FFMPEG)
    MESSAGE (FATAL_ERROR "UltramicroscopyDeployment Module requires Ffmpeg Module")
ENDIF()
IF(NOT VRN_MODULE_HDF5)
    MESSAGE(FATAL_ERROR "UltramicroscopyDeployment Module requires HDF5 Module")
ENDIF()
IF(NOT VRN_MODULE_OPENCL)
    MESSAGE(FATAL_ERROR "UltramicroscopyDeployment Module requires OpenCL Module")
ENDIF()
IF(NOT VRN_MODULE_OPENMP)
    MESSAGE(FATAL_ERROR "UltramicroscopyDeployment Module requires OpenMP Module")
ENDIF()
IF(NOT VRN_MODULE_PLOTTING)
    MESSAGE(FATAL_ERROR "UltramicroscopyDeployment Module requires Plotting Module")
ENDIF()
IF(NOT VRN_MODULE_RANDOMWALKER)
    MESSAGE(FATAL_ERROR "UltramicroscopyDeployment Module requires RandomWalker Module")
ENDIF()
IF(NOT VRN_MODULE_STAGING)
    MESSAGE(FATAL_ERROR "UltramicroscopyDeployment Module requires Staging Module")
ENDIF()
IF(NOT VRN_MODULE_TIFF)
    MESSAGE(FATAL_ERROR "UltramicroscopyDeployment Module requires TIFF Module")
ENDIF()

SET(MOD_CORE_MODULECLASS UltramicroscopyDeploymentModule)

# deployment
SET(MOD_INSTALL_DIRECTORIES
#    ${MOD_DIR}/workspaces
)

SET(MOD_CORE_SOURCES
    #Processors

)

SET(MOD_CORE_HEADERS
    #Processors

)
