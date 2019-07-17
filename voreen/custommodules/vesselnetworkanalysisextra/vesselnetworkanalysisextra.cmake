################################################################################
# Core module resources
################################################################################
IF(NOT VRN_MODULE_VESSELNETWORKANALYSIS)
    MESSAGE(FATAL_ERROR "VesselNetworkAnalysisExtra Module requires VesselNetworkAnalysis Module")
ENDIF()

IF(NOT VRN_MODULE_PLOTTING)
    MESSAGE(FATAL_ERROR "VesselNetworkAnalysisExtra Module requires Plotting Module")
ENDIF()

SET(MOD_CORE_MODULECLASS VesselNetworkAnalysisExtraModule)

################################################################################
# External dependency: YAML library
################################################################################
MESSAGE(STATUS "Trying to find YAML libraries")
FIND_PACKAGE(YAML-CPP)
IF(YAML_CPP_INCLUDE_DIR) # Apparently YAML_CPP does not declare a *_FOUND
    MESSAGE(STATUS "  - Found YAML library")

    LIST(APPEND MOD_INCLUDE_DIRECTORIES ${YAML_CPP_INCLUDE_DIR})
    MESSAGE(STATUS "Include Directories: " ${YAML_CPP_INCLUDE_DIR})
    MESSAGE(STATUS "Libraries: " ${YAML_CPP_LIBRARIES})
    IF(VRN_MSVC)
        ADD_DEFINITIONS("-DYAML_CPP_DLL")
        FOREACH(elem ${YAML_CPP_LIBRARIES})
            # Note that we set the same library for both, relase and debug.
            # This way the user has to decide how the library should be compiled.
            LIST(APPEND MOD_RELEASE_DLLS ${YAML-CPP_DIR}/../bin/${elem}.dll)
            LIST(APPEND MOD_RELEASE_LIBRARIES ${YAML-CPP_DIR}/../lib/${elem}.lib)
            LIST(APPEND MOD_DEBUG_DLLS ${YAML-CPP_DIR}/../bin/${elem}.dll)
            LIST(APPEND MOD_DEBUG_LIBRARIES ${YAML-CPP_DIR}/../lib/${elem}.lib)
        ENDFOREACH()
    ELSE()
        LIST(APPEND MOD_LIBRARIES ${YAML_CPP_LIBRARIES})
    ENDIF()
ELSE()
    MESSAGE(FATAL_ERROR "Could not find YAML Library.")
ENDIF()

SET(MOD_CORE_SOURCES
    ${MOD_DIR}/processors/aortasegmentation.cpp
    ${MOD_DIR}/processors/localandglobalthreshold.cpp
    ${MOD_DIR}/processors/lymphatictestvesselgenerator.cpp
    ${MOD_DIR}/processors/segmentationlistvalidation.cpp
    ${MOD_DIR}/processors/subgraphextractor.cpp
    ${MOD_DIR}/processors/templatesubgraphextractor.cpp
    ${MOD_DIR}/processors/vesselgraphstatplotter.cpp
    ${MOD_DIR}/processors/volumefloodfill.cpp
    ${MOD_DIR}/processors/createTestVolume.cpp
    ${MOD_DIR}/processors/createVesselAroundPoints.cpp
    ${MOD_DIR}/processors/foldVessel.cpp
    ${MOD_DIR}/processors/unfoldVessel.cpp
    ${MOD_DIR}/processors/volumelistloopfinalizer.cpp
    ${MOD_DIR}/processors/volumelistloopinitiator.cpp
    ${MOD_DIR}/processors/volumemultithreshold.cpp
)

SET(MOD_CORE_HEADERS
    ${MOD_DIR}/processors/aortasegmentation.h
    ${MOD_DIR}/processors/localandglobalthreshold.h
    ${MOD_DIR}/processors/lymphatictestvesselgenerator.h
    ${MOD_DIR}/processors/segmentationlistvalidation.h
    ${MOD_DIR}/processors/subgraphextractor.h
    ${MOD_DIR}/processors/templatesubgraphextractor.h
    ${MOD_DIR}/processors/vesselgraphstatplotter.h
    ${MOD_DIR}/processors/volumefloodfill.h
    ${MOD_DIR}/processors/createTestVolume.h
    ${MOD_DIR}/processors/foldVessel.h
    ${MOD_DIR}/processors/unfoldVessel.h
    ${MOD_DIR}/processors/createVesselAroundPoints.h
    ${MOD_DIR}/processors/volumelistloopfinalizer.h
    ${MOD_DIR}/processors/volumelistloopinitiator.h
    ${MOD_DIR}/processors/volumemultithreshold.h
)
################################################################################
# VRAS library by Jose Alejandro Matute Flores
################################################################################
OPTION(VRN_VESSELNETWORKANALYSISEXTRA_BUILD_VRAS "Build Vessel Reconstruction Library?" ON)
IF(${VRN_VESSELNETWORKANALYSISEXTRA_BUILD_VRAS})
    ADD_SUBDIRECTORY(${CMAKE_CURRENT_SOURCE_DIR}/custommodules/${module_dir}/ext/VRAS/ReconstructionLibrary)
    TARGET_LINK_LIBRARIES(${core_app_we} VascResc)
    ADD_DEFINITIONS("-DVESSELNETWORKANALYSISEXTRA_USE_VRAS")
ENDIF()
# Deployment
SET(MOD_INSTALL_DIRECTORIES
    ${MOD_DIR}/scripts
    ${MOD_DIR}/workspaces
)
