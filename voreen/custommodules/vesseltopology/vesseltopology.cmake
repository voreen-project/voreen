################################################################################
# Core module resources
################################################################################
IF(NOT VRN_MODULE_BIGDATAIMAGEPROCESSING)
    MESSAGE(FATAL_ERROR "VolumeTopology Module requires big data image processing Module")
ENDIF()

IF(NOT VRN_MODULE_HDF5)
    MESSAGE(FATAL_ERROR "VolumeTopology Module requires HDF5 Module")
ENDIF()

IF(NOT VRN_MODULE_PLOTTING)
    MESSAGE(FATAL_ERROR "VolumeTopology Module requires Plotting Module")
ENDIF()

SET(MOD_CORE_MODULECLASS VesselTopologyModule)

################################################################################
# External dependency: Lemon (graph library)
################################################################################
MESSAGE(STATUS "Trying to find Lemon libraries")
FIND_PACKAGE(LEMON QUIET) # Disable output if package cannot be found
IF(LEMON_FOUND)
    MESSAGE(STATUS "  - Found Lemon library")

    LIST(APPEND MOD_INCLUDE_DIRECTORIES ${LEMON_INCLUDE_DIRS})
    MESSAGE(STATUS "Include Directories: " ${LEMON_INCLUDE_DIRS})
    LIST(APPEND MOD_LIBRARIES ${LEMON_LIBRARIES})
    MESSAGE(STATUS "Libraries: " ${LEMON_LIBRARIES})
    SET(VRN_LEMON_FOUND 1)
    ADD_DEFINITIONS("-DLEMON_FOUND")
ELSE()
    MESSAGE(STATUS "Could not find Lemon Library, some Processors will not be available.")
ENDIF()

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
            LIST(APPEND MOD_RELEASE_DLLS ${YAML-CPP_DIR}/Release/${YAML_CPP_LIBRARIES}.dll)
            LIST(APPEND MOD_RELEASE_LIBRARIES ${YAML-CPP_DIR}/Release/${YAML_CPP_LIBRARIES}.lib)
            # Don't copy debug binaries, since they are named the same as the release binaries.
            #LIST(APPEND MOD_DEBUG_DLLS ${YAML-CPP_DIR}/Debug/${YAML_CPP_LIBRARIES}.dll)
            #LIST(APPEND MOD_DEBUG_LIBRARIES ${YAML-CPP_DIR}/Debug/${YAML_CPP_LIBRARIES}.lib)
        ENDFOREACH()
    ELSE()
        LIST(APPEND MOD_LIBRARIES ${YAML_CPP_LIBRARIES})
    ENDIF()
ELSE()
    MESSAGE(FATAL_ERROR "Could not find YAML Library.")
ENDIF()

SET(MOD_CORE_SOURCES
    ${MOD_DIR}/algorithm/idvolume.cpp
    ${MOD_DIR}/algorithm/streaminggraphcreation.cpp
    ${MOD_DIR}/algorithm/surface.cpp
    ${MOD_DIR}/algorithm/vesselgraphnormalization.cpp
    ${MOD_DIR}/algorithm/volumemask.cpp
    ${MOD_DIR}/datastructures/vesselgraph.cpp
    ${MOD_DIR}/datastructures/protovesselgraph.cpp
    ${MOD_DIR}/ports/vesselgraphport.cpp
    ${MOD_DIR}/processors/appropriatespacinglinker.cpp
    ${MOD_DIR}/processors/aortasegmentation.cpp
    ${MOD_DIR}/processors/localandglobalthreshold.cpp
    ${MOD_DIR}/processors/segmentationlistvalidation.cpp
    ${MOD_DIR}/processors/subgraphextractor.cpp
    ${MOD_DIR}/processors/templatesubgraphextractor.cpp
    ${MOD_DIR}/processors/vascusynthgraphloader.cpp
    ${MOD_DIR}/processors/vesselgraphcreator.cpp
    ${MOD_DIR}/processors/vesselgraphglobalstats.cpp
    ${MOD_DIR}/processors/vesselgraphnormalizer.cpp
    ${MOD_DIR}/processors/vesselgraphperturbation.cpp
    ${MOD_DIR}/processors/vesselgraphrenderer.cpp
    ${MOD_DIR}/processors/vesselgraphsave.cpp
    ${MOD_DIR}/processors/vesselgraphskeletonextractor.cpp
    ${MOD_DIR}/processors/vesselgraphsource.cpp
    ${MOD_DIR}/processors/vesselgraphstatplotter.cpp
    ${MOD_DIR}/processors/vesselnessextractor.cpp
    ${MOD_DIR}/processors/volumefloodfill.cpp
    ${MOD_DIR}/processors/createTestVolume.cpp
    ${MOD_DIR}/processors/createVesselAroundPoints.cpp
    ${MOD_DIR}/processors/foldVessel.cpp
    ${MOD_DIR}/processors/unfoldVessel.cpp
    ${MOD_DIR}/processors/volumelistloopfinalizer.cpp
    ${MOD_DIR}/processors/volumelistloopinitiator.cpp
    ${MOD_DIR}/processors/volumemultiplier.cpp
    ${MOD_DIR}/processors/volumemultithreshold.cpp
    ${MOD_DIR}/processors/volumeslicepadding.cpp
    ${MOD_DIR}/processors/volumethinning.cpp
    ${MOD_DIR}/util/tasktimelogger.cpp
)

SET(MOD_CORE_HEADERS
    ${MOD_DIR}/algorithm/boundshierarchy.h
    ${MOD_DIR}/algorithm/idvolume.h
    ${MOD_DIR}/algorithm/streaminggraphcreation.h
    ${MOD_DIR}/algorithm/surface.cpp
    ${MOD_DIR}/algorithm/vesselgraphnormalization.h
    ${MOD_DIR}/algorithm/volumemask.h
    ${MOD_DIR}/datastructures/vesselgraph.h
    ${MOD_DIR}/datastructures/protovesselgraph.h
    ${MOD_DIR}/ports/vesselgraphport.h
    ${MOD_DIR}/processors/appropriatespacinglinker.h
    ${MOD_DIR}/processors/aortasegmentation.h
    ${MOD_DIR}/processors/localandglobalthreshold.h
    ${MOD_DIR}/processors/segmentationlistvalidation.h
    ${MOD_DIR}/processors/subgraphextractor.h
    ${MOD_DIR}/processors/templatesubgraphextractor.h
    ${MOD_DIR}/processors/vascusynthgraphloader.h
    ${MOD_DIR}/processors/vesselgraphcreator.h
    ${MOD_DIR}/processors/vesselgraphglobalstats.h
    ${MOD_DIR}/processors/vesselgraphnormalizer.h
    ${MOD_DIR}/processors/vesselgraphperturbation.h
    ${MOD_DIR}/processors/vesselgraphrenderer.h
    ${MOD_DIR}/processors/vesselgraphsave.h
    ${MOD_DIR}/processors/vesselgraphskeletonextractor.h
    ${MOD_DIR}/processors/vesselgraphsource.h
    ${MOD_DIR}/processors/vesselgraphstatplotter.h
    ${MOD_DIR}/processors/vesselnessextractor.h
    ${MOD_DIR}/processors/volumefloodfill.h
    ${MOD_DIR}/processors/createTestVolume.h
    ${MOD_DIR}/processors/foldVessel.h
    ${MOD_DIR}/processors/unfoldVessel.h
    ${MOD_DIR}/processors/createVesselAroundPoints.h
    ${MOD_DIR}/processors/volumelistloopfinalizer.h
    ${MOD_DIR}/processors/volumelistloopinitiator.h
    ${MOD_DIR}/processors/volumemultiplier.h
    ${MOD_DIR}/processors/volumemultithreshold.h
    ${MOD_DIR}/processors/volumeslicepadding.h
    ${MOD_DIR}/processors/volumethinning.h
    ${MOD_DIR}/util/tasktimelogger.h
)
################################################################################
# VRAS library by Jose Alejandro Matute Flores
################################################################################
OPTION(VRN_VESSELTOPOLOGY_BUILD_VRAS "Build Vessel Reconstruction Library?" ON)
IF(${VRN_VESSELTOPOLOGY_BUILD_VRAS})
    ADD_SUBDIRECTORY(${CMAKE_CURRENT_SOURCE_DIR}/custommodules/${module_dir}/ext/VRAS/ReconstructionLibrary)
    TARGET_LINK_LIBRARIES(${core_app_we} VascResc)
    ADD_DEFINITIONS("-DVESSELTOPOLOGY_USE_VRAS")
ENDIF()

IF(VRN_LEMON_FOUND)
SET(MOD_CORE_HEADERS ${MOD_CORE_HEADERS}
    ${MOD_DIR}/processors/vesselgraphcomparison.h
)
SET(MOD_CORE_SOURCES ${MOD_CORE_SOURCES}
    ${MOD_DIR}/processors/vesselgraphcomparison.cpp
)
ENDIF()

# Deployment
SET(MOD_INSTALL_DIRECTORIES
    ${MOD_DIR}/glsl
    ${MOD_DIR}/scripts
    ${MOD_DIR}/workspaces
    ${MOD_DIR}/tfs
)
