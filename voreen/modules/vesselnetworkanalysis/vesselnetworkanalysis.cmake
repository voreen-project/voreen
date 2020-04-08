################################################################################
# Core module resources
################################################################################
IF(NOT VRN_MODULE_BIGDATAIMAGEPROCESSING)
    MESSAGE(FATAL_ERROR "VesselNetworkAnalysis Module requires big data image processing Module")
ENDIF()

IF(NOT VRN_MODULE_HDF5)
    MESSAGE(FATAL_ERROR "VesselNetworkAnalysis Module requires HDF5 Module")
ENDIF()

SET(MOD_CORE_MODULECLASS VesselNetworkAnalysisModule)

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

SET(MOD_CORE_SOURCES
    ${MOD_DIR}/algorithm/idvolume.cpp
    ${MOD_DIR}/algorithm/streaminggraphcreation.cpp
    ${MOD_DIR}/algorithm/surface.cpp
    ${MOD_DIR}/algorithm/vesselgraphrefinement.cpp
    ${MOD_DIR}/algorithm/volumemask.cpp
    ${MOD_DIR}/datastructures/vesselgraph.cpp
    ${MOD_DIR}/datastructures/protovesselgraph.cpp
    ${MOD_DIR}/ports/vesselgraphport.cpp
    ${MOD_DIR}/ports/vesselgraphlistport.cpp
    ${MOD_DIR}/processors/appropriatespacinglinker.cpp
    ${MOD_DIR}/processors/vascusynthgraphloader.cpp
    ${MOD_DIR}/processors/vesselgraphcreator.cpp
    ${MOD_DIR}/processors/vesselgraphglobalstats.cpp
    ${MOD_DIR}/processors/vesselgraphrefiner.cpp
    ${MOD_DIR}/processors/vesselgraphperturbation.cpp
    ${MOD_DIR}/processors/vesselgraphrenderer.cpp
    ${MOD_DIR}/processors/vesselgraphsave.cpp
    ${MOD_DIR}/processors/vesselgraphcenterlineconverter.cpp
    ${MOD_DIR}/processors/vesselgraphsource.cpp
    ${MOD_DIR}/processors/vesselgraphselector.cpp
    ${MOD_DIR}/processors/vesselnessextractor.cpp
    ${MOD_DIR}/processors/volumemultiplier.cpp
    ${MOD_DIR}/processors/volumeslicepadding.cpp
    ${MOD_DIR}/processors/volumesurfacenoise.cpp
    ${MOD_DIR}/processors/volumethinning.cpp
    ${MOD_DIR}/util/tasktimelogger.cpp
)

SET(MOD_CORE_HEADERS
    ${MOD_DIR}/algorithm/idvolume.h
    ${MOD_DIR}/algorithm/intervalwalker.h
    ${MOD_DIR}/algorithm/streaminggraphcreation.h
    ${MOD_DIR}/algorithm/surface.h
    ${MOD_DIR}/algorithm/vesselgraphrefinement.h
    ${MOD_DIR}/algorithm/volumemask.h
    ${MOD_DIR}/datastructures/vesselgraph.h
    ${MOD_DIR}/datastructures/protovesselgraph.h
    ${MOD_DIR}/datastructures/kdtree.h
    ${MOD_DIR}/ports/vesselgraphport.h
    ${MOD_DIR}/ports/vesselgraphlistport.h
    ${MOD_DIR}/processors/appropriatespacinglinker.h
    ${MOD_DIR}/processors/vascusynthgraphloader.h
    ${MOD_DIR}/processors/vesselgraphcreator.h
    ${MOD_DIR}/processors/vesselgraphglobalstats.h
    ${MOD_DIR}/processors/vesselgraphrefiner.h
    ${MOD_DIR}/processors/vesselgraphperturbation.h
    ${MOD_DIR}/processors/vesselgraphrenderer.h
    ${MOD_DIR}/processors/vesselgraphsave.h
    ${MOD_DIR}/processors/vesselgraphcenterlineconverter.h
    ${MOD_DIR}/processors/vesselgraphsource.h
    ${MOD_DIR}/processors/vesselgraphselector.h
    ${MOD_DIR}/processors/vesselnessextractor.h
    ${MOD_DIR}/processors/volumemultiplier.h
    ${MOD_DIR}/processors/volumeslicepadding.h
    ${MOD_DIR}/processors/volumesurfacenoise.h
    ${MOD_DIR}/processors/volumethinning.h
    ${MOD_DIR}/util/tasktimelogger.h
)

IF(VRN_LEMON_FOUND)
    SET(MOD_CORE_HEADERS ${MOD_CORE_HEADERS}
        ${MOD_DIR}/processors/vesselgraphcomparison.h
    )
    SET(MOD_CORE_SOURCES ${MOD_CORE_SOURCES}
        ${MOD_DIR}/processors/vesselgraphcomparison.cpp
    )
ENDIF()
OPTION(VRN_VESSELNETWORKANALYSIS_BUILD_NETMETS "Build Netmets Library? (Requires cuda)" OFF)
IF(${VRN_VESSELNETWORKANALYSIS_BUILD_NETMETS})
    SET(VRN_NETMETS_VRN_ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR})
    SET(VRN_NETMETS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/modules/${module_dir}/ext/netmets)
    ADD_SUBDIRECTORY(${VRN_NETMETS_DIRECTORY})
    LIST(APPEND MOD_LIBRARIES "netmets")
    ADD_DEFINITIONS("-DVESSELNETWORKANALYSIS_USE_NETMETS")
ENDIF()


# Deployment
SET(MOD_INSTALL_DIRECTORIES
    ${MOD_DIR}/glsl
    ${MOD_DIR}/workspaces
)
SET(MOD_INSTALL_FILES
    # munkres-cpp
    ${MOD_DIR}/ext/munkres-cpp/COPYING
    ${MOD_DIR}/ext/munkres-cpp/AUTHORS.txt
)
