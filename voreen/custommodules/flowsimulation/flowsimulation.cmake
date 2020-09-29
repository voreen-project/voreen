################################################################################
# Core module resources
################################################################################

IF(NOT VRN_MODULE_PYTHON) # for converting
    MESSAGE(WARNING "FlowSimulation Module requires Python Module for converter scripts")
ENDIF()
IF(NOT VRN_MODULE_VESSELNETWORKANALYSIS) # for flow indicator detection
    MESSAGE(WARNING "FlowSimulation Module requires VesselNetworkAnalysis Module for flow indicator detection")
ENDIF()
IF(NOT VRN_MODULE_PLOTTING) # for flow indicator analysis
    MESSAGE(WARNING "FlowSimulation Module requires Plotting Module for flow indicator analysis")
ENDIF()

SET(MOD_CORE_MODULECLASS FlowSimulationModule)

SET(MOD_CORE_SOURCES
    ${MOD_DIR}/flowsimulationmodule.cpp

    # datastructures
    ${MOD_DIR}/datastructures/flowparameters.cpp
    
    # ports
    ${MOD_DIR}/ports/flowparametrizationport.cpp

    # processors
    ${MOD_DIR}/processors/features/lambda2criterion.cpp
    ${MOD_DIR}/processors/features/lambdacicriterion.cpp
    ${MOD_DIR}/processors/geometry/geometryclose.cpp
    ${MOD_DIR}/processors/geometry/geometryoffsetremove.cpp
    ${MOD_DIR}/processors/geometry/geometrysmoothnormals.cpp
    ${MOD_DIR}/processors/render/unalignedsliceviewer.cpp
    ${MOD_DIR}/processors/simulation/flowcharacteristics.cpp
    ${MOD_DIR}/processors/simulation/flowensemblecreator.cpp
    #${MOD_DIR}/processors/simulation/flowindicatorselection.cpp
    ${MOD_DIR}/processors/simulation/flowindicatorrenderer.cpp
    ${MOD_DIR}/processors/simulation/flowparametrizationensemble.cpp
    ${MOD_DIR}/processors/simulation/flowparametrizationrun.cpp
    ${MOD_DIR}/processors/simulation/flowsimulationcluster.cpp
    ${MOD_DIR}/processors/simulation/flowsimulationgeometry.cpp
    ${MOD_DIR}/processors/volume/flowtestdatagenerator.cpp
    ${MOD_DIR}/processors/volume/phaseunwrapping.cpp
    ${MOD_DIR}/processors/volume/volumelistadapter.cpp
    ${MOD_DIR}/processors/volume/volumelistaggregate.cpp
    ${MOD_DIR}/processors/volume/volumelistmultichanneladapter.cpp
    ${MOD_DIR}/processors/volume/volumenoise.cpp
    ${MOD_DIR}/processors/volume/volumeselectormultichannel.cpp

    # utils
    ${MOD_DIR}/utils/utils.cpp
)

SET(MOD_CORE_HEADERS
    ${MOD_DIR}/flowsimulationmodule.h

    # datastructures
    ${MOD_DIR}/datastructures/flowparameters.h

    # ports
    ${MOD_DIR}/ports/flowparametrizationport.h

    # processors
    ${MOD_DIR}/processors/features/lambda2criterion.h
    ${MOD_DIR}/processors/features/lambdacicriterion.h
    ${MOD_DIR}/processors/geometry/geometryclose.h
    ${MOD_DIR}/processors/geometry/geometryoffsetremove.h
    ${MOD_DIR}/processors/geometry/geometrysmoothnormals.h
    ${MOD_DIR}/processors/render/unalignedsliceviewer.h
    ${MOD_DIR}/processors/simulation/flowcharacteristics.h
    ${MOD_DIR}/processors/simulation/flowensemblecreator.h
    #${MOD_DIR}/processors/simulation/flowindicatorselection.h
    ${MOD_DIR}/processors/simulation/flowindicatorrenderer.h
    ${MOD_DIR}/processors/simulation/flowparametrizationensemble.h
    ${MOD_DIR}/processors/simulation/flowparametrizationrun.h
    ${MOD_DIR}/processors/simulation/flowsimulationcluster.h
    ${MOD_DIR}/processors/simulation/flowsimulationgeometry.h
    ${MOD_DIR}/processors/volume/flowtestdatagenerator.h
    ${MOD_DIR}/processors/volume/phaseunwrapping.h
    ${MOD_DIR}/processors/volume/volumelistadapter.h
    ${MOD_DIR}/processors/volume/volumelistaggregate.h
    ${MOD_DIR}/processors/volume/volumelistmultichanneladapter.h
    ${MOD_DIR}/processors/volume/volumenoise.h
    ${MOD_DIR}/processors/volume/volumeselectormultichannel.h

    # utils
    ${MOD_DIR}/utils/utils.h
)

IF(VRN_MODULE_VESSELNETWORKANALYSIS)
    SET(MOD_CORE_HEADERS ${MOD_CORE_HEADERS}
        ${MOD_DIR}/processors/simulation/flowindicatordetection.h
    )
    SET(MOD_CORE_SOURCES ${MOD_CORE_SOURCES}
        ${MOD_DIR}/processors/simulation/flowindicatordetection.cpp
    )
ENDIF()

IF(VRN_MODULE_PLOTTING)
    SET(MOD_CORE_HEADERS ${MOD_CORE_HEADERS}
        ${MOD_DIR}/processors/simulation/flowindicatoranalysis.h
    )
    SET(MOD_CORE_SOURCES ${MOD_CORE_SOURCES}
        ${MOD_DIR}/processors/simulation/flowindicatoranalysis.cpp
    )
ENDIF()

################################################################################
# External dependency: OpenLB library
################################################################################

OPTION(VRN_FLOWSIMULATION_BUILD_OPENLB "Build OpenLB?" ON)
IF(VRN_FLOWSIMULATION_BUILD_OPENLB)
    IF(VRN_MSVC)
        # OpenLB was developed on and for POSIX systems, therefore, windows and MSVC are not supported.
        MESSAGE(FATAL_ERROR "OpenLB currently not supported by MSVC")
    ENDIF()

    IF(VRN_MODULE_OPENMP)
        ADD_DEFINITIONS("-DPARALLEL_MODE_OMP")
    ELSE()
        MESSAGE(WARNING "OpenMP module strongly recommended!")
    ENDIF()

    SET(OpenLB_DIR ${MOD_DIR}/ext/openlb)
    SET(OpenLB_INCLUDE_DIR ${OpenLB_DIR}/src)

    SET(OLB_BUILDTYPE "generic" CACHE STRING "OpenLB Build Type")
    SET_PROPERTY(CACHE OLB_BUILDTYPE PROPERTY STRINGS "precompiled" "generic")
    SET(OpenLB_LIBRARY_PATH ${OpenLB_DIR}/build/${OLB_BUILDTYPE}/lib/libolb.a)
    IF(${OLB_BUILDTYPE} MATCHES "precompiled")
        ADD_DEFINITIONS("-DOLB_PRECOMPILED")
    ENDIF()
    LIST(APPEND MOD_INCLUDE_DIRECTORIES ${OpenLB_INCLUDE_DIR})
    LIST(APPEND MOD_LIBRARIES ${OpenLB_LIBRARY_PATH})

    ADD_CUSTOM_TARGET(OpenLB COMMAND BUILDTYPE=${OLB_BUILDTYPE} make WORKING_DIRECTORY ${OpenLB_DIR})
    ADD_DEFINITIONS("-DVRN_FLOWSIMULATION_USE_OPENLB")

    SET(MOD_CORE_HEADERS ${MOD_CORE_HEADERS}
        ${MOD_DIR}/processors/features/wallshearstressextractor.h
        ${MOD_DIR}/processors/simulation/flowsimulation.h
        ${MOD_DIR}/processors/geometry/geometryinsidetest.h
    )
    SET(MOD_CORE_SOURCES ${MOD_CORE_SOURCES}
        ${MOD_DIR}/processors/features/wallshearstressextractor.cpp
        ${MOD_DIR}/processors/simulation/flowsimulation.cpp
        ${MOD_DIR}/processors/geometry/geometryinsidetest.cpp
    )
ENDIF()

################################################################################
# External dependency: halfedge
################################################################################

SET(MOD_CORE_HEADERS ${MOD_CORE_HEADERS}
    ${MOD_DIR}/ext/halfedge/trimesh.h
    ${MOD_DIR}/ext/halfedge/trimesh_types.h
)
SET(MOD_CORE_SOURCES ${MOD_CORE_SOURCES}
    ${MOD_DIR}/ext/halfedge/trimesh.cpp
)

SET(MOD_INSTALL_FILES
    ${MOD_DIR}/ext/halfedge/README
)

################################################################################
# External dependency: Octree
################################################################################

SET(MOD_CORE_HEADERS ${MOD_CORE_HEADERS}
    ${MOD_DIR}/ext/octree/Octree.hpp
    )

SET(MOD_INSTALL_FILES ${MOD_INSTALL_FILES}
    ${MOD_DIR}/ext/octree/LICENSE
)

# Deployment
SET(MOD_INSTALL_DIRECTORIES
    ${MOD_DIR}/scripts
    ${MOD_DIR}/workspaces
)

################################################################################
# External dependency: unwrap3d
################################################################################

SET(MOD_CORE_HEADERS ${MOD_CORE_HEADERS}
    ${MOD_DIR}/ext/unwrap3d/unwrap3d.h
)

SET(MOD_CORE_SOURCES ${MOD_CORE_SOURCES}
    ${MOD_DIR}/ext/unwrap3d/unwrap3d.cpp
)