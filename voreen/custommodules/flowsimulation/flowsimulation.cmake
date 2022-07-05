################################################################################
# Core module resources
################################################################################

IF(NOT VRN_MODULE_FLOWANALYSIS)
    MESSAGE(FATAL_ERROR "FlowSimulation Module requires Flow Analysis Module")
ENDIF()

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
    ${MOD_DIR}/datastructures/flowsimulationconfig.cpp
    ${MOD_DIR}/datastructures/volumeramremappingproxy.cpp
    
    # ports
    ${MOD_DIR}/ports/flowsimulationconfigport.cpp

    # processors
    ${MOD_DIR}/processors/features/wallshearstress.cpp
    ${MOD_DIR}/processors/geometry/geometryclose.cpp
    ${MOD_DIR}/processors/geometry/geometrysmoothnormals.cpp
    ${MOD_DIR}/processors/render/unalignedsliceviewer.cpp
    ${MOD_DIR}/processors/simulation/flowcharacteristics.cpp
    ${MOD_DIR}/processors/simulation/flowensemblecreator.cpp
    ${MOD_DIR}/processors/simulation/flowindicatorrenderer.cpp
    ${MOD_DIR}/processors/simulation/flowparametrizationensemble.cpp
    ${MOD_DIR}/processors/simulation/flowparametrizationrun.cpp
    ${MOD_DIR}/processors/simulation/flowsimulationcluster.cpp
    ${MOD_DIR}/processors/simulation/flowsimulationgeometry.cpp
    ${MOD_DIR}/processors/volume/connectedcomponentselector.cpp
#    ${MOD_DIR}/processors/volume/debugvolumes.cpp
    ${MOD_DIR}/processors/volume/flowtestdatagenerator.cpp
    ${MOD_DIR}/processors/volume/phaseunwrapping.cpp
    ${MOD_DIR}/processors/volume/vectordecompose.cpp
    ${MOD_DIR}/processors/volume/volumeapplyrealworldmapping.cpp
    ${MOD_DIR}/processors/volume/volumelistadapter.cpp
    ${MOD_DIR}/processors/volume/volumelistaggregate.cpp
    ${MOD_DIR}/processors/volume/volumelistmerger.cpp
    ${MOD_DIR}/processors/volume/volumelistoffset.cpp
    ${MOD_DIR}/processors/volume/volumelistrealworldmapping.cpp
    ${MOD_DIR}/processors/volume/volumelisttimestep.cpp
    ${MOD_DIR}/processors/volume/volumemerger.cpp
    ${MOD_DIR}/processors/volume/volumenoise.cpp
    ${MOD_DIR}/processors/volume/volumeresampleproxy.cpp
    ${MOD_DIR}/processors/volume/volumeselectormultichannel.cpp
    ${MOD_DIR}/processors/volume/volumetimestep.cpp

    # utils
    ${MOD_DIR}/utils/serializationhelper.cpp
    ${MOD_DIR}/utils/utils.cpp
    
    # openlb
    ${MOD_DIR}/ext/openlb/voreen/openlb_parameters.cpp
)

SET(MOD_CORE_HEADERS
    ${MOD_DIR}/flowsimulationmodule.h

    # datastructures
    ${MOD_DIR}/datastructures/flowsimulationconfig.h
    ${MOD_DIR}/datastructures/volumeramremappingproxy.h

    # ports
    ${MOD_DIR}/ports/flowsimulationconfigport.h

    # processors
    ${MOD_DIR}/processors/features/wallshearstress.h
    ${MOD_DIR}/processors/geometry/geometryclose.h
    ${MOD_DIR}/processors/geometry/geometrysmoothnormals.h
    ${MOD_DIR}/processors/render/unalignedsliceviewer.h
    ${MOD_DIR}/processors/simulation/flowcharacteristics.h
    ${MOD_DIR}/processors/simulation/flowensemblecreator.h
    ${MOD_DIR}/processors/simulation/flowindicatorrenderer.h
    ${MOD_DIR}/processors/simulation/flowparametrizationensemble.h
    ${MOD_DIR}/processors/simulation/flowparametrizationrun.h
    ${MOD_DIR}/processors/simulation/flowsimulationcluster.h
    ${MOD_DIR}/processors/simulation/flowsimulationgeometry.h
    ${MOD_DIR}/processors/volume/connectedcomponentselector.h
#    ${MOD_DIR}/processors/volume/debugvolumes.h
    ${MOD_DIR}/processors/volume/flowtestdatagenerator.h
    ${MOD_DIR}/processors/volume/phaseunwrapping.h
    ${MOD_DIR}/processors/volume/vectordecompose.h
    ${MOD_DIR}/processors/volume/volumeapplyrealworldmapping.h
    ${MOD_DIR}/processors/volume/volumelistadapter.h
    ${MOD_DIR}/processors/volume/volumelistaggregate.h
    ${MOD_DIR}/processors/volume/volumelistmerger.h
    ${MOD_DIR}/processors/volume/volumelistoffset.h
    ${MOD_DIR}/processors/volume/volumelistrealworldmapping.h
    ${MOD_DIR}/processors/volume/volumelisttimestep.h
    ${MOD_DIR}/processors/volume/volumemerger.h
    ${MOD_DIR}/processors/volume/volumenoise.h
    ${MOD_DIR}/processors/volume/volumeresampleproxy.h
    ${MOD_DIR}/processors/volume/volumeselectormultichannel.h
    ${MOD_DIR}/processors/volume/volumetimestep.h

    # utils
    ${MOD_DIR}/utils/serializationhelper.h
    ${MOD_DIR}/utils/utils.h
    
    # openlb
    ${MOD_DIR}/ext/openlb/voreen/openlb_parameters.h
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
        ${MOD_DIR}/processors/plotting/flowindicatoranalysis.h
        ${MOD_DIR}/processors/plotting/flowprofilestacking.h
        ${MOD_DIR}/processors/plotting/regionofinterestanalysis.h
        ${MOD_DIR}/processors/plotting/roianalysis.h
    )
    SET(MOD_CORE_SOURCES ${MOD_CORE_SOURCES}
        ${MOD_DIR}/processors/plotting/flowindicatoranalysis.cpp
        ${MOD_DIR}/processors/plotting/flowprofilestacking.cpp
        ${MOD_DIR}/processors/plotting/regionofinterestanalysis.cpp
        ${MOD_DIR}/processors/plotting/roianalysis.cpp
    )
ENDIF()

################################################################################
# External dependency: OpenLB library
################################################################################

IF(VRN_MSVC)
    # OpenLB was developed on and for POSIX systems, therefore, windows/MSVC is currently not supported.
    OPTION(VRN_FLOWSIMULATION_BUILD_OPENLB "Build OpenLB?" OFF)
ELSE()
    OPTION(VRN_FLOWSIMULATION_BUILD_OPENLB "Build OpenLB?" ON)
ENDIF()

IF(VRN_FLOWSIMULATION_BUILD_OPENLB)
    IF(VRN_MSVC)
        # We might add support in the future, but for now we have to output an error message.
        MESSAGE(FATAL_ERROR "OpenLB currently not supported by MSVC")
    ENDIF()

    # OpenLB requires c++17 standard.
    SET(CMAKE_CXX_STANDARD 17)
    SET(CMAKE_CXX_STANDARD_REQUIRED ON)

    ADD_DEFINITIONS("-DPLATFORM_CPU_SISD")
    IF(VRN_MODULE_OPENMP)
        ADD_DEFINITIONS("-DPARALLEL_MODE_OMP")
    ELSE()
        MESSAGE(WARNING "OpenMP module strongly recommended!")
    ENDIF()

    SET(OLB_OPTIONS "")

    # Set CXX FLAGS.
    SET(OLB_CXXFLAGS -std=c++17)
    IF(CMAKE_BUILD_TYPE STREQUAL "Debug")
        SET(OLB_CXXFLAGS ${OLB_CXXFLAGS} -g -DOLB_DEBUG)
    ELSEIF(CMAKE_BUILD_TYPE STREQUAL "Release")
        SET(OLB_CXXFLAGS ${OLB_CXXFLAGS} -O3)
    ENDIF()

    # Set parallel mode.
    SET(OLB_PARALLEL_MODE "OFF" CACHE STRING "OpenLB Parallel Mode")
    SET_PROPERTY(CACHE OLB_PARALLEL_MODE PROPERTY STRINGS "OFF" "MPI" "OMP" "HYBRID")
    LIST(APPEND OLB_OPTIONS "PARALLEL_MODE=${OLB_PARALLEL_MODE}")
    IF(${OLB_PARALLEL_MODE} STREQUAL "OMP" OR ${OLB_PARALLEL_MODE} STREQUAL "HYBRID")
        LIST(APPEND OLB_OPTIONS "OMPFLAGS=-fopenmp")
    ENDIF()

    # Set default compiler.
    SET(OLB_CXX "g++")
    SET(OLB_CC "gcc")

    IF(${OLB_PARALLEL_MODE} STREQUAL "MPI" OR ${OLB_PARALLEL_MODE} STREQUAL "HYBRID")
        SET(OLB_CXX "mpic++")
    ENDIF()

    # Set platforms.
    SET(OLB_PLATFORMS "CPU_SISD") # mandatory

    OPTION(OLB_PLATFORM_CPU_SIMD "Enable OpenLB SIMD Platform?" OFF)
    IF(OLB_PLATFORM_CPU_SIMD)
        LIST(APPEND OLB_PLATFORMS "CPU_SIMD")
    ENDIF()

    OPTION(OLB_PLATFORM_GPU_CUDA "Enable OpenLB CUDA Platform?" OFF)
    IF(OLB_PLATFORM_GPU_CUDA)
        SET(OLB_CXX "nvcc") # Overrides mpic++, which is intentional.
        SET(OLB_CC "nvcc")

        LIST(APPEND OLB_PLATFORMS "GPU_CUDA")

#        # Arch is auto-detected by nvcc (assuming we compile on the same architecture that we run the code on).
#        SET(OLB_PLATFORM_CUDA_ARCH 60 CACHE STRING "CUDA ARCH - see rules.mk")
#        SET_PROPERTY(CACHE OLB_PLATFORM_CUDA_ARCH PROPERTY STRINGS 20 30 35 37 50 52 53 60 61 62 70 72 75 80 86 87)
#        LIST(APPEND OLB_OPTIONS "CUDA_ARCH=${OLB_PLATFORM_CUDA_ARCH}")

        IF(${OLB_PARALLEL_MODE} STREQUAL "MPI" OR ${OLB_PARALLEL_MODE} STREQUAL "HYBRID")
            SET(OLB_CXXFLAGS ${OLB_CXXFLAGS} -lmpi_cxx -lmpi)
        ENDIF()
    ELSE()
        SET(OLB_CXXFLAGS ${OLB_CXXFLAGS} -Wall -march=native -mtune=native)
    ENDIF()

    LIST(APPEND OLB_OPTIONS "CXX=${OLB_CXX}" "CC=${OLB_CC}" "CXXFLAGS=\"${OLB_CXXFLAGS}\"" "PLATFORMS=\"${OLB_PLATFORMS}\"")

    SET(OpenLB_DIR ${MOD_DIR}/ext/openlb)
    SET(OpenLB_INCLUDE_DIR ${OpenLB_DIR}/src)
    LIST(APPEND MOD_INCLUDE_DIRECTORIES ${OpenLB_INCLUDE_DIR})

    # Add OpenLB third party libraries.
    LIST(APPEND MOD_INCLUDE_DIRECTORIES ${OpenLB_DIR}/external/zlib)
    LIST(APPEND MOD_INCLUDE_DIRECTORIES ${OpenLB_DIR}/external/tinyxml)

    # It currently seems to be not possible to execute the build before all other builds,
    # so the user needs to manually build the OpenLB target first..
    ADD_CUSTOM_TARGET(OpenLB COMMAND ${OLB_OPTIONS} make WORKING_DIRECTORY ${OpenLB_DIR}/voreen)
    ADD_DEFINITIONS("-DVRN_FLOWSIMULATION_USE_OPENLB")

    SET(MOD_CORE_HEADERS ${MOD_CORE_HEADERS}
        ${MOD_DIR}/processors/simulation/flowsimulation.h
        ${MOD_DIR}/processors/geometry/geometryinsidetest.h
    )
    SET(MOD_CORE_SOURCES ${MOD_CORE_SOURCES}
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