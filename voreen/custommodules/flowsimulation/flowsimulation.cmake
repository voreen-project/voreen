################################################################################
# Core module resources
################################################################################

SET(MOD_CORE_MODULECLASS FlowSimulationModule)

SET(MOD_CORE_SOURCES
    ${MOD_DIR}/flowsimulationmodule.cpp

    # datastructures
    ${MOD_DIR}/datastructures/flowparameters.cpp

    # ports
    ${MOD_DIR}/ports/flowparametrizationport.cpp

    # processors
    ${MOD_DIR}/processors/geometry/geometryclose.cpp
    ${MOD_DIR}/processors/geometry/geometryoffsetremove.cpp
    ${MOD_DIR}/processors/render/unalignedsliceviewer.cpp
    ${MOD_DIR}/processors/simulation/flowcharacteristics.cpp
    ${MOD_DIR}/processors/simulation/flowensemblecreator.cpp
    ${MOD_DIR}/processors/simulation/flowindicatorselection.cpp
    ${MOD_DIR}/processors/simulation/flowindicatorrenderer.cpp
    ${MOD_DIR}/processors/simulation/flowparametrization.cpp
    ${MOD_DIR}/processors/simulation/flowsimulationcluster.cpp
    ${MOD_DIR}/processors/volume/volumelistadapter.cpp

    # utils
    ${MOD_DIR}/utils/geometryconverter.cpp
)

SET(MOD_CORE_HEADERS
    ${MOD_DIR}/flowsimulationmodule.h

    # datastructures
    ${MOD_DIR}/datastructures/flowparameters.h

    # ports
    ${MOD_DIR}/ports/flowparametrizationport.h

    # processors
    ${MOD_DIR}/processors/geometry/geometryclose.h
    ${MOD_DIR}/processors/geometry/geometryoffsetremove.h
    ${MOD_DIR}/processors/render/unalignedsliceviewer.h
    ${MOD_DIR}/processors/simulation/flowcharacteristics.h
    ${MOD_DIR}/processors/simulation/flowensemblecreator.h
    ${MOD_DIR}/processors/simulation/flowindicatorselection.h
    ${MOD_DIR}/processors/simulation/flowindicatorrenderer.h
    ${MOD_DIR}/processors/simulation/flowparametrization.h
    ${MOD_DIR}/processors/simulation/flowsimulationcluster.h
    ${MOD_DIR}/processors/volume/volumelistadapter.h

    # utils
    ${MOD_DIR}/utils/geometryconverter.h
)

IF(VRN_MODULE_VESSELTOPOLOGY)
    SET(MOD_CORE_HEADERS ${MOD_CORE_HEADERS}
        ${MOD_DIR}/processors/simulation/flowindicatordetection.h
    )
    SET(MOD_CORE_SOURCES ${MOD_CORE_SOURCES}
        ${MOD_DIR}/processors/simulation/flowindicatordetection.cpp
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

    SET(OLB_BUILDTYPE "precompiled" CACHE STRING "OpenLB Build Type")
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
        ${MOD_DIR}/processors/simulation/flowsimulation.h
        ${MOD_DIR}/processors/geometry/implicitrepresentation.h
    )
    SET(MOD_CORE_SOURCES ${MOD_CORE_SOURCES}
        ${MOD_DIR}/processors/simulation/flowsimulation.cpp
        ${MOD_DIR}/processors/geometry/implicitrepresentation.cpp
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


# TODO: Add mailio or alternative to parse palma feedback sent per mail.
#SET(VRN_MAILIO_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/modules/${module_dir}/ext/mailio)
#ADD_SUBDIRECTORY(${VRN_MAILIO_DIRECTORY})
#LIST(APPEND MOD_LIBRARIES "mailio")


# Deployment
SET(MOD_INSTALL_DIRECTORIES
    ${MOD_DIR}/scripts
    ${MOD_DIR}/workspaces
)
