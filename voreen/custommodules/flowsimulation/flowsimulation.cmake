################################################################################
# Core module resources
################################################################################

IF(NOT VRN_MODULE_BIGDATAIMAGEPROCESSING)
    MESSAGE(FATAL_ERROR "FlowSimulation Module requires big data image processing Module")
ENDIF()

SET(MOD_CORE_MODULECLASS FlowSimulationModule)

SET(MOD_CORE_SOURCES
    ${MOD_DIR}/flowsimulationmodule.cpp

    # datastructures
    ${MOD_DIR}/datastructures/flowparameters.cpp

    # flow features
    ${MOD_DIR}/flowfeatures/magnitudefeature.cpp

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

    # flow features
    ${MOD_DIR}/flowfeatures/flowfeature.h
    ${MOD_DIR}/flowfeatures/magnitudefeature.h

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
        ${MOD_DIR}/processors/simulation/flowindicatorselection.h
    )
    SET(MOD_CORE_SOURCES ${MOD_CORE_SOURCES}
        ${MOD_DIR}/processors/simulation/flowindicatordetection.cpp
    )
ENDIF()

################################################################################
# External dependency: OpenLB library
################################################################################

OPTION(VRN_FLOWREEN_BUILD_OPENLB "Build OpenLB?" ON)
IF(VRN_FLOWREEN_BUILD_OPENLB)
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
    SET(OpenLB_LIBRARY_PATH ${OpenLB_DIR}/build/precompiled/lib/libolb.a)
    LIST(APPEND MOD_INCLUDE_DIRECTORIES ${OpenLB_INCLUDE_DIR})
    LIST(APPEND MOD_LIBRARIES ${OpenLB_LIBRARY_PATH})

    ADD_CUSTOM_TARGET(OpenLB COMMAND make WORKING_DIRECTORY ${OpenLB_DIR})
    ADD_DEFINITIONS("-DFLOWREEN_USE_OPENLB")

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
)
