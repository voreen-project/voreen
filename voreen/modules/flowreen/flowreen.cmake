################################################################################
# Core module resources
################################################################################
SET(MOD_CORE_MODULECLASS FlowreenModule)

SET(MOD_CORE_SOURCES
    ${MOD_DIR}/flowreenmodule.cpp

    # datastructures
    ${MOD_DIR}/datastructures/streamline.cpp
    ${MOD_DIR}/datastructures/streamlinebundle.cpp
    ${MOD_DIR}/datastructures/streamlinelist.cpp
    ${MOD_DIR}/datastructures/streamlinelistbase.cpp
    ${MOD_DIR}/datastructures/streamlinelistdecorator.cpp
    ${MOD_DIR}/datastructures/streamlinelistobserver.cpp

    # ports
    ${MOD_DIR}/ports/streamlinelistport.cpp

    # processors
    ${MOD_DIR}/processors/flowdirectionoverlay.cpp
    ${MOD_DIR}/processors/streamlinerenderer3d.cpp
    ${MOD_DIR}/processors/geometry/geometryclose.cpp
    ${MOD_DIR}/processors/geometry/geometryoffsetremove.cpp
    ${MOD_DIR}/processors/streamline/streamlinecombine.cpp
    ${MOD_DIR}/processors/streamline/streamlinecreator.cpp
    ${MOD_DIR}/processors/streamline/streamlinerotation.cpp
    ${MOD_DIR}/processors/streamline/streamlinesave.cpp
    ${MOD_DIR}/processors/streamline/streamlineselector.cpp
    ${MOD_DIR}/processors/streamline/streamlinesource.cpp
    ${MOD_DIR}/processors/streamline/streamlinetoboundingbox.cpp
    ${MOD_DIR}/processors/streamline/pathlinecreator.cpp

    # utils
    ${MOD_DIR}/utils/streamlinebundledetectorbackgroundthread.cpp
    ${MOD_DIR}/utils/streamlinecreatorbackgroundthread.cpp
    ${MOD_DIR}/utils/pathlinecreatorbackgroundthread.cpp
)

IF(VRN_OPENGL_COMPATIBILITY_PROFILE)
    #check if core profile compatible or port
    LIST(APPEND MOD_CORE_SOURCES
        # datastructures
        ${MOD_DIR}/datastructures/deprecated/flow2d.cpp
        ${MOD_DIR}/datastructures/deprecated/flow3d.cpp
        ${MOD_DIR}/datastructures/deprecated/simpletexture.cpp
        ${MOD_DIR}/datastructures/deprecated/streamlinetexture.cpp
        ${MOD_DIR}/datastructures/deprecated/volumeflow3d.cpp
        ${MOD_DIR}/datastructures/deprecated/volumeoperatorflowmagnitude.cpp

        # io
        ${MOD_DIR}/io/flowreader.cpp

        # processors
        ${MOD_DIR}/processors/flowarrowrenderer2d.cpp
        ${MOD_DIR}/processors/flowarrowrenderer3d.cpp
        ${MOD_DIR}/processors/flowmagnitudes3d.cpp
        ${MOD_DIR}/processors/floworthogonalslicerenderer.cpp
        ${MOD_DIR}/processors/flowreenadapter.cpp
        ${MOD_DIR}/processors/flowreenprocessor.cpp
        ${MOD_DIR}/processors/flowslicerenderer.cpp
        ${MOD_DIR}/processors/flowslicerenderer2d.cpp
        ${MOD_DIR}/processors/flowslicerenderer3d.cpp
        ${MOD_DIR}/processors/flowstreamlinestexture3d.cpp
        ${MOD_DIR}/processors/pathlinerenderer3d.cpp
        ${MOD_DIR}/processors/pathlinerenderer3d.cpp

        # utils
        ${MOD_DIR}/utils/colorcodingability.cpp
        ${MOD_DIR}/utils/flowmath.cpp
    )
ENDIF()
    
SET(MOD_CORE_HEADERS
    ${MOD_DIR}/flowreenmodule.h

    # datastructures
    ${MOD_DIR}/datastructures/streamline.h
    ${MOD_DIR}/datastructures/streamlinebundle.h
    ${MOD_DIR}/datastructures/streamlinelist.h
    ${MOD_DIR}/datastructures/streamlinelistbase.h
    ${MOD_DIR}/datastructures/streamlinelistdecorator.h
    ${MOD_DIR}/datastructures/streamlinelistobserver.h

    # ports
    ${MOD_DIR}/ports/streamlinelistport.h

    # processors
    ${MOD_DIR}/processors/flowdirectionoverlay.h
    ${MOD_DIR}/processors/streamlinerenderer3d.h
    ${MOD_DIR}/processors/geometry/geometryclose.h
    ${MOD_DIR}/processors/geometry/geometryoffsetremove.h
    ${MOD_DIR}/processors/streamline/streamlinecombine.h
    ${MOD_DIR}/processors/streamline/streamlinecreator.h
    ${MOD_DIR}/processors/streamline/streamlinerotation.h
    ${MOD_DIR}/processors/streamline/streamlinesave.h
    ${MOD_DIR}/processors/streamline/streamlineselector.h
    ${MOD_DIR}/processors/streamline/streamlinesource.h
    ${MOD_DIR}/processors/streamline/streamlinetoboundingbox.h
    ${MOD_DIR}/processors/streamline/pathlinecreator.h

    # utils
    ${MOD_DIR}/utils/streamlinebundledetectorbackgroundthread.h
    ${MOD_DIR}/utils/streamlinecreatorbackgroundthread.h
    ${MOD_DIR}/utils/pathlinecreatorbackgroundthread.h
)
 
IF(VRN_OPENGL_COMPATIBILITY_PROFILE)
    #check if core profile compatible or port
    LIST(APPEND MOD_CORE_HEADERS
        # # datastructures
        ${MOD_DIR}/datastructures/deprecated/flow2d.h
        ${MOD_DIR}/datastructures/deprecated/flow3d.h
        ${MOD_DIR}/datastructures/deprecated/simpletexture.h
        ${MOD_DIR}/datastructures/deprecated/streamlinetexture.h
        ${MOD_DIR}/datastructures/deprecated/volumeflow3d.h
        ${MOD_DIR}/datastructures/deprecated/volumeoperatorflowmagnitude.h
        ${MOD_DIR}/datastructures/deprecated/volumeoperatorintensitymask.h

        # io
        ${MOD_DIR}/io/flowreader.h

        # processors
        ${MOD_DIR}/processors/flowarrowrenderer2d.h
        ${MOD_DIR}/processors/flowarrowrenderer3d.h
        ${MOD_DIR}/processors/flowmagnitudes3d.h
        ${MOD_DIR}/processors/floworthogonalslicerenderer.h
        ${MOD_DIR}/processors/flowreenadapter.h
        ${MOD_DIR}/processors/flowreenprocessor.h
        ${MOD_DIR}/processors/flowslicerenderer.h
        ${MOD_DIR}/processors/flowslicerenderer2d.h
        ${MOD_DIR}/processors/flowslicerenderer3d.h
        ${MOD_DIR}/processors/flowstreamlinestexture3d.h
        ${MOD_DIR}/processors/pathlinerenderer3d.h

        # utils
        ${MOD_DIR}/utils/colorcodingability.h
        ${MOD_DIR}/utils/flowmath.h
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
ENDIF()

IF(VRN_FLOWREEN_BUILD_OPENLB)
    SET(MOD_CORE_HEADERS ${MOD_CORE_HEADERS}
        ${MOD_DIR}/processors/flowsimulation.h
        ${MOD_DIR}/processors/geometry/implicitrepresentation.h

        ${MOD_DIR}/utils/geometryconverter.h
    )
    SET(MOD_CORE_SOURCES ${MOD_CORE_SOURCES}
        ${MOD_DIR}/processors/flowsimulation.cpp
        ${MOD_DIR}/processors/geometry/implicitrepresentation.cpp

        ${MOD_DIR}/utils/geometryconverter.cpp
    )
ENDIF()

################################################################################
# External dependency: halfedge
################################################################################

SET(MOD_CORE_HEADERS ${MOD_CORE_HEADERS}
    ${MOD_DIR}/processors/flowsimulation.h
    ${MOD_DIR}/ext/halfedge/trimesh.h
    ${MOD_DIR}/ext/halfedge/trimesh_types.h
)
SET(MOD_CORE_SOURCES ${MOD_CORE_SOURCES}
    ${MOD_DIR}/ext/halfedge/trimesh.cpp
)

SET(MOD_INSTALL_FILES
    ${MOD_DIR}/ext/halfedge/README
)


# Deployment
SET(MOD_INSTALL_DIRECTORIES
    ${MOD_DIR}/glsl
    ${MOD_DIR}/data
    ${MOD_DIR}/transferfuncs
)
