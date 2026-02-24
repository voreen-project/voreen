
################################################################################
# Staging module resources
################################################################################
SET(MOD_CORE_MODULECLASS StagingModule)

SET(MOD_CORE_SOURCES
    ${MOD_DIR}/processors/alignedsliceproxygeometry.cpp
    ${MOD_DIR}/processors/arbitraryvolumeclipping.cpp
    ${MOD_DIR}/processors/clipregiongeometrycreator.cpp
    ${MOD_DIR}/processors/geometryslicerenderer.cpp
    ${MOD_DIR}/processors/interactiveregistrationwidget.cpp
    ${MOD_DIR}/processors/multislicerenderer.cpp
    ${MOD_DIR}/processors/multisliceviewer.cpp
    ${MOD_DIR}/processors/particles.cpp
    ${MOD_DIR}/processors/planegeometrycreator.cpp
    ${MOD_DIR}/processors/pong.cpp
    ${MOD_DIR}/processors/preintegrationtablerenderer.cpp
    ${MOD_DIR}/processors/registrationinitializer.cpp
    ${MOD_DIR}/processors/samplingpositiontransformation.cpp
    ${MOD_DIR}/processors/screenspaceambientocclusion.cpp
    ${MOD_DIR}/processors/tabbedview.cpp
    ${MOD_DIR}/processors/toucheventsimulator.cpp
    ${MOD_DIR}/processors/transfuncalphachannelanimation.cpp
    ${MOD_DIR}/processors/transfuncoverlay.cpp
    ${MOD_DIR}/processors/volumerealworldmapping.cpp
    ${MOD_DIR}/processors/volumeuncertaintymeasure.cpp

    ${MOD_DIR}/processors/slicepoints/slicepointrenderer2d.cpp
    ${MOD_DIR}/processors/slicepoints/slicepointrenderer3d.cpp
)

SET(MOD_CORE_HEADERS
    ${MOD_DIR}/processors/alignedsliceproxygeometry.h
    ${MOD_DIR}/processors/arbitraryvolumeclipping.h
    ${MOD_DIR}/processors/clipregiongeometrycreator.h
    ${MOD_DIR}/processors/geometryslicerenderer.h
    ${MOD_DIR}/processors/interactiveregistrationwidget.h
    ${MOD_DIR}/processors/multislicerenderer.h
    ${MOD_DIR}/processors/multisliceviewer.h
    ${MOD_DIR}/processors/particles.h
    ${MOD_DIR}/processors/planegeometrycreator.h
    ${MOD_DIR}/processors/pong.h
    ${MOD_DIR}/processors/preintegrationtablerenderer.h
    ${MOD_DIR}/processors/registrationinitializer.h
    ${MOD_DIR}/processors/samplingpositiontransformation.h
    ${MOD_DIR}/processors/screenspaceambientocclusion.h
    ${MOD_DIR}/processors/tabbedview.h
    ${MOD_DIR}/processors/toucheventsimulator.h
    ${MOD_DIR}/processors/transfuncalphachannelanimation.h
    ${MOD_DIR}/processors/transfuncoverlay.h
    ${MOD_DIR}/processors/volumerealworldmapping.h
    ${MOD_DIR}/processors/volumeuncertaintymeasure.h

    ${MOD_DIR}/processors/slicepoints/slicepointrenderer2d.h
    ${MOD_DIR}/processors/slicepoints/slicepointrenderer3d.h
)

SET(VRN_STAGING_SIMDRAYCASTER_AVAILABLE FALSE)
IF(WIN32)
    IF(CMAKE_SYSTEM_PROCESSOR MATCHES "^(AMD64|amd64|X64|x64|x86_64|i[3-6]86|x86)$")
        SET(VRN_STAGING_SIMDRAYCASTER_AVAILABLE TRUE)
    ENDIF()
ELSE()
    INCLUDE(CheckCXXCompilerFlag)
    IF(CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|amd64|i[3-6]86|x86)$")
        CHECK_CXX_COMPILER_FLAG("-msse3" VRN_STAGING_COMPILER_SUPPORTS_SSE3)
        IF(VRN_STAGING_COMPILER_SUPPORTS_SSE3)
            SET(VRN_STAGING_SIMDRAYCASTER_AVAILABLE TRUE)
        ENDIF()
    ENDIF()
ENDIF()

IF(VRN_STAGING_SIMDRAYCASTER_AVAILABLE)
    LIST(APPEND MOD_CORE_SOURCES
        ${MOD_DIR}/processors/simdraycaster/simdraycaster.cpp
        ${MOD_DIR}/utils/simdraycaster/jobqueue.cpp
        ${MOD_DIR}/utils/simdraycaster/memory.cpp
        ${MOD_DIR}/utils/simdraycaster/performancemetric.cpp
        ${MOD_DIR}/processors/simdraycaster/raycaster_sse3.cpp
    )
    LIST(APPEND MOD_CORE_HEADERS
        ${MOD_DIR}/processors/simdraycaster/raycast_generic.h
        ${MOD_DIR}/processors/simdraycaster/simdraycaster.h
        ${MOD_DIR}/utils/simdraycaster/brickedvolume.h
        ${MOD_DIR}/utils/simdraycaster/brickedvolumebase.h
        ${MOD_DIR}/utils/simdraycaster/jobqueue.h
        ${MOD_DIR}/utils/simdraycaster/memory.h
        ${MOD_DIR}/utils/simdraycaster/performancemetric.h
    )

    LIST(APPEND VRN_MODULE_DEFINITIONS "-DVRN_STAGING_HAS_SIMDRAYCASTER")
    LIST(APPEND VRN_MODULE_DEFINITIONS "-DSIMD_SSE3")

    IF(WIN32)
        LIST(APPEND MOD_CORE_SOURCES ${MOD_DIR}/processors/simdraycaster/raycaster_sse41.cpp)
        LIST(APPEND VRN_MODULE_DEFINITIONS "-DSIMD_SSE41")
    ELSE()
        SET_PROPERTY(SOURCE ${MOD_DIR}/processors/simdraycaster/raycaster_sse3.cpp APPEND PROPERTY COMPILE_OPTIONS "-msse3")
        OPTION(VRN_USE_SSE41 "Use the sse 4.1 Instruction set extension" ON)
        IF(VRN_USE_SSE41)
            CHECK_CXX_COMPILER_FLAG("-msse4.1" VRN_STAGING_COMPILER_SUPPORTS_SSE41)
            IF(VRN_STAGING_COMPILER_SUPPORTS_SSE41)
                LIST(APPEND MOD_CORE_SOURCES ${MOD_DIR}/processors/simdraycaster/raycaster_sse41.cpp)
                SET_PROPERTY(SOURCE ${MOD_DIR}/processors/simdraycaster/raycaster_sse41.cpp APPEND PROPERTY COMPILE_OPTIONS "-msse4.1")
                LIST(APPEND VRN_MODULE_DEFINITIONS "-DSIMD_SSE41")
            ENDIF()
        ENDIF()
    ENDIF()
ELSE()
    MESSAGE(STATUS "Staging module: SIMDRaycaster disabled (no SSE support detected)")
ENDIF()

# deployment
SET(MOD_INSTALL_DIRECTORIES
    ${MOD_DIR}/glsl
    ${MOD_DIR}/textures
)
