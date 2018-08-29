
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

    ${MOD_DIR}/processors/simdraycaster/apvtestraycaster.cpp
    ${MOD_DIR}/processors/simdraycaster/simdraycaster.cpp
	
    ${MOD_DIR}/utils/simdraycaster/jobqueue.cpp
    ${MOD_DIR}/utils/simdraycaster/memory.cpp
    ${MOD_DIR}/utils/simdraycaster/performancemetric.cpp
)

if (UNIX)
    SET_SOURCE_FILES_PROPERTIES(${MOD_DIR}/processors/simdraycaster/raycaster_sse41.cpp PROPERTIES COMPILE_FLAGS "${CMAKE_C_FLAGS} -msse4.1")
    SET_SOURCE_FILES_PROPERTIES(${MOD_DIR}/processors/simdraycaster/raycaster_sse3.cpp PROPERTIES COMPILE_FLAGS "${CMAKE_C_FLAGS} -msse3")
endif (UNIX) 

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

    ${MOD_DIR}/processors/simdraycaster/apvtestraycaster.h
    ${MOD_DIR}/processors/simdraycaster/raycast_generic.h
    ${MOD_DIR}/processors/simdraycaster/simdraycaster.h
	
    ${MOD_DIR}/utils/simdraycaster/brickedvolume.h
    ${MOD_DIR}/utils/simdraycaster/brickedvolumebase.h
    ${MOD_DIR}/utils/simdraycaster/jobqueue.h
    ${MOD_DIR}/utils/simdraycaster/memory.h
    ${MOD_DIR}/utils/simdraycaster/performancemetric.h
)

IF(WIN32)
    SET(MOD_CORE_SOURCES
        ${MOD_CORE_SOURCES}
        ${MOD_DIR}/processors/simdraycaster/raycaster_sse41.cpp
        ${MOD_DIR}/processors/simdraycaster/raycaster_sse3.cpp
    )
    LIST(APPEND VRN_MODULE_DEFINITIONS "-DSIMD_SSE3")
    LIST(APPEND VRN_MODULE_DEFINITIONS "-DSIMD_SSE41")
ELSE(WIN32)
    OPTION (VRN_USE_SSE41 
            "Use the sse 4.1 Instruction set extension" ON) 
    IF(VRN_USE_SSE41)
        SET(MOD_CORE_SOURCES
            ${MOD_CORE_SOURCES}
            ${MOD_DIR}/processors/simdraycaster/raycaster_sse3.cpp
            ${MOD_DIR}/processors/simdraycaster/raycaster_sse41.cpp
        )
        LIST(APPEND VRN_MODULE_DEFINITIONS "-DSIMD_SSE3")
        LIST(APPEND VRN_MODULE_DEFINITIONS "-DSIMD_SSE41 -msse4.1")
    ELSE(VRN_USE_SSE41)
        SET(MOD_CORE_SOURCES
            ${MOD_CORE_SOURCES}
            ${MOD_DIR}/processors/simdraycaster/raycaster_sse3.cpp
        )
        LIST(APPEND VRN_MODULE_DEFINITIONS "-DSIMD_SSE3 -msse3")
    ENDIF(VRN_USE_SSE41)
ENDIF(WIN32)

# deployment
SET(MOD_INSTALL_DIRECTORIES
    ${MOD_DIR}/glsl
    ${MOD_DIR}/textures
)
