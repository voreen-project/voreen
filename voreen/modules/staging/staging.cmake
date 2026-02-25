
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

# deployment
SET(MOD_INSTALL_DIRECTORIES
    ${MOD_DIR}/glsl
    ${MOD_DIR}/textures
)
