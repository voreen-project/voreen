################################################################################
# Bovenkamp module resources
################################################################################

IF(NOT VRN_MODULE_BASE)
    MESSAGE(FATAL_ERROR "Bovenkamp Module requires Base Module")
ENDIF()
IF(NOT VRN_MODULE_DEVIL)
    MESSAGE (FATAL_ERROR "Bovenkamp Module requires Devil Module")
ENDIF()
IF(NOT VRN_MODULE_FFMPEG)
    MESSAGE (FATAL_ERROR "Bovenkamp Module requires Ffmpeg Module")
ENDIF()
IF(NOT VRN_MODULE_FLOWANALYSIS)
    MESSAGE(FATAL_ERROR "Bovenkamp Module requires FlowAnalysis Module")
ENDIF()
    
SET(MOD_CORE_MODULECLASS BovenkampModule)

# deployment
SET(MOD_INSTALL_DIRECTORIES
    ${MOD_DIR}/glsl
    ${MOD_DIR}/workspaces
)
    
SET(MOD_CORE_SOURCES
    #IO
    ${MOD_DIR}/io/volumediskpb.cpp
    
    #Processors
    ${MOD_DIR}/processors/pbreader.cpp
    ${MOD_DIR}/processors/pblinkcontrol.cpp
)

SET(MOD_CORE_HEADERS
    #IO
    ${MOD_DIR}/io/volumediskpb.h

    #Processors
    ${MOD_DIR}/processors/pbreader.h
    ${MOD_DIR}/processors/pblinkcontrol.h
)
