################################################################################
# Core module resources 
################################################################################
SET(MOD_CORE_MODULECLASS StereoscopyModule)

SET(MOD_CORE_SOURCES
    #Processors        
    ${MOD_DIR}/processors/stereocanvasrenderer.cpp   
)

SET(MOD_CORE_HEADERS
    #Processors    
    ${MOD_DIR}/processors/stereocanvasrenderer.h  
)

################################################################################
# Qt module resources 
################################################################################
SET(MOD_QT_MODULECLASS StereoscopyModuleQt)

SET(MOD_QT_SOURCES
    ${MOD_DIR}/qt/processor/stereocanvasrendererwidget.cpp
    ${MOD_DIR}/qt/processor/stereoscopyprocessorwidgetfactory.cpp    
)  
    
SET(MOD_QT_HEADERS
    ${MOD_DIR}/qt/processor/stereocanvasrendererwidget.h
)

SET(MOD_QT_HEADERS_NONMOC
    ${MOD_DIR}/qt/processor/stereoscopyprocessorwidgetfactory.h
)

################################################################################
# Deployment
################################################################################
SET(MOD_INSTALL_DIRECTORIES
    ${MOD_DIR}/glsl
    ${MOD_DIR}/workspaces
)
