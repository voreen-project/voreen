################################################################################
# Core module resources 
################################################################################
SET(MOD_CORE_MODULECLASS DynamicGLSLModule)

SET(MOD_CORE_SOURCES
    ${MOD_DIR}/processors/dynamicglslprocessor.cpp
    ${MOD_DIR}/parser/glslparser.cpp
    ${MOD_DIR}/parser/glslprogram.cpp
    ${MOD_DIR}/parser/glsllexer.cpp
)

SET(MOD_CORE_HEADERS
    ${MOD_DIR}/processors/dynamicglslprocessor.h
    ${MOD_DIR}/parser/glslparser.h
    ${MOD_DIR}/parser/glslprogram.h
    ${MOD_DIR}/parser/glsllexer.h
)
   
# deployment
SET(MOD_INSTALL_DIRECTORIES
    ${MOD_DIR}/glsl
)


################################################################################
# Qt module resources 
################################################################################
SET(MOD_QT_MODULECLASS DynamicGLSLModuleQt)

SET(MOD_QT_SOURCES
    ${MOD_DIR}/qt/dynamicglslwidget.cpp
    ${MOD_DIR}/qt/dynamicglslprocessorwidgetfactory.cpp
)
    
SET(MOD_QT_HEADERS
    ${MOD_DIR}/qt/dynamicglslwidget.h
)

SET(MOD_QT_HEADERS_NONMOC
    ${MOD_DIR}/qt/dynamicglslprocessorwidgetfactory.h
)
