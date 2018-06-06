
################################################################################
# External dependency: JAMA lib 
################################################################################
SET(MOD_INCLUDE_DIRECTORIES "${MOD_DIR}/ext/jama/include")


################################################################################
# Core module resources 
################################################################################
SET(MOD_CORE_MODULECLASS VolumeLabelingModule)
SET(MOD_REQUIRE_OPENGL_COMPATIBILITY_PROFILE ON)

SET(MOD_CORE_SOURCES
    ${MOD_DIR}/processors/idraycaster.cpp
    ${MOD_DIR}/processors/volumelabeling.cpp
    ${MOD_DIR}/utils/labelingmath.cpp
)

SET(MOD_CORE_HEADERS
    ${MOD_DIR}/processors/idraycaster.h
    ${MOD_DIR}/processors/volumelabeling.h
    ${MOD_DIR}/utils/labelingmath.h
)
  
# deployment
SET(MOD_INSTALL_DIRECTORIES
    ${MOD_DIR}/glsl
)
SET(MOD_INSTALL_FILES
    ${MOD_DIR}/ext/jama/license.txt
)


################################################################################
# Qt module resources 
################################################################################
SET(MOD_QT_SOURCES
    ${MOD_DIR}/qt/labelingwidgetqt.cpp
)  
    
SET(MOD_QT_HEADERS
    ${MOD_DIR}/qt/labelingwidgetqt.h
)

SET(MOD_QT_HEADERS_NONMOC
)
