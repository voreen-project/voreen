################################################################################
# Core module resources
################################################################################
SET(MOD_CORE_MODULECLASS SurfaceModule)

# deployment
SET(MOD_INSTALL_DIRECTORIES
    ${MOD_DIR}/glsl
    ${MOD_DIR}/workspaces
)

SET(MOD_CORE_SOURCES
    ${MOD_DIR}/processors/isosurfaceextractor.cpp
)

SET(MOD_CORE_HEADERS
    ${MOD_DIR}/processors/isosurfaceextractor.h
)
