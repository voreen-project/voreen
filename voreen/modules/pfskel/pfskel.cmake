################################################################################
# Dependencies 
################################################################################
IF (NOT VRN_MODULE_ROI)
    MESSAGE(FATAL_ERROR "PFSkel requires the ROI module")
ENDIF()

################################################################################
# Core module resources 
################################################################################
SET(MOD_CORE_MODULECLASS PFSkelModule)
SET(MOD_REQUIRE_OPENGL_COMPATIBILITY_PROFILE ON)

SET(MOD_CORE_SOURCES 
    ${MOD_DIR}/processors/roiskeletonize.cpp
    ${MOD_DIR}/utils/pfskelwrapper.cpp
)

SET(MOD_CORE_HEADERS 
    ${MOD_DIR}/processors/roiskeletonize.h
    ${MOD_DIR}/utils/pfskelwrapper.h
)

