
################################################################################
# Core module resources 
################################################################################
SET(MOD_CORE_MODULECLASS GraphLayoutModule)
SET(MOD_REQUIRE_OPENGL_COMPATIBILITY_PROFILE ON)

IF(NOT VRN_MODULE_PLOTTING)
    MESSAGE(FATAL_ERROR "GraphLayout Module requires Plotting Module")
ENDIF()

SET(MOD_CORE_SOURCES
    ${MOD_DIR}/src/forcedirectedlayout.cpp
    ${MOD_DIR}/src/forcedirectednodegraph.cpp
    ${MOD_DIR}/src/graphplot.cpp
    ${MOD_DIR}/src/nodegraph.cpp
    ${MOD_DIR}/src/nodegraphsource.cpp
)

SET(MOD_CORE_HEADERS
    ${MOD_DIR}/include/forcedirectedlayout.h
    ${MOD_DIR}/include/forcedirectednodegraph.h
    ${MOD_DIR}/include/graphplot.h
    ${MOD_DIR}/include/nodegraph.h
    ${MOD_DIR}/include/nodegraphsource.h
)
