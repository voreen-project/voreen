################################################################################
# Dependencies 
################################################################################
IF (NOT VRN_MODULE_PLOTTING)
    MESSAGE(FATAL_ERROR "Deprecated module requires the Plotting module")
ENDIF()

################################################################################
# Core module resources 
################################################################################
SET(MOD_CORE_MODULECLASS DeprecatedModule)
SET(MOD_REQUIRE_OPENGL_COMPATIBILITY_PROFILE ON)

SET(MOD_CORE_SOURCES
    ${MOD_DIR}/ports/textureport.cpp
    ${MOD_DIR}/processors/buttonoverlayprocessor.cpp
    ${MOD_DIR}/processors/cubeproxygeometry.cpp 
    ${MOD_DIR}/processors/canny.cpp 
    ${MOD_DIR}/processors/geometryclippingwidget.cpp
    ${MOD_DIR}/processors/multiview.cpp
    ${MOD_DIR}/processors/rawtexturesave.cpp 
    ${MOD_DIR}/processors/rawtexturesource.cpp
    ${MOD_DIR}/processors/targettotexture.cpp
    ${MOD_DIR}/processors/textseriessource.cpp
    ${MOD_DIR}/processors/texturetotarget.cpp
    ${MOD_DIR}/processors/volumenormalization.cpp
    ${MOD_DIR}/processors/volumeseriessource.cpp
    
    ${MOD_DIR}/io/philipsusvolumereader.cpp
    ${MOD_DIR}/io/visiblehumanreader.cpp
    ${MOD_DIR}/io/vevovolumereader.cpp
    
    ${MOD_DIR}/octree/octreebrickpoolmanagerdisksinglethreaded.cpp
)

SET(MOD_CORE_HEADERS
    ${MOD_DIR}/ports/textureport.h
    ${MOD_DIR}/processors/buttonoverlayprocessor.h
    ${MOD_DIR}/processors/canny.h 
    ${MOD_DIR}/processors/cubeproxygeometry.h
    ${MOD_DIR}/processors/geometryclippingwidget.h
    ${MOD_DIR}/processors/multiview.h
    ${MOD_DIR}/processors/textseriessource.h
    ${MOD_DIR}/processors/targettotexture.h
    ${MOD_DIR}/processors/texturetotarget.h
    ${MOD_DIR}/processors/volumeseriessource.h
    ${MOD_DIR}/processors/rawtexturesource.h
    ${MOD_DIR}/processors/rawtexturesave.h
    ${MOD_DIR}/processors/volumenormalization.h 

    ${MOD_DIR}/io/philipsusvolumereader.h
    ${MOD_DIR}/io/visiblehumanreader.h
    ${MOD_DIR}/io/vevovolumereader.h
    ${MOD_DIR}/operators/volumeoperatornormalize.h

    ${MOD_DIR}/octree/octreebrickpoolmanagerdisksinglethreaded.h

)

# deployment
SET(MOD_INSTALL_DIRECTORIES
    ${MOD_DIR}/glsl
)
