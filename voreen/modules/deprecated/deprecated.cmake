################################################################################
# Core module resources 
################################################################################
SET(MOD_CORE_MODULECLASS DeprecatedModule)

# Core profile ready
SET(MOD_CORE_SOURCES
    ${MOD_DIR}/processors/volume/volumefiltering.cpp
    ${MOD_DIR}/processors/volume/volumemorphology.cpp
    
    ${MOD_DIR}/io/philipsusvolumereader.cpp
    ${MOD_DIR}/io/visiblehumanreader.cpp
    ${MOD_DIR}/io/vevovolumereader.cpp
    
    ${MOD_DIR}/octree/octreebrickpoolmanagerdisksinglethreaded.cpp
)

# Core profile ready
SET(MOD_CORE_HEADERS
    ${MOD_DIR}/processors/volume/volumefiltering.h
    ${MOD_DIR}/processors/volume/volumemorphology.h

    ${MOD_DIR}/io/philipsusvolumereader.h
    ${MOD_DIR}/io/visiblehumanreader.h
    ${MOD_DIR}/io/vevovolumereader.h

    ${MOD_DIR}/octree/octreebrickpoolmanagerdisksinglethreaded.h
)
