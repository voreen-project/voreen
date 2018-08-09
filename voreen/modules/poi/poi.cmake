SET(MOD_CORE_MODULECLASS POIModule)

# module's core source files, path relative to module dir
SET(MOD_CORE_SOURCES
    ${MOD_DIR}/processors/poicsvexport.cpp
    ${MOD_DIR}/processors/poicsvimport.cpp
    ${MOD_DIR}/processors/poisave.cpp
    ${MOD_DIR}/processors/poistorage.cpp
    ${MOD_DIR}/processors/poitextinfo.cpp
    ${MOD_DIR}/processors/poirenderer3d.cpp
    ${MOD_DIR}/processors/poirenderer2d.cpp
    ${MOD_DIR}/processors/poiselectionmanipulation.cpp
    ${MOD_DIR}/processors/poisource.cpp
    ${MOD_DIR}/processors/poipointsegmentgeometryexporter.cpp
    ${MOD_DIR}/processors/poipointsegmentgeometryimporter.cpp

    ${MOD_DIR}/datastructures/poilist.cpp
)

# module's core header files, path relative to module dir
SET(MOD_CORE_HEADERS
    ${MOD_DIR}/processors/poicsvexport.h
    ${MOD_DIR}/processors/poicsvimport.h
    ${MOD_DIR}/processors/poisave.h
    ${MOD_DIR}/processors/poiselectionmanipulation.h
    ${MOD_DIR}/processors/poistorage.h
    ${MOD_DIR}/processors/poisource.h
    ${MOD_DIR}/processors/poirenderer3d.h
    ${MOD_DIR}/processors/poirenderer2d.h
    ${MOD_DIR}/processors/poitextinfo.h
    ${MOD_DIR}/processors/poipointsegmentgeometryexporter.h
    ${MOD_DIR}/processors/poipointsegmentgeometryimporter.h

    ${MOD_DIR}/datastructures/poilist.h
)

# deployment
SET(MOD_INSTALL_DIRECTORIES
    ${MOD_DIR}/documentation
    ${MOD_DIR}/glsl
    ${MOD_DIR}/workspacess
)