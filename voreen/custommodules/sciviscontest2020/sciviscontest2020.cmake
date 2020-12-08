################################################################################
# Core module resources
################################################################################

IF(NOT VRN_MODULE_ENSEMBLEANALYSIS)
    MESSAGE(FATAL_ERROR "SciVis Contest 2020 module requires Ensemble Analysis Module")
ENDIF()


SET(MOD_CORE_MODULECLASS SciVisContest2020Module)

SET(MOD_CORE_SOURCES
    ${MOD_DIR}/sciviscontest2020module.cpp

    # datastructures
    ${MOD_DIR}/datastructures/vortex.cpp

    # ports
    ${MOD_DIR}/ports/vortexport.cpp

    # processors
    ${MOD_DIR}/processors/approximateparallelvectors.cpp
    ${MOD_DIR}/processors/corelinedensityvolumecreator.cpp
    ${MOD_DIR}/processors/curlprocessor.cpp
    ${MOD_DIR}/processors/particlerenderer.cpp
    ${MOD_DIR}/processors/rotationaldirectionprocessor.cpp
    ${MOD_DIR}/processors/uncertainvectorfieldprocessor.cpp
    ${MOD_DIR}/processors/vortexcollectioncreator.cpp
    ${MOD_DIR}/processors/vortexcollectionsource.cpp
    ${MOD_DIR}/processors/vortexlistselector.cpp
    ${MOD_DIR}/processors/vortexmatchselector.cpp
    ${MOD_DIR}/processors/vortexprocessor.cpp
    ${MOD_DIR}/processors/vortexselector.cpp
    ${MOD_DIR}/processors/vortextracking.cpp
)
    
SET(MOD_CORE_HEADERS
    ${MOD_DIR}/sciviscontest2020module.h

    # datastructures
    ${MOD_DIR}/datastructures/vortex.h

    # ports
    ${MOD_DIR}/ports/vortexport.h

    # processors
    ${MOD_DIR}/processors/approximateparallelvectors.h
    ${MOD_DIR}/processors/corelinedensityvolumecreator.h
    ${MOD_DIR}/processors/curlprocessor.h
    ${MOD_DIR}/processors/flowmapprocessor.h
    ${MOD_DIR}/processors/ftvaprocessor.h
    ${MOD_DIR}/processors/particlerenderer.h
    ${MOD_DIR}/processors/rotationaldirectionprocessor.h
    ${MOD_DIR}/processors/uncertainvectorfieldprocessor.h
    ${MOD_DIR}/processors/vortexcollectioncreator.h
    ${MOD_DIR}/processors/vortexcollectionsource.h
    ${MOD_DIR}/processors/vortexlistselector.h
    ${MOD_DIR}/processors/vortexmatchselector.h
    ${MOD_DIR}/processors/vortexprocessor.h
    ${MOD_DIR}/processors/vortexselector.h
    ${MOD_DIR}/processors/vortextracking.h
)

# Deployment
SET(MOD_INSTALL_DIRECTORIES
    ${MOD_DIR}/glsl
)
