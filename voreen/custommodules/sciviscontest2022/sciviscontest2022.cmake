IF(NOT VRN_MODULE_ENSEMBLEANALYSIS)
    MESSAGE(FATAL_ERROR "SciVis Contest 2022 module requires Ensemble Analysis Module")
ENDIF()

IF(NOT VRN_MODULE_FLOWANALYSIS)
    MESSAGE(FATAL_ERROR "SciVis Contest 2022 module requires Flow Analysis Module")
ENDIF()

IF(NOT VRN_MODULE_PYTHON)
    MESSAGE(WARNING "SciVis Contest 2022 module requires Python Module")
ENDIF()

IF(NOT VRN_MODULE_WEBVIEW)
    MESSAGE(WARNING "SciVis Contest 2022 module requires WebView Module")
ENDIF()


SET(MOD_CORE_MODULECLASS SciVisContest2022Module)

SET(MOD_CORE_SOURCES
    ${MOD_DIR}/processors/binary_geometry/binarygeometry.cpp
    ${MOD_DIR}/processors/binary_geometry/binarygeometrysave.cpp
    ${MOD_DIR}/processors/binary_geometry/binarygeometrysequencesource.cpp
    ${MOD_DIR}/processors/binary_geometry/binarygeometrysource.cpp
    ${MOD_DIR}/processors/geometry/geometryselector.cpp
    ${MOD_DIR}/processors/geometry/geometrysequencesource.cpp
    ${MOD_DIR}/processors/geometry/boundingboxsampler.cpp
    ${MOD_DIR}/processors/geometry/curvecreator.cpp
    ${MOD_DIR}/datastructures/curve.cpp
    ${MOD_DIR}/processors/surface/pathsurfacecreator.cpp
    ${MOD_DIR}/processors/surface/pathsurfacerenderer.cpp
    ${MOD_DIR}/processors/surface/integrators.cpp
    ${MOD_DIR}/processors/surface/streamsurfacecreator.cpp
    ${MOD_DIR}/processors/surface/streamsurfaceviewer.cpp
    ${MOD_DIR}/processors/vorticity/vorticityfieldcreator.cpp
    ${MOD_DIR}/processors/vorticity/divergencefieldcreator.cpp
    ${MOD_DIR}/processors/vorticity/pyrogenicvorticitymapper.cpp
    ${MOD_DIR}/processors/vorticity/materialderivative.cpp
    ${MOD_DIR}/processors/corelinelengthvolumecreator.cpp
    ${MOD_DIR}/processors/similarityvolumecreator.cpp
    ${MOD_DIR}/processors/volumeinterpolation.cpp
    ${MOD_DIR}/processors/volumekernel.cpp
    ${MOD_DIR}/processors/seededstreamlinecreator.cpp
    ${MOD_DIR}/processors/seededpathlinecreator.cpp
    ${MOD_DIR}/processors/geometry/curvecreator.cpp
)

SET(MOD_CORE_HEADERS
    ${MOD_DIR}/processors/binary_geometry/binarygeometry.h
    ${MOD_DIR}/processors/binary_geometry/binarygeometrysave.h
    ${MOD_DIR}/processors/binary_geometry/binarygeometrysequencesource.h
    ${MOD_DIR}/processors/binary_geometry/binarygeometrysource.h
    ${MOD_DIR}/processors/geometry/geometryselector.h
    ${MOD_DIR}/processors/geometry/geometrysequencesource.h
    ${MOD_DIR}/processors/geometry/boundingboxsampler.h
    ${MOD_DIR}/processors/geometry/curvecreator.h
    ${MOD_DIR}/datastructures/curve.h
    ${MOD_DIR}/processors/surface/pathsurfacecreator.h
    ${MOD_DIR}/processors/surface/pathsurfacerenderer.h
    ${MOD_DIR}/processors/surface/integrators.h
    ${MOD_DIR}/processors/surface/streamsurfacecreator.h
    ${MOD_DIR}/processors/surface/streamsurfaceviewer.h
    ${MOD_DIR}/processors/vorticity/vorticityfieldcreator.h
    ${MOD_DIR}/processors/vorticity/divergencefieldcreator.h
    ${MOD_DIR}/processors/vorticity/pyrogenicvorticitymapper.h
    ${MOD_DIR}/processors/vorticity/materialderivative.h
    ${MOD_DIR}/processors/corelinelengthvolumecreator.h
    ${MOD_DIR}/processors/similarityvolumecreator.h
    ${MOD_DIR}/processors/volumeinterpolation.h
    ${MOD_DIR}/processors/volumekernel.h
    ${MOD_DIR}/processors/seededstreamlinecreator.h
    ${MOD_DIR}/processors/seededpathlinecreator.h
    ${MOD_DIR}/processors/geometry/curvecreator.h
    
    ${MOD_DIR}/properties/link/linkevaluatorvectorelement.h
)

# Deployment
SET(MOD_INSTALL_DIRECTORIES
    ${MOD_DIR}/workspaces
)
