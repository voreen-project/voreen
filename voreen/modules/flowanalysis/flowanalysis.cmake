################################################################################
# Core module resources
################################################################################

IF(NOT VRN_MODULE_ENSEMBLEANALYSIS)
    MESSAGE(WARNING "EnsembleAnalysis Module not enabled, some features will not be available")
ENDIF()


SET(MOD_CORE_MODULECLASS FlowAnalysisModule)

SET(MOD_CORE_SOURCES
    ${MOD_DIR}/flowanalysismodule.cpp

    # datastructures
    ${MOD_DIR}/datastructures/streamline.cpp
    ${MOD_DIR}/datastructures/streamlinebundle.cpp
    ${MOD_DIR}/datastructures/streamlinelist.cpp
    ${MOD_DIR}/datastructures/streamlinelistbase.cpp
    ${MOD_DIR}/datastructures/streamlinelistdecorator.cpp
    ${MOD_DIR}/datastructures/streamlinelistobserver.cpp
    ${MOD_DIR}/datastructures/vortex.cpp

    # ports
    ${MOD_DIR}/ports/parallelvectorsolutionpointsport.cpp
    ${MOD_DIR}/ports/streamlinelistport.cpp
    ${MOD_DIR}/ports/vortexport.cpp

    # processors
    # corelines
    ${MOD_DIR}/processors/corelines/corelinecreator.cpp
    ${MOD_DIR}/processors/corelines/corelinedensityvolumecreator.cpp
    ${MOD_DIR}/processors/corelines/parallelvectors.cpp

    # geometry
    ${MOD_DIR}/processors/geometry/streamlinetoboundingbox.cpp
    ${MOD_DIR}/processors/geometry/streamlinetogeometry.cpp

    # render
    ${MOD_DIR}/processors/render/flowdirectionoverlay.cpp
    ${MOD_DIR}/processors/render/streamlinerenderer3d.cpp

    # streamline
    ${MOD_DIR}/processors/streamline/pathlinecreator.cpp
    ${MOD_DIR}/processors/streamline/streamlinebundledetector.cpp
    ${MOD_DIR}/processors/streamline/streamlinecombine.cpp
    ${MOD_DIR}/processors/streamline/streamlinecreator.cpp
    ${MOD_DIR}/processors/streamline/streamlinefilter.cpp
    ${MOD_DIR}/processors/streamline/streamlinerotation.cpp
    ${MOD_DIR}/processors/streamline/streamlinesave.cpp
    ${MOD_DIR}/processors/streamline/streamlineselector.cpp
    ${MOD_DIR}/processors/streamline/streamlinesource.cpp

    # volume
    ${MOD_DIR}/processors/volume/accelerationprocessor.cpp
    ${MOD_DIR}/processors/volume/curlprocessor.cpp
    #${MOD_DIR}/processors/volume/flowmapprocessor.cpp
    #${MOD_DIR}/processors/volume/ftvaprocessor.cpp
    ${MOD_DIR}/processors/volume/helicitydensity.cpp
    ${MOD_DIR}/processors/volume/vortexprocessor.cpp

    # vortex
    ${MOD_DIR}/processors/vortex/rotationaldirectionprocessor.cpp
    ${MOD_DIR}/processors/vortex/vortexmatchselector.cpp
    ${MOD_DIR}/processors/vortex/vortexselector.cpp
    ${MOD_DIR}/processors/vortex/vortextracking.cpp

    # utils
    ${MOD_DIR}/utils/flowutils.cpp
)
    
SET(MOD_CORE_HEADERS
    ${MOD_DIR}/flowanalysismodule.h

    # datastructures
    ${MOD_DIR}/datastructures/streamline.h
    ${MOD_DIR}/datastructures/streamlinebundle.h
    ${MOD_DIR}/datastructures/streamlinelist.h
    ${MOD_DIR}/datastructures/streamlinelistbase.h
    ${MOD_DIR}/datastructures/streamlinelistdecorator.h
    ${MOD_DIR}/datastructures/streamlinelistobserver.h
    ${MOD_DIR}/datastructures/vortex.h

    # ports
    ${MOD_DIR}/ports/parallelvectorsolutionpointsport.h
    ${MOD_DIR}/ports/streamlinelistport.h
    ${MOD_DIR}/ports/vortexport.h

    # processors
    # corelines
    ${MOD_DIR}/processors/corelines/corelinecreator.h
    ${MOD_DIR}/processors/corelines/corelinedensityvolumecreator.h
    ${MOD_DIR}/processors/corelines/parallelvectors.h

    # geometry
    ${MOD_DIR}/processors/geometry/streamlinetoboundingbox.h
    ${MOD_DIR}/processors/geometry/streamlinetogeometry.h

    # render
    ${MOD_DIR}/processors/render/flowdirectionoverlay.h
    ${MOD_DIR}/processors/render/streamlinerenderer3d.h

    # streamline
    ${MOD_DIR}/processors/streamline/pathlinecreator.h
    ${MOD_DIR}/processors/streamline/streamlinebundledetector.h
    ${MOD_DIR}/processors/streamline/streamlinecombine.h
    ${MOD_DIR}/processors/streamline/streamlinecreator.h
    ${MOD_DIR}/processors/streamline/streamlinefilter.h
    ${MOD_DIR}/processors/streamline/streamlinerotation.h
    ${MOD_DIR}/processors/streamline/streamlinesave.h
    ${MOD_DIR}/processors/streamline/streamlineselector.h
    ${MOD_DIR}/processors/streamline/streamlinesource.h

    # volume
    ${MOD_DIR}/processors/volume/accelerationprocessor.h
    ${MOD_DIR}/processors/volume/curlprocessor.h
    ${MOD_DIR}/processors/volume/flowmapprocessor.h
    ${MOD_DIR}/processors/volume/ftvaprocessor.h
    ${MOD_DIR}/processors/volume/helicitydensity.h
    ${MOD_DIR}/processors/volume/vortexprocessor.h

    # vortex
    ${MOD_DIR}/processors/vortex/rotationaldirectionprocessor.h
    ${MOD_DIR}/processors/vortex/vortexmatchselector.h
    ${MOD_DIR}/processors/vortex/vortexselector.h
    ${MOD_DIR}/processors/vortex/vortextracking.h

    # utils
    ${MOD_DIR}/utils/flowutils.h
)

IF(VRN_MODULE_ENSEMBLEANALYSIS)
    SET(MOD_CORE_SOURCES ${MOD_CORE_SOURCES}
        ${MOD_DIR}/processors/ensemble/approximateparallelvectors.cpp
        ${MOD_DIR}/processors/ensemble/particlerenderer.cpp
        ${MOD_DIR}/processors/ensemble/uncertainvectorfieldprocessor.cpp
        ${MOD_DIR}/processors/ensemble/vortexcollectioncreator.cpp
        ${MOD_DIR}/processors/ensemble/vortexcollectionsource.cpp
        ${MOD_DIR}/processors/ensemble/vortexlistselector.cpp
    )
    SET(MOD_CORE_HEADERS ${MOD_CORE_HEADERS}
        ${MOD_DIR}/processors/ensemble/approximateparallelvectors.h
        ${MOD_DIR}/processors/ensemble/particlerenderer.h
        ${MOD_DIR}/processors/ensemble/uncertainvectorfieldprocessor.h
        ${MOD_DIR}/processors/ensemble/vortexcollectioncreator.h
        ${MOD_DIR}/processors/ensemble/vortexcollectionsource.h
        ${MOD_DIR}/processors/ensemble/vortexlistselector.h
    )
ENDIF()

# Deployment
SET(MOD_INSTALL_DIRECTORIES
    ${MOD_DIR}/glsl
)
