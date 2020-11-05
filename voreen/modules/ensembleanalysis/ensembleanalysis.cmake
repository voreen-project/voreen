################################################################################
# Core module resources
################################################################################

# Dependencies
IF(NOT VRN_MODULE_PYTHON) # for converting
    MESSAGE(WARNING "EnsembleAnalysis Module requires Python Module for converter scripts")
ENDIF()
IF(NOT VRN_MODULE_HDF5) # for volume reading
    MESSAGE(WARNING "EnsembleAnalysis Module requires HDF5 Module for efficient memory storage")
ENDIF()
IF(NOT VRN_MODULE_VTK) # for data io
    MESSAGE(WARNING "EnsembleAnalysis Module requires VTK Module to load certain data formats")
ENDIF()
IF(NOT VRN_MODULE_PLOTTING) # for plotting
    MESSAGE(FATAL_ERROR "EnsembleAnalysis Module requires Plotting Module")
ENDIF()

SET(MOD_CORE_MODULECLASS EnsembleAnalysisModule)


# deployment
SET(MOD_INSTALL_DIRECTORIES
    ${MOD_DIR}/workspaces
)

SET(MOD_CORE_SOURCES
    #Datastructures
    ${MOD_DIR}/datastructures/ensembledataset.cpp
    ${MOD_DIR}/datastructures/parallelcoordinatesaxes.cpp
    ${MOD_DIR}/datastructures/similaritymatrix.cpp

    #IO
    ${MOD_DIR}/io/similaritymatrixsave.cpp
    ${MOD_DIR}/io/similaritymatrixsource.cpp

    #Processors
    ${MOD_DIR}/processors/connectedcomponentselector.cpp
    ${MOD_DIR}/processors/ensembledatasource.cpp
    ${MOD_DIR}/processors/ensemblefilter.cpp
    ${MOD_DIR}/processors/ensemblevolumeextractor.cpp
    ${MOD_DIR}/processors/localsimilarityanalysis.cpp
    ${MOD_DIR}/processors/referencevolumecreator.cpp
    ${MOD_DIR}/processors/parallelcoordinatesaxescreator.cpp
    ${MOD_DIR}/processors/parallelcoordinatessource.cpp
    ${MOD_DIR}/processors/parallelcoordinatesviewer.cpp
    ${MOD_DIR}/processors/parallelcoordinatesvoxelselection.cpp
    ${MOD_DIR}/processors/similaritymatrixcombine.cpp
    ${MOD_DIR}/processors/similaritymatrixcreator.cpp
    ${MOD_DIR}/processors/similarityplot.cpp
    ${MOD_DIR}/processors/volumelistmerger.cpp
    ${MOD_DIR}/processors/volumemerger.cpp
    ${MOD_DIR}/processors/metadataadder.cpp

    #Properties
    ${MOD_DIR}/properties/parallelcoordinatessectionsproperty.cpp
    ${MOD_DIR}/properties/parallelcoordinatesselectionproperty.cpp

    #Ports
    ${MOD_DIR}/ports/ensembledatasetport.cpp
    ${MOD_DIR}/ports/parallelcoordinatesaxesport.cpp
    ${MOD_DIR}/ports/similaritymatrixport.cpp

    #Conditions
    ${MOD_DIR}/ports/conditions/portconditionensemble.cpp
    
    #Interaction
    
    #Utils
    ${MOD_DIR}/utils/ensemblehash.cpp
    ${MOD_DIR}/utils/utils.cpp
)

SET(MOD_CORE_HEADERS
    #Datastructures
    ${MOD_DIR}/datastructures/ensembledataset.h
    ${MOD_DIR}/datastructures/parallelcoordinatesaxes.h
    ${MOD_DIR}/datastructures/similaritymatrix.h

    #IO
    ${MOD_DIR}/io/similaritymatrixsave.h
    ${MOD_DIR}/io/similaritymatrixsource.h

    #Processors
    ${MOD_DIR}/processors/connectedcomponentselector.h
    ${MOD_DIR}/processors/ensembledatasource.h
    ${MOD_DIR}/processors/ensemblefilter.h
    ${MOD_DIR}/processors/ensemblevolumeextractor.h
    ${MOD_DIR}/processors/localsimilarityanalysis.h
    ${MOD_DIR}/processors/referencevolumecreator.h
    ${MOD_DIR}/processors/parallelcoordinatesaxescreator.h
    ${MOD_DIR}/processors/parallelcoordinatessource.h
    ${MOD_DIR}/processors/parallelcoordinatesviewer.h
    ${MOD_DIR}/processors/parallelcoordinatesvoxelselection.h
    ${MOD_DIR}/processors/similaritymatrixcombine.h
    ${MOD_DIR}/processors/similaritymatrixcreator.h
    ${MOD_DIR}/processors/similarityplot.h
    ${MOD_DIR}/processors/volumelistmerger.h
    ${MOD_DIR}/processors/volumemerger.h
    ${MOD_DIR}/processors/metadataadder.h

    #Properties
    ${MOD_DIR}/properties/parallelcoordinatessectionsproperty.h
    ${MOD_DIR}/properties/parallelcoordinatesselectionproperty.h

    #Ports
    ${MOD_DIR}/ports/ensembledatasetport.h
    ${MOD_DIR}/ports/parallelcoordinatesaxesport.h
    ${MOD_DIR}/ports/similaritymatrixport.h

    #Conditions
    ${MOD_DIR}/ports/conditions/portconditionensemble.h

    #Interaction

    #Utils
    ${MOD_DIR}/utils/ensemblehash.h
    ${MOD_DIR}/utils/utils.h
)

IF(VRN_MODULE_HDF5)
    SET(MOD_CORE_SOURCES ${MOD_CORE_SOURCES}
        ${MOD_DIR}/processors/volumemerger.cpp
    )
    SET(MOD_CORE_HEADERS ${MOD_CORE_HEADERS}
        ${MOD_DIR}/processors/volumemerger.h
    )
ENDIF()
