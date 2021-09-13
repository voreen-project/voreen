IF(NOT VRN_MODULE_POI) # for converting
    MESSAGE(FATAL_ERROR "ProSeMediVis2020 Module requires POI Module")
ENDIF()

SET(MOD_CORE_MODULECLASS ProseMedvis2020Module)
 
# module's core source files, path relative to module dir
SET(MOD_CORE_SOURCES
    ## Processors
    # MRI
    ${MOD_DIR}/processors/concretevesselgraphcreator.cpp
    ${MOD_DIR}/processors/labelvolumecreator.cpp
    ${MOD_DIR}/processors/templatevesselgraphloader.cpp
    ${MOD_DIR}/processors/concretevesselgraphsource.cpp
    ${MOD_DIR}/processors/abstractvisualizationcreator.cpp
    ${MOD_DIR}/processors/cvglabelswapper.cpp
    # PET
    ${MOD_DIR}/processors/timeseriesextraction.cpp
    ${MOD_DIR}/processors/timeseriesfilter.cpp
    ${MOD_DIR}/processors/meanplot.cpp
    ${MOD_DIR}/processors/ensemblevolumepicking.cpp
    ${MOD_DIR}/processors/clickabletextureoverlay.cpp


    ## Datastructures
    # MRI
    ${MOD_DIR}/datastructures/concretevesselgraph.cpp
    ${MOD_DIR}/datastructures/templatevesselgraph.cpp
    # PET
    ${MOD_DIR}/datastructures/ensemblesamplepointmapping.cpp

    ## Ports
    # MRI
    ${MOD_DIR}/ports/concretevesselgraphport.cpp
    ${MOD_DIR}/ports/templatevesselgraphport.cpp
    # PET
    ${MOD_DIR}/ports/ensemblesamplepointmappingport.cpp

    # Utils
    ${MOD_DIR}/utils/samplepointconfigloader.cpp
)
 
# module's core header files, path relative to module dir
SET(MOD_CORE_HEADERS
    ## Processors
    # MRI
    ${MOD_DIR}/processors/concretevesselgraphcreator.h
    ${MOD_DIR}/processors/labelvolumecreator.h
    ${MOD_DIR}/processors/templatevesselgraphloader.h
    ${MOD_DIR}/processors/concretevesselgraphsource.h
    ${MOD_DIR}/processors/abstractvisualizationcreator.h
    ${MOD_DIR}/processors/cvglabelswapper.h
    # PET
    ${MOD_DIR}/processors/timeseriesextraction.h
    ${MOD_DIR}/processors/meanplot.h
    ${MOD_DIR}/processors/timeseriesfilter.h
    ${MOD_DIR}/processors/ensemblevolumepicking.h
    ${MOD_DIR}/processors/clickabletextureoverlay.h

    ## Datastructures
    # MRI
    ${MOD_DIR}/datastructures/concretevesselgraph.h
    ${MOD_DIR}/datastructures/templatevesselgraph.h
    # PET
    ${MOD_DIR}/datastructures/ensemblesamplepointmapping.h

    ## Ports
    # MRI
    ${MOD_DIR}/ports/concretevesselgraphport.h
    ${MOD_DIR}/ports/templatevesselgraphport.h
    # PET
    ${MOD_DIR}/ports/ensemblesamplepointmappingport.h

    # Utils
    ${MOD_DIR}/utils/samplepointconfigloader.h
)
