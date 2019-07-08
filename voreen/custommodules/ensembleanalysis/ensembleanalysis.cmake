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
IF(NOT VRN_MODULE_PLOTTING) # for plotting
    MESSAGE(FATAL_ERROR "EnsembleAnalysis Module requires Plotting Module")
ENDIF()

SET(MOD_CORE_MODULECLASS EnsembleAnalysisModule)
SET(MOD_OPENGL_CORE_PROFILE_COMPATIBLE ON)

# External libraries
OPTION(VRN_USE_VTK "Use VTK for conversion into HDF5" ON)
IF(${VRN_USE_VTK})
    MESSAGE(STATUS "Using VTK library")
    SET(MOD_DEFINITIONS "-DVRN_USE_VTK")
    IF(UNIX)
        FIND_PACKAGE(VTK REQUIRED)
        INCLUDE(${VTK_USE_FILE})
        SET(MOD_LIBRARIES ${VTK_LIBRARIES})
    ELSE()

        SET(VRN_VTK_VERSION 8.1)

        LIST(APPEND VTK_LIB_NAMES #add missing
            "CommonCore" "CommonDataModel" "CommonMisc" "CommonSystem" "CommonTransforms"
            "CommonExecutionModel" "CommonMath" "expat" "sys" "lz4" "zlib" "IOCore" "IOXML" "IOXMLParser"
        )

        FOREACH(elem ${VTK_LIB_NAMES})
            LIST(APPEND MOD_DEBUG_LIBRARIES   "${MOD_DIR}/ext/vtk/lib/debug/vtk${elem}-${VRN_VTK_VERSION}.lib")
            LIST(APPEND MOD_DEBUG_DLLS        "${MOD_DIR}/ext/vtk/lib/debug/vtk${elem}-${VRN_VTK_VERSION}.dll")
            LIST(APPEND MOD_RELEASE_LIBRARIES "${MOD_DIR}/ext/vtk/lib/release/vtk${elem}-${VRN_VTK_VERSION}.lib")
            LIST(APPEND MOD_RELEASE_DLLS      "${MOD_DIR}/ext/vtk/lib/release/vtk${elem}-${VRN_VTK_VERSION}.dll")
        ENDFOREACH()

        SET(MOD_INCLUDE_DIRECTORIES "${MOD_DIR}/ext/vtk/include")

    ENDIF()
ENDIF()

# deployment
SET(MOD_INSTALL_DIRECTORIES
    ${MOD_DIR}/glsl
    ${MOD_DIR}/scripts
    ${MOD_DIR}/workspaces
)

SET(MOD_CORE_SOURCES
    #Datastructures
    ${MOD_DIR}/datastructures/ensembledataset.cpp
    ${MOD_DIR}/datastructures/fieldplotdata.cpp
    ${MOD_DIR}/datastructures/similaritymatrix.cpp

    #IO
    ${MOD_DIR}/io/fieldplotsave.cpp
    ${MOD_DIR}/io/fieldplotsource.cpp
    ${MOD_DIR}/io/similaritymatrixsave.cpp
    ${MOD_DIR}/io/similaritymatrixsource.cpp

    #Processors
    ${MOD_DIR}/processors/ensembledatasource.cpp
    ${MOD_DIR}/processors/ensemblefilter.cpp
    ${MOD_DIR}/processors/ensemblevolumeextractor.cpp
    ${MOD_DIR}/processors/fieldparallelplotcreator.cpp
    ${MOD_DIR}/processors/fieldparallelplothistogram.cpp
    ${MOD_DIR}/processors/fieldparallelplotviewer.cpp
    ${MOD_DIR}/processors/physicalclippinglinker.cpp
    ${MOD_DIR}/processors/similaritydatavolume.cpp
    ${MOD_DIR}/processors/similaritymatrixcombine.cpp
    ${MOD_DIR}/processors/similaritymatrixcreator.cpp
    ${MOD_DIR}/processors/similarityplot.cpp
    ${MOD_DIR}/processors/volumelistmerger.cpp
    ${MOD_DIR}/processors/volumemerger.cpp

    #Properties
    ${MOD_DIR}/properties/stringlistproperty.cpp
    
    #Ports
    ${MOD_DIR}/ports/ensembledatasetport.cpp
    ${MOD_DIR}/ports/fieldplotdataport.cpp
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
    ${MOD_DIR}/datastructures/fieldplotdata.h
    ${MOD_DIR}/datastructures/similaritymatrix.h

    #IO
    ${MOD_DIR}/io/fieldplotsave.h
    ${MOD_DIR}/io/fieldplotsource.h
    ${MOD_DIR}/io/similaritymatrixsave.h
    ${MOD_DIR}/io/similaritymatrixsource.h

    #Processors
    ${MOD_DIR}/processors/ensembledatasource.h
    ${MOD_DIR}/processors/ensemblefilter.h
    ${MOD_DIR}/processors/ensemblevolumeextractor.h
    ${MOD_DIR}/processors/fieldparallelplotcreator.h
    ${MOD_DIR}/processors/fieldparallelplothistogram.h
    ${MOD_DIR}/processors/fieldparallelplotviewer.h
    ${MOD_DIR}/processors/physicalclippinglinker.h
    ${MOD_DIR}/processors/similaritydatavolume.h
    ${MOD_DIR}/processors/similaritymatrixcombine.h
    ${MOD_DIR}/processors/similaritymatrixcreator.h
    ${MOD_DIR}/processors/similarityplot.h
    ${MOD_DIR}/processors/volumelistmerger.h
    ${MOD_DIR}/processors/volumemerger.h

    #Properties
    ${MOD_DIR}/properties/stringlistproperty.h
    ${MOD_DIR}/properties/link/ensembleanalysislinkevaluatorid.h

    #Ports
    ${MOD_DIR}/ports/ensembledatasetport.h
    ${MOD_DIR}/ports/fieldplotdataport.h
    ${MOD_DIR}/ports/similaritymatrixport.h

    #Conditions
    ${MOD_DIR}/ports/conditions/portconditionensemble.h

    #Interaction

    #Utils
    ${MOD_DIR}/utils/ensemblehash.h
    ${MOD_DIR}/utils/utils.h
)

IF(${VRN_USE_VTK})
    LIST(APPEND MOD_CORE_SOURCES
        ${MOD_DIR}/io/vtivolumereader.cpp
        ${MOD_DIR}/io/vtmvolumereader.cpp
    )
    LIST(APPEND MOD_CORE_HEADERS
        ${MOD_DIR}/io/vtivolumereader.h
        ${MOD_DIR}/io/vtmvolumereader.h
    )
ENDIF()

###############################################################################
# Qt module resources 
################################################################################
SET(MOD_QT_MODULECLASS EnsembleAnalysisModuleQt)

SET(MOD_QT_SOURCES
    #Factories
    ${MOD_DIR}/qt/properties/ensembleanalysispropertywidgetfactory.cpp

    #Properties
    ${MOD_DIR}/qt/properties/stringlistpropertywidget.cpp
)  
    
SET(MOD_QT_HEADERS
    #Factories
    ${MOD_DIR}/qt/properties/ensembleanalysispropertywidgetfactory.h
    
    #Properties
    ${MOD_DIR}/qt/properties/stringlistpropertywidget.h
)

SET(MOD_QT_HEADERS_NONMOC
)
