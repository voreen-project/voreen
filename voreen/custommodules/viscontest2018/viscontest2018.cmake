################################################################################
# Core module resources
################################################################################

# Dependencies
IF(NOT VRN_MODULE_PYTHON) # for converting
    MESSAGE(FATAL_ERROR "VisContest2018 Module requires Python Module")
ENDIF()
IF(NOT VRN_MODULE_HDF5) # for volume reading
    MESSAGE(FATAL_ERROR "VisContest2018 Module requires HDF5 Module")
ENDIF()
IF(NOT VRN_MODULE_PLOTTING) # for plotting
    MESSAGE(FATAL_ERROR "VisContest2018 Module requires Plotting Module")
ENDIF()

SET(MOD_CORE_MODULECLASS VisContest2018Module)
SET(MOD_OPENGL_CORE_PROFILE_COMPATIBLE ON)

# External libraries
IF(UNIX)
    FIND_PACKAGE(VTK REQUIRED)
    INCLUDE(${VTK_USE_FILE})
    SET(MOD_LIBRARIES ${VTK_LIBRARIES})
ELSE()
    
    SET(VRN_VTK_VERSION 8.0)

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

# deployment
SET(MOD_INSTALL_DIRECTORIES
    ${MOD_DIR}/glsl
    ${MOD_DIR}/workspaces
)

SET(MOD_CORE_SOURCES
    #Datastructures
    ${MOD_DIR}/datastructures/ensembledataset.cpp
    ${MOD_DIR}/datastructures/fieldplotdata.cpp
    ${MOD_DIR}/datastructures/similaritydata.cpp
    
    #IO
    ${MOD_DIR}/io/fieldplotsave.cpp
    ${MOD_DIR}/io/fieldplotsource.cpp
    ${MOD_DIR}/io/vtivolumereader.cpp
    
    #Processors
    ${MOD_DIR}/processors/ensembledatasource.cpp
    ${MOD_DIR}/processors/ensemblefilter.cpp
    ${MOD_DIR}/processors/ensemblesimilarityplot.cpp
    ${MOD_DIR}/processors/ensemblevolumeextractor.cpp
    ${MOD_DIR}/processors/fieldparallelplotcreator.cpp
    ${MOD_DIR}/processors/fieldparallelplothistogram.cpp
    ${MOD_DIR}/processors/fieldparallelplotviewer.cpp
    ${MOD_DIR}/processors/mdsplot.cpp
    ${MOD_DIR}/processors/volumeintensityfilter.cpp
    ${MOD_DIR}/processors/probabilityvolumecreator.cpp
    ${MOD_DIR}/processors/waveheightextractor.cpp

    #Properties
    ${MOD_DIR}/properties/stringlistproperty.cpp
    
    #Ports
    ${MOD_DIR}/ports/ensembledatasetport.cpp
    ${MOD_DIR}/ports/fieldplotdataport.cpp

    #Conditions
    ${MOD_DIR}/ports/conditions/portconditionensemble.cpp
    
    #Interaction
    
    #Utils
    ${MOD_DIR}/utils/colorpool.cpp
    ${MOD_DIR}/utils/ensemblehash.cpp
    ${MOD_DIR}/utils/utils.cpp
)

SET(MOD_CORE_HEADERS
    #Datastructures
    ${MOD_DIR}/datastructures/ensembledataset.h
    ${MOD_DIR}/datastructures/fieldplotdata.h
    ${MOD_DIR}/datastructures/similaritydata.h

    #IO
    ${MOD_DIR}/io/fieldplotsave.h
    ${MOD_DIR}/io/fieldplotsource.h
    ${MOD_DIR}/io/vtivolumereader.h

    #Processors
    ${MOD_DIR}/processors/ensembledatasource.h
    ${MOD_DIR}/processors/ensemblefilter.h
    ${MOD_DIR}/processors/ensemblesimilarityplot.h
    ${MOD_DIR}/processors/ensemblevolumeextractor.h
    ${MOD_DIR}/processors/fieldparallelplotcreator.h
    ${MOD_DIR}/processors/fieldparallelplothistogram.h
    ${MOD_DIR}/processors/fieldparallelplotviewer.h
    ${MOD_DIR}/processors/mdsplot.h
    ${MOD_DIR}/processors/volumeintensityfilter.h
    ${MOD_DIR}/processors/probabilityvolumecreator.h
    ${MOD_DIR}/processors/waveheightextractor.h

    #Properties
    ${MOD_DIR}/properties/stringlistproperty.h
    ${MOD_DIR}/properties/link/viscontest2018linkevaluatorid.h

    #Ports
    ${MOD_DIR}/ports/ensembledatasetport.h
    ${MOD_DIR}/ports/fieldplotdataport.h
    ${MOD_DIR}/ports/similaritydataport.h

    #Conditions
    ${MOD_DIR}/ports/conditions/portconditionensemble.h

    #Interaction

    #Utils
    ${MOD_DIR}/utils/colorpool.h
    ${MOD_DIR}/utils/ensemblehash.h
    ${MOD_DIR}/utils/utils.h
)

###############################################################################
# Qt module resources 
################################################################################
SET(MOD_QT_MODULECLASS VisContest2018ModuleQt)

SET(MOD_QT_SOURCES
    #Factories
    ${MOD_DIR}/qt/properties/viscontest2018propertywidgetfactory.cpp

    #Properties
    ${MOD_DIR}/qt/properties/stringlistpropertywidget.cpp
)  
    
SET(MOD_QT_HEADERS
    #Factories
    ${MOD_DIR}/qt/properties/viscontest2018propertywidgetfactory.h
    
    #Properties
    ${MOD_DIR}/qt/properties/stringlistpropertywidget.h
)

SET(MOD_QT_HEADERS_NONMOC
)
