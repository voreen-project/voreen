if(VRN_MSVC2012 OR VRN_MSVC2013)
    MESSAGE(FATAL_ERROR "Big Data Image Processing Module does NOT compile with MSVC 2012 and MSVC 2013 - use MSVC 2015 or higher")
ENDIF()

IF(NOT VRN_MODULE_HDF5)
    MESSAGE(FATAL_ERROR "Big Data Image Processing Module requires HDF5 Module")
ENDIF()

################################################################################
# External dependency: lz4 library
################################################################################

IF(VRN_MSVC)

    SET(LZ4_INCLUDE_DIR ${MOD_DIR}/ext/lz4/include)
    SET(LZ4_LIBRARIES liblz4)
    SET(LZ4_FOUND TRUE)
    
    SET(MOD_INSTALL_FILES
        ${MOD_DIR}/ext/lz4/LICENSE
    )

ELSE()

    # Add path for additional cmake find scripts
    LIST(APPEND CMAKE_MODULE_PATH "${MOD_DIR}/ext/")

    MESSAGE(STATUS "Trying to find lz4 libraries")
    FIND_PACKAGE(LZ4)

ENDIF()

IF(LZ4_FOUND)
    MESSAGE(STATUS "  - Found lz4 library")

    MESSAGE(STATUS "Include Directories: " ${LZ4_INCLUDE_DIR})
    MESSAGE(STATUS "Libraries: " ${LZ4_LIBRARIES})
    LIST(APPEND MOD_INCLUDE_DIRECTORIES ${LZ4_INCLUDE_DIR})
    
    IF(VRN_MSVC)
        FOREACH(elem ${LZ4_LIBRARIES})
            # Don't copy both debug and release binaries, since they are named the same. Prefer debug.
            #LIST(APPEND MOD_RELEASE_DLLS ${LZ4_INCLUDE_DIR}/../lib/${elem}.dll)
            #LIST(APPEND MOD_RELEASE_LIBRARIES ${LZ4_INCLUDE_DIR}/../lib/${elem}.lib)
            LIST(APPEND MOD_DEBUG_DLLS ${LZ4_INCLUDE_DIR}/../lib/${elem}.dll)
            LIST(APPEND MOD_DEBUG_LIBRARIES ${LZ4_INCLUDE_DIR}/../lib/${elem}.lib)
        ENDFOREACH()
    ELSE()
        LIST(APPEND MOD_LIBRARIES ${LZ4_LIBRARIES})
    ENDIF()
ELSE()
    MESSAGE(FATAL_ERROR "Could not find lz4 Library.")
ENDIF()

################################################################################
# Core module resources
################################################################################
SET(MOD_CORE_MODULECLASS BigDataImageProcessingModule)

SET(MOD_CORE_SOURCES
    ${MOD_DIR}/datastructures/lz4slicevolume.cpp
    ${MOD_DIR}/io/lz4slicevolumefilereader.cpp
    ${MOD_DIR}/io/volumedisklz4.cpp
    ${MOD_DIR}/processors/binarymedian.cpp
    ${MOD_DIR}/processors/connectedcomponentanalysis.cpp
    ${MOD_DIR}/processors/largevolumeformatconversion.cpp
    ${MOD_DIR}/processors/segmentationquantification.cpp
    ${MOD_DIR}/processors/volumebricksave.cpp
    ${MOD_DIR}/processors/volumebricksource.cpp
    ${MOD_DIR}/processors/volumefilterlist.cpp
    ${MOD_DIR}/processors/volumeresampletransformation.cpp

    # Volumefiltering
    ${MOD_DIR}/volumefiltering/slicereader.cpp
    ${MOD_DIR}/volumefiltering/parallelvolumefilter.cpp
    ${MOD_DIR}/volumefiltering/gaussianfilter.cpp
    ${MOD_DIR}/volumefiltering/medianfilter.cpp
    ${MOD_DIR}/volumefiltering/binarymedianfilter.cpp

    # Properties
    ${MOD_DIR}/properties/interactivelistproperty.cpp
)
IF(VRN_MODULE_PLOTTING)
    LIST(APPEND MOD_CORE_SOURCES
        ${MOD_DIR}/processors/segmentationslicedensity.cpp
    )
ENDIF()
SET(MOD_CORE_HEADERS
    ${MOD_DIR}/datastructures/lz4slicevolume.h
    ${MOD_DIR}/io/lz4slicevolumefilereader.h
    ${MOD_DIR}/io/volumedisklz4.h
    ${MOD_DIR}/processors/binarymedian.cpp
    ${MOD_DIR}/processors/connectedcomponentanalysis.h
    ${MOD_DIR}/processors/largevolumeformatconversion.h
    ${MOD_DIR}/processors/segmentationquantification.h
    ${MOD_DIR}/processors/volumebricksave.h
    ${MOD_DIR}/processors/volumebricksource.h
    ${MOD_DIR}/processors/volumefilterlist.h
    ${MOD_DIR}/processors/volumeresampletransformation.h

    # Volumefiltering
    ${MOD_DIR}/volumefiltering/slicereader.h
    ${MOD_DIR}/volumefiltering/volumefilter.h
    ${MOD_DIR}/volumefiltering/parallelvolumefilter.h
    ${MOD_DIR}/volumefiltering/gaussianfilter.h
    ${MOD_DIR}/volumefiltering/medianfilter.h
    ${MOD_DIR}/volumefiltering/binarymedianfilter.h

    # Properties
    ${MOD_DIR}/properties/interactivelistproperty.h
)
IF(VRN_MODULE_PLOTTING)
    LIST(APPEND MOD_CORE_HEADERS
        ${MOD_DIR}/processors/segmentationslicedensity.h
    )
ENDIF()

###############################################################################
# Qt module resources
################################################################################
SET(MOD_QT_MODULECLASS BigDataImageProcessingModuleQt)

SET(MOD_QT_SOURCES
        #Factories
        ${MOD_DIR}/qt/properties/bigdataimageprocessingpropertywidgetfactory.cpp

        #Properties
        ${MOD_DIR}/qt/properties/interactivelistpropertywidget.cpp
        )

SET(MOD_QT_HEADERS
        #Factories
        ${MOD_DIR}/qt/properties/bigdataimageprocessingpropertywidgetfactory.h

        #Properties
        ${MOD_DIR}/qt/properties/interactivelistpropertywidget.h
        )

SET(MOD_QT_HEADERS_NONMOC
        )