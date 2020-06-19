IF(NOT VRN_MODULE_BASE)
    MESSAGE(FATAL_ERROR "Big Data Image Processing Module requires Base Module")
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
            # Note that we set the same library for both, relase and debug.
            # This way the user has to decide how the library should be compiled.
            LIST(APPEND MOD_RELEASE_DLLS ${LZ4_INCLUDE_DIR}/../lib/${elem}.dll)
            LIST(APPEND MOD_RELEASE_LIBRARIES ${LZ4_INCLUDE_DIR}/../lib/${elem}.lib)
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
    ${MOD_DIR}/processors/connectedcomponentanalysis.cpp
    ${MOD_DIR}/processors/largevolumeformatconversion.cpp
    ${MOD_DIR}/processors/volumefilterlist.cpp
    ${MOD_DIR}/processors/volumeresampletransformation.cpp

    # Volumefiltering
    ${MOD_DIR}/volumefiltering/binarizationfilter.cpp
    ${MOD_DIR}/volumefiltering/binarymedianfilter.cpp
    ${MOD_DIR}/volumefiltering/gaussianfilter.cpp
    ${MOD_DIR}/volumefiltering/gradientfilter.cpp
    ${MOD_DIR}/volumefiltering/medianfilter.cpp
    ${MOD_DIR}/volumefiltering/morphologyfilter.cpp
    ${MOD_DIR}/volumefiltering/parallelvolumefilter.cpp
    ${MOD_DIR}/volumefiltering/resamplefilter.cpp
    ${MOD_DIR}/volumefiltering/rescalefilter.cpp
    ${MOD_DIR}/volumefiltering/slicereader.cpp
    ${MOD_DIR}/volumefiltering/thresholdingfilter.cpp
    ${MOD_DIR}/volumefiltering/valuemapfilter.cpp
    ${MOD_DIR}/volumefiltering/volumefilter.cpp
    ${MOD_DIR}/volumefiltering/vorticityfilter.cpp

    # Filter Properties
    ${MOD_DIR}/volumefilterproperties/binarizationfilterproperties.cpp
    ${MOD_DIR}/volumefilterproperties/binarymedianfilterproperties.cpp
    ${MOD_DIR}/volumefilterproperties/filterproperties.cpp
    ${MOD_DIR}/volumefilterproperties/gaussianfilterproperties.cpp
    ${MOD_DIR}/volumefilterproperties/gradientfilterproperties.cpp
    ${MOD_DIR}/volumefilterproperties/medianfilterproperties.cpp
    ${MOD_DIR}/volumefilterproperties/morphologyfilterproperties.cpp
    ${MOD_DIR}/volumefilterproperties/resamplefilterproperties.cpp
    ${MOD_DIR}/volumefilterproperties/rescalefilterproperties.cpp
    ${MOD_DIR}/volumefilterproperties/templatefilterproperties.cpp
    ${MOD_DIR}/volumefilterproperties/thresholdingfilterproperties.cpp
    ${MOD_DIR}/volumefilterproperties/valuemapfilterproperties.cpp
    ${MOD_DIR}/volumefilterproperties/vorticityfilterproperties.cpp

    # cell nuclei cluster splitting
    ${MOD_DIR}/util/clustersplittingthread.cpp
    ${MOD_DIR}/util/connectedcomponentqueue.cpp
    ${MOD_DIR}/processors/nucleiclustersplitting.cpp
)
SET(MOD_CORE_HEADERS
    ${MOD_DIR}/algorithm/streamingcomponents.h
    ${MOD_DIR}/datastructures/lz4slicevolume.h
    ${MOD_DIR}/io/lz4slicevolumefilereader.h
    ${MOD_DIR}/io/volumedisklz4.h
    ${MOD_DIR}/processors/connectedcomponentanalysis.h
    ${MOD_DIR}/processors/largevolumeformatconversion.h
    ${MOD_DIR}/processors/volumefilterlist.h
    ${MOD_DIR}/processors/volumeresampletransformation.h

    # Volumefiltering
    ${MOD_DIR}/volumefiltering/binarizationfilter.h
    ${MOD_DIR}/volumefiltering/binarymedianfilter.h
    ${MOD_DIR}/volumefiltering/gaussianfilter.h
    ${MOD_DIR}/volumefiltering/gradientfilter.h
    ${MOD_DIR}/volumefiltering/medianfilter.h
    ${MOD_DIR}/volumefiltering/morphologyfilter.h
    ${MOD_DIR}/volumefiltering/parallelvolumefilter.h
    ${MOD_DIR}/volumefiltering/resamplefilter.h
    ${MOD_DIR}/volumefiltering/rescalefilter.h
    ${MOD_DIR}/volumefiltering/slicereader.h
    ${MOD_DIR}/volumefiltering/thresholdingfilter.h
    ${MOD_DIR}/volumefiltering/valuemapfilter.h
    ${MOD_DIR}/volumefiltering/volumefilter.h
    ${MOD_DIR}/volumefiltering/vorticityfilter.h

    # Filter Properties
    ${MOD_DIR}/volumefilterproperties/binarizationfilterproperties.h
    ${MOD_DIR}/volumefilterproperties/binarymedianfilterproperties.h
    ${MOD_DIR}/volumefilterproperties/filterproperties.h
    ${MOD_DIR}/volumefilterproperties/gaussianfilterproperties.h
    ${MOD_DIR}/volumefilterproperties/gradientfilterproperties.h
    ${MOD_DIR}/volumefilterproperties/medianfilterproperties.h
    ${MOD_DIR}/volumefilterproperties/morphologyfilterproperties.h
    ${MOD_DIR}/volumefilterproperties/rescalefilterproperties.h
    ${MOD_DIR}/volumefilterproperties/templatefilterproperties.h
    ${MOD_DIR}/volumefilterproperties/thresholdingfilterproperties.h
    ${MOD_DIR}/volumefilterproperties/valuemapfilterproperties.cpp
    ${MOD_DIR}/volumefilterproperties/vorticityfilterproperties.h

    # cell nuclei cluster splitting
    ${MOD_DIR}/operators/ternaryvolumeoperator.h
    ${MOD_DIR}/operators/volumeoperatordistancetransform.h
    ${MOD_DIR}/operators/volumeoperatorfastvolumecombine.h
    ${MOD_DIR}/operators/volumeoperatorgradientdescent.h
    ${MOD_DIR}/operators/volumeoperatorwatershed.h
    ${MOD_DIR}/util/clustersplittingthread.h
    ${MOD_DIR}/util/connectedcomponentqueue.h
    ${MOD_DIR}/util/csvwriter.h
    ${MOD_DIR}/processors/nucleiclustersplitting.h
)
