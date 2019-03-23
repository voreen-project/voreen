################################################################################
# Core module resources 
################################################################################

IF(NOT VRN_MODULE_HDF5)
    MESSAGE(FATAL_ERROR "Cell Nucleus Detection Module requires HDF5 Module")
ENDIF()

IF(NOT VRN_MODULE_OPENMP)
    MESSAGE(FATAL_ERROR "Cell Nucleus Detection Module requires OpenMP Module")
ENDIF()


SET(MOD_CORE_MODULECLASS NucleusDetectionModule)

# External dependencies: CBLAS and Caffe library for the training data and classification processors
UNSET(VRN_NUCLEUSDETECTION_CAFFE_FOUND)
UNSET(VRN_NUCLEUSDETECTION_CBLAS_FOUND)

IF(UNIX)
    # CBLAS library is required for feature extraction
    MESSAGE(STATUS "Trying to find CBLAS libraries...")

    find_path(CBLAS_INCLUDE_DIR cblas.h
        /usr/include/atlas
        /usr/local/include/atlas
        /usr/include
        /usr/local/include
    )

    set(CBLAS_NAMES ${CBLAS_NAMES} cblas)

    find_library(CBLAS_LIBRARY
        NAMES ${CBLAS_NAMES}
        PATHS
        /usr/lib64/atlas
        /usr/lib/atlas
        /usr/local/lib64/atlas
        /usr/local/lib/atlas
        /usr/lib64
        /usr/lib
        /usr/local/lib64
        /usr/local/lib
    )

    if (CBLAS_LIBRARY AND CBLAS_INCLUDE_DIR)
        SET(VRN_NUCLEUSDETECTION_CBLAS_FOUND 1)
        set(CBLAS_LIBRARIES ${CBLAS_LIBRARY})

        SET(MOD_INCLUDE_DIRECTORIES ${CBLAS_INCLUDE_DIR})
        SET(MOD_LIBRARIES ${CBLAS_LIBRARIES})
        MESSAGE(STATUS "Found CBLAS - building nucleus detection module with feature extraction")
        ADD_DEFINITIONS(-DVRN_NUCLEUSDETECTION_CBLAS_FOUND)

        # if we found CBLAS, we look for CAFFE
        MESSAGE(STATUS "Trying to find Caffe libraries...")

        FIND_PATH(CAFFE_INCLUDE_DIR NAMES caffe/caffe.hpp caffe/common.hpp caffe/net.hpp caffe/proto/caffe.pb.h caffe/util/io.hpp caffe/vision_layers.hpp HINTS /usr/local/include)

        FIND_LIBRARY(CAFFE_LIBRARIES NAMES caffe HINTS /usr/local/lib)

        IF(CAFFE_LIBRARIES AND CAFFE_INCLUDE_DIR)
            SET(VRN_NUCLEUSDETECTION_CAFFE_FOUND 1)
            SET(MOD_INCLUDE_DIRECTORIES ${MOD_INCLUDE_DIRECTORIES} ${CAFFE_INCLUDE_DIR})
            SET(MOD_LIBRARIES ${CAFFE_LIBRARIES})
            MESSAGE(STATUS "Found Caffe - building nucleus detection module with classification")
            ADD_DEFINITIONS(-DVRN_NUCLEUSDETECTION_CAFFE_FOUND)

            # try to find CUDA for include directories (without knowing if CAFFE was built with GPU support)
            FIND_PACKAGE(CUDA)
            IF (CUDA_FOUND)
                SET(MOD_INCLUDE_DIRECTORIES ${MOD_INCLUDE_DIRECTORIES} ${CUDA_INCLUDE_DIRS})
            ENDIF()

        ELSE()
            MESSAGE("Caffe Library not found - Cell Nucleus Detection Module will be built without classification processor")
        ENDIF()

    ELSE()
        MESSAGE("CBLAS Library not found - Cell Nucleus Detection Module will be built without training data extraction and classification processor")
    ENDIF()

ELSE()
    MESSAGE("CBLAS Library not found - Cell Nucleus Detection Module will be built without training data extraction and classification processor")
ENDIF()


SET(MOD_CORE_SOURCES

    # patch export and import for k-means
    ${MOD_DIR}/processors/patchextractor.cpp
    ${MOD_DIR}/processors/patchlistreader.cpp

    # extraction of tiles as training data for caffe
    ${MOD_DIR}/processors/tileextractor.cpp

    # helper processor which allows to locate a cropped out region in a list of potential original volumes
    ${MOD_DIR}/processors/roidetector.cpp
)

SET(MOD_CORE_HEADERS

    # patch export and import for k-means
    ${MOD_DIR}/processors/patchextractor.h
    ${MOD_DIR}/processors/patchlistreader.h

    # extraction of tiles as training data for caffe
    ${MOD_DIR}/processors/tileextractor.h

    # helper processor which allows to locate a cropped out region in a list of potential original volumes
    ${MOD_DIR}/processors/roidetector.h
)  

# CBLAS library is needed for feature extraction
IF (VRN_NUCLEUSDETECTION_CBLAS_FOUND)
    LIST(APPEND MOD_CORE_SOURCES
    # export of training data for caffe
	${MOD_DIR}/processors/patchtrainingdataextractor.cpp
    ${MOD_DIR}/util/cosinetransform.cpp
    ${MOD_DIR}/util/featureextractors.cpp
	${MOD_DIR}/util/trainingdataextractionthread.cpp

    # abstract base class for processors that extract local features
    ${MOD_DIR}/processors/patchfeatureextractor.cpp
    )
    LIST(APPEND MOD_CORE_HEADERS
    ${MOD_DIR}/processors/patchtrainingdataextractor.h
    ${MOD_DIR}/util/cosinetransform.h
    ${MOD_DIR}/util/featureextractors.h
	${MOD_DIR}/util/trainingdataextractionthread.h

    # abstract base class for processors that extract local features
    ${MOD_DIR}/processors/patchfeatureextractor.h
    )
ENDIF()

# caffe library is needed for classification
IF (VRN_NUCLEUSDETECTION_CAFFE_FOUND)
    LIST(APPEND MOD_CORE_SOURCES
        ${MOD_DIR}/util/featureextractionthread.cpp
        ${MOD_DIR}/processors/patchcaffeclassifier.cpp
        ${MOD_DIR}/processors/tilecaffeclassifier.cpp
    )
    LIST(APPEND MOD_CORE_HEADERS
        ${MOD_DIR}/util/featureextractionthread.h
        ${MOD_DIR}/processors/patchcaffeclassifier.h
        ${MOD_DIR}/processors/tilecaffeclassifier.h
    )
ENDIF()


