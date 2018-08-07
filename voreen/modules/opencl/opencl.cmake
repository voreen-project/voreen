
################################################################################
# External dependency: OpenCL SDK
################################################################################

IF(WIN32)
    FIND_PATH(NVIDIA_OPENCL_BASEDIR
                "$ENV{NVSDKCUDA_ROOT}"
                "$ENV{NVSDKCOMPUTE_ROOT}")
              
    FIND_PATH(AMD_OPENCL_BASEDIR Cl/cl.h
              PATHS
              $ENV{AMDAPPSDKROOT}/include
              $ENV{ATISTREAMSDKROOT}/include)

    IF(AMD_OPENCL_BASEDIR)
        # Setup for AMD/ATI Stream-SDK    
        FIND_PATH( OPENCL_INCLUDE_DIR
                   NAMES CL/cl.hpp OpenCL/cl.hpp CL/cl.h OpenCL/cl.h
                   PATHS
                   $ENV{AMDAPPSDKROOT}/include
                   $ENV{ATISTREAMSDKROOT}/include 
                   ${AMD_OPENCL_BASEDIR}/include)

        MESSAGE(STATUS "  - AMD/ATI Stream OpenCL include path: " ${OPENCL_INCLUDE_DIR})

		FIND_LIBRARY(OPENCL_LIBRARY
					 NAMES OpenCL
					 PATHS 
					 $ENV{AMDAPPSDKROOT}/lib/x86_64
					 $ENV{ATISTREAMSDKROOT}/lib/x86_64
					 ${AMD_OPENCL_BASEDIR}/lib/x86_64)
        MESSAGE(STATUS "  - AMD/ATI Stream OpenCL library path: " ${OPENCL_LIBRARY})
 
        #for openCL 1.2
        ADD_DEFINITIONS( -DCL_USE_DEPRECATED_OPENCL_1_0_APIS )
        ADD_DEFINITIONS( -DCL_USE_DEPRECATED_OPENCL_1_1_APIS )
        ADD_DEFINITIONS( -DCL_USE_DEPRECATED_OPENCL_1_2_APIS )
    ELSEIF(NVIDIA_OPENCL_BASEDIR)
        # Setup for NVIDIA CUDA SDK
        FIND_PATH(OPENCL_INCLUDE_DIR
                  NAMES CL/cl.hpp OpenCL/cl.hpp CL/cl.h OpenCL/cl.h
                  PATHS $ENV{NVSDKCOMPUTE_ROOT}/OpenCL/common/inc/
                        ${NVIDIA_OPENCL_BASEDIR}/../include)
        MESSAGE(STATUS "  - NVidia CUDA OpenCL include path: " ${OPENCL_INCLUDE_DIR})
      
		FIND_LIBRARY(OPENCL_LIBRARY
					 NAMES OpenCL
					 PATHS $ENV{NVSDKCOMPUTE_ROOT}/OpenCL/common/lib/x64
					 ${NVIDIA_OPENCL_BASEDIR}/../lib/x64)

        MESSAGE(STATUS "  - NVidia CUDA OpenCL library path: " ${OPENCL_LIBRARY})
        
        #for openCL 1.2
        ADD_DEFINITIONS( -DCL_USE_DEPRECATED_OPENCL_1_2_APIS )
    ENDIF()

ELSEIF(UNIX)

    FIND_PATH(
        OPENCL_INCLUDE_DIR
        NAMES CL/cl.h
        PATHS ${OPENCL_INCLUDE_DIR}
              ${OPENCL_DIR}/include
              ${OPENCL_DIR}/OpenCL/common/inc
              $ENV{OPENCL_INCLUDE_DIR}
              $ENV{OPENCL_DIR}/include
              $ENV{OPENCL_DIR}/OpenCL/common/inc
              /usr/local/cuda/include
              /usr/local/include
              /usr/include
    )

    FIND_LIBRARY(
        OPENCL_LIBRARY
        NAMES OpenCL
        PATHS ${OPENCL_LIBRARY_DIR}
              ${OPENCL_DIR}/lib
              ${OPENCL_DIR}/lib/x86
              $ENV{OPENCL_LIBRARY_DIR}
              $ENV{OPENCL_DIR}/lib
              $ENV{OPENCL_DIR}/lib/x86
              /usr/local/lib64
              /usr/local/lib
              /usr/lib64
              /usr/lib
              /usr/lib/x86_64-linux-gnu
    )

LIST(APPEND MOD_DEFINITIONS "-DCL_USE_DEPRECATED_OPENCL_1_1_APIS")
LIST(APPEND MOD_DEFINITIONS "-DCL_USE_DEPRECATED_OPENCL_1_2_APIS")

ENDIF()

IF(OPENCL_INCLUDE_DIR AND OPENCL_LIBRARY)
    #MESSAGE(STATUS "  - Found OpenCL library")
    SET(MOD_INCLUDE_DIRECTORIES "${OPENCL_INCLUDE_DIR}")
    SET(MOD_LIBRARIES "${OPENCL_LIBRARY}")
    MARK_AS_ADVANCED(OPENCL_INCLUDE_DIR OPENCL_LIBRARY)
ELSE()
    MESSAGE(FATAL_ERROR "OpenCL library not found: Please set either NVIDIA_OPENCL_BASEDIR or AMD_OPENCL_BASEDIR")
    MARK_AS_ADVANCED(CLEAR OPENCL_INCLUDE_DIR OPENCL_LIBRARY)
ENDIF()


################################################################################
# Core module resources 
################################################################################
SET(MOD_CORE_MODULECLASS OpenCLModule)

SET(MOD_CORE_SOURCES 
    ${MOD_DIR}/utils/clwrapper.cpp
    ${MOD_DIR}/processors/openclprocessor.cpp
    ${MOD_DIR}/processors/dynamicclprocessor.cpp
    ${MOD_DIR}/processors/grayscale_cl.cpp
    ${MOD_DIR}/processors/raycaster_cl.cpp
    ${MOD_DIR}/processors/raytracingentryexitpoints.cpp
    ${MOD_DIR}/processors/singleoctreeraycastercl.cpp
    ${MOD_DIR}/processors/volumegradient_cl.cpp
    ${MOD_DIR}/utils/voreenblascl.cpp
    ${MOD_DIR}/properties/openclproperty.cpp
)

SET(MOD_CORE_HEADERS
    ${MOD_DIR}/utils/clwrapper.h
    ${MOD_DIR}/processors/openclprocessor.h
    ${MOD_DIR}/processors/dynamicclprocessor.h
    ${MOD_DIR}/processors/grayscale_cl.h
    ${MOD_DIR}/processors/raycaster_cl.h
    ${MOD_DIR}/processors/raytracingentryexitpoints.h
    ${MOD_DIR}/processors/singleoctreeraycastercl.h
    ${MOD_DIR}/processors/volumegradient_cl.h
    ${MOD_DIR}/utils/voreenblascl.h
    ${MOD_DIR}/properties/openclproperty.h
)

# deployment
SET(MOD_INSTALL_DIRECTORIES
    ${MOD_DIR}/cl
    ${MOD_DIR}/glsl
)

################################################################################
# Qt module resources 
################################################################################
SET(MOD_QT_MODULECLASS OpenCLModuleQt)

SET(MOD_QT_SOURCES 
    ${MOD_DIR}/qt/openclpropertywidgetfactory.cpp
    ${MOD_DIR}/qt/openclhighlighter.cpp
    ${MOD_DIR}/qt/openclpropertywidget.cpp
    ${MOD_DIR}/qt/openclplugin.cpp
    ${MOD_DIR}/qt/openclprocessorwidgetfactory.cpp
    ${MOD_DIR}/qt/dynamicopenclwidget.cpp
)

SET(MOD_QT_HEADERS
    ${MOD_DIR}/qt/openclpropertywidget.h
    ${MOD_DIR}/qt/openclplugin.h
    ${MOD_DIR}/qt/dynamicopenclwidget.h
)

SET(MOD_QT_HEADERS_NONMOC
    ${MOD_DIR}/qt/openclpropertywidgetfactory.h
    ${MOD_DIR}/qt/openclprocessorwidgetfactory.h
    ${MOD_DIR}/qt/openclhighlighter.h
)
