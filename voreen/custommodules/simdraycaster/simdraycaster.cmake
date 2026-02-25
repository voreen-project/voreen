
################################################################################
# SIMDraycaster module resources
################################################################################
SET(MOD_CORE_MODULECLASS SIMDraycasterModule)

SET(VRN_STAGING_SIMDRAYCASTER_AVAILABLE FALSE)
IF(WIN32)
    IF(CMAKE_SYSTEM_PROCESSOR MATCHES "^(AMD64|amd64|X64|x64|x86_64|i[3-6]86|x86)$")
        SET(VRN_STAGING_SIMDRAYCASTER_AVAILABLE TRUE)
    ENDIF()
ELSE()
    INCLUDE(CheckCXXCompilerFlag)
    IF(CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|amd64|i[3-6]86|x86)$")
        CHECK_CXX_COMPILER_FLAG("-msse3" VRN_STAGING_COMPILER_SUPPORTS_SSE3)
        IF(VRN_STAGING_COMPILER_SUPPORTS_SSE3)
            SET(VRN_STAGING_SIMDRAYCASTER_AVAILABLE TRUE)
        ENDIF()
    ENDIF()
ENDIF()

IF(NOT VRN_STAGING_SIMDRAYCASTER_AVAILABLE)
    MESSAGE(FATAL_ERROR "SIMD instructions not available on the current architecture")
ENDIF()

SET(MOD_CORE_SOURCES
    ${MOD_DIR}/processors/simdraycaster.cpp
	
    ${MOD_DIR}/utils/jobqueue.cpp
    ${MOD_DIR}/utils/memory.cpp
    ${MOD_DIR}/utils/performancemetric.cpp
)

if (UNIX)
    SET_SOURCE_FILES_PROPERTIES(${MOD_DIR}/processors/raycaster_sse41.cpp PROPERTIES COMPILE_FLAGS "${CMAKE_C_FLAGS} -msse4.1")
    SET_SOURCE_FILES_PROPERTIES(${MOD_DIR}/processors/raycaster_sse3.cpp PROPERTIES COMPILE_FLAGS "${CMAKE_C_FLAGS} -msse3")
endif (UNIX) 

SET(MOD_CORE_HEADERS
    ${MOD_DIR}/processors/raycast_generic.h
    ${MOD_DIR}/processors/simdraycaster.h
	
    ${MOD_DIR}/utils/brickedvolume.h
    ${MOD_DIR}/utils/brickedvolumebase.h
    ${MOD_DIR}/utils/jobqueue.h
    ${MOD_DIR}/utils/memory.h
    ${MOD_DIR}/utils/performancemetric.h
)

IF(WIN32)
    SET(MOD_CORE_SOURCES
        ${MOD_CORE_SOURCES}
        ${MOD_DIR}/processors/raycaster_sse41.cpp
        ${MOD_DIR}/processors/raycaster_sse3.cpp
    )
    LIST(APPEND VRN_MODULE_DEFINITIONS "-DSIMD_SSE3")
    LIST(APPEND VRN_MODULE_DEFINITIONS "-DSIMD_SSE41")
ELSE(WIN32)
    OPTION (VRN_USE_SSE41 
            "Use the sse 4.1 Instruction set extension" ON) 
    IF(VRN_USE_SSE41)
        SET(MOD_CORE_SOURCES
            ${MOD_CORE_SOURCES}
            ${MOD_DIR}/processors/raycaster_sse3.cpp
            ${MOD_DIR}/processors/raycaster_sse41.cpp
        )
        LIST(APPEND VRN_MODULE_DEFINITIONS "-DSIMD_SSE3")
        LIST(APPEND VRN_MODULE_DEFINITIONS "-DSIMD_SSE41 -msse4.1")
    ELSE(VRN_USE_SSE41)
        SET(MOD_CORE_SOURCES
            ${MOD_CORE_SOURCES}
            ${MOD_DIR}/processors/raycaster_sse3.cpp
        )
        LIST(APPEND VRN_MODULE_DEFINITIONS "-DSIMD_SSE3 -msse3")
    ENDIF(VRN_USE_SSE41)
ENDIF(WIN32)

# deployment
SET(MOD_INSTALL_DIRECTORIES
)
