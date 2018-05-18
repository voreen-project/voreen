
# Enable OpenMP
if(WIN32)
    LIST(APPEND MOD_DEFINITIONS /openmp)
ELSEIF(UNIX)
    FIND_PACKAGE(OpenMP)
    IF(OpenMP_FOUND)
        MESSAGE(STATUS "- Found OpenMP")
        SET(CMAKE_CXX_FLAGS "${OpenMP_CXX_FLAGS} ${CMAKE_CXX_FLAGS}")
        SET(MOD_LIBRARIES ${OpenMP_CXX_LIB_NAMES})
    ELSE()
        MESSAGE(FATAL_ERROR "Could not find OpenMP")
    ENDIF()
ENDIF()

################################################################################
# Core module resources 
################################################################################
SET(MOD_CORE_MODULECLASS OpenMPModule)

SET(MOD_CORE_SOURCES
    ${MOD_DIR}/src/voreenblasmp.cpp
)

SET(MOD_CORE_HEADERS
    ${MOD_DIR}/include/voreenblasmp.h
)
