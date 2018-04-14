
# Enable OpenMP
if(WIN32)
    LIST(APPEND MOD_DEFINITIONS /openmp)
ELSEIF(UNIX)
    set(CMAKE_CXX_FLAGS "-fopenmp ${CMAKE_CXX_FLAGS}")
    #LIST(APPEND CMAKE_CXX_FLAGS -fopenmp)
    SET(MOD_LIBRARIES -lgomp)
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
