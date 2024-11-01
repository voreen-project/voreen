
################################################################################
# Core module resources
################################################################################
SET(MOD_CORE_MODULECLASS RandomWalkerModule)

OPTION(VRN_RW_MAGMA_SOLVER "Build magma solver? (Requires cuda)" OFF)
IF(VRN_RW_MAGMA_SOLVER)
    add_definitions(-DVRN_RW_USE_MAGMA)

    find_package(CUDA)
    IF(BOOST_FOUND)
	LIST(APPEND MOD_INCLUDE_DIRECTORIES ${CUDA_INCLUDE_DIRS})
	LIST(APPEND MOD_LIBRARIES ${CUDA_LIBRARIES})
    ELSE()
        MESSAGE(ERROR "Cuda not found")
    ENDIF()

    find_package(PkgConfig REQUIRED)
    pkg_check_modules(MAGMA REQUIRED magma)
    IF(MAGMA_FOUND)
	LIST(APPEND MOD_INCLUDE_DIRECTORIES ${MAGMA_INCLUDE_DIRS})
	LIST(APPEND MOD_LIBRARIES ${MAGMA_LIBRARIES})
    ELSE()
        MESSAGE(ERROR "Magma not found")
    ENDIF()
ENDIF()

SET(MOD_CORE_SOURCES
    ${MOD_DIR}/processors/octreewalker.cpp
    ${MOD_DIR}/processors/randomwalker.cpp
    ${MOD_DIR}/processors/randomwalkeranalyzer.cpp
    ${MOD_DIR}/processors/rwmultilabelloopfinalizer.cpp
    ${MOD_DIR}/processors/rwmultilabelloopinitializer.cpp
    ${MOD_DIR}/processors/supervoxelwalker.cpp

    ${MOD_DIR}/solver/randomwalkersolver.cpp
    ${MOD_DIR}/solver/randomwalkerseeds.cpp
    ${MOD_DIR}/solver/randomwalkerweights.cpp

    ${MOD_DIR}/util/memprofiler.cpp
    ${MOD_DIR}/util/noisemodel.cpp
    ${MOD_DIR}/util/preprocessing.cpp
    ${MOD_DIR}/util/seeds.cpp
)

SET(MOD_CORE_HEADERS
    ${MOD_DIR}/processors/octreewalker.h
    ${MOD_DIR}/processors/randomwalker.h
    ${MOD_DIR}/processors/randomwalkeranalyzer.h
    ${MOD_DIR}/processors/rwmultilabelloopinitializer.h
    ${MOD_DIR}/processors/rwmultilabelloopfinalizer.h
    ${MOD_DIR}/processors/supervoxelwalker.h

    ${MOD_DIR}/solver/randomwalkersolver.h
    ${MOD_DIR}/solver/randomwalkerseeds.h
    ${MOD_DIR}/solver/randomwalkerweights.h

    ${MOD_DIR}/util/memprofiler.h
    ${MOD_DIR}/util/noisemodel.h
    ${MOD_DIR}/util/preprocessing.h
    ${MOD_DIR}/util/seeds.h
)


################################################################################
# Qt module resources
################################################################################
SET(MOD_QT_MODULECLASS RandomWalkerModuleQt)

SET(MOD_QT_SOURCES
    ${MOD_DIR}/qt/randomwalkeranalyzerwidget.cpp
    ${MOD_DIR}/qt/randomwalkerprocessorwidgetfactory.cpp
)

SET(MOD_QT_HEADERS
    ${MOD_DIR}/qt/randomwalkeranalyzerwidget.h
)

SET(MOD_QT_HEADERS_NONMOC
    ${MOD_DIR}/qt/randomwalkerprocessorwidgetfactory.h
)

# deployment
SET(MOD_INSTALL_DIRECTORIES
    ${MOD_DIR}/workspaces
)
