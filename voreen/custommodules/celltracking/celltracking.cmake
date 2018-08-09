# module class must reside in sample/samplemodule.h + sample/samplemodule.cpp
SET(MOD_CORE_MODULECLASS CellTrackingModule)

IF (NOT VRN_MODULE_TIFF)
    MESSAGE(FATAL_ERROR "Celltracking module requires the Tiff module")
ENDIF()
#IF (NOT VRN_MODULE_HDF5)
#    MESSAGE(FATAL_ERROR "Celltracking module requires the HDF5 module")
#ENDIF()
IF (NOT VRN_MODULE_PLOTTING)
    MESSAGE(FATAL_ERROR "Celltracking module requires the HDF5 module")
ENDIF()

# module's core source files, path relative to module dir
SET(MOD_CORE_SOURCES
    ${MOD_DIR}/processors/cellmigrationscore.cpp
    ${MOD_DIR}/processors/celltrackconverter.cpp
    ${MOD_DIR}/processors/timestepduration.cpp
    ${MOD_DIR}/processors/xmltiffvolumesource.cpp
)

# module's core header files, path relative to module dir
SET(MOD_CORE_HEADERS
    ${MOD_DIR}/processors/cellmigrationscore.h
    ${MOD_DIR}/processors/celltrackconverter.h
    ${MOD_DIR}/processors/timestepduration.h
    ${MOD_DIR}/processors/xmltiffvolumesource.h
)
