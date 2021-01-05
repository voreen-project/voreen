################################################################################
# Core module resources
################################################################################

SET(MOD_CORE_MODULECLASS VTKModule)


IF(UNIX)
    FIND_PACKAGE(VTK REQUIRED)
    INCLUDE(${VTK_USE_FILE})
    SET(MOD_LIBRARIES ${VTK_LIBRARIES})
ELSE()

    SET(VRN_VTK_VERSION 8.1)

    LIST(APPEND VTK_LIB_NAMES #add missing
        "CommonCore" "CommonDataModel" "CommonMisc" "CommonSystem" "CommonTransforms"
        "CommonExecutionModel" "CommonMath" "expat" "hdf5" "hdf5_hl" "jpeg" "metaio" "NetCDF" "netcdfcpp"
        "png" "sys" "tiff" "lz4" "zlib" "DICOMParser" "IOCore" "IOImage" "IONetCDF" "IOXML" "IOXMLParser"
    )

    FOREACH(elem ${VTK_LIB_NAMES})
        LIST(APPEND MOD_DEBUG_LIBRARIES   "${MOD_DIR}/ext/vtk/lib/debug/vtk${elem}-${VRN_VTK_VERSION}.lib")
        LIST(APPEND MOD_DEBUG_DLLS        "${MOD_DIR}/ext/vtk/lib/debug/vtk${elem}-${VRN_VTK_VERSION}.dll")
        LIST(APPEND MOD_RELEASE_LIBRARIES "${MOD_DIR}/ext/vtk/lib/release/vtk${elem}-${VRN_VTK_VERSION}.lib")
        LIST(APPEND MOD_RELEASE_DLLS      "${MOD_DIR}/ext/vtk/lib/release/vtk${elem}-${VRN_VTK_VERSION}.dll")
    ENDFOREACH()

    SET(MOD_INCLUDE_DIRECTORIES "${MOD_DIR}/ext/vtk/include")

    SET(MOD_INSTALL_FILES
        ${MOD_DIR}/ext/vtk/Copyright.txt
    )

ENDIF()

SET(MOD_CORE_SOURCES
    ${MOD_DIR}/io/netcdfvolumereader.cpp
    ${MOD_DIR}/io/niftivolumewriter.cpp
    ${MOD_DIR}/io/vtivolumereader.cpp
    ${MOD_DIR}/io/vtivolumewriter.cpp
    ${MOD_DIR}/io/vtmvolumereader.cpp
)

SET(MOD_CORE_HEADERS
    ${MOD_DIR}/io/netcdfvolumereader.h
    ${MOD_DIR}/io/niftivolumewriter.h
    ${MOD_DIR}/io/vtivolumereader.h
    ${MOD_DIR}/io/vtivolumewriter.h
    ${MOD_DIR}/io/vtmvolumereader.h
)


