################################################################################
# External dependency: TIFF library
################################################################################

IF(WIN32)
    SET(MOD_INCLUDE_DIRECTORIES "${MOD_DIR}/ext/libtiff/include")
    
    # deployment    
    SET(MOD_INSTALL_FILES
        ${MOD_DIR}/ext/libtiff/COPYRIGHT
    )
    
    # Define used libs:
    LIST(APPEND TIFF_LIB_NAMES "libtiff" "libtiff_i")
    
    # Define used dlls:
    LIST(APPEND TIFF_DLL_NAMES "libtiff")
    
    # set debug and release libraries
    IF(VRN_MSVC2017 OR VRN_MSVC2019 OR VRN_MSVC2022)            
        SET(TIFF_LIB_DIR "${MOD_DIR}/ext/libtiff/lib/msvc2017")
    ELSE()
        MESSAGE(FATAL_ERROR "Unsupported MSVC toolchain")
    ENDIF()
    
    # set libraries
    FOREACH(elem ${TIFF_LIB_NAMES})
        LIST(APPEND MOD_LIBRARIES "${TIFF_LIB_DIR}/${elem}.lib")
    ENDFOREACH()

    # set dlls
    FOREACH(elem ${TIFF_DLL_NAMES})
        LIST(APPEND MOD_DLLS "${TIFF_LIB_DIR}/${elem}.dll")
    ENDFOREACH()
    
ELSEIF(UNIX)
    FIND_PACKAGE(TIFF REQUIRED)
    IF(TIFF_FOUND)
        MESSAGE(STATUS "  - Found TIFF library")
        SET(MOD_INCLUDE_DIRECTORIES ${TIFF_INCLUDE_DIR})
        SET(MOD_LIBRARIES ${TIFF_LIBRARIES})
    ELSE()
        MESSAGE(FATAL_ERROR "Tiff library not found!")
    ENDIF()
ENDIF()


################################################################################
# Core module resources 
################################################################################
SET(MOD_CORE_MODULECLASS TiffModule)

SET(MOD_CORE_SOURCES 
    ${MOD_DIR}/io/tiffvolumereader.cpp
    ${MOD_DIR}/io/ometiffvolumereader.cpp
    ${MOD_DIR}/io/volumediskometiff.cpp
)

SET(MOD_CORE_HEADERS 
    ${MOD_DIR}/io/tiffvolumereader.h
    ${MOD_DIR}/io/ometiffvolumereader.h
    ${MOD_DIR}/io/volumediskometiff.h
)
   
