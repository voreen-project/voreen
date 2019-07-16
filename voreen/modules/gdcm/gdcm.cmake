
################################################################################
# External dependency: GDCM library
################################################################################

IF(WIN32)
    
    #set these properties
	SET(VRN_GDCM_VERSION 2.8)
    
    # include path
    SET(GDCM_INCLUDE_DIR ${MOD_DIR}/ext/gdcm-${VRN_GDCM_VERSION}/include)
    IF(NOT EXISTS ${GDCM_INCLUDE_DIR}/gdcmReader.h)
        MESSAGE(FATAL_ERROR "GDCM 2.x (>1) headers not found (${GDCM_INCLUDE_DIR}/gdcmReader.h). "
            "Copy GDCM 2.x (>1)library to modules/gdcm/ext/gdcm-${VRN_GDCM_VERSION} (see http://voreen.uni-muenster.de)")
    ENDIF()
    SET(MOD_INCLUDE_DIRECTORIES ${GDCM_INCLUDE_DIR})
    
    # Define used libs:
    LIST(APPEND GDCM_LIB_NAMES 
        "gdcmcharls" "gdcmCommon" "gdcmDICT" "gdcmDSED" "gdcmexpat" "gdcmgetopt" "gdcmIOD"
        "gdcmjpeg8" "gdcmjpeg12" "gdcmjpeg16" "gdcmMEXD" "gdcmMSFF" "gdcmopenjp2" "gdcmzlib"
        "socketxx" #needed for debug builds
    )
    
    # Define used dlls (equal to lib in this case):
    LIST(APPEND GDCM_DLL_NAMES ${GDCM_LIB_NAMES})
    
    # set debug and release libraries
    IF(VRN_MSVC2012)            
		SET(GDCM_LIB_DIR "${MOD_DIR}/ext/gdcm-${VRN_GDCM_VERSION}/lib/msvc2012")
	ELSEIF(VRN_MSVC2013)            
		SET(GDCM_LIB_DIR "${MOD_DIR}/ext/gdcm-${VRN_GDCM_VERSION}/lib/msvc2013")
	ELSEIF(VRN_MSVC2015)            
		SET(GDCM_LIB_DIR "${MOD_DIR}/ext/gdcm-${VRN_GDCM_VERSION}/lib/msvc2015")
	ELSEIF(VRN_MSVC2017)            
		SET(GDCM_LIB_DIR "${MOD_DIR}/ext/gdcm-${VRN_GDCM_VERSION}/lib/msvc2017")
	ENDIF()

    IF(NOT EXISTS ${GDCM_LIB_DIR}/debug/gdcmCommon.lib)
        MESSAGE(FATAL_ERROR "GDCM 2.x (>1) library not found (${GDCM_LIB_DIR}/debug/gdcmCommon.lib). "
            "Copy GDCM 2.x (>1) library to modules/gdcm/ext/gdcm-${VRN_GDCM_VERSION} (see http://voreen.uni-muenster.de)")
    ENDIF()
    
    # set libraries
    FOREACH(elem ${GDCM_LIB_NAMES})
        LIST(APPEND MOD_DEBUG_LIBRARIES    "${GDCM_LIB_DIR}/debug/${elem}.lib")
        LIST(APPEND MOD_RELEASE_LIBRARIES  "${GDCM_LIB_DIR}/release/${elem}.lib")
    ENDFOREACH()

    # set dlls
    FOREACH(elem ${GDCM_DLL_NAMES})
        LIST(APPEND MOD_DEBUG_DLLS         "${GDCM_LIB_DIR}/debug/${elem}.dll")
        LIST(APPEND MOD_RELEASE_DLLS       "${GDCM_LIB_DIR}/release/${elem}.dll")
    ENDFOREACH()
    
    # deployment
    SET(MOD_INSTALL_FILES
        ${MOD_DIR}/ext/gdcm-${VRN_GDCM_VERSION}/Copyright.txt
    )
    
ELSEIF(UNIX)
    find_package(GDCM)
    IF(GDCM_FOUND)
        SET(GDCM_USE_VTK 0)
        INCLUDE(${GDCM_USE_FILE})

        IF((GDCM_MAJOR_VERSION EQUAL 2) AND (GDCM_MINOR_VERSION GREATER 1) OR (GDCM_MAJOR_VERSION EQUAL 3))
            MESSAGE( STATUS "GDCM 2.2 or newer detected." )
            #LIST(APPEND MOD_LIBRARIES -lgdcmMEXD -lgdcmopenjpeg -lgdcmuuid -lgdcmzlib)
            LIST(APPEND MOD_LIBRARIES -lgdcmCommon -lgdcmDICT -lgdcmDSED -lgdcmIOD -lgdcmMEXD -lgdcmMSFF -lgdcmjpeg12 -lgdcmjpeg16 -lgdcmjpeg8)
        ELSE()
            MESSAGE(FATAL_ERROR "Only GDCM major version 2 with minor version 2 or greater is supported.")
        ENDIF()
    ELSE()
        MESSAGE(FATAL_ERROR "GDCM library not found.")
    ENDIF()

ENDIF()


################################################################################
# Core module resources 
################################################################################
SET(MOD_CORE_MODULECLASS GdcmModule)

SET(MOD_CORE_SOURCES 
    ${MOD_DIR}/customdicomdict.cpp
    ${MOD_DIR}/dicominfo.cpp
    ${MOD_DIR}/dicomdictentry.cpp
    ${MOD_DIR}/dicomdict.cpp
    ${MOD_DIR}/io/dicomdirparser.cpp
    ${MOD_DIR}/io/dicomnetworkconnector.cpp
    ${MOD_DIR}/io/gdcmvolumereader.cpp
    ${MOD_DIR}/io/volumediskdicom.cpp
)

SET(MOD_CORE_HEADERS 
    ${MOD_DIR}/customdicomdict.h
    ${MOD_DIR}/dicominfo.h
    ${MOD_DIR}/dicomdictentry.h
    ${MOD_DIR}/dicomdict.h
    ${MOD_DIR}/io/dicomdirparser.h
    ${MOD_DIR}/io/dicomnetworkconnector.h
    ${MOD_DIR}/io/gdcmvolumereader.h
    ${MOD_DIR}/io/volumediskdicom.h
)

# deployment
LIST(APPEND MOD_INSTALL_DIRECTORIES ${MOD_DIR}/dicts)


################################################################################
# Qt module resources 
################################################################################
#SET(MOD_CQT_MODULECLASS GdcmModule)

SET(MOD_QT_SOURCES 
    ${MOD_DIR}/qt/dicomconnectiondialog.cpp
    ${MOD_DIR}/qt/dicomhierarchymodel.cpp
)

SET(MOD_QT_HEADERS 
    ${MOD_DIR}/qt/dicomconnectiondialog.h
    ${MOD_DIR}/qt/dicomhierarchymodel.h
)
   
SET(MOD_QT_HEADERS_NONMOC
)

