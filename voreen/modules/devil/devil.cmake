################################################################################
# External dependency: DevIL library
################################################################################

SET(MOD_DEFINITIONS -DTGT_HAS_DEVIL)

IF(WIN32)
    SET(MOD_INCLUDE_DIRECTORIES "${MOD_DIR}/ext/il/include")

	SET(MOD_LIBRARIES
		"${MOD_DIR}/ext/il/lib/DevIL.lib"
		"${MOD_DIR}/ext/il/lib/ILU.lib"
		#"${MOD_DIR}/ext/il/lib/ILUT.lib" #< currently not needed
	)

	SET(MOD_DEBUG_DLLS
		"${MOD_DIR}/ext/il/lib/DevIL.dll"
		"${MOD_DIR}/ext/il/lib/ILU.dll"
		#"${MOD_DIR}/ext/il/lib/ILUT.dll" #< currently not needed
	)
	SET(MOD_RELEASE_DLLS ${MOD_DEBUG_DLLS})

    LIST(APPEND MOD_DEBUG_DLLS "${MOD_DIR}/ext/jpeg/jpeg62.dll")
    LIST(APPEND MOD_RELEASE_DLLS "${MOD_DIR}/ext/jpeg/jpeg62.dll")

    # deployment
    SET(MOD_INSTALL_FILES
        ${MOD_DIR}/ext/il/lgpl.txt
        ${MOD_DIR}/ext/il/libpng-LICENSE.txt
        ${MOD_DIR}/ext/jpeg/license.txt
    )

ELSEIF(UNIX)
    FIND_PACKAGE(DevIL REQUIRED)
    IF(IL_INCLUDE_DIR AND IL_LIBRARIES AND ILU_LIBRARIES)
        MESSAGE(STATUS "  - Found DevIL library")
        SET(MOD_INCLUDE_DIRECTORIES ${IL_INCLUDE_DIR})
        SET(MOD_LIBRARIES
            ${IL_LIBRARIES}
            ${ILU_LIBRARIES}
            #${ILUT_LIBRARIES} #< currently not needed
        )
    ELSE()
        MESSAGE(FATAL_ERROR "DevIL library not found!")
    ENDIF()
ENDIF()


################################################################################
# Core module resources
################################################################################
SET(MOD_CORE_MODULECLASS DevILModule)

#SET(MOD_CORE_SOURCES )

#SET(MOD_CORE_HEADERS )

