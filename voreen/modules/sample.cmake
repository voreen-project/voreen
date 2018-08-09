################################################################################
# Sample CMake file for Voreen modules.                                        #
#                                                                              #
# Copy this file to the module directory and rename it to match the dir name,  #
# e.g. mymodule/mymodule.cmake                                                 #
#                                                                              #
# Unnecessary output variables and sections may be removed.                    #
################################################################################


################################################################################
# External dependencies (optional)
################################################################################

# sets the preprocessor definition "MY_DEFINE"
SET(MOD_DEFINITIONS -DMY_DEFINE)

IF(WIN32)
    # include directories of external libraries
    SET(MOD_INCLUDE_DIRECTORIES
        ${MOD_DIR}/ext/lib1/include
        ${MOD_DIR}/ext/lib2/include
    )

    # external libraries to link.
	# Note, that 32 bit builds are no longer officially supported.
	
	# debug and release libraries
	SET(MOD_DEBUG_LIBRARIES 
		${MOD_DIR}/ext/lib1/debug/lib1.lib
		${MOD_DIR}/ext/lib2/debug/lib2.lib
	)
	SET(MOD_RELEASE_LIBRARIES 
		${MOD_DIR}/ext/lib1/release/lib1.lib
		${MOD_DIR}/ext/lib2/release/lib2.lib
	)

	# if no separate debug/release libs are present, use MOD_LIBRARIES
	#SET(MOD_LIBRARIES )
			
	# debug/release DLLs to be copied to the project's binary directory
	SET(MOD_DEBUG_DLLS
		${MOD_DIR}/ext/lib1/debug/lib1.dll
		${MOD_DIR}/ext/lib2/debug/lib2.dll
	)
	SET(MOD_RELEASE_DLLS
		${MOD_DIR}/ext/lib1/release/lib1.dll
		${MOD_DIR}/ext/lib2/release/lib2.dll
	)

ELSEIF(UNIX)
    # use cmake module for detecting library, if available
    FIND_PACKAGE(MyLibrary REQUIRED)
    IF(MYLIBRARY_FOUND) #< look into cmake module documentation for specific output variables
        MESSAGE(STATUS "  - Found MyLibrary library")
        SET(MOD_INCLUDE_DIRECTORIES ${MYLIBRARY_INCLUDE_DIR})
        SET(MOD_LIBRARIES ${MYLIBRARY_LIBRARIES}
    ELSE()
        MESSAGE(FATAL_ERROR "MyLibrary library not found!")
    ENDIF()
    
    # alternative: let compiler/linker find library
    SET(MOD_LIBRARIES -lmylibrary1 -lmylibrary2)
ENDIF()


################################################################################
# Core module resources 
################################################################################

# module class must reside in mymodule/mymodule.h + mymodule/mymodule.cpp
SET(MOD_CORE_MODULECLASS MyModule)

# module's core source files, path relative to module dir
SET(MOD_CORE_SOURCES
    ${MOD_DIR}/processors/myproc1.cpp
    ${MOD_DIR}/processors/myproc2.cpp
)

# module's core header files, path relative to module dir
SET(MOD_CORE_HEADERS
    ${MOD_DIR}/processors/myproc1.h
    ${MOD_DIR}/processors/myproc2.h
)   

# module applications without Qt dependency
SET(MOD_CORE_APPLICATIONS
    ${MOD_DIR}/apps/mycoretool.cpp
)


################################################################################
# Qt module resources (optional)
################################################################################
SET(MOD_QT_MODULECLASS MyModuleQt)

SET(MOD_QT_SOURCES
    ${MOD_DIR}/qt/mywidget.cpp
    ${MOD_DIR}/qt/nonmocwidget.cpp
)  
    
# widget headers that have to be run through Qt's moc
SET(MOD_QT_HEADERS
    ${MOD_DIR}/qt/mywidget.h
)

# widget headers that do NOT have to be run through moc
SET(MOD_QT_HEADERS_NONMOC
    ${MOD_DIR}/qt/nonmocwidget.h
)

SET(MOD_QT_FORMS_HEADERS
    #${MOD_DIR}/myform.ui
)
    
# module applications that are linked against the voreen_qt library
SET(MOD_QT_APPLICATIONS
    ${MOD_DIR}/apps/myqttool.cpp
)


################################################################################
# VoreenVE module resources (optional)
################################################################################
#SET(MOD_VE_MODULECLASS MyModuleVE)

SET(MOD_VE_SOURCES
)  
    
SET(MOD_VE_HEADERS
)

SET(MOD_VE_HEADERS_NONMOC
)
