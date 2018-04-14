
################################################################################
# External dependency: Python library
################################################################################

IF(WIN32)
    SET(MOD_INCLUDE_DIRECTORIES "${MOD_DIR}/ext/python27/include")

	SET(MOD_DEBUG_LIBRARIES 
		"${MOD_DIR}/ext/python27/lib/debug/python27_d.lib"
	)
	SET(MOD_RELEASE_LIBRARIES 
		"${MOD_DIR}/ext/python27/lib/release/python27.lib"
	)
	
	SET(MOD_DEBUG_DLLS
		"${MOD_DIR}/ext/python27/lib/debug/python27_d.dll"
	)
	SET(MOD_RELEASE_DLLS 
		"${MOD_DIR}/ext/python27/lib/release/python27.dll"
	)
    
    # deployment
    SET(MOD_INSTALL_DIRECTORIES
        ${MOD_DIR}/scripts
    )
    SET(MOD_INSTALL_FILES
        ${MOD_DIR}/ext/python27/LICENSE
    )

ELSEIF(UNIX)
    # this is a hack to prevent CMake from finding python 3 instead of 2.7
    MESSAGE(STATUS "Trying to find Python 2.7 version...")
    FOREACH(i RANGE 20)
        SET(python_base "2.7.")
        FIND_PACKAGE(PythonLibs ${python_base}${i} EXACT) #REQUIRED)
        IF(PYTHONLIBS_FOUND)
            BREAK()
        ENDIF()
    ENDFOREACH(i)
    #FIND_PACKAGE(PythonLibs REQUIRED)
    IF(PYTHONLIBS_FOUND)
        MESSAGE(STATUS "  - Found Python library")
        SET(MOD_INCLUDE_DIRECTORIES ${PYTHON_INCLUDE_DIRS})
        SET(MOD_LIBRARIES ${PYTHON_LIBRARIES})
    ELSE()
        MESSAGE(FATAL_ERROR "Python library not found!")
    ENDIF()
ENDIF()


################################################################################
# Core module resources 
################################################################################
SET(MOD_CORE_MODULECLASS PythonModule)

SET(MOD_CORE_SOURCES 
    ${MOD_DIR}/core/pythonscript.cpp
    ${MOD_DIR}/core/pyvoreen.cpp
)

SET(MOD_CORE_HEADERS
    ${MOD_DIR}/core/pythonscript.h
    ${MOD_DIR}/core/pyvoreen.h
)


################################################################################
# Qt module resources 
################################################################################
SET(MOD_QT_MODULECLASS PythonModuleQt)

SET(MOD_QT_SOURCES 
    ${MOD_DIR}/qt/pyvoreenqt.cpp
    ${MOD_DIR}/qt/pythonhighlighter.cpp
    ${MOD_DIR}/qt/menuentity/pythoneditor.cpp
)

SET(MOD_QT_HEADERS
    ${MOD_DIR}/qt/menuentity/pythoneditor.h
)

SET(MOD_QT_HEADERS_NONMOC
    ${MOD_DIR}/qt/pyvoreenqt.h
    ${MOD_DIR}/qt/pythonhighlighter.h
)

SET(MOD_QT_RESOURCES
    ${MOD_DIR}/qt/menuentity/python.qrc
)
