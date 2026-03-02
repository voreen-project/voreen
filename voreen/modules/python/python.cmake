IF(NOT VRN_MODULE_BASE)
    MESSAGE(FATAL_ERROR "Python Module requires Base Module")
ENDIF()

################################################################################
# External dependency: Python library
################################################################################

IF(WIN32)
    SET(VRN_USE_PYTHON_VERSION 312)
    SET(MOD_DEFINITIONS "-DVRN_USE_PYTHON_VERSION=\"Python${VRN_USE_PYTHON_VERSION}\"")

    SET(MOD_INCLUDE_DIRECTORIES "${MOD_DIR}/ext/Python${VRN_USE_PYTHON_VERSION}/include")

    SET(MOD_RELEASE_DLLS
        "${MOD_DIR}/ext/Python${VRN_USE_PYTHON_VERSION}/python${VRN_USE_PYTHON_VERSION}.dll"
    )
    SET(MOD_DEBUG_DLLS
        "${MOD_DIR}/ext/Python${VRN_USE_PYTHON_VERSION}/python${VRN_USE_PYTHON_VERSION}_d.dll"
    )
    SET(MOD_RELEASE_LIBRARIES
        "${MOD_DIR}/ext/Python${VRN_USE_PYTHON_VERSION}/libs/python${VRN_USE_PYTHON_VERSION}.lib"
    )
    SET(MOD_DEBUG_LIBRARIES
        "${MOD_DIR}/ext/Python${VRN_USE_PYTHON_VERSION}/libs/python${VRN_USE_PYTHON_VERSION}_d.lib"
    )
    
    # deployment
    SET(MOD_INSTALL_DIRECTORIES
        ${MOD_DIR}/ext/Python${VRN_USE_PYTHON_VERSION}/lib
        ${MOD_DIR}/scripts
        ${MOD_DIR}/workspaces
    )
    SET(MOD_INSTALL_FILES
        ${MOD_DIR}/ext/Python${VRN_USE_PYTHON_VERSION}/LICENSE.txt
    )

ELSEIF(UNIX)
    MESSAGE(STATUS "Trying to find Python 3 version...")
    
    FIND_PACKAGE(Python3 COMPONENTS Interpreter Development REQUIRED)
    
    SET(MOD_INCLUDE_DIRECTORIES ${Python3_INCLUDE_DIRS})
    SET(MOD_LIBRARIES ${Python3_LIBRARIES})
    MESSAGE(STATUS ${MOD_LIBRARIES})

    # On macOS, Python framework libraries often use an @rpath install name
    # (e.g. @rpath/Python3.framework/Versions/3.9/Python3). Ensure the parent
    # Frameworks directory is part of runtime search paths.
    IF(APPLE)
        FOREACH(_py_lib ${Python3_LIBRARIES})
            IF(_py_lib MATCHES ".*/Frameworks/[^/]+\\.framework/Versions/[^/]+/lib/[^/]+$")
                GET_FILENAME_COMPONENT(_py_lib_dir "${_py_lib}" DIRECTORY)              # .../Versions/X.Y/lib
                GET_FILENAME_COMPONENT(_py_version_dir "${_py_lib_dir}" DIRECTORY)      # .../Versions/X.Y
                GET_FILENAME_COMPONENT(_py_versions_dir "${_py_version_dir}" DIRECTORY) # .../Versions
                GET_FILENAME_COMPONENT(_py_framework_dir "${_py_versions_dir}" DIRECTORY)# .../*.framework
                GET_FILENAME_COMPONENT(_py_frameworks_dir "${_py_framework_dir}" DIRECTORY) # .../Frameworks

                LIST(APPEND CMAKE_BUILD_RPATH "${_py_frameworks_dir}")
                LIST(APPEND CMAKE_INSTALL_RPATH "${_py_frameworks_dir}")
            ENDIF()
        ENDFOREACH()
        LIST(REMOVE_DUPLICATES CMAKE_BUILD_RPATH)
        LIST(REMOVE_DUPLICATES CMAKE_INSTALL_RPATH)
    ENDIF()
    
    # deployment
    SET(MOD_INSTALL_DIRECTORIES
        ${MOD_DIR}/scripts
        ${MOD_DIR}/workspaces
    )
ENDIF()


################################################################################
# Core module resources 
################################################################################
SET(MOD_CORE_MODULECLASS PythonModule)

SET(MOD_CORE_SOURCES
    ${MOD_DIR}/core/pythonscript.cpp
    ${MOD_DIR}/core/pyvoreen.cpp
    ${MOD_DIR}/core/pyvoreenobjects.cpp
    ${MOD_DIR}/properties/pythonproperty.cpp
    ${MOD_DIR}/processors/dynamicpythonprocessor.cpp
)

SET(MOD_CORE_HEADERS
    ${MOD_DIR}/core/pythonoutputlistener.h
    ${MOD_DIR}/core/pythonscript.h
    ${MOD_DIR}/core/pyvoreen.h
    ${MOD_DIR}/core/pyvoreenobjects.h
    ${MOD_DIR}/properties/pythonproperty.h
    ${MOD_DIR}/processors/dynamicpythonprocessor.h
)


################################################################################
# Qt module resources 
################################################################################
SET(MOD_QT_MODULECLASS PythonModuleQt)

SET(MOD_QT_SOURCES
    ${MOD_DIR}/qt/dynamicpythonwidget.cpp
    ${MOD_DIR}/qt/pyvoreenqt.cpp
    ${MOD_DIR}/qt/pythonhighlighter.cpp
    ${MOD_DIR}/qt/pythonplugin.cpp
    ${MOD_DIR}/qt/pythonprocessorwidgetfactory.cpp
    ${MOD_DIR}/qt/pythonpropertywidget.cpp
    ${MOD_DIR}/qt/pythonpropertywidgetfactory.cpp
    ${MOD_DIR}/qt/menuentity/pythoneditor.cpp
)

SET(MOD_QT_HEADERS
    ${MOD_DIR}/qt/dynamicpythonwidget.h
    ${MOD_DIR}/qt/pythonplugin.h
    ${MOD_DIR}/qt/pythonprocessorwidgetfactory.h
    ${MOD_DIR}/qt/pythonpropertywidget.h
    ${MOD_DIR}/qt/pythonpropertywidgetfactory.h
    ${MOD_DIR}/qt/menuentity/pythoneditor.h
)

SET(MOD_QT_HEADERS_NONMOC
    ${MOD_DIR}/qt/pyvoreenqt.h
    ${MOD_DIR}/qt/pythonhighlighter.h
)

SET(MOD_QT_RESOURCES
    ${MOD_DIR}/qt/menuentity/python.qrc
)
