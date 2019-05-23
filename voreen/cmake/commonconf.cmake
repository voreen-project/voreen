IF(NOT COMMONCONF_PROCESSED)

SET(VRN_HOME ${CMAKE_CURRENT_SOURCE_DIR})
MESSAGE(STATUS "Voreen Home: ${VRN_HOME}")

# include macros and config
INCLUDE(${VRN_HOME}/cmake/macros.cmake)
IF(EXISTS ${VRN_HOME}/config.cmake)
    MESSAGE(STATUS "Including custom configuration file 'config.cmake'")
    INCLUDE(${VRN_HOME}/config.cmake)
ELSE()
    INCLUDE(${VRN_HOME}/config-default.cmake)
ENDIF()

# set release mode by default when under unix
IF(UNIX)
	IF(NOT CMAKE_BUILD_TYPE) 
		SET(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel." FORCE) 
	ENDIF()
ENDIF()

#PCH is currently problematic under linux -> display message and disable if set
IF(UNIX AND VRN_PRECOMPILED_HEADER)
	MESSAGE(WARNING "Precompiled Headers currently not supported for linux builds - disabling PCH")
	SET(VRN_PRECOMPILED_HEADER OFF CACHE BOOL "Use pre-compiled headers?" FORCE)
ENDIF()

# set/create binary output path
IF(NOT VRN_BINARY_OUTPUT_DIR)
    SET(VRN_BINARY_OUTPUT_DIR ${VRN_HOME}/bin)
ENDIF()
IF(NOT EXISTS ${VRN_BINARY_OUTPUT_DIR})
    MESSAGE(STATUS "VRN_BINARY_OUTPUT_DIR does not exist: ${VRN_BINARY_OUTPUT_DIR}. Creating it ...")
    FILE(MAKE_DIRECTORY ${VRN_BINARY_OUTPUT_DIR})
    IF(NOT EXISTS ${VRN_BINARY_OUTPUT_DIR})
        MESSAGE(FATAL_ERROR "Failed to create VRN_BINARY_OUTPUT_DIR: ${VRN_BINARY_OUTPUT_DIR}")
    ENDIF()
ENDIF()
MESSAGE(STATUS "Output Path: ${VRN_BINARY_OUTPUT_DIR}")
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${VRN_BINARY_OUTPUT_DIR})
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${VRN_BINARY_OUTPUT_DIR})
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${VRN_BINARY_OUTPUT_DIR})

# detect compiler and architecture
IF(${CMAKE_GENERATOR} STREQUAL "Visual Studio 11" OR ${CMAKE_GENERATOR} STREQUAL "Visual Studio 11 2012" OR
   ${CMAKE_GENERATOR} STREQUAL "Visual Studio 12" OR ${CMAKE_GENERATOR} STREQUAL "Visual Studio 12 2013" OR
   ${CMAKE_GENERATOR} STREQUAL "Visual Studio 14" OR ${CMAKE_GENERATOR} STREQUAL "Visual Studio 14 2015" OR
   ${CMAKE_GENERATOR} STREQUAL "Visual Studio 15" OR ${CMAKE_GENERATOR} STREQUAL "Visual Studio 15 2017")
   MESSAGE(FATAL_ERROR "32 Bit is no longer supported: ${CMAKE_GENERATOR}. Please use a 64 Bit generator.")
ELSEIF(${CMAKE_GENERATOR} STREQUAL "Visual Studio 11 Win64" OR
       ${CMAKE_GENERATOR} STREQUAL "Visual Studio 11 2012 Win64")
    SET(VRN_MSVC2012 TRUE)
    SET(VRN_MSVC TRUE)
    MESSAGE(STATUS "Visual Studio 2012 Build (64 Bit)")
ELSEIF(${CMAKE_GENERATOR} STREQUAL "Visual Studio 12 Win64" OR
       ${CMAKE_GENERATOR} STREQUAL "Visual Studio 12 2013 Win64")
    SET(VRN_MSVC2013 TRUE)
    SET(VRN_MSVC TRUE)
    MESSAGE(STATUS "Visual Studio 2013 Build (64 bit)")
ELSEIF(${CMAKE_GENERATOR} STREQUAL "Visual Studio 14 Win64" OR
       ${CMAKE_GENERATOR} STREQUAL "Visual Studio 14 2015 Win64")
    SET(VRN_MSVC2015 TRUE)
    SET(VRN_MSVC TRUE)
    MESSAGE(STATUS "Visual Studio 2015 Build (64 bit)")
ELSEIF(${CMAKE_GENERATOR} STREQUAL "Visual Studio 15 Win64" OR
       ${CMAKE_GENERATOR} STREQUAL "Visual Studio 15 2017 Win64")
    SET(VRN_MSVC2017 TRUE)
    SET(VRN_MSVC TRUE)
    MESSAGE(STATUS "Visual Studio 2017 Build (64 bit)")
ELSEIF(${CMAKE_GENERATOR} MATCHES "Unix" OR ${CMAKE_GENERATOR} MATCHES "Ninja")
    SET(VRN_UNIX TRUE)
    MESSAGE(STATUS "Unix Build")
ELSEIF(${CMAKE_GENERATOR} MATCHES "Xcode")
    SET(VRN_UNIX TRUE)
    SET(VRN_XCODE TRUE)
    MESSAGE(STATUS "Xcode Build")
ELSE()
    MESSAGE(FATAL_ERROR "Unsupported or unknown generator: ${CMAKE_GENERATOR}. Please use VC11, VC12 VC14, VC15 or Unix or Xcode.")
ENDIF()

# common include directories
LIST(APPEND VRN_COMMON_INCLUDE_DIRECTORIES "${VRN_HOME}" "${VRN_HOME}/include" "${VRN_HOME}/ext")
LIST(APPEND VRN_COMMON_INCLUDE_DIRECTORIES ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}) 

# Configure shared library-build, static builds are no longer supported
SET(BUILD_SHARED_LIBS TRUE)

# print pch status
IF(VRN_PRECOMPILED_HEADER)
    MESSAGE(STATUS "Precompiled Headers: Enabled")
ELSE()
    MESSAGE(STATUS "Precompiled Headers: Disabled")
ENDIF()

# platform-dependent configuration
IF(VRN_MSVC)
    LIST(APPEND VRN_DEFINITIONS -DNOMINMAX -D_CRT_SECURE_NO_DEPRECATE -DPSAPI_VERSION=1)
    
    # Windows API dependencies
    LIST(APPEND VRN_EXTERNAL_LIBRARIES netapi32 version psapi)

    # Disable warnings for Microsoft compiler:
    # C4305: 'identifier' : truncation from 'type1' to 'type2'
    # C4800: 'type' : forcing value to bool 'true' or 'false' (performance warning
    # C4290: C++ exception specification ignored except to indicate a function is
    #        not __declspec(nothrow)
    # C4068: unknown pragma
    # C4251  class needs to have dll interface (used for std classes)
    # C4355: 'this' : used in base member initializer list 
    #        occurs in processors' constructors when initializing event properties, 
    #        but is safe there, since the 'this' pointer is only stored and not accessed.
    # C4390: ';' : empty controlled statement found; is this the intent?
    #        occurs when OpenGL error logging macros are disabled
    LIST(APPEND VRN_DEFINITIONS /wd4305 /wd4800 /wd4290 /wd4068 /wd4251 /wd4355 /wd4390)
    
    # enable parallel builds in Visual Studio
    LIST(APPEND VRN_DEFINITIONS /MP)

    # prevent error: number of sections exceeded object file format limit
    LIST(APPEND VRN_DEFINITIONS /bigobj)
    
    # prevents rarely-used header files from being automatically included by windows.h (esp. winsock) 
    LIST(APPEND VRN_DEFINITIONS -DWIN32_LEAN_AND_MEAN)
    
    # disable warning on std::copy call with unchecked parameters
    LIST(APPEND VRN_DEFINITIONS -D_SCL_SECURE_NO_WARNINGS)
    
    # Linking against Windows DLLs requires explicit instantiation of templates
    LIST(APPEND VRN_DEFINITIONS -DDLL_TEMPLATE_INST)

    IF(NOT VRN_GENERATE_MANIFEST)
        # Do not embed manifest into binaries in debug mode (slows down incremental linking)
        SET(CMAKE_SHARED_LINKER_FLAGS_DEBUG "${CMAKE_SHARED_LINKER_FLAGS_DEBUG} /MANIFEST:NO")
        SET(CMAKE_EXE_LINKER_FLAGS_DEBUG    "${CMAKE_EXE_LINKER_FLAGS_DEBUG} /MANIFEST:NO")
    ENDIF()

    # set RAM usage to 1000% for PCH (only update on win64 or problems with opencl)
    IF(VRN_PRECOMPILED_HEADER)
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Zm1000")
    ENDIF()
    
    # enable/disable incremental linking in debug builds
    If(VRN_INCREMENTAL_LINKING)
        IF(NOT VRN_PRECOMPILED_HEADER)
            MESSAGE("Incremental linking only available when using precompiled headers (VRN_PRECOMPILED_HEADER).")
            SET(CMAKE_SHARED_LINKER_FLAGS_DEBUG "${CMAKE_SHARED_LINKER_FLAGS_DEBUG} /INCREMENTAL:NO")
            SET(CMAKE_EXE_LINKER_FLAGS_DEBUG    "${CMAKE_EXE_LINKER_FLAGS_DEBUG} /INCREMENTAL:NO")
        ELSE()
            SET(CMAKE_SHARED_LINKER_FLAGS_DEBUG "${CMAKE_SHARED_LINKER_FLAGS_DEBUG} /INCREMENTAL")
            SET(CMAKE_EXE_LINKER_FLAGS_DEBUG    "${CMAKE_EXE_LINKER_FLAGS_DEBUG} /INCREMENTAL")
            IF((MSVC_TOOLSET_VERSION EQUAL 140) OR (MSVC_TOOLSET_VERSION GREATER 140))
                MESSAGE(STATUS "MSVC: Edit and Continue supported")
                SET(CMAKE_SHARED_LINKER_FLAGS_DEBUG "${CMAKE_SHARED_LINKER_FLAGS_DEBUG} /EDITANDCONTINUE")
                SET(CMAKE_EXE_LINKER_FLAGS_DEBUG    "${CMAKE_EXE_LINKER_FLAGS_DEBUG} /EDITANDCONTINUE")
            ENDIF()
        ENDIF()
    ELSE()
        SET(CMAKE_SHARED_LINKER_FLAGS_DEBUG "${CMAKE_SHARED_LINKER_FLAGS_DEBUG} /INCREMENTAL:NO")
        SET(CMAKE_EXE_LINKER_FLAGS_DEBUG    "${CMAKE_EXE_LINKER_FLAGS_DEBUG} /INCREMENTAL:NO")
    ENDIF()
    
    # suppress ZERO_CHECK target
    IF(VRN_DISABLE_ZERO_CHECK)
        SET(CMAKE_SUPPRESS_REGENERATION TRUE)
    ENDIF()
    
    # Windows deployment   
    IF(VRN_DEPLOYMENT)
        LIST(APPEND VRN_DEFINITIONS "-DVRN_DEPLOYMENT") 
        MESSAGE(STATUS "Windows deployment build:")

        MESSAGE(STATUS "* Adding install target")
        SET(VRN_ADD_INSTALL_TARGET ON)
        MESSAGE(STATUS "* Install directory (CMAKE_INSTALL_PREFIX): ${CMAKE_INSTALL_PREFIX}")

        MESSAGE(STATUS "* Adding Visual Studio redist libraries to install target")
        IF(VRN_MSVC2012)
            GET_FILENAME_COMPONENT(VS_DIR "[HKEY_LOCAL_MACHINE\\SOFTWARE\\WOW6432Node\\Microsoft\\VisualStudio\\SxS\\VS7;11.0]" REALPATH)
            IF(NOT EXISTS ${VS_DIR})
                MESSAGE(WARNING "Visual Studio directory not found: ${VS_DIR}")
            ELSE()
                MESSAGE(STATUS "  - ${VS_DIR}/VC/redist/x64/Microsoft.VC110.CRT/msvcp110.dll")
                MESSAGE(STATUS "  - ${VS_DIR}/VC/redist/x64/Microsoft.VC110.CRT/msvcr110.dll")
                MESSAGE(STATUS "  - ${VS_DIR}/VC/redist/x64/Microsoft.VC110.OpenMP/vcomp110.dll")
                INSTALL(FILES "${VS_DIR}/VC/redist/x64/Microsoft.VC110.CRT/msvcp110.dll" DESTINATION .)
                INSTALL(FILES "${VS_DIR}/VC/redist/x64/Microsoft.VC110.CRT/msvcr110.dll" DESTINATION .)
                INSTALL(FILES "${VS_DIR}/VC/redist/x64/Microsoft.VC110.OpenMP/vcomp110.dll" DESTINATION .)
            ENDIF()
        ELSEIF(VRN_MSVC2013)
            GET_FILENAME_COMPONENT(VS_DIR "[HKEY_LOCAL_MACHINE\\SOFTWARE\\WOW6432Node\\Microsoft\\VisualStudio\\SxS\\VS7;12.0]" REALPATH)
            IF(NOT EXISTS ${VS_DIR})
                MESSAGE(WARNING "Visual Studio directory not found: ${VS_DIR}")
            ELSE()
                MESSAGE(STATUS "  - ${VS_DIR}/VC/redist/x64/Microsoft.VC120.CRT/msvcp120.dll")
                MESSAGE(STATUS "  - ${VS_DIR}/VC/redist/x64/Microsoft.VC120.CRT/msvcr120.dll")
                MESSAGE(STATUS "  - ${VS_DIR}/VC/redist/x64/Microsoft.VC120.OpenMP/vcomp120.dll")
                INSTALL(FILES "${VS_DIR}/VC/redist/x64/Microsoft.VC120.CRT/msvcp120.dll" DESTINATION .)
                INSTALL(FILES "${VS_DIR}/VC/redist/x64/Microsoft.VC120.CRT/msvcr120.dll" DESTINATION .)
                INSTALL(FILES "${VS_DIR}/VC/redist/x64/Microsoft.VC120.OpenMP/vcomp120.dll" DESTINATION .)
            ENDIF()
		ELSEIF(VRN_MSVC2015)
            GET_FILENAME_COMPONENT(VS_DIR "[HKEY_LOCAL_MACHINE\\SOFTWARE\\WOW6432Node\\Microsoft\\VisualStudio\\SxS\\VS7;14.0]" REALPATH)
            IF(NOT EXISTS ${VS_DIR})
                MESSAGE(WARNING "Visual Studio directory not found: ${VS_DIR}")
            ELSE()
                MESSAGE(STATUS "  - ${VS_DIR}/VC/redist/x64/Microsoft.VC140.CRT/msvcp140.dll")
                MESSAGE(STATUS "  - ${VS_DIR}/VC/redist/x64/Microsoft.VC140.CRT/vcruntime.dll")
                MESSAGE(STATUS "  - ${VS_DIR}/VC/redist/x64/Microsoft.VC140.OpenMP/vcomp140.dll")
                INSTALL(FILES "${VS_DIR}/VC/redist/x64/Microsoft.VC140.CRT/msvcp140.dll" DESTINATION .)
                INSTALL(FILES "${VS_DIR}/VC/redist/x64/Microsoft.VC140.CRT/vcruntime140.dll" DESTINATION .)
                INSTALL(FILES "${VS_DIR}/VC/redist/x64/Microsoft.VC140.OpenMP/vcomp140.dll" DESTINATION .)
            ENDIF()
		ELSEIF(VRN_MSVC2017)
            GET_FILENAME_COMPONENT(VS_DIR "[HKEY_LOCAL_MACHINE\\SOFTWARE\\WOW6432Node\\Microsoft\\VisualStudio\\SxS\\VS7;15.0]" REALPATH)
            GET_FILENAME_COMPONENT(VS_VER "[HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\DevDiv\\vc\\Servicing\\14.0\\RuntimeDebug;Version]" NAME)
            IF(NOT EXISTS "${VS_DIR}/VC/Redist/MSVC/${VS_VER}")
                MESSAGE(WARNING "Visual Studio directory not found: ${VS_DIR}/VC/Redist/MSVC/${VS_VER}")
            ELSE()
                MESSAGE(STATUS "  - ${VS_DIR}/VC/Redist/MSVC/${VS_VER}/x64/Microsoft.VC141.CRT/msvcp140.dll")
                MESSAGE(STATUS "  - ${VS_DIR}/VC/redist/MSVC/${VS_VER}/x64/Microsoft.VC141.CRT/vcruntime140.dll")
                MESSAGE(STATUS "  - ${VS_DIR}/VC/redist/MSVC/${VS_VER}/x64/Microsoft.VC141.OpenMP/vcomp140.dll")
                INSTALL(FILES "${VS_DIR}/VC/Redist/MSVC/${VS_VER}/x64/Microsoft.VC141.CRT/msvcp140.dll" DESTINATION .)
                INSTALL(FILES "${VS_DIR}/VC/Redist/MSVC/${VS_VER}/x64/Microsoft.VC141.CRT/vcruntime140.dll" DESTINATION .)
                INSTALL(FILES "${VS_DIR}/VC/Redist/MSVC/${VS_VER}/x64/Microsoft.VC141.OpenMP/vcomp140.dll" DESTINATION .)
            ENDIF()
        ELSE()
            MESSAGE(WARNING "Deploying redist libraries only supported for Visual Studio 2012, 2013, 2015 or 2017.")
        ENDIF()
    ELSE(VRN_DEPLOYMENT)
        # hardcode Voreen base path, if binary output dir has been modified and we are not in deployment mode
        IF(NOT ${VRN_BINARY_OUTPUT_DIR} STREQUAL "${VRN_HOME}/bin")
            LIST(APPEND VRN_DEFINITIONS "-DVRN_BASE_PATH=\"${VRN_HOME}\"")
        ENDIF()
    ENDIF(VRN_DEPLOYMENT)

ELSEIF(UNIX)

    LIST(APPEND VRN_DEFINITIONS "-DUNIX")
    LIST(APPEND VRN_DEFINITIONS "-D__STDC_CONSTANT_MACROS")

    IF(VRN_DEPLOYMENT)
        LIST(APPEND VRN_DEFINITIONS "-DVRN_DEPLOYMENT")
        MESSAGE(STATUS "Unix deployment build")
        
        MESSAGE(STATUS "* Adding install target")
        SET(VRN_ADD_INSTALL_TARGET ON)
        MESSAGE(STATUS "* Install directory (CMAKE_INSTALL_PREFIX): ${CMAKE_INSTALL_PREFIX}")
    ELSE()
        # hardcode Voreen base path, if binary output dir has been modified and we are not in deployment mode
        IF(NOT ${VRN_BINARY_OUTPUT_DIR} STREQUAL "${VRN_HOME}/bin")
            LIST(APPEND VRN_DEFINITIONS "-DVRN_BASE_PATH=\"${VRN_HOME}\"")
        ENDIF()
    ENDIF()

    include(CheckCXXCompilerFlag)
    CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
    CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
    if(COMPILER_SUPPORTS_CXX11)
        set(CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS}")
    elseif(COMPILER_SUPPORTS_CXX0X)
        set(CMAKE_CXX_FLAGS "-std=c++0x ${CMAKE_CXX_FLAGS}")
    else()
        message(STATUS "Compiler ${CMAKE_CXX_COMPILER} does not support C++11.")
    endif()

    set(CMAKE_CXX_FLAGS_RELEASE "-O3 ${CMAKE_CXX_FLAGS_RELEASE}")
    # enable optimization level 1 for debug build as -Wuninitialized is ignored otherwise
    #set(CMAKE_CXX_FLAGS_DEBUG "-O1 ${CMAKE_CXX_FLAGS_DEBUG}")

    # add linker flags to look for libraries in the executable directory to avoid problems with moving Voreen after compiling
    IF(APPLE)
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-rpath,'$ORIGIN' ")
    ELSE()
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-rpath='$ORIGIN' ")
    ENDIF()

    # Warning switches are compiler dependent:
    IF(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        # disable warnings
        set(CMAKE_CXX_FLAGS "-Wno-unused-function ${CMAKE_CXX_FLAGS}")
        set(CMAKE_CXX_FLAGS "-Wno-unused-variable ${CMAKE_CXX_FLAGS}")
        set(CMAKE_CXX_FLAGS "-Wno-reorder ${CMAKE_CXX_FLAGS}")
        set(CMAKE_CXX_FLAGS "-Wno-parentheses ${CMAKE_CXX_FLAGS}")

        #enable warnings
        set(CMAKE_CXX_FLAGS "-Wunused-but-set-variable ${CMAKE_CXX_FLAGS}")
        set(CMAKE_CXX_FLAGS "-Wuninitialized -Wmaybe-uninitialized ${CMAKE_CXX_FLAGS}")
        set(CMAKE_CXX_FLAGS "-Wall ${CMAKE_CXX_FLAGS}")
    ELSEIF(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        # disable warnings
        set(CMAKE_CXX_FLAGS "-Wno-unused-function ${CMAKE_CXX_FLAGS}")
        set(CMAKE_CXX_FLAGS "-Wno-unused-variable ${CMAKE_CXX_FLAGS}")
        set(CMAKE_CXX_FLAGS "-Wno-reorder ${CMAKE_CXX_FLAGS}")
        set(CMAKE_CXX_FLAGS "-Wno-parentheses ${CMAKE_CXX_FLAGS}")
        set(CMAKE_CXX_FLAGS "-Wno-inconsistent-missing-override ${CMAKE_CXX_FLAGS}") #someone should probably fix this at some point
        set(CMAKE_CXX_FLAGS "-Wno-potentially-evaluated-expression ${CMAKE_CXX_FLAGS}")
        set(CMAKE_CXX_FLAGS "-Wno-missing-braces ${CMAKE_CXX_FLAGS}") #but we enable missing-field-initializers!
        set(CMAKE_CXX_FLAGS "-Wno-invalid-source-encoding ${CMAKE_CXX_FLAGS}") #Strings with umlauts are decoded using qt

        #enable warnings
        set(CMAKE_CXX_FLAGS "-Wunused-value -Wunused-const-variable ${CMAKE_CXX_FLAGS}")
        set(CMAKE_CXX_FLAGS "-Wmissing-field-initializers ${CMAKE_CXX_FLAGS}")
        set(CMAKE_CXX_FLAGS "-Wuninitialized ${CMAKE_CXX_FLAGS}")
        set(CMAKE_CXX_FLAGS "-Wall ${CMAKE_CXX_FLAGS}")
    ELSE()
        message(Status "Unknown compiler ID ${CMAKE_CXX_COMPILER_ID}")
    ENDIF()
ENDIF()

# macosx
IF(APPLE)
    FIND_LIBRARY(COREFOUNDATION_LIBRARY CoreFoundation )
    LIST(APPEND VRN_EXTERNAL_LIBRARIES ${COREFOUNDATION_LIBRARY})
    
    # disable warnings
    set(CMAKE_CXX_FLAGS "-Wno-deprecated-register ${CMAKE_CXX_FLAGS}")
    set(CMAKE_CXX_FLAGS "-Wno-unused-function ${CMAKE_CXX_FLAGS}")
    set(CMAKE_CXX_FLAGS "-Wno-unused-private-field ${CMAKE_CXX_FLAGS}")
    set(CMAKE_CXX_FLAGS "-Wno-unused-variable ${CMAKE_CXX_FLAGS}")
    set(CMAKE_CXX_FLAGS "-Wno-reorder ${CMAKE_CXX_FLAGS}")
    set(CMAKE_CXX_FLAGS "-Wno-invalid-source-encoding ${CMAKE_CXX_FLAGS}")
    set(CMAKE_CXX_FLAGS "-Wall ${CMAKE_CXX_FLAGS}")
    
    # on apple build against the new native libc++
    set(CMAKE_CXX_FLAGS "-std=c++11 -stdlib=libc++ ${CMAKE_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "-stdlib=libc++ ${CMAKE_EXE_LINKER_FLAGS}")
    set(CMAKE_SHARED_LINKER_FLAGS "-stdlib=libc++ ${CMAKE_SHARED_LINKER_FLAGS}")
ENDIF(APPLE)

# use STL in tinyXML
LIST(APPEND VRN_DEFINITIONS "-DTIXML_USE_STL") 

# tgt configuration
LIST(APPEND VRN_DEFINITIONS "-DTGT_WITHOUT_DEFINES") # don't use tgt's build system
IF(VRN_MSVC)
    SET(TGT_WITH_WMI TRUE)  #< enable Windows Management Instrumentation for hardware detection
ENDIF()

# set Voreen debug flags for debug builds
SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DTGT_DEBUG -DVRN_DEBUG")
 
# minimum Qt version
SET(VRN_REQUIRED_QT_VERSION "5.5")

SET(VRN_NON_INTERACTIVE OFF CACHE BOOL "Assume that voreen is used non-interactively, e.g., do not ask user when an assertion fails.")
IF(VRN_NON_INTERACTIVE)
   LIST(APPEND VRN_DEFINITIONS "-DTGT_NON_INTERACTIVE_ASSERT")
ENDIF()


# detect libraries
MESSAGE(STATUS "--------------------------------------------------------------------------------")
MESSAGE(STATUS "Detecting Common Mandatory Libraries:")

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${VRN_HOME}/cmake")

# OpenGL
FIND_PACKAGE(OpenGL REQUIRED)
IF(OPENGL_FOUND)
    MESSAGE(STATUS "* Found OpenGL")
        IF(VRN_OPENGL_COMPATIBILITY_PROFILE)
            MESSAGE(STATUS "    * OpenGL compatibility context build")
            LIST(APPEND VRN_DEFINITIONS "-DVRN_OPENGL_COMPATIBILITY_PROFILE")
        ENDIF()
    LIST(APPEND VRN_COMMON_INCLUDE_DIRECTORIES ${OPENGL_INCLUDE_DIR})
    LIST(APPEND VRN_EXTERNAL_LIBRARIES ${OPENGL_LIBRARIES})
ELSE(OPENGL_FOUND)
    MESSAGE(FATAL_ERROR "OpenGL not found!")
ENDIF(OPENGL_FOUND)
    
# GLEW
IF (VRN_OPENGL_COMPATIBILITY_PROFILE)
FIND_PACKAGE(GlewVRN REQUIRED)
IF(GLEW_FOUND)
    MESSAGE(STATUS "* Found GLEW")
    LIST(APPEND VRN_DEFINITIONS ${GLEW_DEFINITIONS})
    LIST(APPEND VRN_COMMON_INCLUDE_DIRECTORIES ${GLEW_INCLUDE_DIR})
    LIST(APPEND VRN_EXTERNAL_LIBRARIES ${GLEW_LIBRARY})
    LIST(APPEND VRN_EXTERNAL_DEBUG_DLLS ${GLEW_DLL_DEBUG})
    LIST(APPEND VRN_EXTERNAL_RELEASE_DLLS ${GLEW_DLL_RELEASE})
    LIST(APPEND VRN_EXTERNAL_LICENSE_FILES ${GLEW_LICENSE_FILE})
ELSE(GLEW_FOUND)
    MESSAGE(FATAL_ERROR "GLEW not found!")
ENDIF(GLEW_FOUND)
ENDIF()

# Boost	
FIND_PACKAGE(BoostVRN REQUIRED)
IF(Boost_FOUND)
    MESSAGE(STATUS "* Found Boost")
    LIST(APPEND VRN_DEFINITIONS ${Boost_DEFINITIONS})
    LIST(APPEND VRN_COMMON_INCLUDE_DIRECTORIES ${Boost_INCLUDE_DIRS})
    LIST(APPEND VRN_EXTERNAL_LIBRARIES ${Boost_LIBRARIES})
    LIST(APPEND VRN_EXTERNAL_DEBUG_DLLS ${Boost_DEBUG_DLLS})
    LIST(APPEND VRN_EXTERNAL_RELEASE_DLLS ${Boost_RELEASE_DLLS})
ELSE()
    MESSAGE(FATAL_ERROR "Boost not found!")
ENDIF()

# efsw
IF(WIN32) # default under windows is ON, since the win32 implementation of efsw contains some bugs.
    SET(VRN_USE_GENERIC_FILE_WATCHER ON CACHE BOOL "Using generic file watcher prevents racing conditions under windows and allows remote file systems to be watched.")
ELSE()
    SET(VRN_USE_GENERIC_FILE_WATCHER OFF CACHE BOOL "Using generic file watcher prevents racing conditions under windows and allows remote file systems to be watched.")
ENDIF()
IF(VRN_USE_GENERIC_FILE_WATCHER)
    LIST(APPEND VRN_DEFINITIONS "-DVRN_USE_GENERIC_FILE_WATCHER")
ENDIF()
LIST(APPEND VRN_EXTERNAL_LICENSE_FILES "${VRN_HOME}/ext/efsw/LICENSE")

# tinyxml
LIST(APPEND VRN_EXTERNAL_LICENSE_FILES "${VRN_HOME}/ext/tinyxml/license.txt")
    
# eigen
MESSAGE(STATUS "* Found Eigen")
LIST(APPEND VRN_COMMON_INCLUDE_DIRECTORIES "${VRN_HOME}/ext/eigen")
LIST(APPEND VRN_EXTERNAL_LICENSE_FILES "${VRN_HOME}/ext/eigen/COPYING.BSD")
LIST(APPEND VRN_EXTERNAL_LICENSE_FILES "${VRN_HOME}/ext/eigen/COPYING.GPL")
LIST(APPEND VRN_EXTERNAL_LICENSE_FILES "${VRN_HOME}/ext/eigen/COPYING.LGPL")
LIST(APPEND VRN_EXTERNAL_LICENSE_FILES "${VRN_HOME}/ext/eigen/COPYING.MINPACK")
LIST(APPEND VRN_EXTERNAL_LICENSE_FILES "${VRN_HOME}/ext/eigen/COPYING.MPL2")
LIST(APPEND VRN_EXTERNAL_LICENSE_FILES "${VRN_HOME}/ext/eigen/COPYING.README")
LIST(APPEND VRN_DEFINITIONS -DEIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS)

# rapidjson
MESSAGE(STATUS "* Found rapidjson")
LIST(APPEND VRN_COMMON_INCLUDE_DIRECTORIES "${VRN_HOME}/ext/rapidjson")
LIST(APPEND VRN_EXTERNAL_LICENSE_FILES "${VRN_HOME}/ext/rapidjson/license.txt")

message(STATUS "Collected common include directories: ${VRN_COMMON_INCLUDE_DIRECTORIES}")

#######################
# modules
#######################
MESSAGE(STATUS "--------------------------------------------------------------------------------")
    
# collect module dirs
SET(MODULE_BASEDIR_LIST ${VRN_HOME}/modules) #< framework modules
IF(VRN_CUSTOM_MODULEDIR)
    IF(EXISTS ${VRN_CUSTOM_MODULEDIR})
        LIST(APPEND MODULE_BASEDIR_LIST ${VRN_CUSTOM_MODULEDIR})
    ELSE()
        MESSAGE(WARNING "Custom module dir ${VRN_CUSTOM_MODULEDIR} does not exist!")
    ENDIF()
ENDIF()
FOREACH(num RANGE 0 10)
    IF(VRN_CUSTOM_MODULEDIR_${num})
        IF(EXISTS ${VRN_CUSTOM_MODULEDIR_${num}})
            LIST(APPEND MODULE_BASEDIR_LIST ${VRN_CUSTOM_MODULEDIR_${num}})
        ELSE()
            MESSAGE(WARNING "Custom module dir does not exist: ${VRN_CUSTOM_MODULEDIR_${num}}")
        ENDIF()
    ENDIF()
ENDFOREACH()
LIST(REMOVE_DUPLICATES MODULE_BASEDIR_LIST)

# include modules
SET(VRN_MODULE_CORE ON CACHE BOOL "Core module is always included" FORCE)
MARK_AS_ADVANCED(VRN_MODULE_CORE)

IF(VRN_BUILD_VOREENBIOLOGY)
    SET(VRN_MODULE_VOREENBIOLOGY ON CACHE INTERNAL "VoreenBiology module is included if VoreenBiology is built." FORCE)
    MARK_AS_ADVANCED(VRN_MODULE_VOREENBIOLOGY)
ELSE(VRN_BUILD_VOREENBIOLOGY)
    SET(VRN_MODULE_VOREENBIOLOGY OFF CACHE INTERNAL "VoreenBiology module is included if VoreenBiology is built." FORCE)
    MARK_AS_ADVANCED(VRN_MODULE_VOREENBIOLOGY)
ENDIF(VRN_BUILD_VOREENBIOLOGY)

FOREACH(module_basedir ${MODULE_BASEDIR_LIST})

    MESSAGE(STATUS "Including Voreen Modules from ${module_basedir}:")

    IF(EXISTS ${module_basedir}/modulelist.cmake)
        INCLUDE(${module_basedir}/modulelist.cmake)
    ENDIF()
    
    # iterate over subdirectories of module dir and include each enabled module
    LIST_SUBDIRECTORIES(module_dir_list ${module_basedir} false)
    FOREACH(module_dir ${module_dir_list})
        SET(module_file ${module_basedir}/${module_dir}/${module_dir}.cmake)
        STRING(TOLOWER ${module_dir} module_lower)
        STRING(TOUPPER ${module_dir} module_upper)
        IF(EXISTS ${module_file})
            IF(VRN_MODULE_${module_upper})
                # include module .cmake file
                SET(MOD_DIR ${module_basedir}/${module_dir})
                INCLUDE(${module_file})

                FILE(RELATIVE_PATH module_file_rel ${module_basedir} ${module_file})
                MESSAGE(STATUS "* ${module_file_rel}")

                IF(NOT VRN_OPENGL_COMPATIBILITY_PROFILE AND MOD_REQUIRE_OPENGL_COMPATIBILITY_PROFILE)
                    MESSAGE(FATAL_ERROR "  - Module 'VRN_MODULE_${module_upper}' requires OpenGL compatibility profile. Either disable this module or make it core profile ready.")
                ELSE()
                    # add module availability macro
                    LIST(APPEND VRN_MODULE_DEFINITIONS "-DVRN_MODULE_${module_upper}")

                    # add module definitions
                    LIST(APPEND VRN_MODULE_DEFINITIONS          ${MOD_DEFINITIONS})

                    # add external dependencies
                    LIST(APPEND VRN_MODULE_INCLUDE_DIRECTORIES  ${MOD_INCLUDE_DIRECTORIES})
                    LIST(APPEND VRN_EXTERNAL_LIBRARIES          ${MOD_LIBRARIES})
                    FOREACH(lib ${MOD_DEBUG_LIBRARIES})
                        LIST(APPEND VRN_EXTERNAL_LIBRARIES      debug ${lib})
                    ENDFOREACH()
                    FOREACH(lib ${MOD_RELEASE_LIBRARIES})
                        LIST(APPEND VRN_EXTERNAL_LIBRARIES      optimized ${lib})
                    ENDFOREACH()
                    LIST(APPEND VRN_EXTERNAL_DEBUG_DLLS         ${MOD_DEBUG_DLLS})
                    LIST(APPEND VRN_EXTERNAL_RELEASE_DLLS       ${MOD_RELEASE_DLLS})

                    # add install resources
                    LIST(APPEND VRN_MODULE_INSTALL_DIRECTORIES  ${MOD_INSTALL_DIRECTORIES})
                    LIST(APPEND VRN_MODULE_INSTALL_FILES        ${MOD_INSTALL_FILES})

                    # add core resources
                    IF(MOD_CORE_MODULECLASS)
                        LIST(APPEND VRN_MODULE_CORE_MODULECLASSES ${MOD_CORE_MODULECLASS})
                        STRING(TOLOWER ${MOD_CORE_MODULECLASS} moduleclass_lower)
                        LIST(APPEND VRN_MODULE_CORE_MODULECLASSES_INCLUDES "${module_basedir}/${module_dir}/${moduleclass_lower}.h")
                        LIST(APPEND VRN_MODULE_CORE_SOURCES "${module_basedir}/${module_dir}/${moduleclass_lower}.cpp")
                        LIST(APPEND VRN_MODULE_CORE_HEADERS "${module_basedir}/${module_dir}/${moduleclass_lower}.h")
                    ENDIF()
                    LIST(APPEND VRN_MODULE_CORE_SOURCES         ${MOD_CORE_SOURCES})
                    LIST(APPEND VRN_MODULE_CORE_HEADERS         ${MOD_CORE_HEADERS})
                    LIST(APPEND VRN_MODULE_CORE_APPLICATIONS    ${MOD_CORE_APPLICATIONS})

                    # add qt resources
                    IF(MOD_QT_MODULECLASS)
                        LIST(APPEND VRN_MODULE_QT_MODULECLASSES ${MOD_QT_MODULECLASS})
                        STRING(TOLOWER ${MOD_QT_MODULECLASS} moduleclass_lower)
                        LIST(APPEND VRN_MODULE_QT_MODULECLASSES_INCLUDES "${module_basedir}/${module_dir}/${moduleclass_lower}.h")
                        LIST(APPEND VRN_MODULE_QT_SOURCES        "${module_basedir}/${module_dir}/${moduleclass_lower}.cpp")
                        LIST(APPEND VRN_MODULE_QT_HEADERS_NONMOC "${module_basedir}/${module_dir}/${moduleclass_lower}.h")
                    ENDIF()
                    LIST(APPEND VRN_MODULE_QT_SOURCES           ${MOD_QT_SOURCES})
                    LIST(APPEND VRN_MODULE_QT_HEADERS           ${MOD_QT_HEADERS})
                    LIST(APPEND VRN_MODULE_QT_HEADERS_NONMOC    ${MOD_QT_HEADERS_NONMOC})
                    LIST(APPEND VRN_MODULE_QT_FORMS_HEADERS     ${MOD_QT_FORMS_HEADERS})
                    LIST(APPEND VRN_MODULE_QT_APPLICATIONS      ${MOD_QT_APPLICATIONS})
                    LIST(APPEND VRN_MODULE_QT_RESOURCES         ${MOD_QT_RESOURCES})

                ENDIF(NOT VRN_OPENGL_COMPATIBILITY_PROFILE AND MOD_REQUIRE_OPENGL_COMPATIBILITY_PROFILE)
                    
                UNSET(MOD_REQUIRE_OPENGL_COMPATIBILITY_PROFILE)
                UNSET(MOD_DIR)
                UNSET(MOD_DEFINITIONS)
                UNSET(MOD_INCLUDE_DIRECTORIES)
                UNSET(MOD_LIBRARIES)
                UNSET(MOD_DEBUG_LIBRARIES)
                UNSET(MOD_RELEASE_LIBRARIES)
                UNSET(MOD_DEBUG_DLLS)
                UNSET(MOD_RELEASE_DLLS)
                UNSET(MOD_INSTALL_DIRECTORIES)
                UNSET(MOD_INSTALL_FILES)
                
                UNSET(MOD_CORE_MODULECLASS)
                UNSET(MOD_CORE_SOURCES)
                UNSET(MOD_CORE_HEADERS)
                UNSET(MOD_CORE_APPLICATIONS)
                
                UNSET(MOD_QT_MODULECLASS)
                UNSET(MOD_QT_SOURCES)
                UNSET(MOD_QT_HEADERS)
                UNSET(MOD_QT_HEADERS_NONMOC)
                UNSET(MOD_QT_FORMS_HEADERS)
                UNSET(MOD_QT_APPLICATIONS)
                UNSET(MOD_QT_RESOURCES)
            ELSEIF(NOT DEFINED VRN_MODULE_${module_upper})
                # add missing module include option
                SET(VRN_MODULE_${module_upper} OFF CACHE BOOL "Include module \"${module_lower}\"?")
            ENDIF(VRN_MODULE_${module_upper})
        ENDIF(EXISTS ${module_file}) 
    ENDFOREACH(module_dir ${module_dir_list})
    
ENDFOREACH(module_basedir ${MODULE_BASEDIR_LIST})

message(STATUS "Collected module include directories: ${VRN_MODULE_INCLUDE_DIRECTORIES}")

MESSAGE(STATUS "--------------------------------------------------------------------------------")

# Create commonly used list of all include directories
# FIRST include module directories so that commonly system directories do not shadow specialized
# include directories such as that provided by a specific version of a library ffmpeg.
SET(VRN_INCLUDE_DIRECTORIES ${VRN_MODULE_INCLUDE_DIRECTORIES} ${VRN_COMMON_INCLUDE_DIRECTORIES})

# define framework and module install files (note: DLLs are installed by CMakeLists.txt in root dir)
IF (VRN_ADD_INSTALL_TARGET)

    # framework install files
    if(VRN_BUILD_LIB_VOREENCORE)
        INCLUDE(${VRN_HOME}/cmake/installcore.cmake)
    ENDIF()
    if(VRN_BUILD_LIB_VOREENQT)
        INCLUDE(${VRN_HOME}/cmake/installqt.cmake)
    ENDIF()
    if(VRN_BUILD_VOREENVE)
        INCLUDE(${VRN_HOME}/cmake/installvoreenve.cmake)
    ENDIF()
    if(VRN_BUILD_VOREENBIOLOGY)
        INCLUDE(${VRN_HOME}/cmake/installvoreenbiology.cmake)
    ENDIF()
        
    # module install directories
    FOREACH(install_dir ${VRN_MODULE_INSTALL_DIRECTORIES})
        FILE(RELATIVE_PATH install_dir_rel ${VRN_HOME} ${install_dir})
        INSTALL(DIRECTORY ${install_dir}/ DESTINATION ${install_dir_rel})
    ENDFOREACH()
    
    # module install files
    FOREACH(install_file ${VRN_MODULE_INSTALL_FILES})
        FILE(RELATIVE_PATH install_file_rel ${VRN_HOME} ${install_file})
        GET_FILENAME_COMPONENT(install_path_rel ${install_file_rel} PATH)
        INSTALL(FILES ${install_file} DESTINATION ${install_path_rel})
    ENDFOREACH()
        
    # license files
    FOREACH(install_file ${VRN_EXTERNAL_LICENSE_FILES})
        FILE(RELATIVE_PATH install_file_rel ${VRN_HOME} ${install_file})
        GET_FILENAME_COMPONENT(install_path_rel ${install_file_rel} PATH)
        INSTALL(FILES ${install_file} DESTINATION ${install_path_rel})
    ENDFOREACH()
ENDIF()

SET(COMMONCONF_PROCESSED TRUE)
ENDIF(NOT COMMONCONF_PROCESSED)
