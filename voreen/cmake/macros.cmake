
MACRO(LIST_SUBDIRECTORIES Result Directory AbsolutePath)
  FILE(GLOB sub-dirs RELATIVE ${Directory} ${Directory}/*)
  SET(${Result} "")
  FOREACH(d ${sub-dirs})
    IF(IS_DIRECTORY ${Directory}/${d})
        IF(${AbsolutePath})
            LIST(APPEND ${Result} ${Directory}/${d})
        ELSE()
            LIST(APPEND ${Result} ${d})
        ENDIF()
    ENDIF()
  ENDFOREACH()
  LIST(SORT ${Result})
ENDMACRO(LIST_SUBDIRECTORIES)

# copies the passed debug and release DLLs to bin/Debug and bin/Release, resp.
MACRO(COPY_EXTERNAL_DLLS DebugDLLs ReleaseDLLs failOnError)
    MESSAGE(STATUS "Copying external DLLs")
    
    IF(NOT CMAKE_RUNTIME_OUTPUT_DIRECTORY)
        MESSAGE(FATAL_ERROR "CMAKE_RUNTIME_OUTPUT_DIRECTORY not set")
    ENDIF()
    
    SET(debug_dir "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Debug")
    IF(NOT EXISTS ${debug_dir})
        FILE(MAKE_DIRECTORY ${debug_dir})    
    ENDIF()
    FOREACH(dllPath ${${DebugDLLs}})
        IF(EXISTS ${dllPath})
            GET_FILENAME_COMPONENT(dllName ${dllPath} NAME)
            IF(EXISTS ${debug_dir}/${dllName})
                FILE(REMOVE ${debug_dir}/${dllName})
            ENDIF()
            FILE(COPY ${dllPath} DESTINATION ${debug_dir})
        ELSEIF(${failOnError})
            MESSAGE(FATAL_ERROR "Debug DLL not found: ${dllPath}")
        ELSE()
            MESSAGE(WARNING "Debug DLL not found: ${dllPath}")
        ENDIF()
    ENDFOREACH()
    
    SET(release_dir "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Release")
    IF(NOT EXISTS ${release_dir})
        FILE(MAKE_DIRECTORY ${release_dir})    
    ENDIF()
    FOREACH(dllPath ${${ReleaseDLLs}})
        IF(EXISTS ${dllPath})            
            GET_FILENAME_COMPONENT(dllName ${dllPath} NAME)
            IF(EXISTS ${release_dir}/${dllName})
                FILE(REMOVE ${release_dir}/${dllName})
            ENDIF()
            FILE(COPY ${dllPath} DESTINATION ${release_dir})
        ELSEIF(${failOnError})
            MESSAGE(FATAL_ERROR "Release DLL not found: ${dllPath}")
        ELSE()
            MESSAGE(WARNING "Release DLL not found: ${dllPath}")
        ENDIF()
    ENDFOREACH()
ENDMACRO()

# adds custom commands to a target that copy the passed debug and release DLLs 
# to the corresponding output directories (currently not in use)
MACRO(ADD_COPY_DLL_COMMAND Target DebugDLLs ReleaseDLLs)
    #MESSAGE(${${DebugDLLs}})
    MESSAGE(${Target})
    FOREACH(dll ${${DebugDLLs}})
        MESSAGE(${dll})
        ADD_CUSTOM_COMMAND(
            TARGET ${Target}
            COMMAND ${CMAKE_COMMAND} -E copy ${dll} ${VRN_HOME}/bin/Debug
        )
    ENDFOREACH()
    FOREACH(dll ${${ReleaseDLLs}})
        ADD_CUSTOM_COMMAND(
            TARGET ${Target}
            COMMAND ${CMAKE_COMMAND} -E copy ${dll} ${VRN_HOME}/bin/Release
        )
    ENDFOREACH()
ENDMACRO(ADD_COPY_DLL_COMMAND)

# adapted from: http://stackoverflow.com/questions/148570/using-pre-compiled-headers-with-cmake
MACRO(ADD_MSVC_PRECOMPILED_HEADER PrecompiledHeader PrecompiledSource SourcesVar)
  IF(MSVC)
    GET_FILENAME_COMPONENT(PrecompiledBasename ${PrecompiledHeader} NAME_WE)
    SET(PrecompiledBinary "${CMAKE_CURRENT_BINARY_DIR}/${PrecompiledBasename}.pch")
    SET(Sources ${${SourcesVar}})

    SET_SOURCE_FILES_PROPERTIES(${PrecompiledSource}
                                PROPERTIES COMPILE_FLAGS "/Yc\"${PrecompiledHeader}\" /Fp\"${PrecompiledBinary}\""
                                           OBJECT_OUTPUTS "${PrecompiledBinary}")
    SET_SOURCE_FILES_PROPERTIES(${Sources}
                                PROPERTIES COMPILE_FLAGS "/Yu\"${PrecompiledBinary}\" /FI\"${PrecompiledBinary}\" /Fp\"${PrecompiledBinary}\""
                                           OBJECT_DEPENDS "${PrecompiledBinary}")  

    # Add precompiled header to SourcesVar
    LIST(APPEND ${SourcesVar} ${PrecompiledSource})
  ENDIF(MSVC)
ENDMACRO(ADD_MSVC_PRECOMPILED_HEADER)

# adapted from: http://www.mail-archive.com/cmake@cmake.org/msg04394.html
MACRO(ADD_GCC_PRECOMPILED_HEADER _targetName _input )

    GET_FILENAME_COMPONENT(_name ${_input} NAME)
    SET(_source "${CMAKE_CURRENT_SOURCE_DIR}/${_input}")
    SET(_outdir "${CMAKE_CURRENT_BINARY_DIR}/${_name}.gch")
    MAKE_DIRECTORY(${_outdir})
    SET(_output "${_outdir}/${CMAKE_BUILD_TYPE}.c++")
    STRING(TOUPPER "CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}" _flags_var_name)
    SET(_compiler_FLAGS ${${_flags_var_name}})

    LIST(APPEND _compiler_FLAGS ${CMAKE_CXX_FLAGS})

    GET_TARGET_PROPERTY(_type ${_targetName} TYPE)
    IF(${_type} MATCHES SHARED_LIBRARY)
        LIST(APPEND _compiler_FLAGS "-fPIC")
    ENDIF()

    GET_DIRECTORY_PROPERTY(_directory_flags INCLUDE_DIRECTORIES)
    FOREACH(item ${_directory_flags})
    LIST(APPEND _compiler_FLAGS "-I${item}")
    ENDFOREACH(item)

    GET_DIRECTORY_PROPERTY(_directory_flags DEFINITIONS)
    LIST(APPEND _compiler_FLAGS ${_directory_flags})

    # make sure to get potential build-type dependent compile flags (introduced for QT_DEBUG / QT_NO_DEBUG)
    IF(CMAKE_BUILD_TYPE)
        STRING(TOUPPER ${CMAKE_BUILD_TYPE} _upper_build_type)
        GET_DIRECTORY_PROPERTY(_directory_flags "COMPILE_DEFINITIONS_${_upper_build_type}")
        FOREACH(item ${_directory_flags})
            LIST(APPEND _compiler_FLAGS "-D${item}")
        ENDFOREACH(item)
    ELSE()
        GET_DIRECTORY_PROPERTY(_directory_flags COMPILE_DEFINITIONS)
        FOREACH(item ${_directory_flags})
            LIST(APPEND _compiler_FLAGS "-D${item}")
        ENDFOREACH(item)
    ENDIF()

    # copy pch file to build directory, gcc expects it in the same directory
    FILE(COPY ${_input} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

    SEPARATE_ARGUMENTS(_compiler_FLAGS)
    LIST(REMOVE_DUPLICATES _compiler_FLAGS)

    #MESSAGE("_compiler_FLAGS: ${_compiler_FLAGS}")
    #message("${CMAKE_CXX_COMPILER} ${_compiler_FLAGS} -x c++-header -o ${_output} ${_source}")
    ADD_CUSTOM_COMMAND(
        OUTPUT ${_output}
        COMMAND ${CMAKE_CXX_COMPILER}
                                ${_compiler_FLAGS}
                                -x c++-header
                                -o ${_output} ${_source}
        DEPENDS ${_source} )
        ADD_CUSTOM_TARGET(${_targetName}_gch DEPENDS ${_output})
    ADD_DEPENDENCIES(${_targetName} ${_targetName}_gch)
    #SET(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-include ${_name} -Winvalid-pch -H")
    #SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -include ${_name} -Winvalid-pch")

    SET_TARGET_PROPERTIES(${_targetName} PROPERTIES
        COMPILE_FLAGS "-I${CMAKE_CURRENT_BINARY_DIR} -include ${_name} -Winvalid-pch"
    )

ENDMACRO(ADD_GCC_PRECOMPILED_HEADER)

# Assigns the passed Files to source groups according to their file path.
# The path prefixes as well as dot operators are removed from the path.
# E.g.: path/to/my/file.cpp               =>  path\\to\\my
#       <PathPrefix>/path/to/my/file.cpp  =>  path\\to\\my
#       ../../path/to/my/file.cpp         =>  path\\to\\my
MACRO(DEFINE_SOURCE_GROUPS_FROM_SUBDIR Files BasePath RemovePathPrefixes)
    FOREACH(f ${${Files}})
        SET(f_cat "")
        
        # extract relative file path
        IF(IS_ABSOLUTE ${f})
            FILE(RELATIVE_PATH f_rel ${BasePath} ${f})
            GET_FILENAME_COMPONENT(f_rel ${f_rel} PATH)
        ELSE()
            GET_FILENAME_COMPONENT(f_rel ${f} PATH)
        ENDIF()
        
        # create source group specifier from rel path
        IF(f_rel)
            # remove ../
            string(REGEX REPLACE "\\.\\./" "" f_cat ${f_rel})
            # remove specified prefixes
            FOREACH(prefix ${${RemovePathPrefixes}})
                IF(f_cat)
                    string(REGEX REPLACE "${prefix}" "" f_cat ${f_cat})
                ENDIF()
            ENDFOREACH()
            
            # convert path separators into source group separators: 
            # path/to/my => path\\to\\my
            IF(f_cat)
                string(REGEX REPLACE "\\\\|/" "\\\\\\\\" f_cat ${f_cat})
            ENDIF()
        ENDIF()
        
        # set extracted source group
        IF(f_cat)
            SOURCE_GROUP("${f_cat}" FILES ${f})
        ELSE()
            SOURCE_GROUP("" FILES ${f})
        ENDIF()
        
    ENDFOREACH()
ENDMACRO(DEFINE_SOURCE_GROUPS_FROM_SUBDIR)
