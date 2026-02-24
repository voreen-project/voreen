CMAKE_MINIMUM_REQUIRED(VERSION 3.16.3 FATAL_ERROR)
PROJECT(Voreen)

include(cmake/commonconf.cmake)

################################################################################
# handle target dependencies
################################################################################
IF(VRN_BUILD_VOREENBIOLOGY OR VRN_BUILD_VOREENVE OR VRN_BUILD_VOREENTOOL OR VRN_BUILD_SIMPLEQT OR VRN_MODULE_QT_APPLICATIONS)
    SET(VRN_BUILD_LIB_VOREENQT ON)
ENDIF()

IF(VRN_BUILD_LIB_VOREENQT OR VRN_BUILD_SIMPLEGLUT OR VRN_BUILD_TESTAPPS OR
       VRN_BUILD_BLASTEST OR VRN_BUILD_ITKWRAPPER)
    SET(VRN_BUILD_LIB_VOREENCORE ON)
ENDIF()

IF(VRN_BUILD_LIB_VOREENCORE)
    SET(VRN_BUILD_LIB_TGT ON)
    SET(VRN_BUILD_LIB_EFSW ON)
ENDIF()

################################################################################
# include targets
################################################################################
IF(VRN_BUILD_LIB_TGT)
    ADD_SUBDIRECTORY(ext/tgt)
ENDIF()

IF(VRN_BUILD_LIB_EFSW)
    ADD_SUBDIRECTORY(ext/efsw)
ENDIF()

IF(VRN_BUILD_LIB_VOREENCORE)
    ADD_SUBDIRECTORY(src/core)
ENDIF()

IF(VRN_BUILD_LIB_VOREENQT)
    ADD_SUBDIRECTORY(src/qt)
ENDIF()

IF(VRN_BUILD_VOREENVE)
    ADD_SUBDIRECTORY(apps/voreenve)
ENDIF()

IF(VRN_BUILD_VOREENTOOL)
    ADD_SUBDIRECTORY(apps/voreentool)
ENDIF()

IF(VRN_BUILD_SIMPLEQT)
    ADD_SUBDIRECTORY(apps/simple/simple-qt)
ENDIF()

IF(VRN_BUILD_SIMPLEGLUT)
    ADD_SUBDIRECTORY(apps/simple/simple-glut)
ENDIF()

#-------------------------------------------------------------------------------

IF(VRN_BUILD_TESTAPPS AND EXISTS ${VRN_HOME}/apps/tests)
    ADD_SUBDIRECTORY(apps/tests/descriptiontest)
    ADD_SUBDIRECTORY(apps/tests/networkevaluatortest)
    ADD_SUBDIRECTORY(apps/tests/octreetest)
    ADD_SUBDIRECTORY(apps/tests/processorcreatetest)
    ADD_SUBDIRECTORY(apps/tests/processornetworktest)
    ADD_SUBDIRECTORY(apps/tests/processorinittest)
    ADD_SUBDIRECTORY(apps/tests/serializertest)
    ADD_SUBDIRECTORY(apps/tests/volumeorigintest)
    IF(EXISTS ${VRN_HOME}/apps/tests/regressiontest)
        ADD_SUBDIRECTORY(apps/tests/regressiontest)
    ENDIF()
    IF(VRN_MODULE_OPENCL AND EXISTS ${VRN_HOME}/apps/tests/voreenblastest)
        ADD_SUBDIRECTORY(apps/tests/voreenblastest)
    ENDIF()
ENDIF()

IF(VRN_BUILD_ITKWRAPPER AND EXISTS ${VRN_HOME}/apps/itk_wrapper)
    ADD_SUBDIRECTORY(apps/itk_wrapper)
ENDIF()

IF(VRN_BUILD_VOREENBIOLOGY)
    ADD_SUBDIRECTORY(apps/voreenbiology)
ENDIF()

#-------------------------------------------------------------------------------

# add core module applications
IF(VRN_MODULE_CORE_APPLICATIONS)
    MESSAGE(STATUS "Adding core module applications:")
    FOREACH(core_app ${VRN_MODULE_CORE_APPLICATIONS})
        GET_FILENAME_COMPONENT(core_app_we ${core_app} NAME_WE)
        FILE(RELATIVE_PATH core_app_rel ${VRN_HOME} ${core_app})
        MESSAGE(STATUS "* ${core_app_rel}")
        ADD_EXECUTABLE(${core_app_we} ${core_app})
        ADD_DEFINITIONS(${VRN_DEFINITIONS} ${VRN_MODULE_DEFINITIONS})
        INCLUDE_DIRECTORIES(${VRN_INCLUDE_DIRECTORIES})
        TARGET_LINK_LIBRARIES(${core_app_we} tgt voreen_core ${VRN_EXTERNAL_LIBRARIES} )
    ENDFOREACH()
ENDIF()

# add qt module applications
IF(VRN_MODULE_QT_APPLICATIONS)
    MESSAGE(STATUS "Adding Qt module applications:")
    SET(VRN_QT_COMPONENTS Gui Widgets Svg OpenGLWidgets)
    IF(QT_USE_QTNETWORK)
        LIST(APPEND VRN_QT_COMPONENTS Network)
    ENDIF()
    FIND_PACKAGE(Qt6 ${VRN_REQUIRED_QT_VERSION} COMPONENTS ${VRN_QT_COMPONENTS} REQUIRED)
    # create list of libraries
    set(QT_LIBRARIES Qt6::Core Qt6::Gui Qt6::Widgets Qt6::Svg Qt6::OpenGLWidgets)
    IF(QT_USE_QTNETWORK)
        LIST(APPEND QT_LIBRARIES Qt6::Network)
    ENDIF()

    SET(cur_qt_app_name "")
    SET(cur_qt_app_cpp "")
    SET(cur_qt_app_cpp_rel "")
    SET(cur_qt_app_ui_headers)
    FOREACH(qt_app_file ${VRN_MODULE_QT_APPLICATIONS})
        GET_FILENAME_COMPONENT(qt_app_ext ${qt_app_file} EXT)
        IF(${qt_app_ext} STREQUAL ".cpp") # cpp file => new binary
            IF(NOT ${cur_qt_app_name} STREQUAL "")  # finish previous binary, if present
                MESSAGE(STATUS "* ${cur_qt_app_cpp_rel}")
                ADD_EXECUTABLE(${cur_qt_app_name} ${cur_qt_app_cpp} ${cur_qt_app_ui_headers})
                ADD_DEFINITIONS(${VRN_DEFINITIONS} ${VRN_MODULE_DEFINITIONS})
                INCLUDE_DIRECTORIES(${VRN_INCLUDE_DIRECTORIES})
                TARGET_LINK_LIBRARIES(${cur_qt_app_name} tgt voreen_core voreen_qt ${VRN_EXTERNAL_LIBRARIES} )
                # add Qt dependency
                TARGET_LINK_LIBRARIES(${cur_qt_app_name} ${QT_LIBRARIES})
                SET_TARGET_PROPERTIES(${cur_qt_app_name} PROPERTIES AUTOUIC ON)
            ENDIF()
            
            # init new binary 
            SET(cur_qt_app_cpp ${qt_app_file})
            GET_FILENAME_COMPONENT(cur_qt_app_name ${qt_app_file} NAME_WE)
            FILE(RELATIVE_PATH cur_qt_app_cpp_rel ${VRN_HOME} ${qt_app_file})
            SET(cur_qt_app_ui_headers)
        ELSEIF(${qt_app_ext} STREQUAL ".ui") # ui file => add as ui header to current binary
            LIST(APPEND cur_qt_app_ui_headers ${qt_app_file})
        ENDIF()
    ENDFOREACH()
    
    IF(NOT ${cur_qt_app_name} STREQUAL "") # finish last binary, if present
        MESSAGE(STATUS "* ${cur_qt_app_cpp_rel}")
        ADD_EXECUTABLE(${cur_qt_app_name} ${cur_qt_app_cpp} ${cur_qt_app_ui_headers})
        ADD_DEFINITIONS(${VRN_DEFINITIONS} ${VRN_MODULE_DEFINITIONS})
        INCLUDE_DIRECTORIES(${VRN_INCLUDE_DIRECTORIES})
        TARGET_LINK_LIBRARIES(${cur_qt_app_name} tgt voreen_core voreen_qt ${VRN_EXTERNAL_LIBRARIES} )
        # add Qt dependency
        TARGET_LINK_LIBRARIES(${cur_qt_app_name} ${QT_LIBRARIES})
        SET_TARGET_PROPERTIES(${cur_qt_app_name} PROPERTIES AUTOUIC ON)
     ENDIF()
ENDIF()

################################################################################
# copy/install external DLLs
################################################################################
IF(WIN32)
    LIST(REMOVE_DUPLICATES VRN_EXTERNAL_DEBUG_DLLS)
    LIST(REMOVE_DUPLICATES VRN_EXTERNAL_RELEASE_DLLS)

    IF(VRN_COPY_EXTERNAL_DLLS)
        COPY_EXTERNAL_DLLS(VRN_EXTERNAL_DEBUG_DLLS VRN_EXTERNAL_RELEASE_DLLS false)
    ENDIF()
    
    IF(VRN_ADD_INSTALL_TARGET)
        INSTALL(FILES ${VRN_EXTERNAL_DEBUG_DLLS} 
            DESTINATION .
            CONFIGURATIONS Debug
        )
        INSTALL(FILES ${VRN_EXTERNAL_RELEASE_DLLS}
            DESTINATION .
            CONFIGURATIONS Release
        )
    ENDIF()    
ENDIF()

################################################################################
# Deployment
################################################################################
ADD_SUBDIRECTORY(cmake/post_install)
