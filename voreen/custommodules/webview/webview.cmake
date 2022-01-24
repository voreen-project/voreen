IF(NOT VRN_MODULE_BASE)
    MESSAGE(FATAL_ERROR "WebView Module requires Base Module")
ENDIF()

################################################################################
# Core module resources 
################################################################################
SET(MOD_CORE_MODULECLASS WebViewModule)

SET(MOD_CORE_SOURCES
    ${MOD_DIR}/processors/qtsplitter.cpp
    ${MOD_DIR}/processors/webviewprocessor.cpp
)

SET(MOD_CORE_HEADERS
    ${MOD_DIR}/processors/qtsplitter.h
    ${MOD_DIR}/processors/webviewprocessor.h
)


################################################################################
# Qt module resources 
################################################################################
SET(MOD_QT_MODULECLASS WebViewModuleQt)

SET(QT_USE_QTWEBENGINEWIDGETS TRUE)

SET(MOD_QT_SOURCES
    ${MOD_DIR}/qt/qtsplitterwidget.cpp
    ${MOD_DIR}/qt/webviewprocessorwidgetfactory.cpp
    ${MOD_DIR}/qt/webviewwidget.cpp
)

SET(MOD_QT_HEADERS
    ${MOD_DIR}/qt/qtsplitterwidget.h
    ${MOD_DIR}/qt/webviewprocessorwidgetfactory.h
    ${MOD_DIR}/qt/webviewwidget.h
)
