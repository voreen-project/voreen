################################################################################
# Core module resources 
################################################################################
SET(MOD_CORE_MODULECLASS WebViewModule)

SET(MOD_CORE_SOURCES
    ${MOD_DIR}/processors/webviewprocessor.cpp
)

SET(MOD_CORE_HEADERS
    ${MOD_DIR}/processors/webviewprocessor.h
)


################################################################################
# Qt module resources 
################################################################################
SET(MOD_QT_MODULECLASS WebViewModuleQt)

SET(QT_USE_QTWEBENGINEWIDGETS TRUE)

SET(MOD_QT_SOURCES
    ${MOD_DIR}/qt/webviewwidget.cpp
    ${MOD_DIR}/qt/webviewprocessorwidgetfactory.cpp
)

SET(MOD_QT_HEADERS
    ${MOD_DIR}/qt/webviewwidget.h
    ${MOD_DIR}/qt/webviewprocessorwidgetfactory.h
)
