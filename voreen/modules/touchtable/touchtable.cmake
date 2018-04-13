IF(NOT VRN_MODULE_BASE) #for orientationoverlay
    MESSAGE(FATAL_ERROR "TouchTable Module requires Base Module")
ENDIF()
IF(NOT VRN_MODULE_DEVIL) #for texture handling
    MESSAGE (FATAL_ERROR "TouchTable Module requires Devil Module")
ENDIF()
IF(NOT VRN_MODULE_STAGING) #for touch event simulator
    MESSAGE (FATAL_ERROR "TouchTable Module requires Staging Module")
ENDIF()


SET(MOD_CORE_MODULECLASS TouchtableModule)
SET(MOD_REQUIRE_OPENGL_COMPATIBILITY_PROFILE ON)

IF(UNIX)
    IF(NOT APPLE)
        FIND_LIBRARY(LIBXI Xi)
        IF(LIBXI)
            MESSAGE(STATUS "Found libXi at ${LIBXI}, using it to handle touch events.")
            SET(MOD_LIBRARIES ${LIBXI})
            LIST(APPEND VRN_DEFINITIONS "-DUSE_XINPUT2")
            
            # TODO: suitable position?
            LIST(APPEND VRN_QT_COMPONENTS X11Extras)
            LIST(APPEND QT_LIBRARIES Qt5::X11Extras)
            
        ENDIF()
    ENDIF()
ENDIF()

# module's core source files, path relative to module dir
SET(MOD_CORE_SOURCES
    ${MOD_DIR}/processors/bodyparts3d/bodyparts3dsource.cpp
    ${MOD_DIR}/processors/bodyparts3d/bodyparts3drenderer.cpp
    ${MOD_DIR}/processors/bodyparts3d/bodyparts3dmanager.cpp
    ${MOD_DIR}/processors/touchinteraction/touchtableoverlay.cpp
    ${MOD_DIR}/processors/touchinteraction/widgets/touchtablewidget.cpp
    ${MOD_DIR}/processors/touchinteraction/widgets/touchtableboolwidget.cpp
    ${MOD_DIR}/processors/touchinteraction/widgets/touchtablebodypartswidget.cpp
    ${MOD_DIR}/processors/touchinteraction/widgets/touchtableclippingwidget.cpp	
    ${MOD_DIR}/processors/touchinteraction/widgets/touchtablemenuwidget.cpp
    ${MOD_DIR}/processors/touchinteraction/widgets/touchtablemenuwidgetscrollable.cpp
    ${MOD_DIR}/processors/touchinteraction/widgets/touchtabletransfuncwidget.cpp
    ${MOD_DIR}/processors/touchinteraction/widgets/touchtablecamerawidget.cpp
    ${MOD_DIR}/processors/touchinteraction/widgets/touchtablebuttonwidget.cpp
    ${MOD_DIR}/processors/touchinteraction/widgets/touchtablemovewidget.cpp
    ${MOD_DIR}/processors/touchinteraction/widgets/touchtablevolumewidget.cpp
    ${MOD_DIR}/processors/touchinteraction/widgets/touchtablesnapshotwidget.cpp
    ${MOD_DIR}/processors/touchinteraction/widgets/touchtablelightsourcewidget.cpp
    #${MOD_DIR}/processors/touchinteraction/widgets/touchtableanimationwidget.cpp
    ${MOD_DIR}/processors/touchinteraction/touchtablemenuframe.cpp
    ${MOD_DIR}/processors/touchinteraction/touchtablescrollablemenu.cpp
    ${MOD_DIR}/ports/bodyparts3d/bodypartsport.cpp
    ${MOD_DIR}/ports/bodyparts3d/bodypartstextureport.cpp
    ${MOD_DIR}/processors/touchinteraction/touchtablecontrolelement.cpp
    ${MOD_DIR}/processors/touchinteraction/touchtableslider.cpp	
    ${MOD_DIR}/datastructures/geometry/trianglemeshgeometrysinglecolor.cpp
    ${MOD_DIR}/processors/touchinteraction/touchtablepreviewicon.cpp
    ${MOD_DIR}/processors/touchinteraction/touchtabletimeline.cpp
    ${MOD_DIR}/processors/touchinteraction/touchtablescrollabletimelinemenu.cpp
    ${MOD_DIR}/processors/touchpainter.cpp
)

# module's core header files, path relative to module dir
SET(MOD_CORE_HEADERS
    ${MOD_DIR}/processors/bodyparts3d/bodyparts3dsource.h
    ${MOD_DIR}/processors/bodyparts3d/bodyparts3drenderer.h
    ${MOD_DIR}/processors/bodyparts3d/bodyparts3dmanager.h
    ${MOD_DIR}/processors/touchinteraction/touchtableoverlay.h
    ${MOD_DIR}/processors/touchinteraction/widgets/touchtablewidget.h
    ${MOD_DIR}/processors/touchinteraction/widgets/touchtableboolwidget.h
    ${MOD_DIR}/processors/touchinteraction/widgets/touchtablebodypartswidget.h
    ${MOD_DIR}/processors/touchinteraction/widgets/touchtablemenuwidget.h
    ${MOD_DIR}/processors/touchinteraction/widgets/touchtablemenuwidgetscrollable.h
    ${MOD_DIR}/processors/touchinteraction/widgets/touchtableclippingwidget.h
    ${MOD_DIR}/processors/touchinteraction/widgets/touchtabletransfuncwidget.h	
    ${MOD_DIR}/processors/touchinteraction/widgets/touchtablecamerawidget.h
    ${MOD_DIR}/processors/touchinteraction/widgets/touchtablebuttonwidget.h
    ${MOD_DIR}/processors/touchinteraction/widgets/touchtablemovewidget.h
    ${MOD_DIR}/processors/touchinteraction/widgets/touchtablevolumewidget.h
    ${MOD_DIR}/processors/touchinteraction/widgets/touchtablesnapshotwidget.h
    ${MOD_DIR}/processors/touchinteraction/widgets/touchtablelightsourcewidget.h
    #${MOD_DIR}/processors/touchinteraction/widgets/touchtableanimationwidget.h
    ${MOD_DIR}/processors/touchinteraction/touchtablemenuframe.h
    ${MOD_DIR}/processors/touchinteraction/touchtablescrollablemenu.h
    ${MOD_DIR}/ports/bodyparts3d/bodypartsport.h
    ${MOD_DIR}/ports/bodyparts3d/bodypartstextureport.h
    ${MOD_DIR}/processors/touchinteraction/touchtablecontrolelement.cpp
    ${MOD_DIR}/processors/touchinteraction/touchtablecontrolelement.h
    ${MOD_DIR}/processors/touchinteraction/touchtableslider.h	
    ${MOD_DIR}/datastructures/geometry/trianglemeshgeometrysinglecolor.h
    ${MOD_DIR}/processors/touchinteraction/touchtablepreviewicon.h
    ${MOD_DIR}/processors/touchinteraction/touchtabletimeline.h
    ${MOD_DIR}/processors/touchinteraction/touchtablescrollabletimelinemenu.h
    ${MOD_DIR}/processors/touchpainter.h
)



SET(MOD_QT_MODULECLASS TouchtableModuleQt)
SET(MOD_QT_SOURCES
    ${MOD_DIR}/qt/processor/bodyparts3dwidget.cpp
    ${MOD_DIR}/qt/processor/bodyparts3dprocessorwidgetfactory.cpp
    ${MOD_DIR}/qt/datastructures/bodyparts3ditem.cpp
    ${MOD_DIR}/qt/datastructures/bodyparts3dmodel.cpp
    ${MOD_DIR}/qt/datastructures/bodypartscatalogbase.cpp
    ${MOD_DIR}/qt/datastructures/bodypartscatalog30.cpp
    ${MOD_DIR}/qt/datastructures/bodypartscatalog40.cpp
)  
    
SET(MOD_QT_HEADERS
    ${MOD_DIR}/qt/processor/bodyparts3dwidget.h
    ${MOD_DIR}/qt/datastructures/bodyparts3ditem.h
    ${MOD_DIR}/qt/datastructures/bodyparts3dmodel.h
    ${MOD_DIR}/qt/datastructures/bodypartscatalogbase.h
    ${MOD_DIR}/qt/datastructures/bodypartscatalog30.h
    ${MOD_DIR}/qt/datastructures/bodypartscatalog40.h
)

SET(MOD_QT_HEADERS_NONMOC
    ${MOD_DIR}/qt/processor/bodyparts3dprocessorwidgetfactory.h
)
