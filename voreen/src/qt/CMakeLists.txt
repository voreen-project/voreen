################################################################################
# Project file for the Voreen-Qt library
################################################################################
PROJECT(voreen_qt)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.0 FATAL_ERROR)
INCLUDE(../../cmake/commonconf.cmake)

# Supresses CMake warning when using version > 2.8.10
# For details, see 'cmake --help-policy CMP0020'
IF(WIN32)
    CMAKE_POLICY(SET CMP0020 NEW)
ENDIF()

MESSAGE(STATUS "Configuring voreen_qt Library")

# please insert new files in alphabetical order!
SET(VRN_QT_SOURCES 
    aboutboxbase.cpp 
    helpbrowser.cpp 
    progressdialog.cpp
    versionqt.cpp 
    voreenapplicationqt.cpp
    voreenmoduleqt.cpp

    ## main window ##
    mainwindow/voreenqtmainwindow.cpp
    mainwindow/voreenqtsplashscreen.cpp
    mainwindow/voreenqtstartupwizard.cpp
    
    #menuentity
    mainwindow/menuentities/voreenqtmenuentity.cpp
    mainwindow/menuentities/animationmenuentity.cpp
    mainwindow/menuentities/applicationmodeconfigmenuentity.cpp
    mainwindow/menuentities/consolemenuentity.cpp
    mainwindow/menuentities/generalsettingsmenuentity.cpp
    mainwindow/menuentities/inputmappingmenuentity.cpp
    mainwindow/menuentities/networkmodesettingsmenuentity.cpp
    mainwindow/menuentities/newworkspacemenuentity.cpp
    mainwindow/menuentities/openworkspacemenuentity.cpp
    mainwindow/menuentities/processorlistmenuentity.cpp
    mainwindow/menuentities/propertylistmenuentity.cpp
    mainwindow/menuentities/rebuildshadermenuentity.cpp
    mainwindow/menuentities/rendertargetviewermenuentity/rendertargetviewer.cpp
    mainwindow/menuentities/rendertargetviewermenuentity/rendertargetviewermenuentity.cpp
    mainwindow/menuentities/rendertargetviewermenuentity/rendertargetviewerpainter.cpp
    mainwindow/menuentities/saveworkspaceasmenuentity.cpp
    mainwindow/menuentities/saveworkspacemenuentity.cpp
    mainwindow/menuentities/screenshotmenuentity/canvasrenderscreenshotplugin.cpp
    mainwindow/menuentities/screenshotmenuentity/networkscreenshotplugin.cpp
    mainwindow/menuentities/screenshotmenuentity/screenshotmenuentity.cpp
    mainwindow/menuentities/screenshotmenuentity/screenshotpluginbase.cpp
    mainwindow/menuentities/volumeviewermenuentity.cpp
    mainwindow/menuentities/workspacedescriptionmenuentity.cpp
    
    ## utils ##
    utils/voreenqtworkspacehandler.cpp
    
    ## network editor ##

    #editor
    networkeditor/networkeditor.cpp
    
    #dialogs
    networkeditor/dialogs/propertylinkdialog.cpp
    
    #utils
    networkeditor/dialogs/utils/splitgraphicsview.cpp

    #style
    networkeditor/styles/nwestyle_base.cpp
    networkeditor/styles/nwestyle_classic.cpp
    networkeditor/styles/nwestyle_classic_print.cpp
    networkeditor/styles/nwestyle_material.cpp

    #graphlayouts
    networkeditor/graphlayouts/nwegl_base.cpp
    networkeditor/graphlayouts/nwegl_sugiyama.cpp
    
    #gi
    networkeditor/graphicitems/nwebasegraphicsitem.cpp

    #core
    networkeditor/graphicitems/core/portgraphicsitem.cpp
    networkeditor/graphicitems/core/portownergraphicsitem.cpp
    networkeditor/graphicitems/core/processorgraphicsitem.cpp
    networkeditor/graphicitems/core/propertygraphicsitem.cpp
    networkeditor/graphicitems/core/propertylistgraphicsitem.cpp    
    networkeditor/graphicitems/core/propertyownergraphicsitem.cpp

    #connections
    networkeditor/graphicitems/connections/connectionbasegraphicsitem.cpp
    networkeditor/graphicitems/connections/portownerlinkarrowgraphicsitem.cpp
    networkeditor/graphicitems/connections/propertylinkarrowgraphicsitem.cpp
    networkeditor/graphicitems/connections/portarrowgraphicsitem.cpp
    networkeditor/graphicitems/connections/portsizelinkarrowgraphicsitem.cpp

    #textboxes
    networkeditor/graphicitems/textboxes/frameboxgraphicsitem.cpp
    networkeditor/graphicitems/textboxes/textboxbasegraphicsitem.cpp
    networkeditor/graphicitems/textboxes/textboxgraphicsitem.cpp

    #tooltips
    networkeditor/graphicitems/tooltips/tooltipbasegraphicsitem.cpp
    networkeditor/graphicitems/tooltips/tooltipprocessorgraphicsitem.cpp
    networkeditor/graphicitems/tooltips/tooltipportgraphicsitem.cpp

    #utils
    networkeditor/graphicitems/utils/renamabletextgraphicsitem.cpp
    networkeditor/graphicitems/utils/progressbargraphicsitem.cpp
    networkeditor/graphicitems/utils/widgettogglebuttongraphicsitem.cpp
    networkeditor/graphicitems/utils/propertylistbuttongraphicsitem.cpp
    
    #meta
    networkeditor/meta/textboxmetadata.cpp
    networkeditor/meta/zoommetadata.cpp
    
    #animation        
    animation/animationexportwidget.cpp
    animation/keyframegraphicsitem.cpp
    animation/animationeditor.cpp
    animation/currentframegraphicsitem.cpp
    animation/animationinfowidgetbase.cpp
    animation/animationintervalwidget.cpp
    animation/animationkeyframewidget.cpp
    animation/animationpropertylistwidget.cpp
    animation/animationoverviewwidget.cpp
    animation/animationtimelinewidget.cpp
    
    ## widgets ##
    widgets/applicationmodeconfigdialog.cpp
    widgets/choicelistcombobox.cpp
    widgets/colorselectorwidget.cpp
    widgets/clickablelabel.cpp
    widgets/codeedit.cpp
    widgets/consoleplugin.cpp
    widgets/custombutton.cpp
    widgets/customlabel.cpp
    widgets/enterexitpushbutton.cpp
    widgets/eventpropertywidget.cpp
    widgets/expandableheaderbutton.cpp
    widgets/glslhighlighter.cpp
    widgets/inputmappingdialog.cpp
    widgets/keydetectorwidget.cpp
    widgets/lightwidget.cpp
    widgets/lineeditresetwidget.cpp
    widgets/networkmodeconfigdialog.cpp
    widgets/processorlistwidget.cpp
    widgets/propertylistwidget.cpp
    widgets/propertyselectiondialog.cpp
    widgets/rawvolumewidget.cpp
    widgets/shaderplugin.cpp
    widgets/sliderspinboxwidget.cpp
    widgets/stringselectiondialog.cpp
    widgets/syntaxhighlighter.cpp
    widgets/volumeiohelper.cpp
    widgets/volumelistingdialog.cpp
    widgets/volumereaderselectiondialog.cpp
    widgets/volumeviewer.cpp
    widgets/volumeviewhelper.cpp
    widgets/voreensettingsdialog.cpp
    widgets/voreentoolwindow.cpp
    
    widgets/processor/qprocessorwidget.cpp

    widgets/property/appmodepropertygroupwidget.cpp
    widgets/property/boolpropertywidget.cpp
    widgets/property/floatboundingboxpropertywidget.cpp
    widgets/property/intboundingboxpropertywidget.cpp
    widgets/property/buttonpropertywidget.cpp
    widgets/property/camerapropertywidget.cpp
    widgets/property/camerawidget.cpp
    widgets/property/colorpropertywidget.cpp
    widgets/property/corepropertywidgetfactory.cpp
    widgets/property/filedialogpropertywidget.cpp
    widgets/property/floatmat2propertywidget.cpp
    widgets/property/floatmat3propertywidget.cpp
    widgets/property/floatmat4propertywidget.cpp
    widgets/property/floatpropertywidget.cpp
    widgets/property/floatvec2propertywidget.cpp
    widgets/property/floatvec3propertywidget.cpp
    widgets/property/floatvec4propertywidget.cpp
    widgets/property/fontpropertywidget.cpp
    widgets/property/intpropertywidget.cpp
    widgets/property/intvec2propertywidget.cpp
    widgets/property/intvec3propertywidget.cpp
    widgets/property/intvec4propertywidget.cpp
    widgets/property/lightpropertywidget.cpp
    widgets/property/matrixpropertywidget.cpp
    widgets/property/optionpropertywidget.cpp
    widgets/property/planepropertywidget.cpp
    widgets/property/processorpropertieswidget.cpp
    widgets/property/progresspropertywidget.cpp
    widgets/property/propertyownerwidget.cpp
    widgets/property/propertyvectorwidget.cpp
    widgets/property/qpropertywidget.cpp
    widgets/property/qpropertywidgetwithtoolwindow.cpp
    widgets/property/shaderpropertywidget.cpp
    widgets/property/stringexpressionpropertywidget.cpp
    widgets/property/stringpropertywidget.cpp
    widgets/property/temppathpropertywidget.cpp
    widgets/property/grouppropertywidget.cpp
    widgets/property/vecpropertywidget.cpp
    widgets/property/voxeltypepropertywidget.cpp
    widgets/property/volumeinfopropertywidget.cpp
    widgets/property/volumeurlpropertywidget.cpp
    widgets/property/volumeurllistpropertywidget.cpp
    widgets/property/numeric/intintervalpropertywidget.cpp
    widgets/property/numeric/floatintervalpropertywidget.cpp

    #color
    widgets/property/color/colorswitchpropertywidget.cpp

    #string
    widgets/property/string/stringtablepropertywidget.cpp
    
    #trunsfunc
    widgets/property/transfunc/transfuncpropertyeditorbase.cpp
    widgets/property/transfunc/transfuncpropertywidgetbase.cpp
    widgets/property/transfunc/transfuncpropertywidgetpainterbase.cpp

    widgets/property/transfunc/1d/transfunc1dpropertyeditor.cpp
    widgets/property/transfunc/1d/transfunc1dpropertywidget.cpp
    widgets/property/transfunc/1d/transfunc1dpropertywidgetpainter.cpp
	widgets/property/transfunc/1d/1dkeys/transfunc1dkeyspropertyeditor.cpp
    widgets/property/transfunc/1d/1dkeys/transfunc1dkeyspropertyeditorcanvas.cpp
    widgets/property/transfunc/1d/1dkeys/transfunc1dkeyspropertywidget.cpp
    widgets/property/transfunc/1d/1dkeys/transfunc1dkeyspropertywidgetpainter.cpp
    widgets/property/transfunc/1d/1dgaussian/transfunc1dgaussianpropertyeditor.cpp
    widgets/property/transfunc/1d/1dgaussian/transfunc1dgaussianpropertyeditorcanvas.cpp
    widgets/property/transfunc/1d/1dgaussian/transfunc1dgaussianpropertywidget.cpp
    widgets/property/transfunc/1d/1dgaussian/transfunc1dgaussianpropertywidgetpainter.cpp

    widgets/property/transfunc/2d/transfunc2dpropertyeditor.cpp
    widgets/property/transfunc/2d/transfunc2dpropertywidget.cpp
    widgets/property/transfunc/2d/transfunc2dpropertywidgetpainter.cpp
    widgets/property/transfunc/2d/2dprimitives/transfunc2dprimitivespropertyeditor.cpp
    widgets/property/transfunc/2d/2dprimitives/transfunc2dprimitivespropertyeditorcanvas.cpp
    widgets/property/transfunc/2d/2dprimitives/transfunc2dprimitivespropertywidget.cpp
    widgets/property/transfunc/2d/2dprimitives/transfunc2dprimitivespropertywidgetpainter.cpp

    widgets/property/transfunc/utils/alphapicker.cpp
    widgets/property/transfunc/utils/colorluminancepicker.cpp
    widgets/property/transfunc/utils/colorpicker.cpp
    widgets/property/transfunc/utils/doubleslider.cpp
    widgets/property/transfunc/utils/histogrampainter.cpp
    widgets/property/transfunc/utils/quadrupleslider.cpp
    widgets/property/transfunc/utils/transfunciohelperqt.cpp
    widgets/property/transfunc/utils/transfuncpropertywidgetpainterslider.cpp

    #Qt implementations of tgt classes
    ../../ext/tgt/qt/qtcanvas.cpp
    ../../ext/tgt/qt/qtmainglcontext.cpp
    ../../ext/tgt/qt/qttimer.cpp
)

# headers
SET(VRN_QT_HEADERS
    ../../include/voreen/qt/aboutboxbase.h
    ../../include/voreen/qt/helpbrowser.h
    
    ## main window ##
    ../../include/voreen/qt/mainwindow/voreenqtmainwindow.h
    ../../include/voreen/qt/mainwindow/voreenqtsplashscreen.h
    ../../include/voreen/qt/mainwindow/voreenqtstartupwizard.h

    #menuentity
    ../../include/voreen/qt/mainwindow/menuentities/voreenqtmenuentity.h
    ../../include/voreen/qt/mainwindow/menuentities/animationmenuentity.h
    ../../include/voreen/qt/mainwindow/menuentities/applicationmodeconfigmenuentity.h
    ../../include/voreen/qt/mainwindow/menuentities/consolemenuentity.h
    ../../include/voreen/qt/mainwindow/menuentities/generalsettingsmenuentity.h
    ../../include/voreen/qt/mainwindow/menuentities/inputmappingmenuentity.h
    ../../include/voreen/qt/mainwindow/menuentities/networkmodesettingsmenuentity.h
    ../../include/voreen/qt/mainwindow/menuentities/newworkspacemenuentity.h
    ../../include/voreen/qt/mainwindow/menuentities/openworkspacemenuentity.h
    ../../include/voreen/qt/mainwindow/menuentities/processorlistmenuentity.h
    ../../include/voreen/qt/mainwindow/menuentities/propertylistmenuentity.h
    ../../include/voreen/qt/mainwindow/menuentities/rebuildshadermenuentity.h
    ../../include/voreen/qt/mainwindow/menuentities/rendertargetviewermenuentity/rendertargetviewer.h
    ../../include/voreen/qt/mainwindow/menuentities/rendertargetviewermenuentity/rendertargetviewermenuentity.h
    ../../include/voreen/qt/mainwindow/menuentities/rendertargetviewermenuentity/rendertargetviewerpainter.h
     ../../include/voreen/qt/mainwindow/menuentities/saveworkspaceasmenuentity.h
    ../../include/voreen/qt/mainwindow/menuentities/saveworkspacemenuentity.h
    ../../include/voreen/qt/mainwindow/menuentities/screenshotmenuentity/canvasrendererscreenshotplugin.h
    ../../include/voreen/qt/mainwindow/menuentities/screenshotmenuentity/networkscreenshotplugin.h
    ../../include/voreen/qt/mainwindow/menuentities/screenshotmenuentity/screenshotmenuentity.h
    ../../include/voreen/qt/mainwindow/menuentities/screenshotmenuentity/screenshotpluginbase.h
    ../../include/voreen/qt/mainwindow/menuentities/volumeviewermenuentity.h
    ../../include/voreen/qt/mainwindow/menuentities/workspacedescriptionmenuentity.h

    ## utils ##
    ../../include/voreen/qt/utils/voreenqtworkspacehandler.h
    
    ## network editor ##

    #networkeditor
    ../../include/voreen/qt/networkeditor/networkeditor.h
    
    #dialogs
    ../../include/voreen/qt/networkeditor/dialogs/propertylinkdialog.h
    
    #utils
    ../../include/voreen/qt/networkeditor/dialogs/utils/splitgraphicsview.h

    #gi
    ../../include/voreen/qt/networkeditor/graphicitems/nwebasegraphicsitem.h

    #core
    ../../include/voreen/qt/networkeditor/graphicitems/core/portgraphicsitem.h
    ../../include/voreen/qt/networkeditor/graphicitems/core/portownergraphicsitem.h
    ../../include/voreen/qt/networkeditor/graphicitems/core/processorgraphicsitem.h
    ../../include/voreen/qt/networkeditor/graphicitems/core/propertygraphicsitem.h
    ../../include/voreen/qt/networkeditor/graphicitems/core/propertylistgraphicsitem.h
    ../../include/voreen/qt/networkeditor/graphicitems/core/propertyownergraphicsitem.h

    #connections
    ../../include/voreen/qt/networkeditor/graphicitems/connections/propertylinkarrowgraphicsitem.h
    
    #textboxes
    ../../include/voreen/qt/networkeditor/graphicitems/textboxes/frameboxgraphicsitem.h
    ../../include/voreen/qt/networkeditor/graphicitems/textboxes/textboxbasegraphicsitem.h
    ../../include/voreen/qt/networkeditor/graphicitems/textboxes/textboxgraphicsitem.h
 
    #tooltips
    ../../include/voreen/qt/networkeditor/graphicitems/tooltips/tooltipbasegraphicsitem.h

    #utils
    ../../include/voreen/qt/networkeditor/graphicitems/utils/renamabletextgraphicsitem.h
    ../../include/voreen/qt/networkeditor/graphicitems/utils/propertylistbuttongraphicsitem.h
    ../../include/voreen/qt/networkeditor/graphicitems/utils/widgettogglebuttongraphicsitem.h
    
    #meta
    ../../include/voreen/qt/networkeditor/meta/textboxmetadata.h
    ../../include/voreen/qt/networkeditor/meta/zoommetadata.h
    
    #animation       
    ../../include/voreen/qt/animation/animationexportwidget.h
    ../../include/voreen/qt/animation/keyframegraphicsitem.h
    ../../include/voreen/qt/animation/animationeditor.h
    ../../include/voreen/qt/animation/animationinfowidgetbase.h
    ../../include/voreen/qt/animation/animationintervalwidget.h
    ../../include/voreen/qt/animation/animationkeyframewidget.h
    ../../include/voreen/qt/animation/animationpropertylistwidget.h
    ../../include/voreen/qt/animation/animationoverviewwidget.h
    ../../include/voreen/qt/animation/animationtimelinewidget.h

    ## widgets ##
    ../../include/voreen/qt/widgets/applicationmodeconfigdialog.h
    ../../include/voreen/qt/widgets/clickablelabel.h
    ../../include/voreen/qt/widgets/colorselectorwidget.h
    ../../include/voreen/qt/widgets/codeedit.h
    ../../include/voreen/qt/widgets/consoleplugin.h
    ../../include/voreen/qt/widgets/custombutton.h
    ../../include/voreen/qt/widgets/customlabel.h
    ../../include/voreen/qt/widgets/enterexitpushbutton.h
    ../../include/voreen/qt/widgets/eventpropertywidget.h
    ../../include/voreen/qt/widgets/expandableheaderbutton.h
    ../../include/voreen/qt/widgets/inputmappingdialog.h
    ../../include/voreen/qt/widgets/keydetectorwidget.h
    ../../include/voreen/qt/widgets/lightwidget.h
    ../../include/voreen/qt/widgets/lineeditresetwidget.h
    ../../include/voreen/qt/widgets/networkmodeconfigdialog.h
    ../../include/voreen/qt/widgets/processorlistwidget.h
    ../../include/voreen/qt/widgets/propertylistwidget.h
    ../../include/voreen/qt/widgets/propertyselectiondialog.h
    ../../include/voreen/qt/widgets/rawvolumewidget.h
    ../../include/voreen/qt/widgets/shaderplugin.h
    ../../include/voreen/qt/widgets/sliderspinboxwidget.h
    ../../include/voreen/qt/widgets/stringselectiondialog.h
    ../../include/voreen/qt/widgets/volumeiohelper.h
    ../../include/voreen/qt/widgets/volumelistingdialog.h
    ../../include/voreen/qt/widgets/volumereaderselectiondialog.h
    ../../include/voreen/qt/widgets/volumeviewer.h
    ../../include/voreen/qt/widgets/voreensettingsdialog.h
    ../../include/voreen/qt/widgets/voreentoolwindow.h
    
    ../../include/voreen/qt/widgets/processor/qprocessorwidget.h

    ../../include/voreen/qt/widgets/property/appmodepropertygroupwidget.h
    ../../include/voreen/qt/widgets/property/boolpropertywidget.h
    ../../include/voreen/qt/widgets/property/floatboundingboxpropertywidget.h
    ../../include/voreen/qt/widgets/property/intboundingboxpropertywidget.h
    ../../include/voreen/qt/widgets/property/buttonpropertywidget.h
    ../../include/voreen/qt/widgets/property/camerapropertywidget.h
    ../../include/voreen/qt/widgets/property/camerawidget.h
    ../../include/voreen/qt/widgets/property/colorpropertywidget.h
    ../../include/voreen/qt/widgets/property/filedialogpropertywidget.h
    ../../include/voreen/qt/widgets/property/floatmat2propertywidget.h
    ../../include/voreen/qt/widgets/property/floatmat3propertywidget.h
    ../../include/voreen/qt/widgets/property/floatmat4propertywidget.h
    ../../include/voreen/qt/widgets/property/floatpropertywidget.h
    ../../include/voreen/qt/widgets/property/floatvec2propertywidget.h
    ../../include/voreen/qt/widgets/property/floatvec3propertywidget.h
    ../../include/voreen/qt/widgets/property/floatvec4propertywidget.h
    ../../include/voreen/qt/widgets/property/fontpropertywidget.h
    ../../include/voreen/qt/widgets/property/intpropertywidget.h
    ../../include/voreen/qt/widgets/property/intvec2propertywidget.h
    ../../include/voreen/qt/widgets/property/intvec3propertywidget.h
    ../../include/voreen/qt/widgets/property/intvec4propertywidget.h
    ../../include/voreen/qt/widgets/property/lightpropertywidget.h
    ../../include/voreen/qt/widgets/property/matrixpropertywidget.h
    ../../include/voreen/qt/widgets/property/optionpropertywidget.h
    ../../include/voreen/qt/widgets/property/planepropertywidget.h
    ../../include/voreen/qt/widgets/property/processorpropertieswidget.h
    ../../include/voreen/qt/widgets/property/progresspropertywidget.h
    ../../include/voreen/qt/widgets/property/propertyownerwidget.h
    ../../include/voreen/qt/widgets/property/propertyvectorwidget.h
    ../../include/voreen/qt/widgets/property/qpropertywidget.h
    ../../include/voreen/qt/widgets/property/shaderpropertywidget.h
    ../../include/voreen/qt/widgets/property/stringexpressionpropertywidget.h
    ../../include/voreen/qt/widgets/property/stringpropertywidget.h
    ../../include/voreen/qt/widgets/property/temppathpropertywidget.h
    ../../include/voreen/qt/widgets/property/grouppropertywidget.h
    ../../include/voreen/qt/widgets/property/voxeltypepropertywidget.h
    ../../include/voreen/qt/widgets/property/volumeinfopropertywidget.h
    ../../include/voreen/qt/widgets/property/volumeurllistpropertywidget.h
    ../../include/voreen/qt/widgets/property/volumeurlpropertywidget.h
    ../../include/voreen/qt/widgets/property/numeric/intintervalpropertywidget.h
    ../../include/voreen/qt/widgets/property/numeric/floatintervalpropertywidget.h

        #color
    ../../include/voreen/qt/widgets/property/color/colorswitchpropertywidget.h
        #string
    ../../include/voreen/qt/widgets/property/string/stringtablepropertywidget.h
        #trunsfunc
    ../../include/voreen/qt/widgets/property/transfunc/transfuncpropertyeditorbase.h
    ../../include/voreen/qt/widgets/property/transfunc/transfuncpropertywidgetbase.h
    ../../include/voreen/qt/widgets/property/transfunc/transfuncpropertywidgetpainterbase.h

    ../../include/voreen/qt/widgets/property/transfunc/1d/transfunc1dpropertyeditor.h
    ../../include/voreen/qt/widgets/property/transfunc/1d/transfunc1dpropertywidget.h
    ../../include/voreen/qt/widgets/property/transfunc/1d/transfunc1dpropertywidgetpainter.h
    ../../include/voreen/qt/widgets/property/transfunc/1d/1dkeys/transfunc1dkeyspropertyeditor.h
    ../../include/voreen/qt/widgets/property/transfunc/1d/1dkeys/transfunc1dkeyspropertyeditorcanvas.h
    ../../include/voreen/qt/widgets/property/transfunc/1d/1dkeys/transfunc1dkeyspropertywidget.h
    ../../include/voreen/qt/widgets/property/transfunc/1d/1dkeys/transfunc1dkeyspropertywidgetpainter.h
	../../include/voreen/qt/widgets/property/transfunc/1d/1dgaussian/transfunc1dgaussianpropertyeditor.h
    ../../include/voreen/qt/widgets/property/transfunc/1d/1dgaussian/transfunc1dgaussianpropertyeditorcanvas.h
    ../../include/voreen/qt/widgets/property/transfunc/1d/1dgaussian/transfunc1dgaussianpropertywidget.h
    ../../include/voreen/qt/widgets/property/transfunc/1d/1dgaussian/transfunc1dgaussianpropertywidgetpainter.h

    ../../include/voreen/qt/widgets/property/transfunc/2d/transfunc2dpropertyeditor.h
    ../../include/voreen/qt/widgets/property/transfunc/2d/transfunc2dpropertywidget.h
    ../../include/voreen/qt/widgets/property/transfunc/2d/transfunc2dpropertywidgetpainter.h
    ../../include/voreen/qt/widgets/property/transfunc/2d/2dprimitives/transfunc2dprimitivespropertyeditor.h
    ../../include/voreen/qt/widgets/property/transfunc/2d/2dprimitives/transfunc2dprimitivespropertyeditorcanvas.h
    ../../include/voreen/qt/widgets/property/transfunc/2d/2dprimitives/transfunc2dprimitivespropertywidget.h
    ../../include/voreen/qt/widgets/property/transfunc/2d/2dprimitives/transfunc2dprimitivespropertywidgetpainter.h
    
    ../../include/voreen/qt/widgets/property/transfunc/utils/alphapicker.h
    ../../include/voreen/qt/widgets/property/transfunc/utils/colorluminancepicker.h
    ../../include/voreen/qt/widgets/property/transfunc/utils/colorpicker.h
    ../../include/voreen/qt/widgets/property/transfunc/utils/doubleslider.h
    ../../include/voreen/qt/widgets/property/transfunc/utils/histogrampainter.h
    ../../include/voreen/qt/widgets/property/transfunc/utils/quadrupleslider.h
    ../../include/voreen/qt/widgets/property/transfunc/utils/transfunciohelperqt.h
    ../../include/voreen/qt/widgets/property/transfunc/utils/transfuncpropertywidgetpainterslider.h    
)

# headers that do not need to run through moc
SET(VRN_QT_HEADERS_NONMOC
    ../../include/voreen/qt/voreenqtapi.h
    ../../include/voreen/qt/voreenmoduleqt.h
    ../../include/voreen/qt/voreenapplicationqt.h
    ../../include/voreen/qt/progressdialog.h
    ../../include/voreen/qt/versionqt.h
    
    ## network editor ##

    #connections
    ../../include/voreen/qt/networkeditor/graphicitems/connections/connectionbasegraphicsitem.h
    ../../include/voreen/qt/networkeditor/graphicitems/connections/portownerlinkarrowgraphicsitem.h
    ../../include/voreen/qt/networkeditor/graphicitems/connections/portarrowgraphicsitem.h
    ../../include/voreen/qt/networkeditor/graphicitems/connections/portsizelinkarrowgraphicsitem.h 
    
    #tooltips
    ../../include/voreen/qt/networkeditor/graphicitems/tooltips/tooltipportgraphicsitem.h
    ../../include/voreen/qt/networkeditor/graphicitems/tooltips/tooltipprocessorgraphicsitem.h    
    
    #settings
    ../../include/voreen/qt/networkeditor/editor_settings.h

    #utils
    ../../include/voreen/qt/networkeditor/graphicitems/utils/progressbargraphicsitem.h
    
    #style
    ../../include/voreen/qt/networkeditor/styles/nwestyle_base.h
    ../../include/voreen/qt/networkeditor/styles/nwestyle_classic.h
    ../../include/voreen/qt/networkeditor/styles/nwestyle_classic_print.h
    ../../include/voreen/qt/networkeditor/styles/nwestyle_material.h

    #graphlayouts
    ../../include/voreen/qt/networkeditor/graphlayouts/nwegl_base.h
    ../../include/voreen/qt/networkeditor/graphlayouts/nwegl_sugiyama.h
    
    #animation    
    ../../include/voreen/qt/animation/currentframegraphicsitem.h

    ## widgets ##
    
    ../../include/voreen/qt/widgets/choicelistcombobox.h
    ../../include/voreen/qt/widgets/glslhighlighter.h
    ../../include/voreen/qt/widgets/syntaxhighlighter.h
    ../../include/voreen/qt/widgets/volumeviewhelper.h    
    
    ../../include/voreen/qt/widgets/property/corepropertywidgetfactory.h
    ../../include/voreen/qt/widgets/property/qpropertywidgetwithtoolwindow.h
    ../../include/voreen/qt/widgets/property/vecpropertywidget.h

    #Qt implementations of tgt classes
    ../../ext/tgt/qt/qtcanvas.h
    ../../ext/tgt/qt/qtmainglcontext.h
    ../../ext/tgt/qt/qttimer.h
)

# TODO: still necessary?
# we need to link against cocoa
IF(APPLE)
    FIND_LIBRARY(COCOA_LIB Cocoa)
    LIST(APPEND VRN_EXTERNAL_LIBRARIES ${COCOA_LIB})
ENDIF()

# add module sources and headers
LIST(APPEND VRN_QT_SOURCES ${VRN_MODULE_QT_SOURCES})
LIST(APPEND VRN_QT_HEADERS ${VRN_MODULE_QT_HEADERS})
LIST(APPEND VRN_QT_HEADERS_NONMOC ${VRN_MODULE_QT_HEADERS_NONMOC})

# QMAKE_CHANGED must be set to true, to allow changing Qt versions by changing the qmake.exe in QT_QMAKE_EXECUTABLE
SET(QT_QMAKE_CHANGED TRUE)
# find Qt library
SET(VRN_QT_COMPONENTS Core Gui Widgets Svg)
IF(QT_USE_QTNETWORK)
    LIST(APPEND VRN_QT_COMPONENTS Network)
ENDIF()
FIND_PACKAGE(Qt5 COMPONENTS ${VRN_QT_COMPONENTS} REQUIRED)
# create list of libraries
set(QT_LIBRARIES Qt5::Core Qt5::Gui Qt5::Widgets Qt5::Svg ${QT5SVG_LIBRARIES})
set(QT_INCLUDE_DIRECTORIES ${Qt5Core_INCLUDES} ${Qt5Gui_INCLUDES} ${Qt5Widgets_INCLUDES} ${Qt5Svg_INCLUDES})
IF(QT_USE_QTNETWORK)
    LIST(APPEND QT_LIBRARIES Qt5::Network)
    LIST(APPEND QT_INCLUDE_DIRECTORIES ${Qt5Network_INCLUDES})
ENDIF()

################################################################################
# generate module registration header
################################################################################
SET(REGISTRATION_FILE "${CMAKE_BINARY_DIR}/gen_moduleregistration_qt.h" )
MESSAGE(STATUS "* Generating module registration header ${REGISTRATION_FILE}")

SET(REGISTRATION_SOURCE "// WARNING: This file is generated and will be overwritten on next cmake run!\n" )
LIST(APPEND REGISTRATION_SOURCE "#include \"voreen/qt/voreenapplicationqt.h\"\n\n" )
LIST(APPEND REGISTRATION_SOURCE "// include qt module class headers\n" )
FOREACH(inc ${VRN_MODULE_QT_MODULECLASSES_INCLUDES})
    LIST(APPEND REGISTRATION_SOURCE "#include \"${inc}\"\n" )
ENDFOREACH()
LIST(APPEND REGISTRATION_SOURCE "\nnamespace voreen {\n" )
LIST(APPEND REGISTRATION_SOURCE "\n// instantiate qt module classes\n" )
LIST(APPEND REGISTRATION_SOURCE "void registerAllQtModules(VoreenApplicationQt* vappQt) {\n" )
IF(VRN_MODULE_QT_MODULECLASSES)
    LIST(APPEND REGISTRATION_SOURCE "    VoreenModuleQt* module\;\n\n")
    LIST(LENGTH VRN_MODULE_QT_MODULECLASSES num_modules)
    MATH(EXPR max_index "${num_modules} - 1")
    FOREACH(i RANGE ${max_index})
        LIST(GET VRN_MODULE_QT_MODULECLASSES ${i} class)
        LIST(GET VRN_MODULE_QT_MODULECLASSES_INCLUDES ${i} inc)
        FILE(RELATIVE_PATH inc_rel ${VRN_HOME} ${inc})
        GET_FILENAME_COMPONENT(mod_path_rel ${inc_rel} PATH)
        LIST(APPEND REGISTRATION_SOURCE "    module = new ${class}(\"${mod_path_rel}\")\;\n")
        LIST(APPEND REGISTRATION_SOURCE "    vappQt->registerQtModule(module)\;\n\n")
    ENDFOREACH()
ENDIF()
LIST(APPEND REGISTRATION_SOURCE "}\n\n" )
LIST(APPEND REGISTRATION_SOURCE "} // namespace\n" )

FILE(WRITE ${REGISTRATION_FILE} ${REGISTRATION_SOURCE})
LIST(APPEND VRN_QT_HEADERS_NONMOC ${REGISTRATION_FILE})


################################################################################
# define library target
################################################################################

# Qt moc
QT5_WRAP_CPP(VRN_QT_MOC_SOURCES ${VRN_QT_HEADERS})
LIST(APPEND VRN_QT_SOURCES ${VRN_QT_MOC_SOURCES})

# Qt ui and resources
QT5_ADD_RESOURCES(VRN_QT_RESOURCES_RCC vrn_qt.qrc ${VRN_MODULE_QT_RESOURCES})

# add msvc pch
IF (VRN_PRECOMPILED_HEADER AND MSVC)
    ADD_MSVC_PRECOMPILED_HEADER(pch_qt.h pch_qt.cpp VRN_QT_SOURCES)
ENDIF()

# define target
ADD_LIBRARY(voreen_qt 
    ${VRN_QT_SOURCES} ${VRN_QT_HEADERS} ${VRN_QT_HEADERS_NONMOC} 
    ${VRN_QT_RESOURCES_RCC} 
    ${VRN_QT_FORMS_HEADERS}) 
ADD_DEFINITIONS(${VRN_DEFINITIONS} ${VRN_MODULE_DEFINITIONS})
IF(MSVC)
    ADD_DEFINITIONS("-DVRN_QT_BUILD_DLL")
ENDIF()
INCLUDE_DIRECTORIES(${VRN_INCLUDE_DIRECTORIES})
TARGET_LINK_LIBRARIES(voreen_qt tgt voreen_core ${VRN_EXTERNAL_LIBRARIES})

# group source files by subdir
IF(VRN_GROUP_SOURCE_FILES)
    LIST(APPEND RemovePathPrefixes "include/voreen/qt")
    DEFINE_SOURCE_GROUPS_FROM_SUBDIR(VRN_QT_SOURCES ${VRN_HOME} RemovePathPrefixes)
    DEFINE_SOURCE_GROUPS_FROM_SUBDIR(VRN_QT_HEADERS ${VRN_HOME} RemovePathPrefixes)
    DEFINE_SOURCE_GROUPS_FROM_SUBDIR(VRN_QT_HEADERS_NONMOC ${VRN_HOME} RemovePathPrefixes)
ENDIF()

# add Qt dependency
INCLUDE_DIRECTORIES(${QT_INCLUDE_DIRECTORIES})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}) #< location of generated ui headers in out-of-source builds
ADD_DEFINITIONS(${QT_DEFINITIONS})
TARGET_LINK_LIBRARIES(voreen_qt ${QT_LIBRARIES})

# add gcc pch
IF (VRN_PRECOMPILED_HEADER AND CMAKE_COMPILER_IS_GNUCXX)
    ADD_GCC_PRECOMPILED_HEADER(voreen_qt pch_qt.h)
ENDIF()

# Copy Qt5 DLLs using windeployqt.exe tool.
# The tool will be executed right after voreen_qt has been build, since 
# the binary will be scanned for Qt dependencies by the tool.
IF(WIN32)
    UNSET(WINDEPLOYQT_EXE CACHE) #Old path might not longer exist.
    FIND_PROGRAM(WINDEPLOYQT_EXE windeployqt PATHS ${Qt5_DIR}/../../../bin/)
    IF(NOT WINDEPLOYQT_EXE)
        MESSAGE(FATAL_ERROR "Unable to find windeployqt.exe")
    ENDIF()

    # We copy the compiler runtime manually since they are not only needed by voreen_qt.
    SET(WINDEPLOYQT_ARGS "--no-compiler-runtime")
    IF(CMAKE_BUILD_TYPE STREQUAL "Debug")
        LIST(APPEND WINDEPLOYQT_ARGS "--debug")
    ELSEIF(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
        LIST(APPEND WINDEPLOYQT_ARGS "--release-with-debug-info")
    ELSEIF(CMAKE_BUILD_TYPE STREQUAL "Release")
        LIST(APPEND WINDEPLOYQT_ARGS "--release")
    ENDIF()

    # Execute windeployqt tool right after the library has been build.
    ADD_CUSTOM_COMMAND(TARGET voreen_qt POST_BUILD COMMAND ${WINDEPLOYQT_EXE} --dir $<TARGET_FILE_DIR:voreen_qt> ${WINDEPLOYQT_ARGS} $<TARGET_FILE:voreen_qt>)
    
ENDIF()

################################################################################
# deployment
################################################################################

IF(VRN_ADD_INSTALL_TARGET)
    INSTALL(TARGETS voreen_qt
        RUNTIME DESTINATION .
        LIBRARY DESTINATION .
    )
    
    # Execute the windeployqt tool again for the actual deployment, in particular building the INSTALL target.
    # Since CMake's generator expressions can not be used here, we need to specify our paths directly.
    # This approach should suffice because we already know we are in a win32 environment.
    if(WIN32)
        INSTALL(CODE "EXECUTE_PROCESS(COMMAND ${WINDEPLOYQT_EXE} --dir ${CMAKE_INSTALL_PREFIX} ${WINDEPLOYQT_ARGS} ${CMAKE_INSTALL_PREFIX}/voreen_qt.dll)")
    ENDIF()
ENDIF()
