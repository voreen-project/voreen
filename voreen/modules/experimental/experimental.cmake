
################################################################################
# Dependencies 
################################################################################
IF (NOT VRN_MODULE_PLOTTING)
    MESSAGE(FATAL_ERROR "Experimental module requires the Plotting module")
ENDIF()

################################################################################
# Core module resources 
################################################################################
SET(MOD_CORE_MODULECLASS ExperimentalModule)
SET(MOD_REQUIRE_OPENGL_COMPATIBILITY_PROFILE ON)

# needed for clock_gettime
IF(UNIX)
    SET(MOD_LIBRARIES rt)
ENDIF()

SET(MOD_CORE_SOURCES
    #${MOD_DIR}/processors/lineprofile.cpp
    #${MOD_DIR}/processors/lp_plot.cpp
    #${MOD_DIR}/utils/LineProfile/lp_eigenfunctionfit.cpp
    #${MOD_DIR}/utils/LineProfile/lp_graphic.cpp
    #${MOD_DIR}/utils/LineProfile/lp_maxthreshold.cpp
    #${MOD_DIR}/utils/LineProfile/lp_measure.cpp
    #${MOD_DIR}/utils/LineProfile/lp_plotmodifier.cpp
    #${MOD_DIR}/utils/LineProfile/lp_texteditor.cpp
    ${MOD_DIR}/processors/arrowbillboardtest.cpp
    ${MOD_DIR}/processors/crosshairrenderer.cpp
    ${MOD_DIR}/processors/depthoffield.cpp
    ${MOD_DIR}/processors/divergence.cpp
    ${MOD_DIR}/processors/gabor.cpp
    ${MOD_DIR}/processors/geometryboundingbox.cpp
    ${MOD_DIR}/processors/geometrydelay.cpp
    ${MOD_DIR}/processors/geometryeventblocker.cpp
    ${MOD_DIR}/processors/illuminationlineraycaster.cpp
    ${MOD_DIR}/processors/imageabstraction.cpp
    ${MOD_DIR}/processors/jacobian.cpp
    ${MOD_DIR}/processors/manualsegmentation.cpp
    ${MOD_DIR}/processors/manualsegmentationstorage.cpp
    ${MOD_DIR}/processors/markstats.cpp
    ${MOD_DIR}/processors/meshfrustumclipping.cpp 
    ${MOD_DIR}/processors/mousepositionrenderer.cpp
    ${MOD_DIR}/processors/multivolumecrosssectionanalyzer.cpp
    ${MOD_DIR}/processors/multivolumegeometryraycaster.cpp
    ${MOD_DIR}/processors/normalestimation.cpp
    ${MOD_DIR}/processors/pwivolume.cpp
    ${MOD_DIR}/processors/regiongrowing.cpp
    ${MOD_DIR}/processors/radarglyphrendererbase.cpp
    ${MOD_DIR}/processors/radarglyphrenderer2d.cpp
    ${MOD_DIR}/processors/radarglyphrenderer3d.cpp
    ${MOD_DIR}/processors/radarglyphvolumegenerator.cpp
    ${MOD_DIR}/processors/seedpointgenerator.cpp
    ${MOD_DIR}/processors/slicecolormapper.cpp
    ${MOD_DIR}/processors/slicegenerator.cpp
    ${MOD_DIR}/processors/swivolumecompare.cpp
    ${MOD_DIR}/processors/tiledraycastinginitiator.cpp
    ${MOD_DIR}/processors/tiledraycastingfinalizer.cpp
    #${MOD_DIR}/processors/volumebrickloopfinalizer.cpp
    #${MOD_DIR}/processors/volumebrickloopinitiator.cpp
    ${MOD_DIR}/processors/volumecenter.cpp
    ${MOD_DIR}/processors/volumecollector.cpp
    ${MOD_DIR}/processors/volumecurvature.cpp
    ${MOD_DIR}/processors/volumepositionsimplifier.cpp
    ${MOD_DIR}/processors/volumeseginfo.cpp
    ${MOD_DIR}/processors/volumeselectortime.cpp
    #${MOD_DIR}/processors/volumesliceloopfinalizer.cpp
    #${MOD_DIR}/processors/volumesliceloopinitiator.cpp
    ${MOD_DIR}/processors/volumestreamprocessor.cpp
    
    # properties
    ${MOD_DIR}/properties/volumestreamproperty.cpp
)

SET(MOD_CORE_HEADERS
    ${MOD_DIR}/processors/illuminationlineraycaster.h
    #${MOD_DIR}/processors/lineprofile.h
    #${MOD_DIR}/processors/lp_plot.h
    #${MOD_DIR}/utils/LineProfile/lp_eigenfunctionfit.h
    #${MOD_DIR}/utils/LineProfile/lp_graphic.h
    #${MOD_DIR}/utils/LineProfile/lp_maxthreshold.h
    #${MOD_DIR}/utils/LineProfile/lp_measure.h
    #${MOD_DIR}/utils/LineProfile/lp_plotmodifier.h
    #${MOD_DIR}/utils/LineProfile/lp_texteditor.h
    ${MOD_DIR}/processors/arrowbillboardtest.h
    ${MOD_DIR}/processors/crosshairrenderer.h
    ${MOD_DIR}/processors/depthoffield.h
    ${MOD_DIR}/processors/divergence.h
    ${MOD_DIR}/processors/gabor.h
    ${MOD_DIR}/processors/geometryboundingbox.h
    ${MOD_DIR}/processors/geometrydelay.h
    ${MOD_DIR}/processors/geometryeventblocker.h
    ${MOD_DIR}/processors/imageabstraction.h
    ${MOD_DIR}/processors/jacobian.h
    ${MOD_DIR}/processors/manualsegmentation.h
    ${MOD_DIR}/processors/manualsegmentationstorage.h
    ${MOD_DIR}/processors/markstats.h
    ${MOD_DIR}/processors/meshfrustumclipping.h
    ${MOD_DIR}/processors/mousepositionrenderer.h
    ${MOD_DIR}/processors/multivolumecrosssectionanalyzer.h
    ${MOD_DIR}/processors/multivolumegeometryraycaster.h
    ${MOD_DIR}/processors/normalestimation.h
    ${MOD_DIR}/processors/radarglyphrendererbase.h
    ${MOD_DIR}/processors/radarglyphrenderer2d.h
    ${MOD_DIR}/processors/radarglyphrenderer3d.h
    ${MOD_DIR}/processors/radarglyphvolumegenerator.h
    ${MOD_DIR}/processors/pwivolume.h
    ${MOD_DIR}/processors/regiongrowing.h
    ${MOD_DIR}/processors/seedpointgenerator.h
    ${MOD_DIR}/processors/slicecolormapper.h
    ${MOD_DIR}/processors/slicegenerator.h
    ${MOD_DIR}/processors/swivolumecompare.h
    ${MOD_DIR}/processors/tiledraycastinginitiator.h
    ${MOD_DIR}/processors/tiledraycastingfinalizer.h
    #${MOD_DIR}/processors/volumebrickloopfinalizer.h
    #${MOD_DIR}/processors/volumebrickloopinitiator.h
    ${MOD_DIR}/processors/volumecenter.h
    ${MOD_DIR}/processors/volumecollector.h
    ${MOD_DIR}/processors/volumecurvature.h
    ${MOD_DIR}/processors/volumepositionsimplifier.h
    ${MOD_DIR}/processors/volumeseginfo.h
    ${MOD_DIR}/processors/volumeselectortime.h
    #${MOD_DIR}/processors/volumesliceloopfinalizer.h
    #${MOD_DIR}/processors/volumesliceloopinitiator.h
    ${MOD_DIR}/processors/volumestreamprocessor.h
    
    # properties
    ${MOD_DIR}/properties/volumestreamproperty.h
)
   

################################################################################
# Qt module resources 
################################################################################
SET(MOD_QT_MODULECLASS ExperimentalModuleQt)

SET(MOD_QT_SOURCES
    #${MOD_DIR}/qt/segmentationplugin.cpp must be adapted to new SegmentationRaycaster
    ${MOD_DIR}/qt/experimentalprocessorwidgetfactory.cpp
    ${MOD_DIR}/qt/experimentalpropertywidgetfactory.cpp
    ${MOD_DIR}/qt/volumestreampropertywidget.cpp
    ${MOD_DIR}/qt/volumestreamwidget.cpp
)  
    
SET(MOD_QT_HEADERS
    #${MOD_DIR}/qt/segmentationplugin.h  must be adapted to new SegmentationRaycaster
    ${MOD_DIR}/qt/volumestreampropertywidget.h
    ${MOD_DIR}/qt/volumestreamwidget.h
)

SET(MOD_QT_HEADERS_NONMOC
    ${MOD_DIR}/qt/experimentalprocessorwidgetfactory.h
    ${MOD_DIR}/qt/experimentalpropertywidgetfactory.h
)
