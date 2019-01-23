
################################################################################
# Core module resources
################################################################################
SET(MOD_CORE_MODULECLASS BaseModule)

SET(MOD_CORE_SOURCES
    # Processor sources
    ${MOD_DIR}/processors/entryexitpoints/eepgeometryintegrator.cpp
    ${MOD_DIR}/processors/entryexitpoints/meshentryexitpoints.cpp
    
    ${MOD_DIR}/processors/geometry/boundingboxrenderer.cpp
    ${MOD_DIR}/processors/geometry/camerapositionrenderer.cpp
    ${MOD_DIR}/processors/geometry/geometryprocessor.cpp
    ${MOD_DIR}/processors/geometry/geometryrenderer.cpp
    ${MOD_DIR}/processors/geometry/geometryclipping.cpp
    ${MOD_DIR}/processors/geometry/geometrysequencecreator.cpp
    ${MOD_DIR}/processors/geometry/geometryslabclipping.cpp
    ${MOD_DIR}/processors/geometry/geometrytransformation.cpp
    ${MOD_DIR}/processors/geometry/geometrytransformationvolume.cpp
    ${MOD_DIR}/processors/geometry/lightwidgetrenderer.cpp
    ${MOD_DIR}/processors/geometry/planemanipulation.cpp
    ${MOD_DIR}/processors/geometry/planewidgetprocessor.cpp
    ${MOD_DIR}/processors/geometry/pointlistrenderer.cpp
    ${MOD_DIR}/processors/geometry/pointsegmentlistrenderer.cpp
    ${MOD_DIR}/processors/geometry/quadricrenderer.cpp
    ${MOD_DIR}/processors/geometry/slicepositionrenderer.cpp
    ${MOD_DIR}/processors/geometry/trianglemeshconverter.cpp
    
    ${MOD_DIR}/processors/image/background.cpp
    ${MOD_DIR}/processors/image/binaryimageprocessor.cpp
    ${MOD_DIR}/processors/image/colordepth.cpp
    ${MOD_DIR}/processors/image/compositor.cpp
    ${MOD_DIR}/processors/image/convolution.cpp
    ${MOD_DIR}/processors/image/depthdarkening.cpp
    ${MOD_DIR}/processors/image/edgedetect.cpp
    ${MOD_DIR}/processors/image/explosioncompositor.cpp
    ${MOD_DIR}/processors/image/fade.cpp
    ${MOD_DIR}/processors/image/fxaa.cpp
    ${MOD_DIR}/processors/image/gaussian.cpp
    ${MOD_DIR}/processors/image/grayscale.cpp
    ${MOD_DIR}/processors/image/imagemasking.cpp
    ${MOD_DIR}/processors/image/imagemorphology.cpp
    ${MOD_DIR}/processors/image/imageoverlay.cpp
    ${MOD_DIR}/processors/image/imagethreshold.cpp
    ${MOD_DIR}/processors/image/loopcompositor.cpp
    ${MOD_DIR}/processors/image/mean.cpp
    ${MOD_DIR}/processors/image/median.cpp
    ${MOD_DIR}/processors/image/nonminmaxsuppression.cpp
    ${MOD_DIR}/processors/image/orientationoverlay.cpp
    ${MOD_DIR}/processors/image/quadview.cpp
    ${MOD_DIR}/processors/image/regionofinterest2d.cpp
    ${MOD_DIR}/processors/image/splitter.cpp
    ${MOD_DIR}/processors/image/textoverlay.cpp
    ${MOD_DIR}/processors/image/tripleview.cpp
    ${MOD_DIR}/processors/image/unaryimageprocessor.cpp
    ${MOD_DIR}/processors/image/unsharpmasking.cpp
    
    ${MOD_DIR}/processors/proxygeometry/explosionproxygeometry.cpp
    ${MOD_DIR}/processors/proxygeometry/multipassproxygeometry.cpp
    ${MOD_DIR}/processors/proxygeometry/multiplanarproxygeometry.cpp
    ${MOD_DIR}/processors/proxygeometry/multivolumeproxygeometry.cpp
    ${MOD_DIR}/processors/proxygeometry/optimizedproxygeometry.cpp
    
    ${MOD_DIR}/processors/render/cpuraycaster.cpp
    ${MOD_DIR}/processors/render/multiplanarslicerenderer.cpp
    ${MOD_DIR}/processors/render/multivolumeraycaster.cpp
    ${MOD_DIR}/processors/render/rgbraycaster.cpp
    ${MOD_DIR}/processors/render/segmentationraycaster.cpp
    ${MOD_DIR}/processors/render/simpleraycaster.cpp
    ${MOD_DIR}/processors/render/singlevolumeraycaster.cpp
    ${MOD_DIR}/processors/render/singlevolumeslicer.cpp
    ${MOD_DIR}/processors/render/sliceviewer.cpp
    
    ${MOD_DIR}/processors/utility/clockprocessor.cpp
    ${MOD_DIR}/processors/utility/distancemeasure.cpp
    ${MOD_DIR}/processors/utility/imagesequenceloopfinalizer.cpp
    ${MOD_DIR}/processors/utility/imagesequenceloopinitiator.cpp
    ${MOD_DIR}/processors/utility/intensitymeasure.cpp
    ${MOD_DIR}/processors/utility/metadataextractor.cpp
    ${MOD_DIR}/processors/utility/multiscale.cpp
    ${MOD_DIR}/processors/utility/renderloopfinalizer.cpp
    ${MOD_DIR}/processors/utility/renderloopinitiator.cpp
    ${MOD_DIR}/processors/utility/renderstore.cpp
    ${MOD_DIR}/processors/utility/scale.cpp
    ${MOD_DIR}/processors/utility/segmentationvalidation.cpp
    ${MOD_DIR}/processors/utility/volumelistmodalityfilter.cpp
    ${MOD_DIR}/processors/utility/volumeinformation.cpp
    ${MOD_DIR}/processors/utility/volumepicking.cpp
    
    ${MOD_DIR}/processors/volume/vectormagnitude.cpp
    ${MOD_DIR}/processors/volume/volumechannelmerger.cpp
    ${MOD_DIR}/processors/volume/volumechannelseparator.cpp
    ${MOD_DIR}/processors/volume/volumecombine.cpp
    ${MOD_DIR}/processors/volume/volumecomposer.cpp
    ${MOD_DIR}/processors/volume/volumecreate.cpp
    ${MOD_DIR}/processors/volume/volumecreatebase.cpp
    ${MOD_DIR}/processors/volume/volumecrop.cpp
    ${MOD_DIR}/processors/volume/volumecubify.cpp
    ${MOD_DIR}/processors/volume/volumedecomposer.cpp
    ${MOD_DIR}/processors/volume/volumedistancetransform.cpp
    ${MOD_DIR}/processors/volume/volumeformatconversion.cpp
    ${MOD_DIR}/processors/volume/volumegradient.cpp
    ${MOD_DIR}/processors/volume/volumehalfsample.cpp
    ${MOD_DIR}/processors/volume/volumeinversion.cpp
    ${MOD_DIR}/processors/volume/volumeinputswitch.cpp
    ${MOD_DIR}/processors/volume/volumelandmarkregistration.cpp
    ${MOD_DIR}/processors/volume/volumelistspacing.cpp
    ${MOD_DIR}/processors/volume/volumemasking.cpp
    ${MOD_DIR}/processors/volume/volumemirror.cpp
    ${MOD_DIR}/processors/volume/volumeoctreelevelextractor.cpp
    ${MOD_DIR}/processors/volume/volumeresample.cpp
    ${MOD_DIR}/processors/volume/volumespacing.cpp
    ${MOD_DIR}/processors/volume/volumetransformation.cpp
    ${MOD_DIR}/processors/volume/volumeoffset.cpp
    ${MOD_DIR}/processors/volume/volumevorticity.cpp

    # Volume readers/writers sources
    ${MOD_DIR}/io/analyzevolumereader.cpp
    ${MOD_DIR}/io/brukervolumereader.cpp
    ${MOD_DIR}/io/ecat7volumereader.cpp
    ${MOD_DIR}/io/inveonvolumereader.cpp
    ${MOD_DIR}/io/multivolumereader.cpp
    ${MOD_DIR}/io/mhdvolumereader.cpp
    ${MOD_DIR}/io/mhdvolumewriter.cpp
    ${MOD_DIR}/io/nrrdvolumereader.cpp
    ${MOD_DIR}/io/nrrdvolumewriter.cpp
    ${MOD_DIR}/io/quadhidacvolumereader.cpp
    ${MOD_DIR}/io/synth2dreader.cpp
    ${MOD_DIR}/io/rawvoxvolumereader.cpp
    ${MOD_DIR}/io/tuvvolumereader.cpp
)


SET(MOD_CORE_HEADERS
    # Processor headers
    ${MOD_DIR}/processors/entryexitpoints/eepgeometryintegrator.h
    ${MOD_DIR}/processors/entryexitpoints/meshentryexitpoints.h
    
    ${MOD_DIR}/processors/geometry/boundingboxrenderer.h
    ${MOD_DIR}/processors/geometry/camerapositionrenderer.h
    ${MOD_DIR}/processors/geometry/geometryprocessor.h
    ${MOD_DIR}/processors/geometry/geometryrenderer.h
    ${MOD_DIR}/processors/geometry/geometryclipping.h
    ${MOD_DIR}/processors/geometry/geometrysequencecreator.h
    ${MOD_DIR}/processors/geometry/geometryslabclipping.h
    ${MOD_DIR}/processors/geometry/geometrytransformation.h
    ${MOD_DIR}/processors/geometry/geometrytransformationvolume.h
    ${MOD_DIR}/processors/geometry/lightwidgetrenderer.h
    ${MOD_DIR}/processors/geometry/planemanipulation.h
    ${MOD_DIR}/processors/geometry/planewidgetprocessor.h
    ${MOD_DIR}/processors/geometry/pointlistrenderer.h
    ${MOD_DIR}/processors/geometry/pointsegmentlistrenderer.h
    ${MOD_DIR}/processors/geometry/quadricrenderer.h
    ${MOD_DIR}/processors/geometry/slicepositionrenderer.h
    ${MOD_DIR}/processors/geometry/trianglemeshconverter.h
    
    ${MOD_DIR}/processors/image/background.h
    ${MOD_DIR}/processors/image/binaryimageprocessor.h
    ${MOD_DIR}/processors/image/colordepth.h
    ${MOD_DIR}/processors/image/compositor.h
    ${MOD_DIR}/processors/image/convolution.h
    ${MOD_DIR}/processors/image/depthdarkening.h
    ${MOD_DIR}/processors/image/edgedetect.h
    ${MOD_DIR}/processors/image/explosioncompositor.h
    ${MOD_DIR}/processors/image/fade.h
    ${MOD_DIR}/processors/image/fxaa.h
    ${MOD_DIR}/processors/image/gaussian.h
    ${MOD_DIR}/processors/image/grayscale.h
    ${MOD_DIR}/processors/image/imagemasking.h
    ${MOD_DIR}/processors/image/imagemorphology.h
    ${MOD_DIR}/processors/image/imageoverlay.h
    ${MOD_DIR}/processors/image/imagethreshold.h
    ${MOD_DIR}/processors/image/loopcompositor.h
    ${MOD_DIR}/processors/image/mean.h
    ${MOD_DIR}/processors/image/median.h
    ${MOD_DIR}/processors/image/nonminmaxsuppression.h
    ${MOD_DIR}/processors/image/orientationoverlay.h
    ${MOD_DIR}/processors/image/quadview.h
    ${MOD_DIR}/processors/image/regionofinterest2d.h
    ${MOD_DIR}/processors/image/splitter.h
    ${MOD_DIR}/processors/image/textoverlay.h
    ${MOD_DIR}/processors/image/tripleview.h
    ${MOD_DIR}/processors/image/unaryimageprocessor.h
    ${MOD_DIR}/processors/image/unsharpmasking.h
    
    ${MOD_DIR}/processors/proxygeometry/explosionproxygeometry.h
    ${MOD_DIR}/processors/proxygeometry/multipassproxygeometry.h
    ${MOD_DIR}/processors/proxygeometry/multiplanarproxygeometry.h
    ${MOD_DIR}/processors/proxygeometry/multivolumeproxygeometry.h
    ${MOD_DIR}/processors/proxygeometry/optimizedproxygeometry.h
    
    ${MOD_DIR}/processors/render/cpuraycaster.h
    ${MOD_DIR}/processors/render/multiplanarslicerenderer.h
    ${MOD_DIR}/processors/render/multivolumeraycaster.h
    ${MOD_DIR}/processors/render/rgbraycaster.h
    ${MOD_DIR}/processors/render/segmentationraycaster.h
    ${MOD_DIR}/processors/render/simpleraycaster.h
    ${MOD_DIR}/processors/render/singlevolumeraycaster.h
    ${MOD_DIR}/processors/render/singlevolumeslicer.h
    ${MOD_DIR}/processors/render/sliceviewer.h
    
    ${MOD_DIR}/processors/utility/clockprocessor.h
    ${MOD_DIR}/processors/utility/distancemeasure.h
    ${MOD_DIR}/processors/utility/imagesequenceloopfinalizer.h
    ${MOD_DIR}/processors/utility/imagesequenceloopinitiator.h
    ${MOD_DIR}/processors/utility/intensitymeasure.h
    ${MOD_DIR}/processors/utility/metadataextractor.h
    ${MOD_DIR}/processors/utility/multiscale.h
    ${MOD_DIR}/processors/utility/renderloopfinalizer.h
    ${MOD_DIR}/processors/utility/renderloopinitiator.h
    ${MOD_DIR}/processors/utility/renderstore.h
    ${MOD_DIR}/processors/utility/scale.h
    ${MOD_DIR}/processors/utility/segmentationvalidation.h
    ${MOD_DIR}/processors/utility/volumelistmodalityfilter.h
    ${MOD_DIR}/processors/utility/volumeinformation.h
    ${MOD_DIR}/processors/utility/volumepicking.h
    
    ${MOD_DIR}/processors/volume/vectormagnitude.h
    ${MOD_DIR}/processors/volume/volumechannelmerger.h
    ${MOD_DIR}/processors/volume/volumechannelseparator.h
    ${MOD_DIR}/processors/volume/volumecombine.h
    ${MOD_DIR}/processors/volume/volumecomposer.h
    ${MOD_DIR}/processors/volume/volumecreate.h
    ${MOD_DIR}/processors/volume/volumecreatebase.h
    ${MOD_DIR}/processors/volume/volumecrop.h
    ${MOD_DIR}/processors/volume/volumecubify.h
    ${MOD_DIR}/processors/volume/volumedecomposer.h
    ${MOD_DIR}/processors/volume/volumedistancetransform.h
    ${MOD_DIR}/processors/volume/volumeformatconversion.h
    ${MOD_DIR}/processors/volume/volumegradient.h
    ${MOD_DIR}/processors/volume/volumehalfsample.h
    ${MOD_DIR}/processors/volume/volumeinversion.h
    ${MOD_DIR}/processors/volume/volumeinputswitch.h
    ${MOD_DIR}/processors/volume/volumelandmarkregistration.h
    ${MOD_DIR}/processors/volume/volumelistspacing.h
    ${MOD_DIR}/processors/volume/volumemasking.h
    ${MOD_DIR}/processors/volume/volumemirror.h
    ${MOD_DIR}/processors/volume/volumeoctreelevelextractor.h
    ${MOD_DIR}/processors/volume/volumeresample.h
    ${MOD_DIR}/processors/volume/volumespacing.h
    ${MOD_DIR}/processors/volume/volumetransformation.h
    ${MOD_DIR}/processors/volume/volumeoffset.h
    ${MOD_DIR}/processors/volume/volumevorticity.h

    # Volume readers/writers headers
    ${MOD_DIR}/io/analyzevolumereader.h
    ${MOD_DIR}/io/brukervolumereader.h
    ${MOD_DIR}/io/ecat7volumereader.h
    ${MOD_DIR}/io/inveonvolumereader.h
    ${MOD_DIR}/io/multivolumereader.h
    ${MOD_DIR}/io/mhdvolumereader.h
    ${MOD_DIR}/io/mhdvolumewriter.h
    ${MOD_DIR}/io/nrrdvolumereader.h
    ${MOD_DIR}/io/nrrdvolumewriter.h
    ${MOD_DIR}/io/quadhidacvolumereader.h
    ${MOD_DIR}/io/synth2dreader.h
    ${MOD_DIR}/io/rawvoxvolumereader.h
    ${MOD_DIR}/io/tuvvolumereader.h
)


################################################################################
# Qt module resources
################################################################################
SET(MOD_QT_MODULECLASS BaseModuleQt)

SET(MOD_QT_SOURCES
    ${MOD_DIR}/qt/processorwidgets/baseprocessorwidgetfactory.cpp
)

SET(MOD_QT_HEADERS
)

SET(MOD_QT_HEADERS_NONMOC
    ${MOD_DIR}/qt/processorwidgets/baseprocessorwidgetfactory.h
)


################################################################################
# VoreenVE module resources (currently not present)
################################################################################
# SET(MOD_VE_MODULECLASS BaseModuleVE)

# SET(MOD_VE_SOURCES
# )  
    
# SET(MOD_VE_HEADERS
# )

# SET(MOD_VE_HEADERS_NONMOC
# )


################################################################################
# Deployment
################################################################################
SET(MOD_INSTALL_DIRECTORIES
    ${MOD_DIR}/glsl
)
