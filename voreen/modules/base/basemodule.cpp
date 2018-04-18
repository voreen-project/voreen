/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
 * For a list of authors please refer to the file "CREDITS.txt".                   *
 *                                                                                 *
 * This file is part of the Voreen software package. Voreen is free software:      *
 * you can redistribute it and/or modify it under the terms of the GNU General     *
 * Public License version 2 as published by the Free Software Foundation.          *
 *                                                                                 *
 * Voreen is distributed in the hope that it will be useful, but WITHOUT ANY       *
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR   *
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.      *
 *                                                                                 *
 * You should have received a copy of the GNU General Public License in the file   *
 * "LICENSE.txt" along with this file. If not, see <http://www.gnu.org/licenses/>. *
 *                                                                                 *
 * For non-commercial academic use see the license exception specified in the file *
 * "LICENSE-academic.txt". To get information about commercial licensing please    *
 * contact the authors.                                                            *
 *                                                                                 *
 ***********************************************************************************/

#include "modules/base/basemodule.h"

///
/// Processors
///

// entry-exit points
#include "processors/entryexitpoints/eepgeometryintegrator.h"
#include "processors/entryexitpoints/meshentryexitpoints.h"

// geometry
#include "processors/geometry/boundingboxrenderer.h"
#include "processors/geometry/camerapositionrenderer.h"
#include "processors/geometry/geometryprocessor.h"
#include "processors/geometry/geometryrenderer.h"
#include "processors/geometry/geometryclipping.h"
#include "processors/geometry/geometrysequencecreator.h"
#include "processors/geometry/geometryslabclipping.h"
#include "processors/geometry/geometrytransformation.h"
#include "processors/geometry/geometrytransformationvolume.h"
#include "processors/geometry/lightwidgetrenderer.h"
#include "processors/geometry/pointlistrenderer.h"
#include "processors/geometry/pointsegmentlistrenderer.h"
#include "processors/geometry/planemanipulation.h"
#include "processors/geometry/planewidgetprocessor.h"
#include "processors/geometry/quadricrenderer.h"
#include "processors/geometry/slicepositionrenderer.h"
#include "processors/geometry/trianglemeshconverter.h"

// image
#include "processors/image/background.h"
#include "processors/image/binaryimageprocessor.h"
#include "processors/image/unaryimageprocessor.h"
#include "processors/image/gaussian.h"
#include "processors/image/median.h"
#include "processors/image/mean.h"
#include "processors/image/nonminmaxsuppression.h"
#include "processors/image/colordepth.h"
#include "processors/image/compositor.h"
#include "processors/image/convolution.h"
#include "processors/image/depthdarkening.h"
#include "processors/image/edgedetect.h"
#include "processors/image/explosioncompositor.h"
#include "processors/image/fade.h"
#include "processors/image/fxaa.h"
#include "processors/image/grayscale.h"
#include "processors/image/imagemasking.h"
#include "processors/image/imagemorphology.h"
#include "processors/image/imageoverlay.h"
#include "processors/image/imagethreshold.h"
#include "processors/image/loopcompositor.h"
#include "processors/image/multiview.h"
#include "processors/image/orientationoverlay.h"
#include "processors/image/quadview.h"
#include "processors/image/regionofinterest2d.h"
#include "processors/image/splitter.h"
#include "processors/image/textoverlay.h"
#include "processors/image/tripleview.h"
#include "processors/image/unsharpmasking.h"

// proxy geometry
#include "processors/proxygeometry/explosionproxygeometry.h"
#include "processors/proxygeometry/multipassproxygeometry.h"
#include "processors/proxygeometry/multiplanarproxygeometry.h"
#include "processors/proxygeometry/multivolumeproxygeometry.h"
#include "processors/proxygeometry/optimizedproxygeometry.h"

// render
#include "processors/render/cpuraycaster.h"
#include "processors/render/multiplanarslicerenderer.h"
#include "processors/render/rgbraycaster.h"
#include "processors/render/segmentationraycaster.h"
#include "processors/render/simpleraycaster.h"
#include "processors/render/singlevolumeraycaster.h"
#include "processors/render/multivolumeraycaster.h"
#include "processors/render/singlevolumeslicer.h"
#include "processors/render/sliceviewer.h"

// utility
#include "processors/utility/clockprocessor.h"
#include "processors/utility/distancemeasure.h"
#include "processors/utility/intensitymeasure.h"
#include "processors/utility/imagesequenceloopinitiator.h"
#include "processors/utility/imagesequenceloopfinalizer.h"
#include "processors/utility/metadataextractor.h"
#include "processors/utility/multiscale.h"
#include "processors/utility/renderloopinitiator.h"
#include "processors/utility/renderloopfinalizer.h"
#include "processors/utility/renderstore.h"
#include "processors/utility/scale.h"
#include "processors/utility/segmentationvalidation.h"
#include "processors/utility/volumelistmodalityfilter.h"
#include "processors/utility/volumeinformation.h"
#include "processors/utility/volumepicking.h"

// volume
#include "processors/volume/vectormagnitude.h"
#include "processors/volume/volumechannelmerger.h"
#include "processors/volume/volumechannelseparator.h"
#include "processors/volume/volumecrop.h"
#include "processors/volume/volumedistancetransform.h"
#include "processors/volume/volumeinversion.h"
#include "processors/volume/volumecombine.h"
#include "processors/volume/volumecomposer.h"
#include "processors/volume/volumedecomposer.h"
#include "processors/volume/volumefiltering.h"
#include "processors/volume/volumeformatconversion.h"
#include "processors/volume/volumegradient.h"
#include "processors/volume/volumeinputswitch.h"
#include "processors/volume/volumelandmarkregistration.h"
#include "processors/volume/volumelistspacing.h"
#include "processors/volume/volumespacing.h"
#include "processors/volume/volumeoffset.h"
#include "processors/volume/volumetransformation.h"
#include "processors/volume/volumemasking.h"
#include "processors/volume/volumemorphology.h"
#include "processors/volume/volumeresample.h"
#include "processors/volume/volumecreate.h"
#include "processors/volume/volumecubify.h"
#include "processors/volume/volumehalfsample.h"
#include "processors/volume/volumemirror.h"
#include "processors/volume/volumeoctreelevelextractor.h"
#include "processors/volume/volumevorticity.h"


///
/// volume readers/writers
///
#include "io/analyzevolumereader.h"
#include "io/brukervolumereader.h"
#include "io/ecat7volumereader.h"
#include "io/inveonvolumereader.h"
#include "io/multivolumereader.h"
#include "io/mhdvolumereader.h"
#include "io/mhdvolumewriter.h"
#include "io/nrrdvolumereader.h"
#include "io/nrrdvolumewriter.h"
#include "io/quadhidacvolumereader.h"
#include "io/synth2dreader.h"
#include "io/rawvoxvolumereader.h"
#include "io/tuvvolumereader.h"


namespace voreen {

const std::string BaseModule::loggerCat_("voreen.BaseModule");

BaseModule::BaseModule(const std::string& modulePath) : VoreenModule(modulePath) {
    setID("Base");
    setGuiName("Base");

    // base incompatibel classes because of GLU usage
    // entry-exit points
    registerSerializableType(new EEPGeometryIntegrator());
    registerSerializableType(new MeshEntryExitPoints());

    // geometry
    registerSerializableType(new BoundingBoxRenderer());
    registerSerializableType(new CameraPositionRenderer());
    registerSerializableType(new GeometryProcessor());
    registerSerializableType(new GeometryRenderer());
    registerSerializableType(new GeometryClipping());
    registerSerializableType(new GeometrySequenceCreator());
    registerSerializableType(new GeometrySlabClipping());
    registerSerializableType(new GeometryTransformation());
    registerSerializableType(new GeometryTransformationVolume());
    registerSerializableType(new LightWidgetRenderer());
    registerSerializableType(new PlaneManipulation());
    registerSerializableType(new PlaneWidgetProcessor());
    registerSerializableType(new PointListRenderer());
    registerSerializableType(new PointSegmentListRenderer());
    registerSerializableType(new QuadricRenderer());
    registerSerializableType(new RegionOfInterest2D::RegionOfInterestGeometry());
    registerSerializableType(new SlicePositionRenderer());
    registerSerializableType(new TriangleMeshConverter());

    // image
    registerSerializableType(new Background());
    registerSerializableType(new BinaryImageProcessor());
    registerSerializableType(new ColorDepth());
    registerSerializableType(new Compositor());
    registerSerializableType(new Convolution());
    registerSerializableType(new DepthDarkening());
    registerSerializableType(new EdgeDetect());
    registerSerializableType(new ExplosionCompositor());
    registerSerializableType(new Fade());
    registerSerializableType(new FXAA());
    registerSerializableType(new Gaussian());
    registerSerializableType(new Grayscale());
    registerSerializableType(new ImageMasking());
    registerSerializableType(new ImageMorphology());
    registerSerializableType(new ImageOverlay());
    registerSerializableType(new ImageThreshold());
    registerSerializableType(new Mean());
    registerSerializableType(new Median());
    registerSerializableType(new MultiView());
    registerSerializableType(new NonMinMaxSuppression());
    registerSerializableType(new OrientationOverlay());
    registerSerializableType(new QuadView());
    registerSerializableType(new RegionOfInterest2D());
    registerSerializableType(new Splitter());
    registerSerializableType(new TextOverlay());
    registerSerializableType(new TripleView());
    registerSerializableType(new UnaryImageProcessor());
    registerSerializableType(new UnsharpMasking());

    // proxy geometry
    registerSerializableType(new ExplosionProxyGeometry());
    registerSerializableType(new MultiPassProxyGeometry());
    registerSerializableType(new MultiPlanarProxyGeometry());
    registerSerializableType(new MultiVolumeProxyGeometry());
    registerSerializableType(new OptimizedProxyGeometry());

    // render
    registerSerializableType(new CPURaycaster());
    registerSerializableType(new MultiplanarSliceRenderer());
    registerSerializableType(new MultiVolumeRaycaster());
    registerSerializableType(new SegmentationRaycaster());
    registerSerializableType(new SimpleRaycaster());
    registerSerializableType(new SingleVolumeRaycaster());
    registerSerializableType(new SingleVolumeSlicer());
    registerSerializableType(new SliceViewer());
    registerSerializableType(new RGBRaycaster());

    // utility
    registerSerializableType(new ClockProcessor());
    registerSerializableType(new DistanceMeasure());
    registerSerializableType(new IntensityMeasure());
    registerSerializableType(new ImageSequenceLoopInitiator());
    registerSerializableType(new ImageSequenceLoopFinalizer());
    registerSerializableType(new MetaDataExtractor());
    registerSerializableType(new MultiScale());
    registerSerializableType(new RenderLoopInitiator());
    registerSerializableType(new RenderLoopFinalizer());
    registerSerializableType(new RenderStore());
    registerSerializableType(new SegmentationValidation());
    registerSerializableType(new SingleScale());
    registerSerializableType(new VolumeListModalityFilter());
    registerSerializableType(new VolumeInformation());
    registerSerializableType(new VolumePicking());

    // volume
    registerSerializableType(new VectorMagnitude());
    registerSerializableType(new VolumeChannelMerger());
    registerSerializableType(new VolumeChannelSeparator());
    registerSerializableType(new VolumeCombine());
    registerSerializableType(new VolumeComposer());
    registerSerializableType(new VolumeCreate());
    registerSerializableType(new VolumeCrop());
    registerSerializableType(new VolumeCubify());
    registerSerializableType(new VolumeDecomposer());
    registerSerializableType(new VolumeDistanceTransform());
    registerSerializableType(new VolumeFiltering());
    registerSerializableType(new VolumeFormatConversion());
    registerSerializableType(new VolumeGradient());
    registerSerializableType(new VolumeHalfsample());
    registerSerializableType(new VolumeInputSwitch());
    registerSerializableType(new VolumeInversion());
    registerSerializableType(new VolumeLandmarkRegistration());
    registerSerializableType(new VolumeListSpacing());
    registerSerializableType(new VolumeMasking());
    registerSerializableType(new VolumeMirror());
    registerSerializableType(new VolumeMorphology());
    registerSerializableType(new VolumeOctreeLevelExtractor());
    registerSerializableType(new VolumeResample());
    registerSerializableType(new VolumeSpacing());
    registerSerializableType(new VolumeTransformation());
    registerSerializableType(new VolumeOffset());
    registerSerializableType(new VolumeVorticity());

    // volume readers/writers
    registerVolumeReader(new AnalyzeVolumeReader());
    registerVolumeReader(new BrukerVolumeReader());
    registerVolumeReader(new ECAT7VolumeReader());
    registerVolumeReader(new InveonVolumeReader());
    registerVolumeReader(new MultiVolumeReader(0));
    registerVolumeReader(new MhdVolumeReader());
    registerVolumeWriter(new MhdVolumeWriter());
    registerVolumeReader(new NrrdVolumeReader());
    registerVolumeWriter(new NrrdVolumeWriter());
    registerVolumeReader(new QuadHidacVolumeReader());
    registerVolumeReader(new Synth2DReader());
    registerVolumeReader(new RawVoxVolumeReader());
    registerVolumeReader(new TUVVolumeReader());

    // shader paths
    addShaderPath(getModulePath("glsl"));
    addShaderPath(VoreenApplication::app()->getBasePath() + "/ext/tgt/glsl/immediatemode"); //to use immediatemode.vert shader
}

void BaseModule::initialize() {
    VoreenModule::initialize();
}

void BaseModule::deinitialize() {
    VoreenModule::deinitialize();
}

} // namespace
