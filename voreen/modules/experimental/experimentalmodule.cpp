/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
 * Department of Computer Science.                                                 *
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

#include "experimentalmodule.h"

#include "processors/arrowbillboardtest.h"
#include "processors/crosshairrenderer.h"
#include "processors/depthoffield.h"
#include "processors/divergence.h"
#include "processors/gabor.h"
#include "processors/geometrydelay.h"
#include "processors/geometryeventblocker.h"
#include "processors/geometryboundingbox.h"
#include "processors/illuminationlineraycaster.h"
#include "processors/imageabstraction.h"
#include "processors/manualsegmentation.h"
#include "processors/manualsegmentationstorage.h"
#include "processors/markstats.h"
//#include "processors/lineprofile.h"
//#include "processors/lp_plot.h"
#include "processors/meshfrustumclipping.h"
#include "processors/mousepositionrenderer.h"
#include "processors/multivolumecrosssectionanalyzer.h"
#include "processors/multivolumegeometryraycaster.h"
#include "processors/normalestimation.h"
#include "processors/radarglyphrenderer2d.h"
#include "processors/radarglyphrenderer3d.h"
#include "processors/radarglyphvolumegenerator.h"
#include "processors/pwivolume.h"
#include "processors/regiongrowing.h"
#include "processors/seedpointgenerator.h"
#include "processors/slicecolormapper.h"
#include "processors/slicegenerator.h"
#include "processors/swivolumecompare.h"
#include "processors/tiledraycastinginitiator.h"
#include "processors/tiledraycastingfinalizer.h"
//#include "processors/volumebrickloopfinalizer.h"
//#include "processors/volumebrickloopinitiator.h"
#include "processors/volumecenter.h"
#include "processors/volumecollector.h"
#include "processors/volumecurvature.h"
#include "processors/volumepositionsimplifier.h"
#include "processors/volumeseginfo.h"
#include "processors/volumeselectortime.h"
//#include "processors/volumesliceloopfinalizer.h"
//#include "processors/volumesliceloopinitiator.h"
#include "processors/volumestreamprocessor.h"

namespace voreen {

ExperimentalModule::ExperimentalModule(const std::string& modulePath)
    : VoreenModule(modulePath)
{
    setID("Experimental");
    setGuiName("Experimental");
    registerProcessor(new ArrowBillboardTest);
    registerSerializableType(new IlluminationLineRaycaster());
    registerSerializableType(new MeshFrustumClipping());
    registerSerializableType(new CrosshairRenderer());
    registerSerializableType(new DepthOfField());
    registerSerializableType(new Divergence());
    registerSerializableType(new Gabor());
    registerSerializableType(new GeometryDelay());
    registerSerializableType(new GeometryEventBlocker());
    registerSerializableType(new GeometryBoundingBox());
    registerSerializableType(new ImageAbstraction());
    registerSerializableType(new ManualSegmentation());
    registerSerializableType(new ManualSegmentationStorage());
    registerSerializableType(new MarkStats());
    //registerSerializableType(new LineProfile());
    //registerSerializableType(new LP_Plot());
    registerSerializableType(new MousePositionRenderer());
    registerSerializableType(new MultiVolumeCrossSectionAnalyzer());
#ifdef GL_ATOMIC_COUNTER_BUFFER //disable compilation for old gl headers
    registerSerializableType(new MultiVolumeGeometryRaycaster());
#endif
    registerSerializableType(new NormalEstimation());
    registerSerializableType(new RadarGlyphRenderer2D());
    registerSerializableType(new RadarGlyphRenderer3D());
    registerSerializableType(new RadarGlyphVolumeGenerator());
    registerSerializableType(new PWIVolume());
    registerSerializableType(new RegionGrowingProcessor());
    registerSerializableType(new SeedpointGenerator());
    registerSerializableType(new SliceColorMapper());
    registerSerializableType(new SliceGenerator());
    registerSerializableType(new SWIVolumeCompare());
    registerSerializableType(new TiledRaycastingInitiator());
    registerSerializableType(new TiledRaycastingFinalizer());
    //registerSerializableType(new VolumeBrickLoopFinalizer());
    //registerSerializableType(new VolumeBrickLoopInitiator());
    registerSerializableType(new VolumeCenter());
    registerSerializableType(new VolumeCollector());
    registerSerializableType(new VolumeCurvature());
    registerSerializableType(new VolumePositionSimplifier());
    registerSerializableType(new VolumeSegInfo());
    registerSerializableType(new VolumeSelectorTime());
    //registerSerializableType(new VolumeSliceLoopFinalizer());
    //registerSerializableType(new VolumeSliceLoopInitiator());
    registerSerializableType(new VolumeStreamProcessor());

    addShaderPath(getModulePath("glsl"));
}

} // namespace

