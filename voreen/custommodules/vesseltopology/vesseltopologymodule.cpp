/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#include "vesseltopologymodule.h"
#include "processors/appropriatespacinglinker.h"
#include "processors/aortasegmentation.h"
#include "processors/localandglobalthreshold.h"
#include "processors/segmentationlistvalidation.h"
#include "processors/subgraphextractor.h"
#include "processors/templatesubgraphextractor.h"
#include "processors/vascusynthgraphloader.h"
#ifdef LEMON_FOUND
#include "processors/vesselgraphcomparison.h"
#endif
#include "processors/vesselgraphcreator.h"
#include "processors/vesselgraphglobalstats.h"
#include "processors/vesselgraphrefiner.h"
#include "processors/vesselgraphperturbation.h"
#include "processors/vesselgraphrenderer.h"
#include "processors/vesselgraphsave.h"
#include "processors/vesselgraphskeletonextractor.h"
#include "processors/vesselgraphsource.h"
#include "processors/vesselgraphselector.h"
#include "processors/vesselgraphstatplotter.h"
#include "processors/vesselnessextractor.h"
#include "processors/volumefloodfill.h"
#include "processors/createTestVolume.h"
#include "processors/createVesselAroundPoints.h"
#include "processors/foldVessel.h"
#include "processors/unfoldVessel.h"
#include "processors/volumemultiplier.h"
#include "processors/volumemultithreshold.h"
#include "processors/volumeslicepadding.h"
#include "processors/volumethinning.h"

#include "processors/volumelistloopinitiator.h"
#include "processors/volumelistloopfinalizer.h"

namespace voreen {

VesselTopologyModule::VesselTopologyModule(const std::string& modulePath)
    : VoreenModule(modulePath)
{
    setID("VesselTopology");
    setGuiName("VesselTopology");

    addShaderPath(getModulePath("glsl"));

    registerProcessor(new AppropriateSpacingLinker());
    registerProcessor(new AortaSegmentation());
    registerProcessor(new LocalAndGlobalThreshold());
    registerProcessor(new SegmentationListValidation());
    registerProcessor(new SubGraphExtractor());
    registerProcessor(new TemplateSubgraphExtractor());
    registerProcessor(new VascuSynthGraphLoader());
#ifdef LEMON_FOUND
    registerProcessor(new VesselGraphComparison());
#endif
    registerProcessor(new VesselGraphCreator());
    registerProcessor(new VesselGraphGlobalStats());
    registerProcessor(new VesselGraphRefiner());
    registerProcessor(new VesselGraphPerturbation());
    registerProcessor(new VesselGraphRenderer());
    registerProcessor(new VesselGraphSave());
    registerProcessor(new VesselGraphSkeletonExtractor());
    registerProcessor(new VesselGraphSource());
    registerProcessor(new VesselGraphSelector());
    registerProcessor(new VesselGraphStatPlotter());
    registerProcessor(new VesselnessExtractor());
    registerProcessor(new VolumeFloodFill());
    registerProcessor(new createTestVolume());
    registerProcessor(new foldVessel());
    registerProcessor(new unfoldVessel());
    registerProcessor(new createVesselAroundPoints());
    registerProcessor(new VolumeMultiThreshold());
    registerProcessor(new VolumeMultiplier());
    registerProcessor(new VolumeSlicePadding());
    registerProcessor(new VolumeThinning());

    registerProcessor(new VolumeListLoopInitiator());
    registerProcessor(new VolumeListLoopFinalizer());
}

} // namespace
