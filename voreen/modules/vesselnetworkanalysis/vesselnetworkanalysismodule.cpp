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

#include "vesselnetworkanalysismodule.h"

#include"processors/appropriatespacinglinker.h"
#include"processors/vascusynthgraphloader.h"
#ifdef LEMON_FOUND
#include "processors/vesselgraphcomparison.h"
#endif
#include"processors/vesselgraphcreator.h"
#include"processors/vesselgraphglobalstats.h"
#include"processors/vesselgraphrefiner.h"
#include"processors/vesselgraphperturbation.h"
#include"processors/vesselgraphrenderer.h"
#include"processors/vesselgraphsave.h"
#include"processors/vesselgraphskeletonextractor.h"
#include"processors/vesselgraphsource.h"
#include"processors/vesselgraphselector.h"
#include"processors/vesselnessextractor.h"
#include"processors/volumemultiplier.h"
#include"processors/volumeslicepadding.h"
#include"processors/volumesurfacenoise.h"
#include"processors/volumethinning.h"

namespace voreen {

VesselNetworkAnalysisModule::VesselNetworkAnalysisModule(const std::string& modulePath)
    : VoreenModule(modulePath)
{
    setID("VesselNetworkAnalysis");
    setGuiName("VesselNetworkAnalysis");

    addShaderPath(getModulePath("glsl"));

    registerProcessor(new AppropriateSpacingLinker());
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
    registerProcessor(new VesselnessExtractor());
    registerProcessor(new VolumeMultiplier());
    registerProcessor(new VolumeSlicePadding());
    registerProcessor(new VolumeSurfaceNoise());
    registerProcessor(new VolumeThinning());
}

} // namespace
