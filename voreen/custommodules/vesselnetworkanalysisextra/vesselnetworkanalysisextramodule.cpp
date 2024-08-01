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

#include "vesselnetworkanalysisextramodule.h"

#include "processors/aortasegmentation.h"
#include "processors/createTestVolume.h"
#include "processors/createVesselAroundPoints.h"
#include "processors/foldVessel.h"
#include "processors/interactiveprojectionlabeling.h"
#include "processors/localandglobalthreshold.h"
#include "processors/lymphatictestvesselgenerator.h"
#include "processors/segmentationlistvalidation.h"
#include "processors/subgraphextractor.h"
#include "processors/templatesubgraphextractor.h"
#include "processors/unfoldVessel.h"
#include "processors/vesselgraphstatplotter.h"
#include "processors/volumefloodfill.h"
#include "processors/volumelistloopfinalizer.h"
#include "processors/volumelistloopinitiator.h"
#include "processors/volumemultithreshold.h"

namespace voreen {

VesselNetworkAnalysisExtraModule::VesselNetworkAnalysisExtraModule(const std::string& modulePath)
    : VoreenModule(modulePath)
{
    setID("VesselNetworkAnalysisExtra");
    setGuiName("VesselNetworkAnalysisExtra");

    addShaderPath(getModulePath("glsl"));

    registerProcessor(new AortaSegmentation());
    registerProcessor(new InteractiveProjectionLabeling());
    registerProcessor(new LocalAndGlobalThreshold());
    registerProcessor(new LymphaticTestVesselGenerator());
    registerProcessor(new SegmentationListValidation());
    registerProcessor(new SubGraphExtractor());
    registerProcessor(new TemplateSubgraphExtractor());
    registerProcessor(new VesselGraphStatPlotter());
    registerProcessor(new VolumeFloodFill());
    registerProcessor(new VolumeMultiThreshold());
    registerProcessor(new createTestVolume());
    registerProcessor(new createVesselAroundPoints());
    registerProcessor(new foldVessel());
    registerProcessor(new unfoldVessel());

    registerProcessor(new VolumeListLoopInitiator());
    registerProcessor(new VolumeListLoopFinalizer());
}

} // namespace
