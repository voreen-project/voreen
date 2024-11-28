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

#include "vesselgraphtoflowindicator.h"

namespace voreen {

const std::string VesselGraphToFlowIndicator::loggerCat_("voreen.flowsimulation.VesselGraphToFlowIndicator");

VesselGraphToFlowIndicator::VesselGraphToFlowIndicator()
    : Processor()
    , flowParametrizationInport_(Port::INPORT, "flowparametrization.inport", "Flow Parametrization Input")
    , flowParametrizationOutport_(Port::OUTPORT, "flowparametrization.outport", "Flow Parametrization Output")
    , vesselGraphInport_(Port::INPORT, "vesselgraph.inport", "Vessel Graph Input")
    , maskPort_(Port::INPORT, "mask.inport", "(Optional) Mask Input")
{
    addPort(flowParametrizationInport_);
    addPort(flowParametrizationOutport_);
    addPort(vesselGraphInport_);
    addPort(maskPort_);
}

void VesselGraphToFlowIndicator::process() {

    // TODO: set by properties.
    float relativeRadiusCorrection = 1.0f;
    int numberOfReferenceNodes = 2;
    int padding = 1;

    // Add indicators.
    auto* config = new FlowSimulationConfig(*flowParametrizationInport_.getData());

    // Setup mask.
    std::function<bool(tgt::vec3)> insideMask;
    auto* maskVolume = maskPort_.getData();
    if(maskVolume) {
        tgt::mat4 worldToVoxel = maskVolume->getWorldToVoxelMatrix();
        VolumeRAMRepresentationLock lock(maskVolume);
        insideMask = [=](const tgt::vec3& position) {
            return lock->getVoxelNormalized(worldToVoxel * position) != 0.0f;
        };
    }
    else {
        insideMask = [](tgt::vec3) { return true; };
    }

    // Add flow indicators.
    auto vesselGraph = vesselGraphInport_.getData();
    for(auto& edge : vesselGraph->getEdges()) {

        size_t numVoxels = edge.getVoxels().size();
        if(numVoxels == 0) {
            continue;
        }

        for(size_t i=padding; i<numVoxels-padding; i++) {

            size_t mid = i;
            size_t num = numberOfReferenceNodes; // Number of reference nodes in both directions.

            size_t frontIdx = mid > num ? (mid - num) : 0;
            size_t backIdx  = std::min(mid + num, numVoxels - 1);

            const VesselSkeletonVoxel* ref   = &edge.getVoxels().at(mid);
            const VesselSkeletonVoxel* front = &edge.getVoxels().at(frontIdx);
            const VesselSkeletonVoxel* back  = &edge.getVoxels().at(backIdx);

            // Calculate average radius.
            float radius = 0.0f;
            for (size_t j = frontIdx; j <= backIdx; j++) {
                radius += edge.getVoxels().at(j).avgDistToSurface_;
            }
            radius /= (backIdx - frontIdx + 1);

            // Default is candidate. However, we set initial values as if it was a velocity inlet,
            // since it might be classified as one later on.
            FlowIndicator indicator;
            // indicator.id_ = -1; // Will be set when adding to the config.
            indicator.name_ = "Indicator " + std::to_string(config->getFlowIndicators().size());
            indicator.type_ = FlowIndicatorType::FIT_MEASURE;
            indicator.flowProfile_ = FlowProfile::FP_NONE;

            indicator.center_ = ref->pos_;
            indicator.normal_ = tgt::normalize(back->pos_ - front->pos_);
            indicator.length_ = tgt::distance(front->pos_, back->pos_);
            indicator.radius_ = relativeRadiusCorrection * radius;

            if(insideMask(indicator.center_)) {
                config->addFlowIndicator(indicator);
            }
        }
    }

    flowParametrizationOutport_.setData(config);
}

bool VesselGraphToFlowIndicator::isReady() const {

    if (!isInitialized()) {
        setNotReadyErrorMessage("Not initialized");
        return false;
    }

    if (!vesselGraphInport_.isReady()) {
        setNotReadyErrorMessage("No vessel graph");
        return false;
    }

    if (!flowParametrizationInport_.isReady()) {
        setNotReadyErrorMessage("No flow parametrization");
        return false;
    }

    // Mask is optional.

    return true;
}

}   // namespace
