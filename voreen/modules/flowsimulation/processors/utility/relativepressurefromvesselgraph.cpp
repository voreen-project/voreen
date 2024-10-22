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

#include "relativepressurefromvesselgraph.h"

#include "voreen/core/ports/conditions/portconditionvolumetype.h"

#include "../../utils/utils.h"

namespace voreen {

const std::string RelativePressureFromVesselGraph::loggerCat_("voreen.flowsimulation.RelativePressureFromVesselGraph");

RelativePressureFromVesselGraph::RelativePressureFromVesselGraph()
    : Processor()
    , vesselGraphInport_(Port::INPORT, "vesselgraph.inport", "Vessel Graph Input")
    , pressureVolumeInport_(Port::INPORT, "pressurevolume.inport", "Pressure Volume Input")
    , relativePressureOutport_(Port::OUTPORT, "relativepressure.outport", "Relative Pressure Output", false, Processor::VALID)
    , relativePressure_("relativePressure", "Relative Pressure", 0.0f, -10000.0f, 10000.0f, Processor::VALID)
{
    addPort(vesselGraphInport_);
    addPort(pressureVolumeInport_);
    pressureVolumeInport_.addCondition(new PortConditionVolumeChannelCount(1));

    addPort(relativePressureOutport_);

    addProperty(relativePressure_);
    relativePressure_.setNumDecimals(5);
}

bool RelativePressureFromVesselGraph::isReady() const {
    bool ready = Processor::isReady();

    if(ready && vesselGraphInport_.getData()->getEdges().size() != 1) {
        setNotReadyErrorMessage("VesselGraph needs exactly one edge");
        return false;
    }

    return ready;
}

void RelativePressureFromVesselGraph::process() {

    // TODO: set by properties.
    float relativeRadiusCorrection = 1.0f;


    auto* volume = pressureVolumeInport_.getData();

    auto vesselGraph = vesselGraphInport_.getData();
    for(auto& edge : vesselGraph->getEdges()) {

        size_t numVoxels = edge.getVoxels().size();
        if(numVoxels == 0) {
            continue;
        }

        auto samplesP0 = utils::sampleSphere(volume, edge.getNode1().pos_, edge.getNode1().radius_ * relativeRadiusCorrection);
        auto samplesP1 = utils::sampleSphere(volume, edge.getNode2().pos_, edge.getNode2().radius_ * relativeRadiusCorrection);

        auto calcAveragePressure = [](const std::vector<tgt::vec3>& samples) {
            float meanMagnitude = 0.0f;
            for (const auto& sample : samples) {
                meanMagnitude += tgt::length(sample) / static_cast<float>(samples.size());
            }
            return std::vector<float>(1, meanMagnitude);
        };

        auto p0 = calcAveragePressure(samplesP0);
        auto p1 = calcAveragePressure(samplesP1);

        float relativePressure = p1[0] - p0[0];

        relativePressure_.set(relativePressure);

        auto* relativePressureVolumeData = new VolumeRAM_Float(tgt::svec3::one);
        auto* relativePressureVolume = new Volume(relativePressureVolumeData, tgt::vec3::one, tgt::vec3::zero);
        relativePressureOutport_.setData(relativePressureVolume);
    }
}

}   // namespace
