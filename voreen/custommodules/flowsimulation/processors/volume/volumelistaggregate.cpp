/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#include "volumelistaggregate.h"

#include "voreen/core/ports/conditions/portconditionvolumelist.h"

namespace voreen {

VolumeListAggregate::VolumeListAggregate()
    : Processor()
    , inport_(Port::INPORT, "volumelist.input", "Volume Input", false)
    , outport_(Port::OUTPORT, "volume.output", "Volume Output", false)
    , aggregationFunction_("aggregationFunction", "Aggregation Function")
{
    addPort(inport_);
    inport_.addCondition(new PortConditionVolumeListEnsemble());
    addPort(outport_);

    addProperty(aggregationFunction_);
    aggregationFunction_.addOption("mean", "Mean", MEAN);
    aggregationFunction_.addOption("min", "Min", MIN);
    aggregationFunction_.addOption("max", "Max", MAX);
    aggregationFunction_.addOption("variance", "Variance", VARIANCE);
    aggregationFunction_.addOption("sumOfSquares", "Sum of Squares", SUM_OF_SQUARES);
}

VolumeListAggregate::~VolumeListAggregate() {}

Processor* VolumeListAggregate::create() const {
    return new VolumeListAggregate();
}

void VolumeListAggregate::process() {
    const VolumeList* data = inport_.getData();
    if(!data || data->empty()) {
        outport_.setData(nullptr);
    } else if(data->size() == 1) {
        outport_.setData(data->first(), false);
    } else {
        // Take the first as a reference (we already know we got an ensemble).
        const VolumeBase* reference = data->first();

        Volume* combined = reference->clone();
        VolumeRAM* combinedRepresentation = combined->getWritableRepresentation<VolumeRAM>();

        // For calculating std. dev., we also need a running mean.
        std::unique_ptr<VolumeRAM> tmp;

        AggregationFunction aggregationFunction = aggregationFunction_.getValue();
        switch(aggregationFunction) {
        case MEAN:
        case SUM_OF_SQUARES:
            combinedRepresentation->clear();
            break;
        case MIN:
        case MAX:
            // Use reference values as initialization.
            break;
        case VARIANCE:
            combinedRepresentation->clear();
            tmp.reset(combinedRepresentation->clone());
            break;
        default:
            tgtAssert(false, "unhandled aggregation function");
        }

        tgt::svec3 dim = combined->getDimensions();
        for(size_t i=0; i<data->size(); i++) {
            VolumeRAMRepresentationLock representation(data->at(i));
            for(size_t z=0; z < dim.z; z++) {
                for(size_t y=0; y < dim.y; y++) {
                    for(size_t x=0; x < dim.x; x++) {
                        for(size_t channel=0; channel<combined->getNumChannels(); channel++) {

                            float currentValue = representation->getVoxelNormalized(x, y, z, channel);
                            float aggregatedValue = combinedRepresentation->getVoxelNormalized(x, y, z, channel);

                            switch(aggregationFunction) {
                            case MEAN:
                                aggregatedValue += (currentValue - aggregatedValue) / i;
                                break;
                            case MIN:
                                aggregatedValue = std::min(aggregatedValue, currentValue);
                                break;
                            case MAX:
                                aggregatedValue = std::max(aggregatedValue, currentValue);
                                break;
                            case VARIANCE: {
                                // Running mean.
                                float runningMean = tmp->getVoxelNormalized(x, y, z, channel);
                                float delta = currentValue - runningMean;
                                runningMean += delta / (i+1.0f);
                                tmp->setVoxelNormalized(runningMean, x, y, z, channel);

                                if(i > 2) {
                                    aggregatedValue += delta * (currentValue - runningMean);
                                }
                                break;
                            }
                            case SUM_OF_SQUARES:
                                aggregatedValue += currentValue * currentValue;
                                break;
                            default:
                                tgtAssert(false, "unhandled aggregation function");
                            }

                            combinedRepresentation->setVoxelNormalized(aggregatedValue, x, y, z, channel);
                        }
                    }
                }
            }
        }

        outport_.setData(combined, true);
    }
}

}   // namespace
