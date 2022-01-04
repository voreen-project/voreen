/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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
    : AsyncComputeProcessor<VolumeListAggregateInput, VolumeListAggregateOutput>()
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
    aggregationFunction_.addOption("l2norm", "L2 Norm", L2_NORM);
}

VolumeListAggregate::~VolumeListAggregate() {}

Processor* VolumeListAggregate::create() const {
    return new VolumeListAggregate();
}

VolumeListAggregateInput VolumeListAggregate::prepareComputeInput() {

    auto list = inport_.getThreadSafeData();

    if(!list || list->empty()) {
        throw InvalidInputException("No volume input", InvalidInputException::S_WARNING);
    }

    return ComputeInput { std::move(list) };
}

VolumeListAggregateOutput VolumeListAggregate::compute(ComputeInput input, ProgressReporter& progressReporter) const {

    auto data = std::move(input.inputData);

    // Take the first as a reference (we already know we got an ensemble).
    const VolumeBase* reference = data->first();

    std::unique_ptr<Volume> combined(reference->clone());
    VolumeRAM* combinedRepresentation = combined->getWritableRepresentation<VolumeRAM>();

    // For calculating std. dev., we also need a running mean.
    std::unique_ptr<VolumeRAM> tmp;

    AggregationFunction aggregationFunction = aggregationFunction_.getValue();
    switch(aggregationFunction) {
        case MEAN:
        case L2_NORM:
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
        progressReporter.setProgress(1.0f * i / data->size());
        VolumeRAMRepresentationLock representation(data->at(i));
#ifdef VRN_MODULE_OPENMP
#pragma omp parallel for
#endif
        for(long z=0; z < static_cast<long>(dim.z); z++) {
            for(size_t y=0; y < dim.y; y++) {
                for(size_t x=0; x < dim.x; x++) {
                    for(size_t channel=0; channel<combined->getNumChannels(); channel++) {

                        float currentValue = representation->getVoxelNormalized(x, y, z, channel);
                        float aggregatedValue = combinedRepresentation->getVoxelNormalized(x, y, z, channel);

                        switch(aggregationFunction) {
                            case MEAN:
                                aggregatedValue += (currentValue - aggregatedValue) / (i+1);
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
                            case L2_NORM:
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

    if(aggregationFunction == L2_NORM) {
        for(size_t i=0; i<combinedRepresentation->getNumVoxels(); i++) {
            for(size_t channel = 0; channel < combinedRepresentation->getNumChannels(); channel++) {
                float value = combinedRepresentation->getVoxelNormalized(i, channel);
                combinedRepresentation->setVoxelNormalized(std::sqrt(value), i);
            }
        }
    }

    progressReporter.setProgress(1.0f);
    return ComputeOutput { std::move(combined) };
}

void VolumeListAggregate::processComputeOutput(ComputeOutput output) {
    outport_.setData(output.outputData.release(), true);
}

}   // namespace
