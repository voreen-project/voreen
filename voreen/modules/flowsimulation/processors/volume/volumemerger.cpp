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

#include <modules/core/io/vvdvolumewriter.h>
#include "volumemerger.h"

#include "voreen/core/datastructures/volume/volumefactory.h"

namespace voreen {

const std::string VolumeMerger::loggerCat_ = "voreen.ensembleanalysis.VolumeMerger";

VolumeMerger::VolumeMerger()
    : AsyncComputeProcessor<ComputeInput, ComputeOutput>()
    , inport_(Port::INPORT, "volumelist.input", "Volume List Input", false)
    , outport_(Port::OUTPORT, "volume.output", "Volume Output", false)
    , intersectionResolutionStrategy_("intersectionResolution", "Intersections Resolution")
    , padding_("padding", "Padding", 0, 0, 10)
{
    addPort(inport_);
    addPort(outport_);

    addProperty(intersectionResolutionStrategy_);
    intersectionResolutionStrategy_.addOption("none", "None", IntersectionResolutionStrategy::IRS_NONE);
    intersectionResolutionStrategy_.addOption("last", "Last", IntersectionResolutionStrategy::IRS_LAST);
    intersectionResolutionStrategy_.addOption("max", "Max", IntersectionResolutionStrategy::IRS_MAX);
    intersectionResolutionStrategy_.addOption("min", "Min", IntersectionResolutionStrategy::IRS_MIN);
    intersectionResolutionStrategy_.addOption("avg", "Avg", IntersectionResolutionStrategy::IRS_AVG);
    addProperty(padding_);
}

VolumeMerger::~VolumeMerger() {
}

Processor* VolumeMerger::create() const {
    return new VolumeMerger();
}

void VolumeMerger::setPadding(int padding) {
    padding_.set(padding);
}

int VolumeMerger::getPadding() const {
    return padding_.get();
}

void VolumeMerger::setIntersectionResolutionStrategy(VolumeMerger::IntersectionResolutionStrategy strategy) {
    intersectionResolutionStrategy_.selectByValue(strategy);
}

VolumeMerger::IntersectionResolutionStrategy VolumeMerger::getIntersectionResolutionStrategy() const {
    return intersectionResolutionStrategy_.getValue();
}

VolumeMergerComputeInput VolumeMerger::prepareComputeInput() {
    if(!inport_.hasData()) {
        throw InvalidInputException("No input", InvalidInputException::S_WARNING);
    }

    auto inputList = inport_.getThreadSafeData();
    if(inputList->empty()) {
        throw InvalidInputException("No volume", InvalidInputException::S_ERROR);
    }

    auto* referenceVolume = inputList->first();

    // Gather reference values.
    const tgt::vec3 spacing = referenceVolume->getSpacing();
    const size_t numChannels = referenceVolume->getNumChannels();
    const RealWorldMapping rwm = referenceVolume->getRealWorldMapping();
    const tgt::vec3 rwPadding(spacing * static_cast<float>(padding_.get()));

    const bool intersectionsAllowed = getIntersectionResolutionStrategy() != IRS_NONE;

    tgt::Bounds globalBounds;
    std::vector<tgt::Bounds> processedBounds;
    for(size_t i=0; i < inputList->size(); i++) {

        if (inputList->at(i)->getSpacing() != spacing) {
            throw InvalidInputException("Spacing must match", InvalidInputException::S_ERROR);
        }

        if(inputList->at(i)->getNumChannels() != numChannels) {
            throw InvalidInputException("Num Channels must match", InvalidInputException::S_ERROR);
        }

        if(inputList->at(i)->getRealWorldMapping() != rwm) {
            throw InvalidInputException("Real World Mapping must match", InvalidInputException::S_ERROR);
        }

        // TODO: add more checks

        tgt::Bounds localBounds(inputList->at(i)->getLLF() + rwPadding, inputList->at(i)->getURB() - rwPadding);
        if(!intersectionsAllowed) {
            for (size_t j = 0; j < processedBounds.size(); j++) {
                if (processedBounds[j].intersects(localBounds)) {
                    throw InvalidInputException("Volumes must not intersect", InvalidInputException::S_ERROR);
                }
            }
            processedBounds.push_back(localBounds);
        }

        globalBounds.addVolume(localBounds);
    }

    outport_.setData(nullptr);

    const std::string format = referenceVolume->getFormat();
    const tgt::svec3 dim = tgt::iround(globalBounds.diagonal() / spacing);

    VolumeRAM* outputVolumeData = nullptr;
    try {
        outputVolumeData = VolumeFactory().create(format, dim);
    }
    catch (const std::bad_alloc& ) {
        throw InvalidInputException("Could not create output volume.", InvalidInputException::S_ERROR);
    }

    // We take all meta data values from the reference volume, but need to take care that the transformation is reset.
    auto offset = globalBounds.getLLF() + spacing * 0.5f; // We have to account for the extra voxel.
    std::unique_ptr<Volume> outputVolume(new Volume(outputVolumeData, referenceVolume));
    outputVolume->setSpacing(spacing);
    outputVolume->setOffset(offset);
    outputVolume->setPhysicalToWorldMatrix(tgt::mat4::identity);
    outputVolume->setRealWorldMapping(rwm);

    std::function<float(float, float)> collisionFunction;
    switch (getIntersectionResolutionStrategy()) {
    case IRS_LAST:
        collisionFunction = [] (float lhs, float rhs) {
            return rhs;
        };
        break;

    case IRS_MAX:
        collisionFunction = [] (float lhs, float rhs) {
            return std::max(lhs, rhs);
        };
        break;

    case IRS_MIN:
        collisionFunction = [] (float lhs, float rhs) {
            return std::min(lhs, rhs);
        };
        break;

    case IRS_AVG:
        collisionFunction = [] (float lhs, float rhs) {
            return (lhs + rhs ) / 2;
        };
        break;

    case IRS_NONE:
    default:
        collisionFunction = [] (float lhs, float rhs) {
            return 0.0f; // Should not happen!
        };
        break;
    }

    return VolumeMergerComputeInput{
            std::move(inputList),
            padding_.get(),
            collisionFunction,
            std::move(outputVolume)
    };
}
VolumeMergerComputeOutput VolumeMerger::compute(VolumeMergerComputeInput input, ProgressReporter& progressReporter) const {

    const VolumeList* inputVolumes = input.inputVolumes;
    std::unique_ptr<Volume> outputVolume = std::move(input.outputVolume);
    tgt::mat4 worldToVoxelMatrix = outputVolume->getWorldToVoxelMatrix();
    const int padding = input.padding_;

    VolumeRAM* outputVolumeData = outputVolume->getWritableRepresentation<VolumeRAM>();
    outputVolumeData->clear();

    float progressPerVolume = 1.0f / inputVolumes->size();
    for(size_t i=0; i<inputVolumes->size(); i++) {
        progressReporter.setProgress(i * progressPerVolume);

        const VolumeBase* volume = inputVolumes->at(i);
        VolumeRAMRepresentationLock volumeData(volume);

        const tgt::Vector3<long> dim = volumeData->getDimensions();
        tgt::mat4 mat = worldToVoxelMatrix * volume->getVoxelToWorldMatrix();

#ifdef VRN_MODULE_OPENMP
        #pragma omp parallel for
#endif
        for(long z=padding; z<dim.z-padding; z++) {
            for(long y=padding; y<dim.y-padding; y++) {
                for(long x=padding; x<dim.x-padding; x++) {
                    tgt::vec3 pos = mat * tgt::vec3(x, y, z);
                    for(size_t channel=0; channel < volumeData->getNumChannels(); channel++) {
                        float value = volumeData->getVoxelNormalized(x, y, z, channel);
                        float outputValue = outputVolumeData->getVoxelNormalized(pos, channel);
                        outputValue = input.collisionFunction_(outputValue, value);
                        outputVolumeData->setVoxelNormalized(outputValue, pos, channel);
                    }
                }
            }
        }
    }

    progressReporter.setProgress(1.0f);

    return {
        std::move(outputVolume)
    };
}
void VolumeMerger::processComputeOutput(VolumeMergerComputeOutput output) {
    outport_.setData(output.outputVolume.release());
}

}   // namespace
