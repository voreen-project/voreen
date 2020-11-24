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

#include <modules/core/io/vvdvolumewriter.h>
#include "volumemerger.h"

#include "voreen/core/datastructures/volume/volumefactory.h"

#include "modules/hdf5/io/hdf5volumereader.h"
#include "modules/hdf5/io/hdf5volumewriter.h"
#include "modules/hdf5/io/hdf5filevolume.h"

namespace voreen {

const std::string VolumeMerger::loggerCat_ = "voreen.ensembleanalysis.VolumeMerger";

VolumeMerger::VolumeMerger()
    : AsyncComputeProcessor<ComputeInput, ComputeOutput>()
    , inport_(Port::INPORT, "volumelist.input", "Volume List Input", false)
    , outport_(Port::OUTPORT, "volume.output", "Volume Output", false)
    , allowIntersections_("allowIntersections", "Allow Intersections?", false)
    , padding_("padding", "Padding", 0, 0, 10)
{
    addPort(inport_);
    addPort(outport_);

    addProperty(allowIntersections_);
    addProperty(padding_);
}

VolumeMerger::~VolumeMerger() {
}

Processor* VolumeMerger::create() const {
    return new VolumeMerger();
}

VolumeMergerComputeInput VolumeMerger::prepareComputeInput() {
    if(!inport_.hasData()) {
        throw InvalidInputException("No input", InvalidInputException::S_WARNING);
    }

    auto inputList = inport_.getThreadSafeData();
    if(inputList->empty()) {
        throw InvalidInputException("No volume", InvalidInputException::S_ERROR);
    }

    // Gather reference values.
    const tgt::vec3 spacing = inputList->first()->getSpacing();
    const size_t numChannels = inputList->first()->getNumChannels();
    const RealWorldMapping rwm = inputList->first()->getRealWorldMapping();
    const tgt::vec3 rwPadding(spacing * static_cast<float>(padding_.get()));

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
        if(!allowIntersections_.get()) {
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

    const std::string format = inputList->first()->getFormat();
    const tgt::svec3 dim = tgt::iround(globalBounds.diagonal() / spacing);

    VolumeRAM* outputVolumeData = nullptr;
    try {
        outputVolumeData = VolumeFactory().create(format, dim);
    }
    catch (const std::bad_alloc& ) {
        throw InvalidInputException("Could not create output volume.", InvalidInputException::S_ERROR);
    }

    std::unique_ptr<Volume> outputVolume(new Volume(outputVolumeData, spacing, globalBounds.getLLF()));
    outputVolume->setRealWorldMapping(rwm);

    return VolumeMergerComputeInput{
            std::move(inputList),
            padding_.get(),
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
        tgt::mat4 voxelToWorldMatrix = volume->getVoxelToWorldMatrix();

#ifdef VRN_MODULE_OPENMP
        #pragma omp parallel for
#endif
        for(long z=padding; z<dim.z-padding; z++) {
            for(long y=padding; y<dim.y-padding; y++) {
                for(long x=padding; x<dim.x-padding; x++) {
                    for(size_t channel=0; channel < volumeData->getNumChannels(); channel++) {
                        float value = volumeData->getVoxelNormalized(x, y, z, channel);
                        tgt::vec3 pos = worldToVoxelMatrix * (voxelToWorldMatrix * tgt::vec3(x, y, z));
                        outputVolumeData->setVoxelNormalized(value, pos, channel);
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
