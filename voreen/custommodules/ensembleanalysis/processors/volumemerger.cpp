/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

VolumeMerger::VolumeMerger()
    : AsyncComputeProcessor<ComputeInput, ComputeOutput>()
    , inport_(Port::INPORT, "volumelist.input", "Volume List Input", false)
    , outport_(Port::OUTPORT, "volume.output", "Volume Output", false)
    , allowIntersections_("allowIntersections", "Allow Intersections?", false)
{
    addPort(inport_);
    addPort(outport_);

    addProperty(allowIntersections_);
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

    auto inputVolPtr = inport_.getThreadSafeData();
    const VolumeList& inputVolumes = *inputVolPtr;

    if(inputVolumes.empty()) {
        throw InvalidInputException("No volume", InvalidInputException::S_ERROR);
    }

    // Gather reference values.
    const tgt::vec3& spacing = inputVolumes.first()->getSpacing();
    const size_t numChannels = inputVolumes.first()->getNumChannels();

    tgt::Bounds globalBounds;
    std::vector<tgt::Bounds> processedBounds;
    for(size_t i=0; i<inputVolumes.size(); i++) {

        if (inputVolumes.at(i)->getSpacing() != spacing) {
            throw InvalidInputException("Spacing must match", InvalidInputException::S_ERROR);
        }

        if(inputVolumes.at(i)->getNumChannels() != numChannels) {
            throw InvalidInputException("Num Channels must match", InvalidInputException::S_ERROR);
        }

        // TODO: add more checks.

        tgt::Bounds localBounds(inputVolumes.at(i)->getLLF(), inputVolumes.at(i)->getURB());

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

    const std::string format = inputVolumes.first()->getFormat();
    const tgt::svec3 dim = globalBounds.diagonal() / spacing;

    std::unique_ptr<VolumeRAM> outputVolume = nullptr;
    try {
        outputVolume.reset(VolumeFactory().create(format, dim));
    }
    catch (const std::bad_alloc& ) {
        throw InvalidInputException("Could not create output volume.", InvalidInputException::S_ERROR);
    }

    return VolumeMergerComputeInput{
            inputVolPtr,
            spacing,
            globalBounds.getLLF(),
            std::move(outputVolume)
    };
}
VolumeMergerComputeOutput VolumeMerger::compute(VolumeMergerComputeInput input, ProgressReporter& progressReporter) const {
    tgtAssert(input.inputVolumes, "No inputVolumes");
    tgtAssert(input.outputVolume, "No outputVolume");

    const VolumeList* inputVolumes = input.inputVolumes;
    std::unique_ptr<VolumeRAM> outputVolume = std::move(input.outputVolume);
    outputVolume->clear();

    float progressPerVolume = 1.0f / inputVolumes->size();
    for(size_t i=0; i<inputVolumes->size(); i++) {
        progressReporter.setProgress(i * progressPerVolume);

        const VolumeBase* volume = inputVolumes->at(i);
        const VolumeRAM* volumeData = volume->getRepresentation<VolumeRAM>();

        const tgt::svec3& dim = volumeData->getDimensions();
        const tgt::svec3 offset((volume->getLLF() - input.offset) / input.spacing);

#ifdef VRN_MODULE_OPENMP
        long dimZ = static_cast<long>(dim.z);
        #pragma omp parallel for
        for(long pz = 0; pz < dimZ; pz++) {
            size_t z = static_cast<size_t>(pz);
#else
        for(size_t z=0; z<dim.z; z++) {
#endif
            for(size_t y=0; y<dim.y; y++) {
                for(size_t x=0; x<dim.x; x++) {
                    for(size_t channel=0; channel < volumeData->getNumChannels(); channel++) {
                        float value = volumeData->getVoxelNormalized(x, y, z, channel);
                        outputVolume->setVoxelNormalized(value, x+offset.x, y+offset.y, z+offset.z, channel);
                    }
                }
            }
        }
    }

    progressReporter.setProgress(1.0f);

    return {
        input.offset,
        input.spacing,
        std::move(outputVolume)
    };
}
void VolumeMerger::processComputeOutput(VolumeMergerComputeOutput output) {
    Volume* vol =  new Volume(output.outputVolume.release(), output.spacing, output.offset);
    outport_.setData(vol);
}

}   // namespace
