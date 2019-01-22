/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2017 University of Muenster, Germany.                        *
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

#include "probabilityvolumecreator.h"

#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"

#include "../ports/conditions/portconditionensemble.h"

namespace voreen {

ProbabilityVolumeCreator::ProbabilityVolumeCreator()
    : AsyncComputeProcessor<ComputeInput, ComputeOutput>()
    , ensembleInport_(Port::INPORT, "ensembledatastructurein", "Ensemble Datastructure Input", false)
    , volumeOutport_(Port::OUTPORT, "volumeout", "Volume Output", false)
    , resampleFactor_("resampleFactor", "Resample Factor", 0.5f, 0.01f, 1.0f)
{
    ensembleInport_.addCondition(new PortConditionEnsembleSingleTimeStep());
    ensembleInport_.addCondition(new PortConditionEnsembleChannelCount(1));
    addPort(ensembleInport_);
    addPort(volumeOutport_);

    addProperty(resampleFactor_);
}

ProbabilityVolumeCreator::~ProbabilityVolumeCreator()
{
}

Processor* ProbabilityVolumeCreator::create() const {
    return new ProbabilityVolumeCreator();
}

ProbabilityVolumeCreatorInput ProbabilityVolumeCreator::prepareComputeInput() {
    const EnsembleDataset* inputPtr = ensembleInport_.getData();
    if (!inputPtr)
        throw InvalidInputException("No input", InvalidInputException::S_WARNING);

    const EnsembleDataset& input = *inputPtr;

    tgt::ivec3 newDims = tgt::vec3(input.getDimensions()) * resampleFactor_.get();
    VolumeRAM_Float* volumeData = new VolumeRAM_Float(newDims, true);

    return ProbabilityVolumeCreatorInput{
        input,
        volumeData,
        resampleFactor_.get()
    };
}
ProbabilityVolumeCreatorOutput ProbabilityVolumeCreator::compute(ProbabilityVolumeCreatorInput input, ProgressReporter& progress) const {

    progress.setProgress(0.0f);

    const std::string& channel = input.dataset.getCommonChannels()[0];

    tgt::vec3 ratio(1.0f / input.resampleFactor);
    tgt::ivec3 newDims = input.volumeData->getDimensions();

    float progressIncrement = 0.95f / (input.dataset.getRuns().size() * newDims.z);

    float max = 0.0f;
    for(const EnsembleDataset::Run& run : input.dataset.getRuns()) {
        progress.setProgressMessage("Calculating Run " + run.name_);

        const VolumeBase* volume = run.timeSteps_[0].channels_.at(channel);
        const VolumeRAM_Float* volumeData = dynamic_cast<const VolumeRAM_Float*>(volume->getRepresentation<VolumeRAM>());

        tgt::vec3 d_a = tgt::vec3(newDims - tgt::ivec3::one) / 2.0f;
        tgt::vec3 d_b = tgt::vec3(volume->getDimensions() - tgt::svec3::one) / 2.0f;

        tgt::ivec3 pos = tgt::ivec3::zero; // iteration variable
        tgt::vec3 nearest; // stores the new position of the target volume

        for (pos.z = 0; pos.z < newDims.z; ++pos.z) {
            nearest.z = (static_cast<float>(pos.z) - d_a.z) * ratio.z + d_b.z;

            for (pos.y = 0; pos.y < newDims.y; ++pos.y) {
                nearest.y = (static_cast<float>(pos.y) - d_a.y) * ratio.y + d_b.y;

                for (pos.x = 0; pos.x < newDims.x; ++pos.x) {
                    nearest.x = (static_cast<float>(pos.x) - d_a.x) * ratio.x + d_b.x;

                    float& voxel = input.volumeData->voxel(pos);
                    voxel += input.dataset.pickSample(volumeData, volume->getSpacing(), nearest);
                    max = std::max(max, voxel);
                }
            }

            // Update progress.
            progress.setProgress(std::min(progress.getProgress() + progressIncrement, 1.0f));
        }
    }

    // Normalize to range [0, 1].
    for(long i=0; i<static_cast<long>(input.volumeData->getNumVoxels()); i++) {
        input.volumeData->voxel(i) /= max;
    }

    progress.setProgress(1.0f);

    Volume* volume = new Volume(input.volumeData, input.dataset.getSpacing(), tgt::vec3::zero);
    return ProbabilityVolumeCreatorOutput{volume};
}

void ProbabilityVolumeCreator::processComputeOutput(ProbabilityVolumeCreatorOutput output) {
    volumeOutport_.setData(output.volume, true);
}

} // namespace
