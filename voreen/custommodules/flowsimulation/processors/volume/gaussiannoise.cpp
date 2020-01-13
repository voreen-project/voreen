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

#include "gaussiannoise.h"

#include "voreen/core/datastructures/volume/histogram.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/utils/statistics.h"

namespace voreen {

const std::string GaussianNoise::loggerCat_("voreen.flowsimulation.gaussiannoise");

GaussianNoise::GaussianNoise()
    : VolumeProcessor()
    , inport_(Port::INPORT, "volumenoise.inport", "Volume Input")
    , outport_(Port::OUTPORT, "volumenoise.outport", "Volume Output")
    , enabledProp_("enabledProp", "Enabled", true)
    , mean_("mean", "Mean", 0.0f, -1.0f, 1.0f)
    , stdDeviation_("stdDeviation", "Std. Deviation", 0.2f, 0.0f, 1.0f)
    , usePredeterminedSeed_("usePredeterminedSeed", "Use Predetermined Seed", false)
    , predeterminedSeed_("predeterminedSeed", "Seed", 0, std::numeric_limits<int>::min(), std::numeric_limits<int>::max())
    , randomDevice_()
{
    addPort(inport_);
    addPort(outport_);
    addProperty(enabledProp_);
    addProperty(mean_);
    addProperty(stdDeviation_);
    addProperty(usePredeterminedSeed_);
    addProperty(predeterminedSeed_);
}

GaussianNoise::~GaussianNoise() {
}

Processor* GaussianNoise::create() const {
    return new GaussianNoise();
}

void GaussianNoise::process() {
    if(!enabledProp_.get()) {
        outport_.setData(inport_.getData(), false);
        return;
    }

    const VolumeBase* inputVolume = inport_.getData();
    VolumeRAM* output = inputVolume->getRepresentation<VolumeRAM>()->clone();

    tgt::vec4 means(mean_.get());
    tgt::vec4 stdDevs(stdDeviation_.get());

    std::mt19937 randomEngine(randomDevice_());
    if(usePredeterminedSeed_.get()) {
        randomEngine.seed(predeterminedSeed_.get());
    }

    for(size_t i=0; i<output->getNumVoxels(); i++) {
        for(size_t channel = 0; channel < output->getNumChannels(); channel++) {
            float value = output->getVoxelNormalized(i, channel);
            float noise = std::normal_distribution<float>(means[channel], stdDevs[channel])(randomEngine);
            output->setVoxelNormalized(value + noise, i, channel);
        }
    }

    VolumeBase* outputVolume = new Volume(output, inputVolume);
    outport_.setData(outputVolume);
}

} // namespace voreen
