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

#include "volumenoise.h"

#include "voreen/core/datastructures/volume/volumeram.h"

namespace voreen {

const std::string VolumeNoise::loggerCat_("voreen.flowsimulation.VolumeNoise");

VolumeNoise::VolumeNoise()
    : VolumeProcessor()
    , inport_(Port::INPORT, "volumenoise.inport", "Volume Input")
    , outport_(Port::OUTPORT, "volumenoise.outport", "Volume Output")
    , enabledProp_("enabledProp", "Enabled", true)
    , noiseClass_("noiseClass", "Noise Class")
    , noiseAmount_("noiseAmount", "noiseAmount", 0.5f, 0.0f, 1.0f)
    , usePredeterminedSeed_("usePredeterminedSeed", "Use Predetermined Seed", false)
    , predeterminedSeed_("predeterminedSeed", "Seed", 0, std::numeric_limits<int>::min(), std::numeric_limits<int>::max())
    , randomDevice_()
{
    addPort(inport_);
    addPort(outport_);
    addProperty(enabledProp_);
    addProperty(noiseClass_);
    noiseClass_.addOption("white", "White", WHITE);
    noiseClass_.addOption("gaussian", "Gaussian", GAUSSIAN);
    noiseClass_.addOption("saltAndPepper", "Salt and Pepper", SALT_AND_PEPPER);
    addProperty(noiseAmount_);
    addProperty(usePredeterminedSeed_);
    addProperty(predeterminedSeed_);
}

VolumeNoise::~VolumeNoise() {
}

Processor* VolumeNoise::create() const {
    return new VolumeNoise();
}

void VolumeNoise::process() {
    if(!enabledProp_.get()) {
        outport_.setData(inport_.getData(), false);
        return;
    }

    const VolumeBase* inputVolume = inport_.getData();
    VolumeRAM* output = inputVolume->getRepresentation<VolumeRAM>()->clone();
    NoiseClass noiseClass = noiseClass_.getValue();

    std::mt19937 randomEngine(randomDevice_());
    if(usePredeterminedSeed_.get()) {
        randomEngine.seed(predeterminedSeed_.get());
    }

    auto applyNoise = [&] (float value) {
        if(noiseClass == WHITE) {
            value += std::uniform_real_distribution<float>(0.0f, noiseAmount_.get())(randomEngine);
        }
        else if(noiseClass == GAUSSIAN) {
            value += std::normal_distribution<float>(0.0f, noiseAmount_.get()*10)(randomEngine);
        }
        else if(noiseClass == SALT_AND_PEPPER) {
            if(std::uniform_real_distribution<float>()(randomEngine) <= noiseAmount_.get()) {
                // Flip a coin to determine salt or pepper.
                value = std::round(std::uniform_real_distribution<float>()(randomEngine));
            }
        }
        return value;
    };

    for(size_t i=0; i<output->getNumVoxels(); i++) {
        for(size_t channel = 0; channel < output->getNumChannels(); channel++) {
            float value = output->getVoxelNormalized(i, channel);
            value = applyNoise(value);
            output->setVoxelNormalized(value, i, channel);
        }
    }

    VolumeBase* outputVolume = new Volume(output, inputVolume);
    outport_.setData(outputVolume);
}

} // namespace voreen
