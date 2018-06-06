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

#include "volumegradient.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorgradient.h"

namespace voreen {

const std::string VolumeGradient::loggerCat_("voreen.base.VolumeGradient");

VolumeGradient::VolumeGradient()
    : CachingVolumeProcessor(),
    inport_(Port::INPORT, "volumehandle.input", "Volume Input"),
    outport_(Port::OUTPORT, "volumehandle.output", "Volume Output", false),
    enableProcessing_("enableProcessing", "Enable"),
    technique_("technique", "Technique")
{
    addPort(inport_);
    addPort(outport_);

    addProperty(enableProcessing_);
    technique_.addOption("central-differences", "Central differences");
    technique_.addOption("sobel", "Sobel");
    technique_.addOption("linear-regression", "Linear regression");
    addProperty(technique_);
}

VolumeGradient::~VolumeGradient() {}

Processor* VolumeGradient::create() const {
    return new VolumeGradient();
}

void VolumeGradient::process() {

    const VolumeBase* inputHandle = inport_.getData();

    if (!enableProcessing_.get()) {
        outport_.setData(inputHandle, false);
        return;
    }

    const VolumeRAM* inputVolume = inputHandle->getRepresentation<VolumeRAM>();
    Volume* outputVolume = 0;

    // expecting a single-channel volume
    if (inputVolume->getNumChannels() == 1) {

        bool bit16 = inputVolume->getBitsAllocated() > 8;
        VolumeOperatorGradient voOpGr;
        if (technique_.get() == "central-differences") {
            if (bit16)
                outputVolume = voOpGr.apply<uint16_t>(inputHandle,VolumeOperatorGradient::VOG_CENTRAL_DIFFERENCE);
            else if (!bit16)
                outputVolume = voOpGr.apply<uint8_t>(inputHandle,VolumeOperatorGradient::VOG_CENTRAL_DIFFERENCE);
            else {
                tgtAssert(false, "Should not get here");
            }
        }
        else if (technique_.get() == "sobel") {
            if (bit16)
                outputVolume = voOpGr.apply<uint16_t>(inputHandle,VolumeOperatorGradient::VOG_SOBEL);
            else if (!bit16)
                outputVolume = voOpGr.apply<uint8_t>(inputHandle,VolumeOperatorGradient::VOG_SOBEL);
            else {
                tgtAssert(false, "Should not get here");
            }
        }
        else if (technique_.get() == "linear-regression") {
            if (bit16)
                outputVolume = voOpGr.apply<uint16_t>(inputHandle,VolumeOperatorGradient::VOG_LINEAR_REGRESSION);
            else if (!bit16)
                outputVolume = voOpGr.apply<uint8_t>(inputHandle,VolumeOperatorGradient::VOG_LINEAR_REGRESSION);
            else {
                tgtAssert(false, "Should not get here");
            }
        }
        else {
            LERROR("Unknown technique");
        }
    }
    else {
        LWARNING("Intensity volume expected, but passed volume consists of " << inputVolume->getNumChannels() << " channels.");
    }

    outport_.setData(outputVolume);
}

}   // namespace
