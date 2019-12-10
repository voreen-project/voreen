/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#include "volumevorticity.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorvorticity.h"

namespace voreen {

const std::string VolumeVorticity::loggerCat_("voreen.base.VolumeVorticity");

VolumeVorticity::VolumeVorticity()
    : CachingVolumeProcessor(),
    inport_(Port::INPORT, "volumehandle.input", "Volume Input"),
    outport_(Port::OUTPORT, "volumehandle.output", "Volume Output", false),
    enableProcessing_("enableProcessing", "Enable"),
    technique_("technique", "Gradient Technique")
{
    addPort(inport_);
    addPort(outport_);

    addProperty(enableProcessing_);
    technique_.addOption("central-differences", "Central differences");
    technique_.addOption("sobel", "Sobel");
    technique_.addOption("linear-regression", "Linear regression");
    addProperty(technique_);
}

VolumeVorticity::~VolumeVorticity() {}

Processor* VolumeVorticity::create() const {
    return new VolumeVorticity();
}

void VolumeVorticity::process() {

    const VolumeBase* inputHandle = inport_.getData();

    if (!enableProcessing_.get()) {
        outport_.setData(inputHandle, false);
        return;
    }

    const VolumeRAM* inputVolume = inputHandle->getRepresentation<VolumeRAM>();
    Volume* outputVolume = 0;

    // expecting a three-channel volume
    if (inputVolume->getNumChannels() == 3) {

        VolumeOperatorVorticity voOpVr;
        if (technique_.get() == "central-differences") {
            if (dynamic_cast<const VolumeRAM_3xUInt16*>(inputHandle->getRepresentation<VolumeRAM>()))
                outputVolume = voOpVr.apply<uint16_t>(inputHandle,VolumeOperatorGradient::VOG_CENTRAL_DIFFERENCE);
            else if (dynamic_cast<const VolumeRAM_3xUInt8*>(inputHandle->getRepresentation<VolumeRAM>()))
                outputVolume = voOpVr.apply<uint8_t>(inputHandle,VolumeOperatorGradient::VOG_CENTRAL_DIFFERENCE);
            else if (dynamic_cast<const VolumeRAM_3xFloat*>(inputHandle->getRepresentation<VolumeRAM>()))
                outputVolume = voOpVr.apply<float>(inputHandle,VolumeOperatorGradient::VOG_CENTRAL_DIFFERENCE);
            else {
                LERROR("Input VolumeRAM type not supported (for this gradient type are supported: 3xUInt16, 3xUInt8, 3xFloat).");
                return;
            }
        } else if (technique_.get() == "sobel") {
            if (dynamic_cast<const VolumeRAM_3xUInt16*>(inputHandle->getRepresentation<VolumeRAM>()))
                outputVolume = voOpVr.apply<uint16_t>(inputHandle,VolumeOperatorGradient::VOG_SOBEL);
            else if (dynamic_cast<const VolumeRAM_3xUInt8*>(inputHandle->getRepresentation<VolumeRAM>()))
                outputVolume = voOpVr.apply<uint8_t>(inputHandle,VolumeOperatorGradient::VOG_SOBEL);
            else {
                LERROR("Input VolumeRAM type not supported (for this gradient type are supported: 3xUInt16, 3xUInt8).");
                return;
            }
        } else if (technique_.get() == "linear-regression") {
            if (dynamic_cast<const VolumeRAM_3xUInt16*>(inputHandle->getRepresentation<VolumeRAM>()))
                outputVolume = voOpVr.apply<uint16_t>(inputHandle,VolumeOperatorGradient::VOG_LINEAR_REGRESSION);
            else if (dynamic_cast<const VolumeRAM_3xUInt8*>(inputHandle->getRepresentation<VolumeRAM>()))
                outputVolume = voOpVr.apply<uint8_t>(inputHandle,VolumeOperatorGradient::VOG_LINEAR_REGRESSION);
            else {
                LERROR("Input VolumeRAM type not supported (for this gradient type are supported: 3xUInt16, 3xUInt8).");
                return;
            }
        } else {
            LERROR("Unknown technique");
            return;
        }
    } else {
        LWARNING("Intensity volume expected, but passed volume consists of " << inputVolume->getNumChannels() << " channels.");
        return;
    }

    outport_.setData(outputVolume);
}

}   // namespace
