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

#include "phaseunwrapping.h"

#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/ports/conditions/portconditionvolumetype.h"

#include "../../ext/unwrap3d/unwrap3d.h"

namespace voreen {


void calculatePhaseMap(const VolumeRAM* input, VolumeRAM* output, const unsigned char* mask, size_t channel, const RealWorldMapping& rwm, float VENC) {
    std::unique_ptr<float []> wrappedPhaseMap(new float [input->getNumVoxels()]);
    std::unique_ptr<float []> unwrappedPhaseMap(new float [input->getNumVoxels()]);

    for(size_t i=0; i<input->getNumVoxels(); i++) {
        float value = rwm.normalizedToRealWorld(input->getVoxelNormalized(i, channel));
        wrappedPhaseMap[i] = value / VENC * tgt::PIf;
    }

    tgt::ivec3 dim = input->getDimensions();

    unwrap3d(wrappedPhaseMap.get(), unwrappedPhaseMap.get(), mask, dim.x, dim.y, dim.z);

    for(size_t i=0; i<input->getNumVoxels(); i++) {
        float value = rwm.realWorldToNormalized(unwrappedPhaseMap[i] / tgt::PIf * VENC);
        output->setVoxelNormalized(value, i, channel);
    }
}


const std::string PhaseUnwrapping::loggerCat_("voreen.flowsimulation.PhaseUnwrapping");

PhaseUnwrapping::PhaseUnwrapping()
    : VolumeProcessor()
    , inport_(Port::INPORT, "phaseunwrapping.inport", "Volume Input")
    , mask_(Port::INPORT, "phaseunwrapping.mask", "Mask")
    , outport_(Port::OUTPORT, "phaseunwrapping.outport", "Volume Output")
    , enabledProp_("enabledProp", "Enabled", true)
    , venc_("venc", "VENC", 1.0f, 1.0f, 600.0f)
{
    addPort(inport_);
    addPort(mask_);
    mask_.addCondition(new PortConditionVolumeTypeUInt8());
    addPort(outport_);
    addProperty(enabledProp_);
    addProperty(venc_);
}

PhaseUnwrapping::~PhaseUnwrapping() {
}

Processor* PhaseUnwrapping::create() const {
    return new PhaseUnwrapping();
}

bool PhaseUnwrapping::isReady() const {
    if(!inport_.isReady()) {
        setNotReadyErrorMessage("No input");
        return false;
    }

    // Mask is optional
    if(mask_.isReady() && mask_.getData()->getDimensions() != inport_.getData()->getDimensions()) {
        setNotReadyErrorMessage("Mask has different dimensions than input volume");
        return false;
    }

    if(!outport_.isReady()) {
        setNotReadyErrorMessage("outport not connected");
        return false;
    }

    return true;
}

void PhaseUnwrapping::process() {
    if(!enabledProp_.get()) {
        outport_.setData(inport_.getData(), false);
        return;
    }

    const VolumeBase* inputVolume = inport_.getData();
    RealWorldMapping rwm = inputVolume->getRealWorldMapping();
    VolumeRAMRepresentationLock input(inputVolume);
    VolumeRAM* output = input->clone();

    const unsigned char* mask = nullptr;
    std::unique_ptr<VolumeRAMRepresentationLock> maskLock;
    if(mask_.hasData()) {
        maskLock.reset(new VolumeRAMRepresentationLock(mask_.getData()));
        mask = reinterpret_cast<const unsigned char*>((*maskLock)->getData());
    }
    else {
        unsigned char* no_mask = new unsigned char[input->getNumVoxels()];
        std::fill(no_mask, no_mask+input->getNumVoxels(), 1);
        mask = no_mask;
    }

    for(size_t channel=0; channel<input->getNumChannels(); channel++) {
        calculatePhaseMap(*input, output, mask, channel, rwm, venc_.get());
    }

    if(!mask_.hasData()) {
        delete [] mask;
    }

    VolumeBase* outputVolume = new Volume(output, inputVolume);
    outport_.setData(outputVolume);
}

} // namespace voreen
