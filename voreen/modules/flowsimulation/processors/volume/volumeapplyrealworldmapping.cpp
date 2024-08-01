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

#include "volumeapplyrealworldmapping.h"

#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/volumefactory.h"
#include "voreen/core/datastructures/volume/volumeminmax.h"

namespace voreen {

const std::string VolumeApplyRealWorldMapping::loggerCat_("voreen.base.VolumeApplyRealWorldMapping");

VolumeApplyRealWorldMapping::VolumeApplyRealWorldMapping()
    : CachingVolumeProcessor()
    , inport_(Port::INPORT, "volumehandle.input", "Volume Input")
    , outport_(Port::OUTPORT, "volumehandle.output", "Volume Output", false)
    , enableProcessing_("enableProcessing", "Enable")
{
    addPort(inport_);
    addPort(outport_);

    addProperty(enableProcessing_);
}

VolumeApplyRealWorldMapping::~VolumeApplyRealWorldMapping() {}

Processor* VolumeApplyRealWorldMapping::create() const {
    return new VolumeApplyRealWorldMapping();
}

void VolumeApplyRealWorldMapping::process() {

    const VolumeBase* inputHandle = inport_.getData();

    if (!enableProcessing_.get()) {
        outport_.setData(inputHandle, false);
        return;
    }

    VolumeRAMRepresentationLock inputVolume(inputHandle);
    RealWorldMapping rwm = inputHandle->getRealWorldMapping();

    std::string format = VolumeFactory().getFormat("float", inputVolume->getNumChannels());
    VolumeRAM* output = VolumeFactory().create(format, inputVolume->getDimensions());

    std::vector<float> min(inputVolume->getNumChannels(), std::numeric_limits<float>::max());
    std::vector<float> max(inputVolume->getNumChannels(), std::numeric_limits<float>::lowest());

    for(size_t i=0; i<inputVolume->getNumVoxels(); i++) {
        for(size_t channel=0; channel<inputVolume->getNumChannels(); channel++) {
            float rwmValue = rwm.normalizedToRealWorld(inputVolume->getVoxelNormalized(i, channel));
            min[channel] = std::min(min[channel], rwmValue);
            max[channel] = std::max(max[channel], rwmValue);
            output->setVoxelNormalized(rwmValue, i, channel);
        }
    }

    Volume* outputVolume = new Volume(output, inputHandle);
    outputVolume->setRealWorldMapping(RealWorldMapping());
    outputVolume->addDerivedData(new VolumeMinMax(min, max, min, max));
    outport_.setData(outputVolume);
}

}   // namespace
