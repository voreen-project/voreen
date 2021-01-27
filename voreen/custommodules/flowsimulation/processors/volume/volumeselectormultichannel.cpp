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

#include "volumeselectormultichannel.h"

namespace voreen {

VolumeSelectorMultiChannel::VolumeSelectorMultiChannel()
    : Processor()
    , inport_(Port::INPORT, "volumelist.input", "Volume Input", false)
    , volumeOutport_(Port::OUTPORT, "volume.output", "Volume Output", false)
    , volumeOutport2_(Port::OUTPORT, "volume.output2", "Volume Output 2", false)
    , volumeOutport3_(Port::OUTPORT, "volume.output3", "Volume Output 3", false)
    , volumeOutport4_(Port::OUTPORT, "volume.output4", "Volume Output 4", false)
    , numChannels_("numChannels", "Num. Channels", 3, 1, 4)
    , selectedVolume_("selectedVolume", "Selected Volume", 0, 0, std::numeric_limits<int>::max())
    , layout_("layout", "Layout", Processor::INVALID_RESULT, false, Property::LOD_ADVANCED)
{
    addPort(inport_);
    addPort(volumeOutport_);
    addPort(volumeOutport2_);
    addPort(volumeOutport3_);
    addPort(volumeOutport4_);

    addProperty(numChannels_);
    ON_CHANGE(numChannels_, VolumeSelectorMultiChannel, adjustPropertiesToInput);
    addProperty(selectedVolume_);
    addProperty(layout_);
    layout_.addOption("xyzxyz", "xyzxyz");
    layout_.addOption("xxyyzz", "xxyyzz");
}

VolumeSelectorMultiChannel::~VolumeSelectorMultiChannel() {}

Processor* VolumeSelectorMultiChannel::create() const {
    return new VolumeSelectorMultiChannel();
}

void VolumeSelectorMultiChannel::adjustPropertiesToInput() {
    const VolumeList* input = inport_.getData();
    if(!input) {
        return;
    }

    selectedVolume_.setMinValue(0);
    selectedVolume_.setMaxValue(input->size() / numChannels_.get() - 1);
}

void VolumeSelectorMultiChannel::process() {

    const VolumeList* input = inport_.getData();
    tgtAssert(input, "no input");

    VolumePort* channelPorts[] = { &volumeOutport_, &volumeOutport2_, &volumeOutport3_, &volumeOutport4_ };
    for (VolumePort* port : channelPorts) {
        port->clear();
    }

    int numVolumes = static_cast<int>(input->size());
    if(numVolumes < numChannels_.get()) {
        return;
    }

    if(layout_.get() == "xyzxyz") {
        for(int channel = 0; channel < numChannels_.get(); channel++) {
            size_t index = selectedVolume_.get() * numChannels_.get() + channel;
            channelPorts[channel]->setData(input->at(index), false);
        }
    }
    else if(layout_.get() == "xxyyzz") {
        for(int channel = 0; channel < numChannels_.get(); channel++) {
            size_t index = channel * selectedVolume_.get() + selectedVolume_.getMaxValue();
            channelPorts[channel]->setData(input->at(index), false);
        }
    }
    else {
        tgtAssert(false, "unknown layout");
    }
}

}   // namespace
