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

#include "volumelistsubset.h"

namespace voreen {

VolumeListSubset::VolumeListSubset()
    : Processor()
    , inport_(Port::INPORT, "volumelist.input", "Volume Input", false)
    , outport_(Port::OUTPORT, "volumelist.output", "Volume Output", false)
    , numChannels_("numChannels", "Num. Channels", 1, 1, 4)
    , timeStep_("timeStep", "Time Step", 0, 0, std::numeric_limits<int>::max())
    , layout_("layout", "Layout")
{
    addPort(inport_);
    addPort(outport_);
    addProperty(numChannels_);
    ON_CHANGE(numChannels_, VolumeListSubset, adjustPropertiesToInput);

    addProperty(timeStep_);
    addProperty(layout_);
    layout_.addOption("xyzxyz", "xyzxyz");
    layout_.addOption("xxyyzz", "xxyyzz");
}

VolumeListSubset::~VolumeListSubset() {}

Processor* VolumeListSubset::create() const {
    return new VolumeListSubset();
}

void VolumeListSubset::adjustPropertiesToInput() {
    const VolumeList* input = inport_.getData();
    if(!input) {
        return;
    }

    timeStep_.setMaxValue(input->size() / numChannels_.get() - 1);
}

void VolumeListSubset::process() {

    const VolumeList* input = inport_.getData();
    tgtAssert(input, "no input");

    if(input->size() < numChannels_.get()) {
        outport_.clear();
        return;
    }

    VolumeList* output = new VolumeList();

    if(layout_.get() == "xyzxyz") {
        for(int channel = 0; channel < numChannels_.get(); channel++) {
            size_t index = timeStep_.get() * numChannels_.get() + channel;
            output->add(input->at(index));
        }
    }
    else if(layout_.get() == "xxyyzz") {
        for(int channel = 0; channel < numChannels_.get(); channel++) {
            size_t index = channel * timeStep_.get() + timeStep_.getMaxValue();
            output->add(input->at(index));
        }
    }
    else {
        tgtAssert(false, "unknown layout");
    }

    outport_.setData(output, true);
}

}   // namespace
