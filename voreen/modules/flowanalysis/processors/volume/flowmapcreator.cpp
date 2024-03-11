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

#include "flowmapcreator.h"

#include "voreen/core/datastructures/volume/volumediskmultichanneladapter.h"
#include "voreen/core/ports/conditions/portconditionvolumelist.h"

namespace voreen {

FlowMapCreator::FlowMapCreator()
    : Processor()
    , inport_(Port::INPORT, "volumelist.input", "Volume List Input", false)
    , outport_(Port::OUTPORT, "volumelist.output", "Volume List Output ", false)
    , bypass_("bypass", "Bypass input", false)
    , numChannels_("numChannels", "Num. Channels", 3, 1, 4)
    , layout_("layout", "Layout", Processor::INVALID_RESULT, false, Property::LOD_ADVANCED)
    , mirrorX_("mirrorX", "Mirror X", false)
    , mirrorY_("mirrorY", "Mirror Y", false)
    , mirrorZ_("mirrorZ", "Mirror Z", false)
    , swizzleChannel1_("swizzleChannel1", "Swizzle Channel 1")
    , swizzleChannel2_("swizzleChannel2", "Swizzle Channel 2")
    , swizzleChannel3_("swizzleChannel3", "Swizzle Channel 3")
    , swizzleChannel4_("swizzleChannel4", "Swizzle Channel 4")
    , negateChannel1_("negateChannel1", "Negate Channel 1", false)
    , negateChannel2_("negateChannel2", "Negate Channel 2", false)
    , negateChannel3_("negateChannel3", "Negate Channel 3", false)
    , negateChannel4_("negateChannel4", "Negate Channel 4", false)
{
    addPort(inport_);
    inport_.addCondition(new PortConditionVolumeListEnsemble());
    inport_.Observable<PortObserver>::addObserver(this);
    addPort(outport_);

    addProperty(bypass_);

    addProperty(numChannels_);
    ON_CHANGE(numChannels_, FlowMapCreator, onChannelCountChanged);
    addProperty(layout_);
    layout_.addOption("xyzxyz", "xyzxyz");
    layout_.addOption("xxyyzz", "xxyyzz");

    addProperty(mirrorX_);
    addProperty(mirrorY_);
    addProperty(mirrorZ_);

    addProperty(swizzleChannel1_);
    addProperty(swizzleChannel2_);
    addProperty(swizzleChannel3_);
    addProperty(swizzleChannel4_);

    addProperty(negateChannel1_);
    addProperty(negateChannel2_);
    addProperty(negateChannel3_);
    addProperty(negateChannel4_);

    // Update GUI according to initial state.
    onChannelCountChanged();
}

FlowMapCreator::~FlowMapCreator() {
    // We need to remove the observer, since outport is destructed before inport (stack hierarchy)
    // and the destruction of inport will clear the (already destructed) outport otherwise.
    inport_.Observable<PortObserver>::removeObserver(this);
}

Processor* FlowMapCreator::create() const {
    return new FlowMapCreator();
}

void FlowMapCreator::onChannelCountChanged() {
    //swizzleChannel1_.setVisibleFlag(numChannels_.get() > 0); // always true.
    swizzleChannel2_.setVisibleFlag(numChannels_.get() > 1);
    swizzleChannel3_.setVisibleFlag(numChannels_.get() > 2);
    swizzleChannel4_.setVisibleFlag(numChannels_.get() > 3);

    //negateChannel1_.setVisibleFlag(numChannels_.get() > 0); // always true.
    negateChannel2_.setVisibleFlag(numChannels_.get() > 1);
    negateChannel3_.setVisibleFlag(numChannels_.get() > 2);
    negateChannel4_.setVisibleFlag(numChannels_.get() > 3);

    std::deque<Option<size_t>> options;
    options.push_back(Option<size_t>("x", "x", 0));
    if(numChannels_.get() > 1) {
        options.push_back(Option<size_t>("y", "y", 1));
    }
    if(numChannels_.get() > 2) {
        options.push_back(Option<size_t>("z", "z", 2));
    }
    if(numChannels_.get() > 3) {
        options.push_back(Option<size_t>("w", "w", 3));
    }

    OptionProperty<size_t>* swizzleProperties[] = {&swizzleChannel1_, &swizzleChannel2_, &swizzleChannel3_, &swizzleChannel4_};
    for(size_t propId = 0; propId < static_cast<size_t>(numChannels_.get()); propId++) {
        bool wasSetBefore = !swizzleProperties[propId]->getOptions().empty();
        swizzleProperties[propId]->setOptions(options);
        if(!wasSetBefore) {
            swizzleProperties[propId]->selectByValue(propId);
        }
    }
}

void FlowMapCreator::process() {

    const VolumeList* input = inport_.getData();
    tgtAssert(input, "no input");

    // Clear old data (order matters!).
    outport_.clear();
    volumes_.clear();

    // Bypass input, in case the input is multi-channel volumes.
    if (bypass_.get()) {
        outport_.setData(input, false);
        return;
    }

    size_t numChannels = static_cast<size_t>(numChannels_.get());
    size_t numVolumes = input->size() / numChannels; // floor(x).

    tgt::bvec3 mirror;
    mirror.x = mirrorX_.get();
    mirror.y = mirrorY_.get();
    mirror.z = mirrorZ_.get();

    std::vector<size_t> swizzle;
    std::vector<bool> negate;
    swizzle.push_back(swizzleChannel1_.getValue());
    negate.push_back(negateChannel1_.get());
    if(numChannels > 1) {
        swizzle.push_back(swizzleChannel2_.getValue());
        negate.push_back(negateChannel2_.get());
    }
    if(numChannels > 2) {
        swizzle.push_back(swizzleChannel3_.getValue());
        negate.push_back(negateChannel3_.get());
    }
    if(numChannels > 3) {
        swizzle.push_back(swizzleChannel4_.getValue());
        negate.push_back(negateChannel4_.get());
    }

    VolumeList* output = new VolumeList();

    for (size_t i = 0; i < numVolumes; i++) {

        std::vector<const VolumeBase*> channels;
        if (layout_.get() == "xyzxyz") {
            for (size_t channel = 0; channel < numChannels; channel++) {
                size_t index = i * numChannels + channel;
                channels.push_back(input->at(index));
            }
        }
        else if (layout_.get() == "xxyyzz") {
            for (size_t channel = 0; channel < numChannels; channel++) {
                size_t index = channel * i + numVolumes;
                channels.push_back(input->at(index));
            }
        }
        else {
            tgtAssert(false, "unknown layout");
        }

        VolumeDisk* vd = new VolumeDiskMultiChannelAdapter(channels, mirror, swizzle, negate);
        VolumeBase* volume = new Volume(vd, input->first());
        output->add(volume);

        // Transfer ownership.
        volumes_.push_back(std::unique_ptr<VolumeBase>(volume));
    }

    outport_.setData(output, true);
}

void FlowMapCreator::dataWillChange(const Port* source) {
    tgtAssert(source == &inport_, "unexpected source");

    // Clear old data (order matters!).
    outport_.clear();
    volumes_.clear();
}

void FlowMapCreator::dataHasChanged(const voreen::Port *source) {
    tgtAssert(source == &inport_, "unexpected source");

    // Restore default first.
    bypass_.setReadOnlyFlag(false);

    if (!inport_.isReady()) {
        return;
    }

    const VolumeList* input = inport_.getData();
    if (!input || input->empty()) {
        return;
    }

    if (input->first()->getNumChannels() > 1) {
        bypass_.set(true);
        bypass_.setReadOnlyFlag(true);
        LWARNING("Input volumes have more than a single channel, bypassing...");
    }
}

}   // namespace
