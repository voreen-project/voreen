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

#include "volumelistmultichanneladapter.h"

#include "voreen/core/datastructures/volume/volumedisk.h"
#include "voreen/core/datastructures/volume/volumefactory.h"
#include "voreen/core/datastructures/volume/volumehash.h"
#include "voreen/core/ports/conditions/portconditionvolumelist.h"

namespace voreen {

class VolumeDisk_MultiChannelAdapter : public VolumeDisk {
public:

    VolumeDisk_MultiChannelAdapter(const std::vector<const VolumeBase*>& channels)
        : VolumeDisk(VolumeFactory().getFormat(channels.front()->getBaseType(), channels.size()), channels.front()->getDimensions())
        , channels_(channels)
    {
        const VolumeBase* ref = channels_.front();
        for (const VolumeBase* channel : channels_) {
            tgtAssert(ref->getFormat() == channel->getFormat(), "Base Type mismatch");
            tgtAssert(ref->getDimensions() == channel->getDimensions(), "Base Type mismatch");
        }
    }

    virtual std::string getHash() const {
        std::string hash;

        for (const VolumeBase* channel : channels_) {
            hash += channel->getHash();
        }

        return VolumeHash(hash).getHash();
    }

    virtual VolumeRAM* loadVolume() const {
        return loadBrick(tgt::svec3::zero, dimensions_);
    }

    virtual VolumeRAM* loadSlices(const size_t firstZSlice, const size_t lastZSlice) const {
        if (firstZSlice > lastZSlice)
            throw VoreenException("last slice must be behind first slice");

        return loadBrick(tgt::svec3(0, 0, firstZSlice), tgt::svec3(dimensions_.x, dimensions_.y, lastZSlice - firstZSlice + 1));
    }

    virtual VolumeRAM* loadBrick(const tgt::svec3& offset, const tgt::svec3& dimensions) const {
        // check parameters
        if (tgt::hmul(dimensions) == 0)
            throw VoreenException("requested brick dimensions are zero");
        if (!tgt::hand(tgt::lessThanEqual(offset + dimensions, dimensions_)))
            throw VoreenException("requested brick (at least partially) outside volume dimensions");

        // Create the output volume.
        VolumeRAM* output = VolumeFactory().create(getFormat(), dimensions);

        for (size_t channel = 0; channel < channels_.size(); channel++) {
            // Check if we have a ram representation already.
            if (channels_[channel]->hasRepresentation<VolumeRAM>()) {
                VolumeRAMRepresentationLock lock(channels_[channel]);
                tgt::svec3 pos;
                for (pos.z = 0; pos.z < dimensions.z; pos.z++) {
                    for (pos.y = 0; pos.y < dimensions.y; pos.y++) {
                        for (pos.x = 0; pos.x < dimensions.x; pos.x++) {
                            float value = lock->getVoxelNormalized(offset + pos);
                            output->setVoxelNormalized(value, pos, channel);
                        }
                    }
                }
            }
            else if(const VolumeDisk* vd = channels_[channel]->getRepresentation<VolumeDisk>()) {
                std::unique_ptr<VolumeRAM> brick(vd->loadBrick(offset, dimensions));
                tgt::svec3 pos;
                for (pos.z = 0; pos.z < dimensions.z; pos.z++) {
                    for (pos.y = 0; pos.y < dimensions.y; pos.y++) {
                        for (pos.x = 0; pos.x < dimensions.x; pos.x++) {
                            float value = brick->getVoxelNormalized(pos);
                            output->setVoxelNormalized(value, pos, channel);
                        }
                    }
                }
            }
            else {
                tgtAssert(false, "Could not get representation for channel");
            }
        }

        return output;
    }

private:

    std::vector<const VolumeBase*> channels_;
};



VolumeListMultiChannelAdapter::VolumeListMultiChannelAdapter()
    : Processor()
    , inport_(Port::INPORT, "volumelist.input", "Volume List Input", false)
    , outport_(Port::OUTPORT, "volumelist.output", "Volume List Output ", false)
    , numChannels_("numChannels", "Num. Channels", 3, 1, 4)
    , layout_("layout", "Layout", Processor::INVALID_RESULT, false, Property::LOD_ADVANCED)
{
    addPort(inport_);
    inport_.addCondition(new PortConditionVolumeListEnsemble());
    inport_.addCondition(new PortConditionVolumeListAdapter(new PortConditionVolumeChannelCount(1)));
    addPort(outport_);

    addProperty(numChannels_);
    addProperty(layout_);
    layout_.addOption("xyzxyz", "xyzxyz");
    layout_.addOption("xxyyzz", "xxyyzz");
}

VolumeListMultiChannelAdapter::~VolumeListMultiChannelAdapter() {}

Processor* VolumeListMultiChannelAdapter::create() const {
    return new VolumeListMultiChannelAdapter();
}

void VolumeListMultiChannelAdapter::process() {

    const VolumeList* input = inport_.getData();
    tgtAssert(input, "no input");

    // Clear old data (order matters!).
    outport_.clear();
    volumes_.clear();

    size_t numChannels = static_cast<size_t>(numChannels_.get());
    size_t numVolumes = input->size() / numChannels; // floor(x).

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

        VolumeDisk* vd = new VolumeDisk_MultiChannelAdapter(channels);
        VolumeBase* volume = new Volume(vd, input->first());
        output->add(volume);

        // Transfer ownership.
        volumes_.push_back(std::unique_ptr<const VolumeBase>(volume));
    }

    outport_.setData(output, true);
}

}   // namespace
