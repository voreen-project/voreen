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

#include "voreen/core/datastructures/volume/volumediskmultichanneladapter.h"

#include "voreen/core/datastructures/volume/volumebase.h"
#include "voreen/core/datastructures/volume/volumedisk.h"
#include "voreen/core/datastructures/volume/volumefactory.h"
#include "voreen/core/utils/hashing.h"

#include <algorithm>

namespace voreen {

VolumeDiskMultiChannelAdapter::VolumeDiskMultiChannelAdapter(const std::vector<const VolumeBase*>& channels,
                                                             const tgt::bvec3& mirror,
                                                             const std::vector<size_t>& swizzle,
                                                             const std::vector<bool>& negate)
    : VolumeDisk(VolumeFactory().getFormat(channels.front()->getBaseType(), channels.size()), channels.front()->getDimensions())
    , channels_(channels)
    , mirror_(mirror)
    , swizzle_(swizzle)
    , negate_(negate)
{
    if(swizzle_.empty()) {
        swizzle_.resize(channels_.size());
        std::iota(swizzle_.begin(), swizzle_.end(), size_t(0));
    }

    if(negate_.empty()) {
        negate_ =  std::vector<bool>(channels_.size(), false);
    }

    tgtAssert(channels_.size() == swizzle_.size(), "size mismatch");
    tgtAssert(channels_.size() == negate_.size(), "size mismatch");
    const VolumeBase* ref = channels_.front();
    for (const VolumeBase* channel : channels_) {
        tgtAssert(ref->getFormat() == channel->getFormat(), "Base Type mismatch");
        tgtAssert(ref->getDimensions() == channel->getDimensions(), "Base Type mismatch");
    }
}

std::string VolumeDiskMultiChannelAdapter::getHash() const {
    std::string hash;

    for (const VolumeBase* channel : channels_) {
        hash += channel->getHash();
    }

    std::stringstream stream;
    stream << mirror_;
    stream << std::accumulate(swizzle_.begin(), swizzle_.end(), "");
    stream << std::accumulate(negate_.begin(), negate_.end(), "");

    return VoreenHash::getHash(hash + stream.str());
}

VolumeRAM* VolumeDiskMultiChannelAdapter::loadVolume() const {
    return loadBrick(tgt::svec3::zero, dimensions_);
}

VolumeRAM* VolumeDiskMultiChannelAdapter::loadSlices(const size_t firstZSlice, const size_t lastZSlice) const {
    if (firstZSlice > lastZSlice)
        throw VoreenException("last slice must be behind first slice");

    return loadBrick(tgt::svec3(0, 0, firstZSlice),
                     tgt::svec3(dimensions_.x, dimensions_.y, lastZSlice - firstZSlice + 1));
}

VolumeRAM* VolumeDiskMultiChannelAdapter::loadBrick(const tgt::svec3& offset, const tgt::svec3& dimensions) const {
    // check parameters
    if (tgt::hmul(dimensions) == 0)
        throw VoreenException("requested brick dimensions are zero");
    if (!tgt::hand(tgt::lessThanEqual(offset + dimensions, dimensions_)))
        throw VoreenException("requested brick (at least partially) outside volume dimensions");

    // Create the output volume.
    VolumeRAM* output = VolumeFactory().create(getFormat(), dimensions);

    for (size_t channel = 0; channel < channels_.size(); channel++) {

        size_t swizzledChannel = swizzle_[channel];

        // Check if we have a ram representation already.
        if (channels_[swizzledChannel]->hasRepresentation<VolumeRAM>()) {
            VolumeRAMRepresentationLock lock(channels_[swizzledChannel]);

            tgt::svec3 pos;
            for (pos.z = 0; pos.z < dimensions.z; pos.z++) {
                size_t z = mirror_.z ? dimensions.z - offset.z - pos.z - 1 : pos.z;
                for (pos.y = 0; pos.y < dimensions.y; pos.y++) {
                    size_t y = mirror_.y ? dimensions.y - offset.y - pos.y - 1 : pos.y;
                    for (pos.x = 0; pos.x < dimensions.x; pos.x++) {
                        size_t x = mirror_.x ? dimensions.x - offset.x - pos.x - 1 : pos.x;
                        float value = lock->getVoxelNormalized(x, y, z);
                        if (negate_[swizzledChannel]) {
                            value = -value;
                        }
                        output->setVoxelNormalized(value, pos, channel);
                    }
                }
            }
        } else if (const VolumeDisk* vd = channels_[swizzledChannel]->getRepresentation<VolumeDisk>()) {
            tgt::svec3 effOffset = offset;
            effOffset.x = mirror_.x ? dimensions_.x - dimensions.x - offset.x : offset.x;
            effOffset.y = mirror_.y ? dimensions_.y - dimensions.y - offset.y : offset.y;
            effOffset.z = mirror_.z ? dimensions_.z - dimensions.z - offset.z : offset.z;

            std::unique_ptr<VolumeRAM> brick(vd->loadBrick(effOffset, dimensions));
            tgt::svec3 pos;
            for (pos.z = 0; pos.z < dimensions.z; pos.z++) {
                size_t z = mirror_.z ? dimensions.z - pos.z - 1 : pos.z;
                for (pos.y = 0; pos.y < dimensions.y; pos.y++) {
                    size_t y = mirror_.y ? dimensions.y - pos.y - 1 : pos.y;
                    for (pos.x = 0; pos.x < dimensions.x; pos.x++) {
                        size_t x = mirror_.x ? dimensions.x - pos.x - 1 : pos.x;
                        float value = brick->getVoxelNormalized(x, y, z);
                        if (negate_[swizzledChannel]) {
                            value = -value;
                        }
                        output->setVoxelNormalized(value, pos, channel);
                    }
                }
            }
        } else {
            tgtAssert(false, "Could not get representation for channel");
        }
    }

    return output;
}

}