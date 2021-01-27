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

#ifndef VRN_VOLUMEOCTREENODEGENERIC_H
#define VRN_VOLUMEOCTREENODEGENERIC_H

#include "volumeoctreebase.h"

namespace voreen {

//-----------------------------------------------------------------------------
// VolumeOctreeNodeGeneric (internal)

template<size_t C = 1>
class VRN_CORE_API VolumeOctreeNodeGeneric : public virtual VolumeOctreeNode {
    friend class VolumeOctreeBase;

public:
    VolumeOctreeNodeGeneric(uint64_t brickAddress, bool inVolume, uint16_t defaultIntensityValue = 0);
    VolumeOctreeNodeGeneric(uint16_t defaultIntensityValue = 0);

    virtual uint16_t getAvgValue(size_t channel = 0) const;
    virtual const uint16_t* getAvgValues() const;
    virtual uint16_t getMinValue(size_t channel = 0) const;
    virtual uint16_t getMaxValue(size_t channel = 0) const;

    virtual size_t getNumChannels() const;

    // serialization
    virtual void serializeContentToBinaryBuffer(char* buffer) const;
    virtual void deserializeContentFromBinaryBuffer(const char* buffer);
    virtual size_t getContentSize() const { return C*3*sizeof(uint16_t) + VolumeOctreeNode::getContentSize(); } ///< numChannels*sizeof(min/max/avg)

    uint16_t avgValues_[C];    ///< The node's average voxel values (one per channel)
    uint16_t minValues_[C];    ///< The node's min voxel values (one per channel)
    uint16_t maxValues_[C];    ///< The node's max voxel values (one per channel)
};

/*****************************************************************************************************************
 *      Implementation of VolumeOctreeNodeGeneric                                                                *
 *****************************************************************************************************************/
template<size_t C>
VolumeOctreeNodeGeneric<C>::VolumeOctreeNodeGeneric(uint64_t brickAddress, bool inVolume, uint16_t defaultIntensityValue)
    : VolumeOctreeNode(brickAddress, inVolume)
{
    for (size_t i=0; i<C; i++) {
        avgValues_[i] = defaultIntensityValue;
        minValues_[i] = defaultIntensityValue;
        maxValues_[i] = defaultIntensityValue;
    }
}

template<size_t C>
VolumeOctreeNodeGeneric<C>::VolumeOctreeNodeGeneric(uint16_t defaultIntensityValue)
    : VolumeOctreeNode()
{
    for (size_t i=0; i<C; i++) {
        avgValues_[i] = defaultIntensityValue;
        minValues_[i] = defaultIntensityValue;
        maxValues_[i] = defaultIntensityValue;
    }
}

template<size_t C>
size_t VolumeOctreeNodeGeneric<C>::getNumChannels() const {
    return C;
}

template<size_t C>
uint16_t VolumeOctreeNodeGeneric<C>::getAvgValue(size_t channel /*= 0*/) const {
    tgtAssert(channel < C, "invalid channel id");
    return avgValues_[channel];
}

template<size_t C>
const uint16_t* VolumeOctreeNodeGeneric<C>::getAvgValues() const {
    tgtAssert(C > 0, "zero or less channels");
    return &avgValues_[0];
}

template<size_t C>
uint16_t VolumeOctreeNodeGeneric<C>::getMinValue(size_t channel /*= 0*/) const {
    tgtAssert(channel < C, "invalid channel id");
    return minValues_[channel];
}

template<size_t C>
uint16_t VolumeOctreeNodeGeneric<C>::getMaxValue(size_t channel /*= 0*/) const {
    tgtAssert(channel < C, "invalid channel id");
    return maxValues_[channel];
}

template<size_t C>
void VolumeOctreeNodeGeneric<C>::serializeContentToBinaryBuffer(char* buffer) const {
    tgtAssert(buffer, "null pointer passed");

    // serialize base node contents
    VolumeOctreeNode::serializeContentToBinaryBuffer(buffer);
    buffer += VolumeOctreeNode::getContentSize();

    // serialize avg/min/max values
    for (size_t channel=0; channel < getNumChannels(); channel++) {
        memcpy(buffer, &avgValues_[channel], sizeof(uint16_t));
        buffer += sizeof(uint16_t);
        memcpy(buffer, &minValues_[channel], sizeof(uint16_t));
        buffer += sizeof(uint16_t);
        memcpy(buffer, &maxValues_[channel], sizeof(uint16_t));
        buffer += sizeof(uint16_t);
    }
}

template<size_t C>
void VolumeOctreeNodeGeneric<C>::deserializeContentFromBinaryBuffer(const char* buffer) {
    tgtAssert(buffer, "null pointer passed");

    // deserialize base node contents
    VolumeOctreeNode::deserializeContentFromBinaryBuffer(buffer);
    buffer += VolumeOctreeNode::getContentSize();

    // deserialize avg/min/max values
    for (size_t channel=0; channel<getNumChannels(); channel++) {
        memcpy(const_cast<uint16_t*>(&avgValues_[channel]), buffer, sizeof(uint16_t));
        buffer += sizeof(uint16_t);
        memcpy(const_cast<uint16_t*>(&minValues_[channel]), buffer, sizeof(uint16_t));
        buffer += sizeof(uint16_t);
        memcpy(const_cast<uint16_t*>(&maxValues_[channel]), buffer, sizeof(uint16_t));
        buffer += sizeof(uint16_t);
    }
}

} //end namespace

#endif
