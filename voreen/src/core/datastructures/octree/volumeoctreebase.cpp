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

#include "voreen/core/datastructures/octree/volumeoctreebase.h"
#include "voreen/core/datastructures/octree/volumeoctreenodegeneric.h"

#include "voreen/core/datastructures/octree/octreeutils.h"

#include "voreen/core/utils/stringutils.h"
#include "voreen/core/datastructures/volume/volume.h"
#include <stdint.h>

#include "tgt/assert.h"
#include "tgt/logmanager.h"
#include "tgt/tgt_math.h"
#include "tgt/stopwatch.h"
#include "tgt/types.h"

using tgt::svec3;
using tgt::vec3;

namespace voreen {

//-------------------------------------------------------------------------------------------------
// VolumeOctreeNode

VolumeOctreeNode::VolumeOctreeNode()
    : brickAddress_(std::numeric_limits<uint64_t>::max())
    , inVolume_(true)
{
    children_[0] = 0;
    children_[1] = 0;
    children_[2] = 0;
    children_[3] = 0;
    children_[4] = 0;
    children_[5] = 0;
    children_[6] = 0;
    children_[7] = 0;
}

VolumeOctreeNode::~VolumeOctreeNode() {
}

bool VolumeOctreeNode::hasBrick() const {
    return brickAddress_ != std::numeric_limits<uint64_t>::max();
}

uint64_t VolumeOctreeNode::getBrickAddress() const {
    return brickAddress_;
}

bool VolumeOctreeNode::isHomogeneous() const {
    return !hasBrick();
}

bool VolumeOctreeNode::inVolume() const {
    return inVolume_;
}

bool VolumeOctreeNode::isLeaf() const {
    return (!children_[0] &&
            !children_[1] &&
            !children_[2] &&
            !children_[3] &&
            !children_[4] &&
            !children_[5] &&
            !children_[6] &&
            !children_[7]    );
}

size_t VolumeOctreeNode::getDepth() const {
    size_t depth = 0;

    if (children_[0]) depth = std::max(depth, children_[0]->getDepth());
    if (children_[1]) depth = std::max(depth, children_[1]->getDepth());
    if (children_[2]) depth = std::max(depth, children_[2]->getDepth());
    if (children_[3]) depth = std::max(depth, children_[3]->getDepth());
    if (children_[4]) depth = std::max(depth, children_[4]->getDepth());
    if (children_[5]) depth = std::max(depth, children_[5]->getDepth());
    if (children_[6]) depth = std::max(depth, children_[6]->getDepth());
    if (children_[7]) depth = std::max(depth, children_[7]->getDepth());

    return depth+1;
}

size_t VolumeOctreeNode::getNodeCount() const {
    size_t nodeCount = 1;

    if (children_[0]) nodeCount += children_[0]->getNodeCount();
    if (children_[1]) nodeCount += children_[1]->getNodeCount();
    if (children_[2]) nodeCount += children_[2]->getNodeCount();
    if (children_[3]) nodeCount += children_[3]->getNodeCount();
    if (children_[4]) nodeCount += children_[4]->getNodeCount();
    if (children_[5]) nodeCount += children_[5]->getNodeCount();
    if (children_[6]) nodeCount += children_[6]->getNodeCount();
    if (children_[7]) nodeCount += children_[7]->getNodeCount();

    return nodeCount;
}

size_t VolumeOctreeNode::getNumBricks() const {
    size_t numBricks = (brickAddress_ < std::numeric_limits<uint64_t>::max() ? 1 : 0);

    if(children_[0]) numBricks += children_[0]->getNumBricks();
    if(children_[1]) numBricks += children_[1]->getNumBricks();
    if(children_[2]) numBricks += children_[2]->getNumBricks();
    if(children_[3]) numBricks += children_[3]->getNumBricks();
    if(children_[4]) numBricks += children_[4]->getNumBricks();
    if(children_[5]) numBricks += children_[5]->getNumBricks();
    if(children_[6]) numBricks += children_[6]->getNumBricks();
    if(children_[7]) numBricks += children_[7]->getNumBricks();

    return numBricks;
}

void VolumeOctreeNode::serializeContentToBinaryBuffer(char* buffer) const {
    tgtAssert(buffer, "null pointer passed");

    memcpy(buffer, &brickAddress_, sizeof(uint64_t));
    buffer += sizeof(uint64_t);
    memcpy(buffer, &inVolume_, sizeof(bool));
}

void VolumeOctreeNode::deserializeContentFromBinaryBuffer(const char* buffer) {
    tgtAssert(buffer, "null pointer passed");

    memcpy(const_cast<uint64_t*>(&brickAddress_), buffer, sizeof(uint64_t));
    buffer += sizeof(uint64_t);
    memcpy(const_cast<bool*>(&inVolume_), buffer, sizeof(bool));
}

//-----------------------------------------------------------------------------
// VolumeOctreeBase

const std::string VolumeOctreeBase::loggerCat_("voreen.VolumeOctreeBase");

VolumeOctreeBase::VolumeOctreeBase(const tgt::svec3& brickDim, const tgt::svec3& volumeDim, size_t numChannels) 
    : VolumeRepresentation(volumeDim)
    , numChannels_(numChannels)
    , brickDim_(brickDim)
    , bytesPerVoxel_(2 * numChannels) //< uint16_t
{
    if (numChannels == 0 || numChannels > 4)
        throw VoreenException("Number of channels (volumes) must be between 1 and 4");

    // check brick dimensions
    if (!tgt::isPowerOfTwo((int)brickDim_.x))
        throw VoreenException("Brick dimensions must be power of two: " + genericToString(brickDim_));

    // check volume dimensions
    if (!tgt::hand(tgt::greaterThan(volumeDim, tgt::svec3::one)))
        throw VoreenException("Volume dimensions must be greater than [1,1,1]: " + genericToString(brickDim_));

    // compute octree dimensions (cubic, power-of-two) from volume dimensions
    octreeDim_ = tgt::svec3(tgt::nextLargerPowerOfTwo((int)tgt::max(volumeDim)));
    tgtAssert(isCubicAndPot(octreeDim_) && tgt::hand(tgt::greaterThanEqual(octreeDim_, getVolumeDim())), "invalid octree dimensions");

    // check brick dimensions against octree dimensions
    if (tgt::hor(tgt::greaterThan(brickDim_, octreeDim_))) {
        throw VoreenException("Brick dimensions " + genericToString(brickDim_) +
            " larger than octree dimensions " + genericToString(octreeDim_));
    }

    // determine (theoretical) tree depth
    numLevels_ = tgt::ilog2((int)(octreeDim_.x / brickDim.x)) + 1;
    tgtAssert(numLevels_ > 0, "invalid level count");

}

VolumeOctreeBase::VolumeOctreeBase()
{}

std::string VolumeOctreeBase::getFormat() const {
    tgtAssert(getNumChannels() > 0, "invalid number of channels");
    if (numChannels_ == 1)
        return getBaseType();
    else
        return "Vector" + itos(getNumChannels()) + "(" + getBaseType() + ")";

}

std::string VolumeOctreeBase::getBaseType() const {
    return "uint16";
}

tgt::svec3 VolumeOctreeBase::getVolumeDim() const {
    return getDimensions();
}

size_t VolumeOctreeBase::getNumChannels() const {
    return numChannels_;
}

size_t VolumeOctreeBase::getBytesPerVoxel() const {
    return bytesPerVoxel_;
}

tgt::svec3 VolumeOctreeBase::getOctreeDim() const {
    return octreeDim_;
}

tgt::svec3 VolumeOctreeBase::getBrickDim() const {
    return brickDim_;
}

size_t VolumeOctreeBase::getNumVoxelsPerBrick() const {
    return tgt::hmul(brickDim_);
}

size_t VolumeOctreeBase::getBrickMemorySize() const {
    return getNumVoxelsPerBrick()*getBytesPerVoxel();
}

size_t VolumeOctreeBase::getNumLevels() const {
    return numLevels_;
}

const std::string& VolumeOctreeBase::getOctreeConfigurationHash() const {
    return octreeConfigurationHash_;
}

void VolumeOctreeBase::setOctreeConfigurationHash(const std::string& hash) {
   octreeConfigurationHash_ = hash;
}

void VolumeOctreeBase::serialize(Serializer& s) const {
    s.serialize("volumeDim", static_cast<tgt::ivec3>(getVolumeDim()));
    s.serialize("octreeDim", static_cast<tgt::ivec3>(octreeDim_));
    s.serialize("bytesPerVoxel", bytesPerVoxel_);

    s.serialize("numLevels", numLevels_);
    s.serialize("brickDim", static_cast<tgt::ivec3>(brickDim_));

    s.serialize("numChannels", numChannels_);

    s.serialize("octreeConfigurationHash", octreeConfigurationHash_);
}

void VolumeOctreeBase::deserialize(Deserializer& s) {
    tgt::ivec3 iVolumeDim;
    s.deserialize("volumeDim", iVolumeDim);
    dimensions_ = static_cast<tgt::svec3>(iVolumeDim);
    numVoxels_ = tgt::hmul(dimensions_);

    tgt::ivec3 iOctreeDim;
    s.deserialize("octreeDim", iOctreeDim);
    octreeDim_ = static_cast<tgt::svec3>(iOctreeDim);
    if (!isCubicAndPot(octreeDim_) || !tgt::hand(tgt::greaterThanEqual(octreeDim_, getVolumeDim())))
        throw SerializationException("Invalid octree dimensions: " + genericToString(octreeDim_));

    s.deserialize("bytesPerVoxel", bytesPerVoxel_);

    s.deserialize("numLevels", numLevels_);
    tgt::ivec3 iBrickDim;
    s.deserialize("brickDim", iBrickDim);
    brickDim_ = static_cast<tgt::svec3>(iBrickDim);

    s.deserialize("numChannels", numChannels_);

    s.deserialize("octreeConfigurationHash", octreeConfigurationHash_);
}

void VolumeOctreeBase::logDescription() const {
    const std::string description = getDescription();
    std::vector<std::string> lines = strSplit(description, "\r\n");
    if (lines.size() == 1)
        lines = strSplit(description, "\n");

    for (size_t i=0; i<lines.size(); i++)
        LINFO(lines.at(i));
}

size_t VolumeOctreeBase::getCompleteTreeNodeCount(size_t treeDepth) {
    if (treeDepth == 0)
        return 0;

    size_t numNodes = 0;
    for (size_t l=0; l<treeDepth; l++)
        numNodes += tgt::iround(powf(8.f, (float)(l)));

    return numNodes;
}

// -----------------------
// node creation functions

VolumeOctreeNode* VolumeOctreeBase::createNode(size_t numChannels) {
    tgtAssert(numChannels > 0 && numChannels <= 4, "number of channels must be between 1 and 4");
    uint16_t* avgValues = new uint16_t[numChannels];
    uint16_t* minValues = new uint16_t[numChannels];
    uint16_t* maxValues = new uint16_t[numChannels];
    for (size_t c=0; c<numChannels; c++) {
        avgValues[c] = 0;
        minValues[c] = 0;
        maxValues[c] = 0;
    }

    VolumeOctreeNode* node = createNode(numChannels, avgValues, minValues, maxValues);
    node->inVolume_ = false;

    delete[] avgValues;
    delete[] minValues;
    delete[] maxValues;

    return node;
}

VolumeOctreeNode* VolumeOctreeBase::createNode(size_t numChannels, uint16_t* avgValues, uint16_t* minValues, uint16_t* maxValues) {
    tgtAssert(numChannels > 0 && numChannels <= 4, "number of channels must be between 1 and 4");
    VolumeOctreeNode* node = createNode(numChannels, avgValues, minValues, maxValues, std::numeric_limits<uint64_t>::max());
    tgtAssert(node->inVolume(), "node not inside volume");
    return node;
}

VolumeOctreeNode* VolumeOctreeBase::createNode(size_t numChannels, uint16_t* avgValues, uint16_t* minValues, uint16_t* maxValues, uint64_t brickAddress) {
    tgtAssert(numChannels > 0 && numChannels <= 4, "number of channels must be between 1 and 4");

    VolumeOctreeNode* children[8];
    for (size_t i=0; i<8; i++)
        children[i] = 0;

    VolumeOctreeNode* node = createNode(numChannels, avgValues, minValues, maxValues, brickAddress, children);
    tgtAssert(node->inVolume(), "node not inside volume");
    return node;
}

VolumeOctreeNode* VolumeOctreeBase::createNode(size_t numChannels, uint16_t* avgValues, uint16_t* minValues, uint16_t* maxValues, uint64_t brickAddress, VolumeOctreeNode* children[8]) {
    tgtAssert(numChannels > 0 && numChannels <= 4, "number of channels must be between 1 and 4");
    tgtAssert(avgValues, "null pointer passed as avgValues");
    tgtAssert(minValues, "null pointer passed as minValues");
    tgtAssert(maxValues, "null pointer passed as maxValues");

    if (numChannels == 1) {
        VolumeOctreeNodeGeneric<1>* node = new VolumeOctreeNodeGeneric<1>();
        node->brickAddress_ = brickAddress;
        node->inVolume_ = true;
        for (size_t i=0; i<8; i++)
            node->children_[i] = children[i];

        tgtAssert(minValues[0] <= avgValues[0] && avgValues[0] <= maxValues[0], "invalid min/avg/max value combination");
        node->avgValues_[0] = avgValues[0];
        node->minValues_[0] = minValues[0];
        node->maxValues_[0] = maxValues[0];

        return node;
    }
    else if (numChannels == 2) {
        VolumeOctreeNodeGeneric<2>* node = new VolumeOctreeNodeGeneric<2>();
        node->brickAddress_ = brickAddress;
        node->inVolume_ = true;
        for (size_t i=0; i<8; i++)
            node->children_[i] = children[i];

        for (size_t c=0; c<2; c++) {
            tgtAssert(minValues[c] <= avgValues[c] && avgValues[c] <= maxValues[c], "invalid min/avg/max value combination");
            node->avgValues_[c] = avgValues[c];
            node->minValues_[c] = minValues[c];
            node->maxValues_[c] = maxValues[c];
        }

        return node;
    }
    else if (numChannels == 3) {
        VolumeOctreeNodeGeneric<3>* node = new VolumeOctreeNodeGeneric<3>();
        node->brickAddress_ = brickAddress;
        node->inVolume_ = true;
        for (size_t i=0; i<8; i++)
            node->children_[i] = children[i];

        for (size_t c=0; c<3; c++) {
            tgtAssert(minValues[c] <= avgValues[c] && avgValues[c] <= maxValues[c], "invalid min/avg/max value combination");
            node->avgValues_[c] = avgValues[c];
            node->minValues_[c] = minValues[c];
            node->maxValues_[c] = maxValues[c];
        }

        return node;
    }
    else if (numChannels == 4) {
        VolumeOctreeNodeGeneric<4>* node = new VolumeOctreeNodeGeneric<4>();
        node->brickAddress_ = brickAddress;
        node->inVolume_ = true;
        for (size_t i=0; i<8; i++)
            node->children_[i] = children[i];

        for (size_t c=0; c<4; c++) {
            tgtAssert(minValues[c] <= avgValues[c] && avgValues[c] <= maxValues[c], "invalid min/avg/max value combination");
            node->avgValues_[c] = avgValues[c];
            node->minValues_[c] = minValues[c];
            node->maxValues_[c] = maxValues[c];
        }

        return node;
    }
    else {
        tgtAssert(false, "invalid channel");
        return 0;
    }

}

//-------------------------------------------------------------------------------------------------
// converter

bool RepresentationConverterOctreeToRAM::canConvert(const VolumeRepresentation* source) const {
    return (dynamic_cast<const VolumeOctreeBase*>(source) != 0);
}

VolumeRepresentation* RepresentationConverterOctreeToRAM::convert(const VolumeRepresentation* source) const {
    tgtAssert(source, "null pointer passed");
    const VolumeOctreeBase* octree = dynamic_cast<const VolumeOctreeBase*>(source);
    if (!octree) {
        LERRORC("voreen.RepresentationConverterOctreeToRAM", "Unable to create VolumeRAM: passed representation is not an octree");
        return 0;
    }

    VolumeRAM* volumeRAM = 0;
    try {
        volumeRAM = octree->createVolume(0);
    }
    catch (VoreenException& e) {
        LERRORC("voreen.RepresentationConverterOctreeToRAM", "Failed to create RAM volume: " << e.what());
    }

    return volumeRAM;
}

} // namespace
