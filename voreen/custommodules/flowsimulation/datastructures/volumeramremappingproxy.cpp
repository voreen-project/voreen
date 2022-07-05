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

#include "volumeramremappingproxy.h"

#include "voreen/core/datastructures/volume/volumefactory.h"

namespace voreen {

const std::string VolumeRAMRemappingProxy::loggerCat_("voreen.VolumeRAMRemappingProxy");


void VolumeRAMRemappingProxy::unimplemented() const {
    LERROR("Function currently not implemented");
}

void VolumeRAMRemappingProxy::readOnly() const {
    LERROR("This is a read-only representation");
}


VolumeRAMRemappingProxy::VolumeRAMRemappingProxy(tgt::svec3 resolution,
                                                 const VolumeBase* originalVolume,
                                                 RemappingFunction remappingFunction,
                                                 float missingValue
                                                 )
    : VolumeRAM(resolution)
    , originalVolume_(originalVolume)
    , data_(originalVolume_)
    , remappingFunction_(std::move(remappingFunction))
    , missingValue_(missingValue)
{
    tgtAssert(originalVolume_, "volume null");
}

size_t VolumeRAMRemappingProxy::getNumBytes() const {
    return 0; // Should be fine? We don't store any additional memory...
}

size_t VolumeRAMRemappingProxy::getNumChannels() const {
    return data_->getNumChannels();
}

size_t VolumeRAMRemappingProxy::getBitsAllocated() const {
    return 0; // Should be fine? We don't store any additional memory...
}

/// Returns the number of bytes that are allocated for each voxel.
size_t VolumeRAMRemappingProxy::getBytesPerVoxel() const {
    return 0; //data_->getBytesPerVoxel(); //TODO: required for the memory manager...
}

/// Returns whether the volume's data type is a signed type.
bool VolumeRAMRemappingProxy::isSigned() const {
    return data_->isSigned();
}

/// Returns whether the volume's data type is an integer type.
bool VolumeRAMRemappingProxy::isInteger() const {
    return data_->isInteger();
}

std::string VolumeRAMRemappingProxy::getFormat() const {
    return data_->getFormat();
}
std::string VolumeRAMRemappingProxy::getBaseType() const {
    return data_->getBaseType();
}


void VolumeRAMRemappingProxy::invalidate() const {
    rawData_.reset();
}


float VolumeRAMRemappingProxy::minNormalizedValue(size_t channel) const {
    return data_->minNormalizedValue(channel);
}

float VolumeRAMRemappingProxy::maxNormalizedValue(size_t channel) const {
    return data_->maxNormalizedValue(channel);
}

float VolumeRAMRemappingProxy::minNormalizedMagnitude() const {
    return data_->minNormalizedMagnitude();
}

float VolumeRAMRemappingProxy::maxNormalizedMagnitude() const {
    return data_->maxNormalizedMagnitude();
}

tgt::vec2 VolumeRAMRemappingProxy::elementRange() const {
    return data_->elementRange();
}


VolumeRAM* VolumeRAMRemappingProxy::clone() const {
    unimplemented();
    return nullptr;
}

VolumeRAM* VolumeRAMRemappingProxy::clone(void* data) const {
    unimplemented();
    return nullptr;
}

VolumeRAM* VolumeRAMRemappingProxy::createNew(const tgt::svec3& dimensions, bool allocMem) const {
    unimplemented();
    return nullptr;
}

void VolumeRAMRemappingProxy::clear() {
    readOnly();
}

const void* VolumeRAMRemappingProxy::getData() const {
    return getRawData()->getData();
}

void* VolumeRAMRemappingProxy::getData() {
    readOnly();
    return nullptr;
}

void* VolumeRAMRemappingProxy::getBrickData(const tgt::svec3& offset, const tgt::svec3& dimension) const {
    return getRawData()->getBrickData(offset, dimension);
}
void* VolumeRAMRemappingProxy::getSliceData(const size_t firstSlice, const size_t lastSlice) const {
    return getRawData()->getSliceData(firstSlice, lastSlice);
}

VolumeRAM* VolumeRAMRemappingProxy::getSubVolume(tgt::svec3 dimensions, tgt::svec3 offset) const {
    VolumeRAM* vol = VolumeFactory().create(getFormat(), dimensions);

#ifdef VRN_MODULE_OPENMP
#pragma omp parallel for
#endif
    for(long z=0; z<static_cast<long>(dimensions.z); z++) {
        for(size_t y=0; y<dimensions.y; y++) {
            for(size_t x=0; x<dimensions.x; x++) {
                tgt::svec3 pos(x, y, z);
                for (size_t channel = 0; channel < getNumChannels(); channel++) {
                    float val = getVoxelNormalized(pos + offset, channel);
                    vol->setVoxelNormalized(val, pos, channel);
                }
            }
        }
    }

    return vol;
}

float VolumeRAMRemappingProxy::getVoxel(tgt::vec3 pos, size_t channel) const {
    if(remappingFunction_(pos)) {
        return data_->getVoxelNormalizedLinear(pos, channel);
    }
    return missingValue_;
}

float VolumeRAMRemappingProxy::getVoxelNormalized(const tgt::svec3& pos, size_t channel) const {
    return getVoxel(pos, channel);
}
float VolumeRAMRemappingProxy::getVoxelNormalized(size_t x, size_t y, size_t z, size_t channel) const {
    return getVoxel(tgt::svec3(x, y, z), channel);
}
float VolumeRAMRemappingProxy::VolumeRAMRemappingProxy::getVoxelNormalized(size_t index, size_t channel) const {
    return getVoxel(indexToPos(index), channel);
}

void VolumeRAMRemappingProxy::setVoxelNormalized(float value, const tgt::svec3& pos, size_t channel) {
    readOnly();
}
void VolumeRAMRemappingProxy::setVoxelNormalized(float value, size_t x, size_t y, size_t z, size_t channel) {
    readOnly();
}
void VolumeRAMRemappingProxy::setVoxelNormalized(float value, size_t index, size_t channel) {
    readOnly();
}

tgt::svec3 VolumeRAMRemappingProxy::indexToPos(size_t index) const {
    tgt::svec3 dim = getDimensions();
    size_t z = index / (dim.x * dim.y);
    index -= (z * dim.x * dim.y);
    size_t y = index / dim.x;
    size_t x = index % dim.x;
    return tgt::svec3(x, y, z);
}

VolumeRAM* VolumeRAMRemappingProxy::getRawData() const {
    if(!rawData_) {
        LWARNING("Gridded RAM representation requested, this may take a while...");
        rawData_.reset(getSubVolume(getDimensions(), tgt::svec3::zero));
    }
    return rawData_.get();
}

} // namespace voreen
