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

#include "sphericalvolumeramproxy.h"

#include "voreen/core/datastructures/volume/volumefactory.h"

namespace voreen {

const float SphericalVolumeRAMProxy::MissingValue(0.0f);

const std::string SphericalVolumeRAMProxy::loggerCat_("voreen.sciviscontest2021.SphericalVolumeRAMProxy");


void SphericalVolumeRAMProxy::unimplemented() const {
    LERROR("Function currently not implemented");
}

void SphericalVolumeRAMProxy::readOnly() const {
    LERROR("This is a read-only representation");
}


SphericalVolumeRAMProxy::SphericalVolumeRAMProxy(size_t resolution,
                                                 const VolumeBase* originalVolume,
                                                 float minRadius, float maxRadius,
                                                 float minClampRadius, float maxClampRadius)
    : VolumeRAM(tgt::svec3(resolution))
    , originalVolume_(originalVolume)
    , data_(originalVolume_)
    , minRadius_(minRadius)
    , maxRadius_(maxRadius)
    , minClampRadius_(minClampRadius)
    , maxClampRadius_(maxClampRadius)
{
    tgtAssert(originalVolume_, "volume null");
}

size_t SphericalVolumeRAMProxy::getNumBytes() const {
    return 0; // Should be fine? We don't store any additional memory...
}

size_t SphericalVolumeRAMProxy::getNumChannels() const {
    return data_->getNumChannels();
}

size_t SphericalVolumeRAMProxy::getBitsAllocated() const {
    return 0; // Should be fine? We don't store any additional memory...
}

/// Returns the number of bytes that are allocated for each voxel.
size_t SphericalVolumeRAMProxy::getBytesPerVoxel() const {
    return 0; //data_->getBytesPerVoxel(); //TODO: required for the memory manager...
}

/// Returns whether the volume's data type is a signed type.
bool SphericalVolumeRAMProxy::isSigned() const {
    return data_->isSigned();
}

/// Returns whether the volume's data type is an integer type.
bool SphericalVolumeRAMProxy::isInteger() const {
    return data_->isInteger();
}

std::string SphericalVolumeRAMProxy::getFormat() const {
    return data_->getFormat();
}
std::string SphericalVolumeRAMProxy::getBaseType() const {
    return data_->getBaseType();
}


void SphericalVolumeRAMProxy::invalidate() const {
    rawData_.reset();
}


float SphericalVolumeRAMProxy::minNormalizedValue(size_t channel) const {
    return data_->minNormalizedValue(channel);
}

float SphericalVolumeRAMProxy::maxNormalizedValue(size_t channel) const {
    return data_->maxNormalizedValue(channel);
}

float SphericalVolumeRAMProxy::minNormalizedMagnitude() const {
    return data_->minNormalizedMagnitude();
}

float SphericalVolumeRAMProxy::maxNormalizedMagnitude() const {
    return data_->maxNormalizedMagnitude();
}

tgt::vec2 SphericalVolumeRAMProxy::elementRange() const {
    return data_->elementRange();
}


VolumeRAM* SphericalVolumeRAMProxy::clone() const {
    unimplemented();
    return nullptr;
}

VolumeRAM* SphericalVolumeRAMProxy::clone(void* data) const {
    unimplemented();
    return nullptr;
}

VolumeRAM* SphericalVolumeRAMProxy::createNew(const tgt::svec3& dimensions, bool allocMem) const {
    unimplemented();
    return nullptr;
}

void SphericalVolumeRAMProxy::clear() {
    readOnly();
}

const void* SphericalVolumeRAMProxy::getData() const {
    return getRawData()->getData();
}

void* SphericalVolumeRAMProxy::getData() {
    readOnly();
    return nullptr;
}

void* SphericalVolumeRAMProxy::getBrickData(const tgt::svec3& offset, const tgt::svec3& dimension) const {
    return getRawData()->getBrickData(offset, dimension);
}
void* SphericalVolumeRAMProxy::getSliceData(const size_t firstSlice, const size_t lastSlice) const {
    return getRawData()->getSliceData(firstSlice, lastSlice);
}

VolumeRAM* SphericalVolumeRAMProxy::getSubVolume(tgt::svec3 dimensions, tgt::svec3 offset) const {
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

float SphericalVolumeRAMProxy::getVoxel(tgt::svec3 pos, size_t channel) const {
    if(convCCtoSC(pos)) {
        return data_->getVoxelNormalized(pos, channel);
    }
    return MissingValue;
}

float SphericalVolumeRAMProxy::getVoxelNormalized(const tgt::svec3& pos, size_t channel) const {
    return getVoxel(pos, channel);
}
float SphericalVolumeRAMProxy::getVoxelNormalized(size_t x, size_t y, size_t z, size_t channel) const {
    return getVoxel(tgt::svec3(x, y, z), channel);
}
float SphericalVolumeRAMProxy::SphericalVolumeRAMProxy::getVoxelNormalized(size_t index, size_t channel) const {
    return getVoxel(indexToPos(index), channel);
}

void SphericalVolumeRAMProxy::setVoxelNormalized(float value, const tgt::svec3& pos, size_t channel) {
    readOnly();
}
void SphericalVolumeRAMProxy::setVoxelNormalized(float value, size_t x, size_t y, size_t z, size_t channel) {
    readOnly();
}
void SphericalVolumeRAMProxy::setVoxelNormalized(float value, size_t index, size_t channel) {
    readOnly();
}

tgt::svec3 SphericalVolumeRAMProxy::indexToPos(size_t index) const {
    tgt::svec3 dim = getDimensions();
    size_t z = index / (dim.x * dim.y);
    index -= (z * dim.x * dim.y);
    size_t y = index / dim.x;
    size_t x = index % dim.x;
    return tgt::svec3(x, y, z);
}

bool SphericalVolumeRAMProxy::convCCtoSC(tgt::svec3& pos) const {

    const tgt::vec3 p = tgt::vec3(pos) - tgt::vec3(getDimensions()) * 0.5f;

    float r = std::sqrt(std::pow(p.x, 2) + std::pow(p.y, 2) + std::pow(p.z, 2));
    float theta = std::acos(p.z/r);
    float phi = std::atan2(p.y, p.x);

    float minClampR = minClampRadius_ / maxRadius_ * getDimensions().x * 0.5f;
    float maxClampR = maxClampRadius_ / maxRadius_ * getDimensions().x * 0.5f;

    if(r < minClampR || r > maxClampR) {
        return false;
    }

    float minR = minRadius_ / maxRadius_ * getDimensions().x * 0.5f;
    float maxR = maxRadius_ / maxRadius_ * getDimensions().x * 0.5f;

    tgt::vec3 sps(phi, r, theta);
    sps.x = (sps.x + tgt::PIf) / (2 * tgt::PIf) * (data_->getDimensions().x-1);    //phi (2PI, 360)
    sps.y = (sps.y - minR)/(maxR - minR) * (data_->getDimensions().y-1); //radius
    sps.z = (sps.z / tgt::PIf) * (data_->getDimensions().z-1);                   //theta (PI, 180)

    pos = tgt::round(sps);

    return true;
}

VolumeRAM* SphericalVolumeRAMProxy::getRawData() const {
    if(!rawData_) {
        LWARNING("Gridded RAM representation requested, this may take a while...");
        rawData_.reset(getSubVolume(getDimensions(), tgt::svec3::zero));
    }
    return rawData_.get();
}

} // namespace voreen
