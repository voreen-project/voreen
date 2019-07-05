/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#include "slicereader.h"

namespace voreen {

// SliceReader -----------------------------------------------------------------------

SliceReader::SliceReader(const tgt::ivec3& signedDim)
    : dim_(signedDim)
{
}

const tgt::ivec3& SliceReader::getSignedDimensions() const {
    return dim_;
}

// CachingSliceReader --------------------------------------------------------------------------

CachingSliceReader::CachingSliceReader(std::unique_ptr<SliceReader>&& base, int neighborhoodSize)
    : SliceReader(tgt::ivec3(base->getDimensions()))
    , base_(std::move(base))
    , slices_(2*neighborhoodSize+1, nullptr)
    , neighborhoodSize_(neighborhoodSize)
{
    tgtAssert(neighborhoodSize_>=0, "neighborhoodSize must be >= 0");
    /* for(VolumeRAM*& slice : slices_) {
        slice = nullptr;
    } */
}

CachingSliceReader::~CachingSliceReader() {
    for(VolumeRAM* slice : slices_) {
        delete slice;
    }
}

int CachingSliceReader::getCurrentZPos() const {
    return base_->getCurrentZPos() - neighborhoodSize_;
}

const tgt::svec3& CachingSliceReader::getDimensions() const {
    return base_->getDimensions();
}

VolumeRAM*& CachingSliceReader::getSlice(int dz) {
    tgtAssert(std::abs(dz) <= neighborhoodSize_, "Invalid dz");
    return slices_[neighborhoodSize_ + dz];
}

VolumeRAM* const& CachingSliceReader::getSlice(int dz) const {
    tgtAssert(std::abs(dz) <= neighborhoodSize_, "Invalid dz");
    return slices_[neighborhoodSize_ + dz];
}

std::string CachingSliceReader::getBaseType() const {
    return base_->getBaseType();
}

size_t CachingSliceReader::getNumChannels() const {
    return base_->getNumChannels();
}

float CachingSliceReader::getVoxelNormalized(const tgt::ivec3& xyz, size_t channel) const {
    int dz=xyz.z-getCurrentZPos();
    VolumeRAM* slice = getSlice(dz);
    tgtAssert(slice, "no slice");
    tgtAssert(channel < getNumChannels(), "Invalid channel");
    return slice->getVoxelNormalized(tgt::svec3(xyz.x, xyz.y, 0), channel);
}
const VolumeRAM* CachingSliceReader::getCurrentSlice() const {
    return getSlice(0);
}

void CachingSliceReader::advance() {
    base_->advance();
    //tgtAssert(getCurrentZPos() < getSignedDimensions().z, "Advanced past volume");

    int dz = -neighborhoodSize_;
    delete getSlice(dz);
    for(; dz < neighborhoodSize_; ++dz) {
        getSlice(dz) = getSlice(dz+1);
    }
    if(getCurrentZPos() + dz < getSignedDimensions().z) {
        getSlice(dz) = base_->getCurrentSlice()->clone(); //TODO make slices unique_ptrs
    } else {
        getSlice(dz) = nullptr;
    }
}

void CachingSliceReader::seek(int z) {
    if(getCurrentZPos() == z) {
        return;
    }

    int dz = -neighborhoodSize_;

    // this loop will only be executed while z+dz < 0 => outside the volume
    for(; z+dz < 0 && dz <= neighborhoodSize_; ++dz) {
        VolumeRAM*& slice = getSlice(dz);
        delete slice;
        slice = nullptr;
    }

    // now z+dz >= 0 => inside the volume!
    base_->seek(z+dz-1);
    for(; z+dz < getSignedDimensions().z && dz <= neighborhoodSize_; ++dz) {
        base_->advance();
        VolumeRAM*& slice = getSlice(dz);
        delete slice;
        slice = base_->getCurrentSlice()->clone();
    }

    // this loop will only be executed while z+dz >= getSignedDimensions().z => outside the volume
    for(; dz <= neighborhoodSize_; ++dz) {
        VolumeRAM*& slice = getSlice(dz);
        delete slice;
        slice = nullptr;
    }
}

int CachingSliceReader::getZExtent() const {
    return neighborhoodSize_;
}

// VolumeSliceReader --------------------------------------------------------------------------

VolumeSliceReader::VolumeSliceReader(const VolumeBase& volume)
    : SliceReader(tgt::ivec3(volume.getDimensions()))
    , volume_(volume)
    , dimensions_(volume_.getDimensions())
    , currentZPos_(std::numeric_limits<int>::max()) // Not initialized, yet
    , currentSlice_(nullptr)
    , numChannels_(volume.getNumChannels())
{
}

VolumeSliceReader::~VolumeSliceReader() {
}

void VolumeSliceReader::advance() {
    seek(currentZPos_+1);
}

void VolumeSliceReader::seek(int z) {
    if(currentZPos_ == z) {
        return;
    }
    currentZPos_ = z;
    currentSlice_.reset(currentZPos_ >= 0 && currentZPos_ < getSignedDimensions().z ? volume_.getSlice(currentZPos_) : nullptr);
}

int VolumeSliceReader::getCurrentZPos() const {
    return currentZPos_;
}

const tgt::svec3& VolumeSliceReader::getDimensions() const {
    return dimensions_;
}

float VolumeSliceReader::getVoxelNormalized(const tgt::ivec3& xyz, size_t channel) const {
    tgtAssert(currentSlice_, "No slice");
    tgtAssert(xyz.z == currentZPos_, "invalid z pos");
    tgtAssert(channel < getNumChannels(), "Invalid channel");
    return currentSlice_->getVoxelNormalized(tgt::svec3(xyz.x, xyz.y, 0), channel);
}

const VolumeRAM* VolumeSliceReader::getCurrentSlice() const {
    tgtAssert(currentSlice_, "no slice");
    return currentSlice_.get();
}

std::string VolumeSliceReader::getBaseType() const {
    return volume_.getBaseType();
}
size_t VolumeSliceReader::getNumChannels() const {
    return numChannels_;
}

// HDF5VolumeSliceReader --------------------------------------------------------------------------

HDF5VolumeSliceReader::HDF5VolumeSliceReader(const HDF5FileVolume& volume)
    : SliceReader(tgt::ivec3(volume.getDimensions()))
    , volume_(volume)
    , dimensions_(volume_.getDimensions())
    , currentZPos_(std::numeric_limits<int>::max()) // Not initialized, yet
    , currentSlice_(nullptr)
    , numChannels_(volume.getNumberOfChannels())
{
    tgtAssert(numChannels_ == 1, "Multichannel volume is not supported by HDF5VolumeSliceReader");
}

HDF5VolumeSliceReader::~HDF5VolumeSliceReader() {
}

void HDF5VolumeSliceReader::advance() {
    seek(currentZPos_+1);
}

void HDF5VolumeSliceReader::seek(int z) {
    if(currentZPos_ == z) {
        return;
    }
    currentZPos_ = z;
    currentSlice_.reset(currentZPos_ >= 0 && currentZPos_ < getSignedDimensions().z ? volume_.loadSlices(currentZPos_, currentZPos_) : nullptr);
}

int HDF5VolumeSliceReader::getCurrentZPos() const {
    return currentZPos_;
}

const tgt::svec3& HDF5VolumeSliceReader::getDimensions() const {
    return dimensions_;
}

float HDF5VolumeSliceReader::getVoxelNormalized(const tgt::ivec3& xyz, size_t channel) const {
    tgtAssert(currentSlice_, "No slice");
    tgtAssert(xyz.z == currentZPos_, "invalid z pos");
    tgtAssert(channel < getNumChannels(), "Invalid channel");
    return currentSlice_->getVoxelNormalized(tgt::svec3(xyz.x, xyz.y, 0), channel);
}

const VolumeRAM* HDF5VolumeSliceReader::getCurrentSlice() const {
    tgtAssert(currentSlice_, "no slice");
    return currentSlice_.get();
}

std::string HDF5VolumeSliceReader::getBaseType() const {
    return volume_.getBaseType();
}
size_t HDF5VolumeSliceReader::getNumChannels() const {
    return numChannels_;
}

// FilteringSliceReader --------------------------------------------------------------------------

FilteringSliceReader::FilteringSliceReader(std::unique_ptr<CachingSliceReader> base, std::unique_ptr<VolumeFilter> filter)
    : SliceReader(tgt::ivec3(base->getDimensions()))
    , base_(std::move(base))
    , filter_(std::move(filter))
{
    tgtAssert(base_->getNumChannels() == filter_->getNumInputChannels(), "Number of channels mismatch");
}

void FilteringSliceReader::advance() {
    base_->advance();
    updateCurrentSlice();
}

void FilteringSliceReader::seek(int z) {
    base_->seek(z);
    updateCurrentSlice();
}

int FilteringSliceReader::getCurrentZPos() const {
    return base_->getCurrentZPos();
}

const tgt::svec3& FilteringSliceReader::getDimensions() const {
    return base_->getDimensions();
}

float FilteringSliceReader::getVoxelNormalized(const tgt::ivec3& xyz, size_t channel) const {
    tgtAssert(xyz.z == getCurrentZPos(), "Invalid z");
    tgtAssert(channel < getNumChannels(), "Invalid channel");
    return currentSlice_->getVoxelNormalized(tgt::svec3(xyz.xy(), 0), channel);
}

const VolumeRAM* FilteringSliceReader::getCurrentSlice() const {
    return currentSlice_.get();
}

std::string FilteringSliceReader::getBaseType() const {
    return filter_->getSliceBaseType();
}

size_t FilteringSliceReader::getNumChannels() const {
    return filter_->getNumOutputChannels();
}

void FilteringSliceReader::updateCurrentSlice() {
    int z = base_->getCurrentZPos();
    if(z < 0 || z >= getSignedDimensions().z) {
        currentSlice_.reset(nullptr);
        return;
    }
    currentSlice_ = filter_->getFilteredSlice(base_.get(), z);
    tgtAssert(currentSlice_->getNumChannels() == getNumChannels(), "filter produced slice with wrong number of channels");
}

// VolumeFilterStackBuilder --------------------------------------------------------------------------

VolumeFilterStackBuilder::VolumeFilterStackBuilder(const VolumeBase& volume)
    : top_(new VolumeSliceReader(volume))
{
}

VolumeFilterStackBuilder& VolumeFilterStackBuilder::addLayer(std::unique_ptr<VolumeFilter> conv) {
    tgtAssert(top_.get(), "No top. Did you .build() already?");
    top_ = std::unique_ptr<SliceReader>(new FilteringSliceReader(
                std::unique_ptr<CachingSliceReader>(new CachingSliceReader(std::move(top_), conv->zExtent())),
                std::move(conv)
                ));
    return *this;
}

std::unique_ptr<SliceReader> VolumeFilterStackBuilder::build(int initZPos) {
    tgtAssert(top_.get(), "No top. Did you .build() already?");
    top_->seek(initZPos);
    //tgtAssert(top_->getNumChannels() == 1, "Built multichannel volume"); // Maybe someone wants that? Then we'd need accessors first
    return std::move(top_);
}

std::unique_ptr<CachingSliceReader> VolumeFilterStackBuilder::buildCaching(int initZPos, int neighborhoodSize) {
    tgtAssert(top_.get(), "No top. Did you .build() already?");
    top_->seek(initZPos);
    //tgtAssert(top_->getNumChannels() == 1, "Built multichannel volume"); // Maybe someone wants that? Then we'd need accessors first
    return std::unique_ptr<CachingSliceReader>(new CachingSliceReader(std::move(top_), neighborhoodSize));
}

void writeSlicesToHDF5File(SliceReader& reader, HDF5FileVolume& file, ProgressReporter* progress) {
    tgt::svec3 dim = reader.getDimensions();
    size_t channels = reader.getNumChannels();
    tgtAssert(dim == file.getDimensions(), "Dimension mismatch");
    tgtAssert(channels == file.getNumberOfChannels(), "numChannels mismatch");

    reader.seek(0);
    for(size_t z = 0; z < dim.z; ++z) {
        if(progress) {
            progress->setProgress(static_cast<float>(z)/dim.z);
        }
        if(z > 0) {
            reader.advance();
        }
        file.writeSlices(reader.getCurrentSlice(), z, 0, channels);
    }
    if(progress) {
        progress->setProgress(1.0f);
    }
}

} // namespace voreen
