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

#include "slicereader.h"

namespace voreen {

// SliceReaderMetaData ---------------------------------------------------------------
SliceReaderMetaData SliceReaderMetaData::fromBase(const SliceReaderMetaData& base) {
    SliceReaderMetaData srmd(base.rwm_, base.dimensions_, base.numChannels_, base.baseType_, base.spacing_);
    return srmd;
}
SliceReaderMetaData SliceReaderMetaData::fromBaseAccurate(const SliceReaderMetaData& base) {
    SliceReaderMetaData srmd = SliceReaderMetaData::fromBase(base);
    const auto& mm = base.getMinMax();
    if(mm) {
        srmd.setMinMax(*mm);
    }
    return srmd;
}

SliceReaderMetaData SliceReaderMetaData::fromVolume(const VolumeBase& vol) {
    SliceReaderMetaData srmd(vol.getRealWorldMapping(), vol.getDimensions(), vol.getNumChannels(), vol.getBaseType(), vol.getSpacing());
    const auto vmm = vol.getDerivedData<VolumeMinMax>();

    std::vector<tgt::vec2> minmax;
    for(int c=0; c<vol.getNumChannels(); ++c) {
        minmax.emplace_back(vmm->getMin(c), vmm->getMax(c));
    }
    srmd.setMinMax(minmax);
    return srmd;
}

SliceReaderMetaData SliceReaderMetaData::fromHDF5Volume(const HDF5FileVolume& volume) {
    std::unique_ptr<RealWorldMapping> rwmPtr(volume.tryReadRealWorldMapping());

    RealWorldMapping rwm;
    if(rwmPtr) {
        rwm = *rwmPtr;
    }

    std::unique_ptr<tgt::vec3> spacingPtr(volume.tryReadSpacing());
    tgt::vec3 spacing;
    if(spacingPtr) {
        spacing = *spacingPtr;
    } else {
        spacing = tgt::vec3::one;
    }
    SliceReaderMetaData srmd(rwm, volume.getDimensions(), volume.getNumberOfChannels(), volume.getBaseType(), spacing);

    std::vector<tgt::vec2> minmax;
    for(int c=0; c<volume.getNumberOfChannels(); ++c) {
        std::unique_ptr<VolumeMinMax> vmm(volume.tryReadVolumeMinMax(c));
        if (vmm) {
            minmax.emplace_back(vmm->getMin(0), vmm->getMax(0));
        } else {
            return srmd;
        }
    }
    srmd.setMinMax(minmax);
    return srmd;
}

SliceReaderMetaData::SliceReaderMetaData(RealWorldMapping rwm, tgt::svec3 dimensions, size_t numChannels, std::string baseType, tgt::vec3 spacing)
    : rwm_(rwm)
    , minmax_(boost::none)
    , dimensions_(dimensions)
    , numChannels_(numChannels)
    , baseType_(baseType)
    , spacing_(spacing)
{
}

void SliceReaderMetaData::setMinMax(std::vector<tgt::vec2> minmax) {
    minmax_ = minmax;
}

void SliceReaderMetaData::setMinMaxNormalized(std::vector<tgt::vec2> minmaxNorm) {
    minmax_ = std::vector<tgt::vec2>();
    for(auto& mm : minmaxNorm) {
        minmax_->emplace_back(
                rwm_.normalizedToRealWorld(mm.x),
                rwm_.normalizedToRealWorld(mm.y)
                );
    }
}

void SliceReaderMetaData::setDimensions(tgt::svec3 dimensions) {
    dimensions_ = dimensions;
}

void SliceReaderMetaData::setNumChannels(size_t numChannels) {
    numChannels_ = numChannels;
}

void SliceReaderMetaData::setBaseType(std::string baseType) {
    baseType_ = baseType;
}

void SliceReaderMetaData::setRealWorldMapping(RealWorldMapping rwm) {
    rwm_ = rwm;
}

void SliceReaderMetaData::setSpacing(tgt::vec3 spacing) {
    spacing_ = spacing;
}

const RealWorldMapping& SliceReaderMetaData::getRealworldMapping() const {
    return rwm_;
}
tgt::vec3 SliceReaderMetaData::getSpacing() const {
    return spacing_;
}

tgt::vec2 SliceReaderMetaData::estimateMinMax() const {
    float min, max;
    if(minmax_) {
        min = minmax_->at(0).x;
        max = minmax_->at(0).y;

        for(int c=1; c<minmax_->size(); ++c) {
            min = std::min(min, minmax_->at(c).x);
            max = std::max(max, minmax_->at(c).x);
        }
    } else {
        min = getRealworldMapping().getOffset();
        max = getRealworldMapping().getScale() + min;
    }
    return tgt::vec2(min,max);
}

std::unique_ptr<VolumeMinMax> SliceReaderMetaData::getVolumeMinMax() const {
    if(!minmax_) {
        return std::unique_ptr<VolumeMinMax>(nullptr);
    }
    std::vector<float> min;
    std::vector<float> max;
    std::vector<float> minNorm;
    std::vector<float> maxNorm;
    for(auto mm : *minmax_) {
        min.push_back(mm.x);
        max.push_back(mm.y);
        minNorm.push_back(rwm_.realWorldToNormalized(mm.x));
        maxNorm.push_back(rwm_.realWorldToNormalized(mm.y));
    }
    return std::unique_ptr<VolumeMinMax>(new VolumeMinMax(min, max, minNorm, maxNorm));
}

const boost::optional<std::vector<tgt::vec2>>& SliceReaderMetaData::getMinMax() const {
    return minmax_;
}
tgt::svec3 SliceReaderMetaData::getDimensions() const {
    return dimensions_;
}
size_t SliceReaderMetaData::getNumChannels() const {
    return numChannels_;
}
std::string SliceReaderMetaData::getBaseType() const {
    return baseType_;
}
RealWorldMapping SliceReaderMetaData::getRealWorldMapping() const {
    return rwm_;
}

// SliceReader -----------------------------------------------------------------------

SliceReader::SliceReader(const tgt::ivec3& signedDim, SliceReaderMetaData&& metadata)
    : dim_(signedDim)
    , metadata_(std::move(metadata))
{
}

const tgt::ivec3& SliceReader::getSignedDimensions() const {
    return dim_;
}
tgt::svec3 SliceReader::getDimensions() const {
    return dim_;
}
const SliceReaderMetaData& SliceReader::getMetaData() const {
    return metadata_;
}

// CachingSliceReader --------------------------------------------------------------------------

CachingSliceReader::CachingSliceReader(std::unique_ptr<SliceReader>&& base, int neighborhoodSize)
    : SliceReader(tgt::ivec3(base->getDimensions()), SliceReaderMetaData::fromBaseAccurate(base->getMetaData()))
    , base_(std::move(base))
    , slices_(2*neighborhoodSize+1, nullptr)
    , neighborhoodSize_(neighborhoodSize)
{
    tgtAssert(neighborhoodSize_>=0, "neighborhoodSize must be >= 0");
}

CachingSliceReader::~CachingSliceReader() {
    for(VolumeRAM* slice : slices_) {
        delete slice;
    }
}

int CachingSliceReader::getCurrentZPos() const {
    return base_->getCurrentZPos() - neighborhoodSize_;
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
    for(;dz <= neighborhoodSize_; ++dz) {
        base_->advance();
        VolumeRAM*& slice = getSlice(dz);
        delete slice;

        const VolumeRAM* baseSlice = base_->getCurrentSlice();
        if(baseSlice) {
            slice = baseSlice->clone();
        } else {
            // This happens if base is outside its volume bounds. That's okay!
            // In that case we will never sample from this slice anyway.
            slice = nullptr;
        }
    }
}

int CachingSliceReader::getZExtent() const {
    return neighborhoodSize_;
}

// VolumeSliceReader --------------------------------------------------------------------------

VolumeSliceReader::VolumeSliceReader(const VolumeBase& volume)
    : SliceReader(tgt::ivec3(volume.getDimensions()), SliceReaderMetaData::fromVolume(volume))
    , volume_(volume)
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

float VolumeSliceReader::getVoxelNormalized(const tgt::ivec3& xyz, size_t channel) const {
    tgtAssert(currentSlice_, "No slice");
    tgtAssert(xyz.z == currentZPos_, "invalid z pos");
    tgtAssert(channel < getNumChannels(), "Invalid channel");
    return currentSlice_->getVoxelNormalized(tgt::svec3(xyz.x, xyz.y, 0), channel);
}

const VolumeRAM* VolumeSliceReader::getCurrentSlice() const {
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
    : SliceReader(tgt::ivec3(volume.getDimensions()), SliceReaderMetaData::fromHDF5Volume(volume))
    , volume_(volume)
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
    currentSlice_.reset(currentZPos_ >= 0 && currentZPos_ < getSignedDimensions().z ? volume_.loadSlices(currentZPos_, currentZPos_, numChannels_) : nullptr);
}

int HDF5VolumeSliceReader::getCurrentZPos() const {
    return currentZPos_;
}

float HDF5VolumeSliceReader::getVoxelNormalized(const tgt::ivec3& xyz, size_t channel) const {
    tgtAssert(currentSlice_, "No slice");
    tgtAssert(xyz.z == currentZPos_, "invalid z pos");
    tgtAssert(channel < getNumChannels(), "Invalid channel");
    return currentSlice_->getVoxelNormalized(tgt::svec3(xyz.x, xyz.y, 0), channel);
}

const VolumeRAM* HDF5VolumeSliceReader::getCurrentSlice() const {
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
    : SliceReader(
            tgt::ivec3(filter->getOverwrittenDimensions() ? *filter->getOverwrittenDimensions() : base->getDimensions()),
            filter->getMetaData(base->getMetaData()))
    , base_(std::move(base))
    , filter_(std::move(filter))
    , z_(std::numeric_limits<int>::max())
    , thisToBaseScale_(float(base_->getDimensions().z) / float(getDimensions().z))
    , thisToBaseOffset_(thisToBaseScale_ * 0.5 - 0.5)
{
    //tgt::vec3 base_begin(-0.5);
    //tgt::vec3 base_end = tgt::vec3(base->getDimensions()) - tgt::vec3(0.5);

    //tgt::vec3 this_begin(-0.5);
    //tgt::vec3 this_end = tgt::vec3(getDimensions()) - tgt::vec3(0.5);

    //tgt::vec3 thisToBaseScale = (base_end - base_begin) / (this_end - this_begin);
    //tgt::vec3 thisToBaseOffset = base_begin - thisToBaseScale * this_begin;

    tgtAssert(base_->getNumChannels() == filter_->getNumInputChannels(), "Number of channels mismatch");
}

void FilteringSliceReader::advance() {
    int currentBaseZ = nearestBaseZ(z_);
    int newBaseZ = nearestBaseZ(z_+1);
    tgtAssert(newBaseZ >= currentBaseZ, "Invalid z change in advance");
    int diff = newBaseZ - currentBaseZ;
    if(diff > 2*filter_->zExtent()+1) {
        // Seek directly to avoid computing unneeded intermediate slices
        seek(z_+1);
    } else {
        // New pos is close so we have an advantage when seeking
        while(base_->getCurrentZPos() < newBaseZ) {
            base_->advance();
        }
        tgtAssert(base_->getCurrentZPos() == newBaseZ, "advance missed target z");
        ++z_;
        updateCurrentSlice();
    }
}

void FilteringSliceReader::seek(int z) {
    if(getCurrentZPos() == z) {
        return;
    }
    base_->seek(nearestBaseZ(z));
    z_ = z;
    updateCurrentSlice();

}

int FilteringSliceReader::getCurrentZPos() const {
    return z_;
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

int FilteringSliceReader::nearestBaseZ(int thisZ) {
    return tgt::round(thisToBaseScale_ * thisZ + thisToBaseOffset_);
}
void FilteringSliceReader::updateCurrentSlice() {
    if(z_ < 0 || z_ >= getSignedDimensions().z) {
        currentSlice_.reset(nullptr);
        return;
    }
    currentSlice_ = filter_->getFilteredSlice(base_.get(), z_);
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
