/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#pragma once

#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/volumebase.h"

#include "voreen/core/io/serialization/serializable.h"
#include "voreen/core/io/serialization/xmldeserializer.h"
#include <vector>
#include <string>

#include "tgt/vector.h"

#include <lz4.h>
#include <boost/optional.hpp>

namespace voreen {

class LZ4SliceVolumeMetadata : public Serializable {
public:
    LZ4SliceVolumeMetadata(tgt::svec3 dimensions);
    LZ4SliceVolumeMetadata(const std::string& xmlfile);
    LZ4SliceVolumeMetadata(const LZ4SliceVolumeMetadata&) = default;

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);

    void save(const std::string& xmlfile) const;

    tgt::svec3 dimensions_;
};

template<typename Voxel>
class LZ4SliceVolumeBuilder;

template<typename Voxel>
class LZ4SliceVolume;

template<typename Voxel>
class LZ4WriteableSlice {
public:
    VolumeAtomic<Voxel>& operator*() { return slice_; }
    VolumeAtomic<Voxel>* operator->() { return &slice_; }

    LZ4WriteableSlice(LZ4WriteableSlice&& other)
        : volume_(other.volume_)
        , sliceNum_(other.sliceNum_)
        , slice_(std::move(other.slice_))
    {
    }
    ~LZ4WriteableSlice();

private:
    LZ4WriteableSlice(LZ4SliceVolume<Voxel>& volume, size_t sliceNum, VolumeAtomic<Voxel>&& slice);

    friend class LZ4SliceVolume<Voxel>;
    friend class LZ4SliceVolumeBuilder<Voxel>;
    LZ4SliceVolume<Voxel>& volume_; //MUST live longer than this object
    size_t sliceNum_;
    VolumeAtomic<Voxel> slice_;
};

template<typename Voxel>
class LZ4SliceVolume {
public:
    static LZ4SliceVolume open(std::string filePath);

    LZ4SliceVolume(LZ4SliceVolume&& other);

    LZ4SliceVolume<Voxel>& operator=(LZ4SliceVolume<Voxel>&& other);
    //LZ4SliceVolume(const LZ4SliceVolume& other) = delete; //Disable

    VolumeAtomic<Voxel> loadSlice(size_t sliceNumber) const;
    void writeSlice(const VolumeAtomic<Voxel>& slice, size_t sliceNumber);
    LZ4WriteableSlice<Voxel> getWritableSlice(size_t sliceNumber);

    tgt::svec3 getDimensions() const;
    size_t getNumSlices() const;

private:
    friend class LZ4SliceVolumeBuilder<Voxel>;
    LZ4SliceVolume(std::string filePath, LZ4SliceVolumeMetadata metadata);

    std::string getSliceFilePath(size_t sliceNum) const;
    tgt::svec3 getSliceDimensions() const;
    size_t getSliceMemorySize() const;

    LZ4SliceVolumeMetadata metadata_;
    std::string filePath_;
};

/*
template<typename Voxel>
class LZ4SliceVolumeSliceCacher {
public:
    LZ4SliceVolumeSliceCacher(LZ4SliceVolume<Voxel>& volume);

    const VolumeAtomic<Voxel>& loadSlice(size_t sliceNumber) const;
private:
    LZ4SliceVolume<Voxel>& volume_;
    std::unique_ptr<VolumeAtomic<Voxel>> slice_;
    size_t sliceNum_;
};
*/

template<typename Voxel, uint64_t neighborhoodExtent>
class LZ4SliceVolumeReader {
    const static uint64_t neighborhoodSize = 2*neighborhoodExtent+1;
    static int sliceStorageIndextoSlicePosOffset(int sliceIndex) {
        return sliceIndex - neighborhoodExtent;
    }
    static int slicePosOffsetToSliceStorageIndex(int slicePosOffset) {
        return slicePosOffset + neighborhoodExtent;
    }
public:
    LZ4SliceVolumeReader(const LZ4SliceVolume<Voxel>& volume);

    const boost::optional<VolumeAtomic<Voxel>>& getSlice(int sliceNumber) const;
    boost::optional<Voxel> getVoxel(tgt::ivec3 pos) const;

    void seek(size_t sliceNumber);
    void advance();
private:
    boost::optional<VolumeAtomic<Voxel>> loadSliceFromVolume(size_t sliceNumber) const;

    const LZ4SliceVolume<Voxel>& volume_;
    std::array<boost::optional<VolumeAtomic<Voxel>>, neighborhoodSize> slices_;
    size_t pos_;
};

template<typename Voxel>
class LZ4SliceVolumeBuilder {
public:
    LZ4SliceVolumeBuilder(std::string filePath, LZ4SliceVolumeMetadata metadata);
    LZ4SliceVolumeBuilder(LZ4SliceVolumeBuilder&& other);

    LZ4WriteableSlice<Voxel> getNextWritableSlice();

    static LZ4SliceVolume<Voxel> finalize(LZ4SliceVolumeBuilder&&);

private:
    LZ4SliceVolume<Voxel> volumeInConstruction_;
    size_t numSlicesPushed_;
};

/// LZ4WritableSlice -----------------------------------------------------------

template<typename Voxel>
LZ4WriteableSlice<Voxel>::LZ4WriteableSlice(LZ4SliceVolume<Voxel>& volume, size_t sliceNum, VolumeAtomic<Voxel>&& slice)
    : volume_(volume)
    , sliceNum_(sliceNum)
    , slice_(std::move(slice))
{
}

template<typename Voxel>
LZ4WriteableSlice<Voxel>::~LZ4WriteableSlice() {
    volume_.writeSlice(slice_, sliceNum_);
}

/// LZ4SliceVolume -------------------------------------------------------------
template<typename Voxel>
LZ4SliceVolume<Voxel> LZ4SliceVolume<Voxel>::open(std::string filePath) {
    return LZ4SliceVolume(filePath, LZ4SliceVolumeMetadata(filePath));
}

template<typename Voxel>
std::string LZ4SliceVolume<Voxel>::getSliceFilePath(size_t sliceNum) const {
    return filePath_ + "_slice" + std::to_string(sliceNum);
}

template<typename Voxel>
tgt::svec3 LZ4SliceVolume<Voxel>::getSliceDimensions() const {
    return tgt::svec3(metadata_.dimensions_.xy(), 1);
}

template<typename Voxel>
size_t LZ4SliceVolume<Voxel>::getSliceMemorySize() const {
    return sizeof(Voxel) * tgt::hmul(getSliceDimensions());
}

template<typename Voxel>
size_t LZ4SliceVolume<Voxel>::getNumSlices() const {
    return metadata_.dimensions_.z;
}

template<typename Voxel>
LZ4SliceVolume<Voxel>::LZ4SliceVolume(LZ4SliceVolume<Voxel>&& other)
    : filePath_(other.filePath_)
    , metadata_(std::move(other.metadata_))
{
    other.filePath_ = "";
}

template<typename Voxel>
LZ4SliceVolume<Voxel>& LZ4SliceVolume<Voxel>::operator=(LZ4SliceVolume<Voxel>&& other) {
    // Destruct the current object, but keep the memory.
    this->~LZ4SliceVolume();
    // Call the move constructor on the memory region of the current object.
    new(this) LZ4SliceVolume(std::move(other));

    return *this;
}

template<typename Voxel>
LZ4SliceVolume<Voxel>::LZ4SliceVolume(std::string filePath, LZ4SliceVolumeMetadata metadata)
    : filePath_(filePath)
    , metadata_(metadata)
{
}

template<typename Voxel>
VolumeAtomic<Voxel> LZ4SliceVolume<Voxel>::loadSlice(size_t sliceNumber) const {
    std::ifstream compressedFile(getSliceFilePath(sliceNumber), std::ifstream::binary);

    const size_t sliceMemorySize = getSliceMemorySize();

    compressedFile.seekg(0,compressedFile.end);
    size_t compressedFileSize = compressedFile.tellg();
    compressedFile.clear(); //?
    compressedFile.seekg(0, std::ios::beg);

    tgtAssert(compressedFileSize < sliceMemorySize, "Invalid slice memory size");

    //Read file into buffer (sliceMemorySize is larger than it needs to be, but this way we do not need to check the actual file size)
    std::unique_ptr<char[]> compressedBuffer(new char[compressedFileSize]);
    compressedFile.read(compressedBuffer.get(), sliceMemorySize); //TODO check if this is valid (because sliceMemorySize may be larger than the file
    compressedFile.close();

    // Decompress buffer
    std::unique_ptr<char[]> decompressedBuffer(new char[sliceMemorySize]);
    size_t bytesDecompressed = LZ4_decompress_safe(compressedBuffer.get(), decompressedBuffer.get(), compressedFileSize, sliceMemorySize);
    tgtAssert(bytesDecompressed == sliceMemorySize, "Invalid memory size (resulting in memory corruption!)");

    return VolumeAtomic<Voxel>(reinterpret_cast<unsigned char*>(decompressedBuffer.release()), getSliceDimensions());
}

template<typename Voxel>
void LZ4SliceVolume<Voxel>::writeSlice(const VolumeAtomic<Voxel>& slice, size_t sliceNumber) {
    tgtAssert(slice.getDimensions() == getSliceDimensions(), "Invalid slice dimensions");
    tgtAssert(sliceNumber < getNumSlices(), "Invalid slice number");

    const size_t sliceMemorySize = getSliceMemorySize();

    std::unique_ptr<char[]> compressedBuffer(new char[sliceMemorySize]);
    size_t compressedSize = LZ4_compress_default((const char *)(slice.getData()), compressedBuffer.get(), sliceMemorySize, sliceMemorySize);
    tgtAssert(compressedSize > 0, "Compression failed");

    std::ofstream outStream(getSliceFilePath(sliceNumber), std::ofstream::binary);
    outStream.write(compressedBuffer.get(), compressedSize);
}

template<typename Voxel>
LZ4WriteableSlice<Voxel> LZ4SliceVolume<Voxel>::getWritableSlice(size_t sliceNumber) {
    return LZ4WriteableSlice<Voxel>(*this, sliceNumber, loadSlice(sliceNumber));
}

template<typename Voxel>
tgt::svec3 LZ4SliceVolume<Voxel>::getDimensions() const {
    return metadata_.dimensions_;
}

LZ4SliceVolume<uint8_t> binarizeVolume(const VolumeBase& volume, float binarizationThresholdSegmentationNormalized);

/// LZ4SliceVolumeReader --------------------------------------------------

template<typename Voxel, uint64_t neighborhoodExtent>
LZ4SliceVolumeReader<Voxel, neighborhoodExtent>::LZ4SliceVolumeReader(const LZ4SliceVolume<Voxel>& volume)
    : volume_(volume)
    , slices_()
    , pos_(-1)
{
    std::fill(slices_.begin(), slices_.end(), boost::none);
}

template<typename Voxel, uint64_t neighborhoodExtent>
void LZ4SliceVolumeReader<Voxel, neighborhoodExtent>::seek(size_t newPos) {
    if(pos_ != newPos) {
        pos_ = newPos;

        for(size_t i=0; i<neighborhoodSize; ++i) {
            auto pos = pos_+slicePosOffsetToSliceStorageIndex(i);
            slices_[i] = std::move(loadSliceFromVolume(pos));
        }
    }
}

template<typename Voxel, uint64_t neighborhoodExtent>
void LZ4SliceVolumeReader<Voxel, neighborhoodExtent>::advance() {
    ++pos_;

    for(size_t i=1; i<neighborhoodSize; ++i) {
        std::swap(slices_[i-1], slices_[i]);
    }
    slices_[neighborhoodSize-1] = loadSliceFromVolume(pos_ + neighborhoodExtent);
}

template<typename Voxel, uint64_t neighborhoodExtent>
const boost::optional<VolumeAtomic<Voxel>>& LZ4SliceVolumeReader<Voxel, neighborhoodExtent>::getSlice(int sliceNumber) const {
    int sliceStorageIndex = slicePosOffsetToSliceStorageIndex(sliceNumber - pos_);
    tgtAssert(0 <= sliceStorageIndex && sliceStorageIndex < neighborhoodSize, "Invalid slice number");
    return slices_[sliceStorageIndex];
}

template<typename Voxel, uint64_t neighborhoodExtent>
boost::optional<Voxel> LZ4SliceVolumeReader<Voxel, neighborhoodExtent>::getVoxel(tgt::ivec3 pos) const {
    const auto& slice = getSlice(pos.z);
    if(slice) {
        return slice->voxel(pos.x, pos.y, 0);
    } else {
        return boost::none;
    }
}

template<typename Voxel, uint64_t neighborhoodExtent>
boost::optional<VolumeAtomic<Voxel>> LZ4SliceVolumeReader<Voxel, neighborhoodExtent>::loadSliceFromVolume(size_t sliceNumber) const {
    if(0 <= sliceNumber && sliceNumber < volume_.getNumSlices()) {
        return volume_.loadSlice(sliceNumber);
    } else {
        return boost::none;
    }
}

/// LZ4SliceVolumeBuilder ------------------------------------------------------
template<typename Voxel>
LZ4SliceVolumeBuilder<Voxel>::LZ4SliceVolumeBuilder(std::string filePath, LZ4SliceVolumeMetadata metadata)
    : volumeInConstruction_(filePath, metadata)
    , numSlicesPushed_(0)
{
}

template<typename Voxel>
LZ4SliceVolumeBuilder<Voxel>::LZ4SliceVolumeBuilder(LZ4SliceVolumeBuilder&& other)
    : volumeInConstruction_(std::move(other.volumeInConstruction_))
    , numSlicesPushed_(other.numSlicesPushed_)
{
    other.numSlicesPushed_ = -1;
}

template<typename Voxel>
LZ4WriteableSlice<Voxel> LZ4SliceVolumeBuilder<Voxel>::getNextWritableSlice() {
    tgtAssert(numSlicesPushed_ < volumeInConstruction_.getNumSlices(), "Cannot push more slices");
    ++numSlicesPushed_;
    return LZ4WriteableSlice<Voxel>(volumeInConstruction_, numSlicesPushed_-1, VolumeAtomic<Voxel>(tgt::svec3(volumeInConstruction_.getDimensions().xy(), 1)));
}

template<typename Voxel>
LZ4SliceVolume<Voxel> LZ4SliceVolumeBuilder<Voxel>::finalize(LZ4SliceVolumeBuilder<Voxel>&& builder) {
    auto tmp = std::move(builder);
    tmp.volumeInConstruction_.metadata_.save(tmp.volumeInConstruction_.filePath_);
    return std::move(tmp.volumeInConstruction_);
}


}
