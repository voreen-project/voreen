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

#pragma once

#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/volumebase.h"

#include "voreen/core/io/serialization/serializable.h"
#include "voreen/core/io/serialization/xmldeserializer.h"
#include "voreen/core/utils/stringutils.h"
#include <vector>
#include <string>

#include "tgt/vector.h"

#include <lz4.h>
#include <boost/optional.hpp>

namespace voreen {

class LZ4SliceVolumeMetadata : public Serializable {
public:
    explicit LZ4SliceVolumeMetadata(tgt::svec3 dimensions);
    LZ4SliceVolumeMetadata(const LZ4SliceVolumeMetadata&) = default;

    static LZ4SliceVolumeMetadata fromVolume(const VolumeBase& vol);

    LZ4SliceVolumeMetadata withOffset(tgt::vec3 offset) const;
    LZ4SliceVolumeMetadata withSpacing(tgt::vec3 spacing) const;
    LZ4SliceVolumeMetadata withPhysicalToWorldTransformation(tgt::mat4 physicalToWorldTransformation) const;
    LZ4SliceVolumeMetadata withRealWorldMapping(RealWorldMapping realWorldMapping) const;

    const tgt::svec3& getDimensions() const;
    const tgt::vec3& getOffset() const;
    const tgt::vec3& getSpacing() const;
    const tgt::mat4& getPhysicalToWorldMatrix() const;
    tgt::mat4 getVoxelToPhysicalMatrix() const;
    tgt::mat4 getVoxelToWorldMatrix() const;
    const RealWorldMapping& getRealWorldMapping() const;

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);

private:
    tgt::svec3 dimensions_;
    tgt::vec3 spacing_;
    tgt::vec3 offset_;
    tgt::mat4 physicalToWorldTransformation_;
    RealWorldMapping realWorldMapping_;
};

class LZ4SliceVolumeMetadataFull : public LZ4SliceVolumeMetadata {
public:
    LZ4SliceVolumeMetadataFull(LZ4SliceVolumeMetadata, std::string format, std::string baseType);
    LZ4SliceVolumeMetadataFull(const LZ4SliceVolumeMetadataFull&) = default;

    static LZ4SliceVolumeMetadataFull load(const std::string& xmlfile);
    void save(const std::string& xmlfile) const;

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);

    const std::string& getFormat() const;
    const std::string& getBaseType() const;

private:
    std::string format_;
    std::string baseType_;
};

template<typename Voxel>
class LZ4SliceVolumeBuilder;

template<typename Voxel>
class LZ4SliceVolume;


template<typename Voxel>
class LZ4WriteableSlab {
public:
    VolumeAtomic<Voxel>& operator*() { return slab_; }
    VolumeAtomic<Voxel>* operator->() { return &slab_; }

    LZ4WriteableSlab(LZ4WriteableSlab&& other)
        : volume_(other.volume_)
        , zBegin_(other.zBegin_)
        , slab_(std::move(other.slab_))
    {
    }
    LZ4WriteableSlab& operator=(LZ4WriteableSlab&& other)
    {
        if(&other != this) {
            this->~LZ4WriteableSlab();
            new(this) LZ4WriteableSlab(std::move(other));
        }
        return *this;
    }
    ~LZ4WriteableSlab();

private:
    LZ4WriteableSlab(LZ4SliceVolume<Voxel>& volume, size_t sliceNum, VolumeAtomic<Voxel>&& slice);

    friend class LZ4SliceVolume<Voxel>;
    friend class LZ4SliceVolumeBuilder<Voxel>;
    LZ4SliceVolume<Voxel>& volume_; //MUST live longer than this object
    size_t zBegin_;
    VolumeAtomic<Voxel> slab_;
};

class LZ4SliceVolumeBase {
public:
    const static std::string FILE_EXTENSION;
    static std::unique_ptr<LZ4SliceVolumeBase> open(std::string filePath);

    LZ4SliceVolumeBase(std::string filePath, LZ4SliceVolumeMetadataFull metadata);
    virtual ~LZ4SliceVolumeBase() { }

    virtual std::unique_ptr<VolumeRAM> loadBaseSlab(size_t beginZ, size_t endZ /*exclusive*/) const = 0;
    virtual std::unique_ptr<LZ4SliceVolumeBase> moveToHeap() && = 0;

    std::unique_ptr<Volume> toVolume() &&;

    const LZ4SliceVolumeMetadataFull& getMetaData() const;
    const tgt::svec3& getDimensions() const;
    size_t getNumSlices() const;
    const std::string& getFilePath() const;

protected:
    LZ4SliceVolumeMetadataFull metadata_;
    std::string filePath_;
};

template<typename Voxel>
class LZ4SliceVolume : public LZ4SliceVolumeBase {
public:
    static LZ4SliceVolume open(std::string filePath);
    void deleteFromDisk() &&;

    LZ4SliceVolume(LZ4SliceVolume&& other);

    LZ4SliceVolume<Voxel>& operator=(LZ4SliceVolume<Voxel>&& other);
    //LZ4SliceVolume(const LZ4SliceVolume& other) = delete; //Disable

    std::unique_ptr<VolumeRAM> loadBaseSlab(size_t beginZ, size_t endZ /*exclusive*/) const;
    virtual std::unique_ptr<LZ4SliceVolumeBase> moveToHeap() &&;

    VolumeAtomic<Voxel> loadSlab(size_t beginZ, size_t endZ /*exclusive*/) const;
    VolumeAtomic<Voxel> loadSlice(size_t sliceNumber) const;
    void writeSlice(const VolumeAtomic<Voxel>& slice, size_t sliceNumber);
    void writeSlab(const VolumeAtomic<Voxel>& slice, size_t sliceNumber);
    LZ4WriteableSlab<Voxel> getWriteableSlice(size_t sliceNumber);

private:
    friend class LZ4SliceVolumeBuilder<Voxel>;
    LZ4SliceVolume(std::string filePath, LZ4SliceVolumeMetadata metadata);

    std::string getSliceFilePath(size_t sliceNum) const;
    tgt::svec3 getSliceDimensions() const;
    size_t getSliceMemorySize() const;
};


template<typename Voxel>
class LZ4SliceVolumeSliceCacher {
public:
    LZ4SliceVolumeSliceCacher(LZ4SliceVolume<Voxel>& volume);

    const VolumeAtomic<Voxel>& getSlice(size_t sliceNumber) const;
private:
    LZ4SliceVolume<Voxel>& volume_;
    mutable VolumeAtomic<Voxel> slice_;
    mutable size_t sliceNum_;
};


template<typename Voxel, uint64_t neighborhoodExtent>
class LZ4SliceVolumeReader {

    const static boost::optional<VolumeAtomic<Voxel>> NO_SLICE;
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
    boost::optional<Voxel> getVoxelRelative(tgt::ivec2 slicePos, int sliceOffset) const;

    void seek(int sliceNumber);
    void advance();

    int getCurrentZPos() const;
    const LZ4SliceVolume<Voxel>& getVolume() const;
private:
    boost::optional<VolumeAtomic<Voxel>> loadSliceFromVolume(int sliceNumber) const;

    const LZ4SliceVolume<Voxel>& volume_;
    std::array<boost::optional<VolumeAtomic<Voxel>>, neighborhoodSize> slices_;
    int pos_;
};

template<typename Voxel>
class LZ4SliceVolumeBuilder {
public:
    LZ4SliceVolumeBuilder(std::string filePath, LZ4SliceVolumeMetadata metadata);
    LZ4SliceVolumeBuilder(LZ4SliceVolumeBuilder&& other);

    LZ4WriteableSlab<Voxel> getNextWriteableSlice();
    LZ4WriteableSlab<Voxel> getNextWriteableSlab(size_t zSize);
    void pushSlice(const VolumeAtomic<Voxel>& slice);
    void fill(Voxel value);

    LZ4SliceVolume<Voxel> finalize() &&;
    tgt::svec3 getDimensions() const;

private:
    LZ4SliceVolume<Voxel> volumeInConstruction_;
    size_t numSlicesPushed_;
};

template<typename Voxel>
class LZ4SliceVolumeVoxelBuilder {
public:
    LZ4SliceVolumeVoxelBuilder(std::string filePath, LZ4SliceVolumeMetadata metadata);
    LZ4SliceVolumeVoxelBuilder(LZ4SliceVolumeVoxelBuilder&& other);

    void pushVoxel(Voxel v);
    void finalizeCurrentSlice();

    LZ4SliceVolume<Voxel> finalize() &&;
    tgt::svec3 getDimensions() const;

private:
    LZ4SliceVolumeBuilder<Voxel> builder_;
    VolumeAtomic<Voxel> currentSlice_;
    size_t numVoxelsPushed_;
};

template<typename Voxel, int SLAB_SIZE, int OVERLAP>
struct OverlappingSlabReader {
    OverlappingSlabReader(const LZ4SliceVolume<Voxel>& volume);
    boost::optional<Voxel> getVoxel(const tgt::ivec3& globalPos) const;
    void advance();

    const LZ4SliceVolume<Voxel>& volume_;
    VolumeAtomic<Voxel> currentSlab_;
    int currentSlabStart_;
};


LZ4SliceVolume<uint8_t> binarizeVolume(const VolumeBase& volume, float binarizationThresholdSegmentationNormalized, ProgressReporter* progress = nullptr);
LZ4SliceVolume<uint8_t> binarizeVolume(const VolumeBase& volume, float binarizationThresholdSegmentationNormalized, ProgressReporter& progress);
LZ4SliceVolume<uint8_t> binarizeVolume(const VolumeBase& volume, float binarizationThresholdSegmentationNormalized, ProgressReporter&& progress);

/// LZ4WriteableSlab -----------------------------------------------------------

template<typename Voxel>
LZ4WriteableSlab<Voxel>::LZ4WriteableSlab(LZ4SliceVolume<Voxel>& volume, size_t zBegin, VolumeAtomic<Voxel>&& slab)
    : volume_(volume)
    , zBegin_(zBegin)
    , slab_(std::move(slab))
{
}

template<typename Voxel>
LZ4WriteableSlab<Voxel>::~LZ4WriteableSlab() {
    if(slab_.getData() != nullptr) {
        volume_.writeSlab(slab_, zBegin_);
    }
}

/// LZ4SliceVolume -------------------------------------------------------------
template<typename Voxel>
LZ4SliceVolume<Voxel> LZ4SliceVolume<Voxel>::open(std::string filePath) {
    auto metadata = LZ4SliceVolumeMetadataFull::load(filePath);
    tgtAssert(metadata.getFormat() == getFormatFromType<Voxel>(), "Opened file with invalid format");
    return LZ4SliceVolume(filePath, metadata);
}
template<typename Voxel>
void LZ4SliceVolume<Voxel>::deleteFromDisk() && {
    LZ4SliceVolume dump = std::move(*this);
    for(size_t z = 0; z < dump.getNumSlices(); ++z) {
        tgt::FileSystem::deleteFile(dump.getSliceFilePath(z));
    }
    tgt::FileSystem::deleteFile(dump.filePath_);
}

template<typename Voxel>
std::string LZ4SliceVolume<Voxel>::getSliceFilePath(size_t sliceNum) const {
    return filePath_ + "_slice" + std::to_string(sliceNum);
}

template<typename Voxel>
tgt::svec3 LZ4SliceVolume<Voxel>::getSliceDimensions() const {
    return tgt::svec3(metadata_.getDimensions().xy(), 1);
}

template<typename Voxel>
size_t LZ4SliceVolume<Voxel>::getSliceMemorySize() const {
    return sizeof(Voxel) * tgt::hmul(getSliceDimensions());
}

template<typename Voxel>
LZ4SliceVolume<Voxel>::LZ4SliceVolume(LZ4SliceVolume<Voxel>&& other)
    : LZ4SliceVolumeBase(std::move(other.filePath_), std::move(other.metadata_))
{
}

template<typename Voxel>
LZ4SliceVolume<Voxel>& LZ4SliceVolume<Voxel>::operator=(LZ4SliceVolume<Voxel>&& other) {
    if(this != &other) {
        // Destruct the current object, but keep the memory.
        this->~LZ4SliceVolume();
        // Call the move constructor on the memory region of the current object.
        new(this) LZ4SliceVolume(std::move(other));
    }

    return *this;
}


template<typename Voxel>
LZ4SliceVolume<Voxel>::LZ4SliceVolume(std::string filePath, LZ4SliceVolumeMetadata metadata)
    : LZ4SliceVolumeBase(
            filePath,
            LZ4SliceVolumeMetadataFull(metadata, getFormatFromType<Voxel>(), getBaseTypeFromType<Voxel>())
            )
{
}

template<typename Voxel>
std::unique_ptr<VolumeRAM> LZ4SliceVolume<Voxel>::loadBaseSlab(size_t beginZ, size_t endZ) const {
    return std::unique_ptr<VolumeRAM>(new VolumeAtomic<Voxel>(loadSlab(beginZ, endZ)));
}

template<typename Voxel>
std::unique_ptr<LZ4SliceVolumeBase> LZ4SliceVolume<Voxel>::moveToHeap() && {
    return std::unique_ptr<LZ4SliceVolumeBase>(new LZ4SliceVolume<Voxel>(std::move(*this)));
}

template<typename Voxel>
VolumeAtomic<Voxel> LZ4SliceVolume<Voxel>::loadSlab(size_t beginZ, size_t endZ) const {
    tgtAssert(beginZ < endZ, "Invalid slab range");

    size_t dim_x = getDimensions().x;
    size_t dim_y = getDimensions().y;
    VolumeAtomic<Voxel> output(tgt::svec3(dim_x, dim_y, endZ - beginZ));

    for(size_t z=beginZ; z<endZ; ++z) {
        auto slice = loadSlice(z);
        auto sliceStart = &slice.voxel(0,0,0);
        std::copy(sliceStart, sliceStart+slice.getNumVoxels(), &output.voxel(0,0,z-beginZ));
    }

    return output;
}

template<typename Voxel>
VolumeAtomic<Voxel> LZ4SliceVolume<Voxel>::loadSlice(size_t sliceNumber) const {
    tgtAssert(sliceNumber < getDimensions().z, "Invalid slice number");

    std::string sliceFileName = getSliceFilePath(sliceNumber);
    std::ifstream compressedFile(sliceFileName, std::ifstream::binary);
    if(compressedFile.fail()) {
        throw std::system_error(errno, std::system_category(), "Failed to open lz4 slice file "+sliceFileName);
    }

    compressedFile.seekg(0,compressedFile.end); // Seek to end
    size_t compressedFileSize = compressedFile.tellg(); // end file position == file size
    compressedFile.clear(); // Clear all error state flags
    compressedFile.seekg(0, std::ios::beg); // Seek to beginning

    //Read file into buffer (buffer size is exactly the size of the file)
    std::unique_ptr<char[]> compressedBuffer(new char[compressedFileSize]);
    compressedFile.read(compressedBuffer.get(), compressedFileSize);
    compressedFile.close();

    // Decompress buffer
    const size_t sliceMemorySize = getSliceMemorySize();
    std::unique_ptr<char[]> decompressedBuffer(new char[sliceMemorySize]);
    size_t bytesDecompressed = LZ4_decompress_safe(compressedBuffer.get(), decompressedBuffer.get(), compressedFileSize, sliceMemorySize);
    tgtAssert(bytesDecompressed == sliceMemorySize, "Invalid memory size (resulting in memory corruption!)");

    return VolumeAtomic<Voxel>(reinterpret_cast<Voxel*>(decompressedBuffer.release()), getSliceDimensions());
}

template<typename Voxel>
void LZ4SliceVolume<Voxel>::writeSlice(const VolumeAtomic<Voxel>& slice, size_t sliceNumber) {
    tgtAssert(slice.getDimensions() == getSliceDimensions(), "Invalid slice dimensions");
    tgtAssert(sliceNumber < getNumSlices(), "Invalid slice number");

    const size_t sliceMemorySize = getSliceMemorySize();

    const size_t dstSize = LZ4_compressBound(sliceMemorySize);
    std::unique_ptr<char[]> compressedBuffer(new char[dstSize]);
    size_t compressedSize = LZ4_compress_default((const char *)(slice.getData()), compressedBuffer.get(), sliceMemorySize, dstSize);
    tgtAssert(compressedSize > 0, "Compression failed");

    std::string sliceFileName = getSliceFilePath(sliceNumber);
    std::ofstream outStream(sliceFileName, std::ofstream::binary | std::ofstream::trunc);
    outStream.write(compressedBuffer.get(), compressedSize);

    if(outStream.fail()) {
        throw std::system_error(errno, std::system_category(), "Failed writing lz4 slice file "+sliceFileName);
    }
}

template<typename Voxel>
void LZ4SliceVolume<Voxel>::writeSlab(const VolumeAtomic<Voxel>& slab, size_t zBegin) {
    tgtAssert(slab.getDimensions().xy() == getSliceDimensions().xy(), "Invalid slab dimensions");
    tgtAssert(zBegin + slab.getDimensions().z <= getNumSlices(), "Invalid slab begin");

    for(size_t dz=0; dz<slab.getDimensions().z; ++dz) {
        Voxel* sliceStart = const_cast<Voxel*>(&slab.voxel(0,0,dz)); //This is fine because we create a const volumeatomic from it
        const VolumeAtomic<Voxel> slice(sliceStart, getSliceDimensions(), false /*do not take ownership*/);
        writeSlice(slice, zBegin + dz);
    }
}

template<typename Voxel>
LZ4WriteableSlab<Voxel> LZ4SliceVolume<Voxel>::getWriteableSlice(size_t sliceNumber) {
    return LZ4WriteableSlab<Voxel>(*this, sliceNumber, loadSlice(sliceNumber));
}

/// LZ4SliceVolumeSliceCacher --------------------------------------------------

template<typename Voxel>
LZ4SliceVolumeSliceCacher<Voxel>::LZ4SliceVolumeSliceCacher(LZ4SliceVolume<Voxel>& volume)
    : volume_(volume)
    , slice_(volume_.loadSlice(0))
    , sliceNum_(0)
{
}

template<typename Voxel>
const VolumeAtomic<Voxel>& LZ4SliceVolumeSliceCacher<Voxel>::getSlice(size_t sliceNumber) const {
    if (sliceNum_ != sliceNumber) {
        sliceNum_ = sliceNumber;
        slice_ = volume_.loadSlice(sliceNum_);
    }

    return slice_;
}

/// LZ4SliceVolumeReader --------------------------------------------------

template<typename Voxel, uint64_t neighborhoodExtent>
const boost::optional<VolumeAtomic<Voxel>> LZ4SliceVolumeReader<Voxel, neighborhoodExtent>::NO_SLICE;

template<typename Voxel, uint64_t neighborhoodExtent>
LZ4SliceVolumeReader<Voxel, neighborhoodExtent>::LZ4SliceVolumeReader(const LZ4SliceVolume<Voxel>& volume)
    : volume_(volume)
    , slices_()
    , pos_(-(int)neighborhoodSize)
{
    std::fill(slices_.begin(), slices_.end(), boost::none);
}

template<typename Voxel, uint64_t neighborhoodExtent>
void LZ4SliceVolumeReader<Voxel, neighborhoodExtent>::seek(int newPos) {
    if(pos_ != newPos) {
        pos_ = newPos;

        for(size_t i=0; i<neighborhoodSize; ++i) {
            auto pos = pos_ + sliceStorageIndextoSlicePosOffset(i);
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
int LZ4SliceVolumeReader<Voxel, neighborhoodExtent>::getCurrentZPos() const {
    return pos_;
}

template<typename Voxel, uint64_t neighborhoodExtent>
const LZ4SliceVolume<Voxel>& LZ4SliceVolumeReader<Voxel, neighborhoodExtent>::getVolume() const {
    return volume_;
}

template<typename Voxel, uint64_t neighborhoodExtent>
const boost::optional<VolumeAtomic<Voxel>>& LZ4SliceVolumeReader<Voxel, neighborhoodExtent>::getSlice(int sliceNumber) const {
    int sliceStorageIndex = slicePosOffsetToSliceStorageIndex(sliceNumber - pos_);
    if(0 <= sliceStorageIndex && sliceStorageIndex < static_cast<int>(neighborhoodSize)) {
        return slices_[sliceStorageIndex];
    } else {
        return NO_SLICE;
    }
}

template<typename Voxel, uint64_t neighborhoodExtent>
boost::optional<Voxel> LZ4SliceVolumeReader<Voxel, neighborhoodExtent>::getVoxel(tgt::ivec3 pos) const {
    const auto& slice = getSlice(pos.z);
    tgt::ivec3 dim(volume_.getDimensions());
    if(slice) {
        if(pos.x < 0 || pos.x >= dim.x || pos.y < 0 || pos.y >= dim.y) {
            return boost::none;
        } else {
            return slice->voxel(pos.x, pos.y, 0);
        }
    } else {
        return boost::none;
    }
}
template<typename Voxel, uint64_t neighborhoodExtent>
boost::optional<Voxel> LZ4SliceVolumeReader<Voxel, neighborhoodExtent>::getVoxelRelative(tgt::ivec2 slicePos, int sliceOffset) const {
    tgtAssert(-static_cast<int>(neighborhoodExtent) <= sliceOffset && sliceOffset <= static_cast<int>(neighborhoodExtent), "Invalid slice offset");
    const auto& slice = getSlice(pos_ + sliceOffset);
    if(slice) {
        if(slicePos.x < 0 || slicePos.x >= (int)volume_.getDimensions().x || slicePos.y < 0 || slicePos.y >= (int)volume_.getDimensions().y) {
            return boost::none;
        } else {
            return slice->voxel(slicePos.x, slicePos.y, 0);
        }
    } else {
        return boost::none;
    }
}

template<typename Voxel, uint64_t neighborhoodExtent>
boost::optional<VolumeAtomic<Voxel>> LZ4SliceVolumeReader<Voxel, neighborhoodExtent>::loadSliceFromVolume(int sliceNumber) const {
    if(0 <= sliceNumber && sliceNumber < static_cast<int>(volume_.getNumSlices())) {
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
LZ4WriteableSlab<Voxel> LZ4SliceVolumeBuilder<Voxel>::getNextWriteableSlice() {
    return getNextWriteableSlab(1);
}

template<typename Voxel>
LZ4WriteableSlab<Voxel> LZ4SliceVolumeBuilder<Voxel>::getNextWriteableSlab(size_t zSize) {
    tgtAssert(numSlicesPushed_ < volumeInConstruction_.getNumSlices(), "Cannot push more slices");
    size_t begin = numSlicesPushed_;
    numSlicesPushed_ += zSize;
    return LZ4WriteableSlab<Voxel>(volumeInConstruction_, begin, VolumeAtomic<Voxel>(tgt::svec3(volumeInConstruction_.getDimensions().xy(), zSize)));
}

template<typename Voxel>
void LZ4SliceVolumeBuilder<Voxel>::pushSlice(const VolumeAtomic<Voxel>& slice) {
    tgtAssert(numSlicesPushed_ < volumeInConstruction_.getNumSlices(), "Cannot push more slices");
    ++numSlicesPushed_;
    return volumeInConstruction_.writeSlice(slice, numSlicesPushed_-1);
}
template<typename Voxel>
void LZ4SliceVolumeBuilder<Voxel>::fill(Voxel value) {
    while(numSlicesPushed_ != volumeInConstruction_.getNumSlices()) {
        auto slice = getNextWriteableSlice();
        slice->fill(value);
    }
}

template<typename Voxel>
LZ4SliceVolume<Voxel> LZ4SliceVolumeBuilder<Voxel>::finalize() && {
    tgtAssert(numSlicesPushed_ == volumeInConstruction_.getNumSlices(), "Invalid number of slices pushed");
    auto tmp = std::move(*this);
    tmp.volumeInConstruction_.metadata_.save(tmp.volumeInConstruction_.filePath_);
    return std::move(tmp.volumeInConstruction_);
}
template<typename Voxel>
tgt::svec3 LZ4SliceVolumeBuilder<Voxel>::getDimensions() const {
    return volumeInConstruction_.getDimensions();
}


/// LZ4SliceVolumeVoxelBuilder ------------------------------------------------------
template<typename Voxel>
LZ4SliceVolumeVoxelBuilder<Voxel>::LZ4SliceVolumeVoxelBuilder(std::string filePath, LZ4SliceVolumeMetadata metadata)
    : builder_(filePath, metadata)
    , currentSlice_(tgt::svec3(metadata.getDimensions().xy(), 1))
    , numVoxelsPushed_(0)
{
}

template<typename Voxel>
LZ4SliceVolumeVoxelBuilder<Voxel>::LZ4SliceVolumeVoxelBuilder(LZ4SliceVolumeVoxelBuilder&& other)
    : builder_(std::move(other.builder_))
    , currentSlice_(std::move(other.currentSlice_))
    , numVoxelsPushed_(other.numVoxelsPushed_)
{
    other.numVoxelsPushed_ = -1;
}

template<typename Voxel>
void LZ4SliceVolumeVoxelBuilder<Voxel>::finalizeCurrentSlice() {
    tgtAssert(numVoxelsPushed_ == currentSlice_.getNumVoxels(), "too many voxels pushed to slice");
    builder_.pushSlice(currentSlice_);
    currentSlice_.clear();
    numVoxelsPushed_ = 0;
}

template<typename Voxel>
void LZ4SliceVolumeVoxelBuilder<Voxel>::pushVoxel(Voxel voxel) {
    tgtAssert(numVoxelsPushed_ < currentSlice_.getNumVoxels(), "too many voxels pushed to slice");
    currentSlice_.voxel(numVoxelsPushed_) = voxel;
    ++numVoxelsPushed_;
    if(numVoxelsPushed_ == currentSlice_.getNumVoxels()) {
        finalizeCurrentSlice();
    }
}

template<typename Voxel>
LZ4SliceVolume<Voxel> LZ4SliceVolumeVoxelBuilder<Voxel>::finalize() && {
    tgtAssert(numVoxelsPushed_ == 0, "Unfinished slice");
    auto tmp = std::move(*this);
    return std::move(tmp.builder_).finalize();
}
template<typename Voxel>
tgt::svec3 LZ4SliceVolumeVoxelBuilder<Voxel>::getDimensions() const {
    return builder_.getDimensions();
}


/// OverlappingSlabReader ------------------------------------------------------

template<typename Voxel, int SLAB_SIZE, int OVERLAP>
OverlappingSlabReader<Voxel, SLAB_SIZE, OVERLAP>::OverlappingSlabReader(const LZ4SliceVolume<Voxel>& volume)
    : volume_(volume)
    , currentSlab_(tgt::svec3(0), false)
    , currentSlabStart_(-SLAB_SIZE)
{
}
template<typename Voxel, int SLAB_SIZE, int OVERLAP>
boost::optional<Voxel> OverlappingSlabReader<Voxel, SLAB_SIZE, OVERLAP>::getVoxel(const tgt::ivec3& globalPos) const {
    int slabZ = globalPos.z - currentSlabStart_ +
        (currentSlabStart_ == 0 ?
        0 :
        /*we have one additional slice! */ OVERLAP);
    if(     0 <= slabZ && slabZ < static_cast<int>(currentSlab_.getDimensions().z)
         && 0 <= globalPos.x && globalPos.x < static_cast<int>(currentSlab_.getDimensions().x)
         && 0 <= globalPos.y && globalPos.y < static_cast<int>(currentSlab_.getDimensions().y))
    {
        tgt::ivec3 pos(globalPos.xy(), slabZ);
        return currentSlab_.voxel(pos);
    } else {
        return boost::none;
    }
}

template<typename Voxel, int SLAB_SIZE, int OVERLAP>
void OverlappingSlabReader<Voxel, SLAB_SIZE, OVERLAP>::advance() {
    currentSlabStart_ += SLAB_SIZE;
    int begin = std::max(currentSlabStart_ - OVERLAP, 0);
    int end = std::min(
            currentSlabStart_+SLAB_SIZE+ 2*OVERLAP,
            static_cast<int>(volume_.getDimensions().z));

    // First free the old slab, then allocate new one (to reduce peak memory usage)
    currentSlab_ = VolumeAtomic<Voxel>(tgt::svec3(0), false);
    currentSlab_ = volume_.loadSlab(begin, end);
}

}
