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

#include "voreen/core/io/serialization/serializable.h"
#include "voreen/core/io/serialization/xmldeserializer.h"
#include <vector>
#include <string>

#include "tgt/vector.h"

#include <lz4.h>

namespace voreen {

class LZ4SliceVolumeMetadata : public Serializable {
public:
    LZ4SliceVolumeMetadata(tgt::svec3 dimensions);
    LZ4SliceVolumeMetadata(const std::string& xmlfile);

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);

    void save(const std::string& xmlfile) const;

    tgt::svec4 dimensions_;
};

template<typename Voxel>
class LZ4SliceVolumeBuilder;

template<typename Voxel>
class LZ4SliceVolume;

template<typename Voxel>
class LZ4WriteableSlice {
public:
    VolumeAtomic<Voxel>& operator*() { return *slice_; }
    VolumeAtomic<Voxel>* operator->() { return slice_; }

    LZ4WriteableSlice(LZ4WriteableSlice&& other)
        : volume_(other.volume_)
        , sliceNum_(other.sliceNum)
        , slice_(std::move(other.slice_))
    {
    }
    ~LZ4WriteableSlice();

private:
    LZ4WriteableSlice(LZ4SliceVolume<Voxel>& volume, size_t sliceNum);

    friend class LZ4SliceVolume<Voxel>;
    LZ4SliceVolume<Voxel>& volume_; //MUST live longer than this object
    size_t sliceNum_;
    VolumeAtomic<Voxel> slice_;
};

template<typename Voxel>
class LZ4SliceVolume {
public:
    static LZ4SliceVolume open(std::string filePath);

    LZ4SliceVolume(LZ4SliceVolume&& other);

    VolumeAtomic<Voxel> loadSlice(size_t sliceNumber) const;
    void writeSlice(const VolumeAtomic<Voxel>& slice, size_t sliceNumber);
    LZ4WriteableSlice<Voxel> getWritableSlice(size_t sliceNumber);

private:
    friend class LZ4SliceVolumeBuilder<Voxel>;
    LZ4SliceVolume(std::string filePath, LZ4SliceVolumeMetadata metadata);
    LZ4SliceVolume(const LZ4SliceVolume& other); //Disable

    std::string getSliceFilePath(size_t sliceNum) const;
    tgt::svec3 getSliceDimensions() const;
    size_t getSliceMemorySize() const;
    size_t getNumSlices() const;

    LZ4SliceVolumeMetadata metadata_;
    std::string filePath_;
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
LZ4WriteableSlice<Voxel>::LZ4WriteableSlice(LZ4SliceVolume<Voxel>& volume, size_t sliceNum)
    : volume_(volume)
    , sliceNum_(sliceNum)
    , slice_(volume_.loadSlice(sliceNum_))
{
}

template<typename Voxel>
LZ4WriteableSlice<Voxel>::~LZ4WriteableSlice() {
    volume_.writeSlice(slice_);
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
    other.filePath = "";
}

template<typename Voxel>
LZ4SliceVolume<Voxel>::LZ4SliceVolume(std::string filePath, LZ4SliceVolumeMetadata metadata)
    : filePath_(filePath)
    , metadata_(metadata)
{
}

template<typename Voxel>
LZ4SliceVolume<Voxel>::LZ4SliceVolume(const LZ4SliceVolume&)
{
    tgtAssert(false, "Tried to copy LZ4SliceVolume");
}

template<typename Voxel>
VolumeAtomic<Voxel> LZ4SliceVolume<Voxel>::loadSlice(size_t sliceNumber) const {
    std::ifstream compressedFile(getSliceName(filePath_, sliceNumber), std::ifstream::binary);

    const size_t sliceMemorySize = getSliceMemorySize();

#ifdef VRN_DEBUG
    //Get compressed file size for buffer
    compressedFile.seekg(0,compressedFile.end);
    size_t compressedFileSize = compressedFile.tellg();
    compressedFile.clear(); //?
    compressedFile.seekg(0, std::ios::beg);

    tgtAssert(compressedFileSize < sliceMemorySize, "Invalid slice memory size");
#endif

    //Read file into buffer (sliceMemorySize is larger than it needs to be, but this way we do not need to check the actual file size)
    std::unique_ptr<char[]> compressedBuffer(new char[sliceMemorySize]);
    compressedFile.read(compressedBuffer.get(), sliceMemorySize); //TODO check if this is valid (because sliceMemorySize may be larger than the file
    compressedFile.close();

    // Decompress buffer
    std::unique_ptr<char[]> decompressedBuffer(new char[sliceMemorySize]);
    size_t bytesDecompressed = LZ4_decompress_fast(compressedBuffer.get(), decompressedBuffer.get(), sliceMemorySize);
    tgtAssert(bytesDecompressed == sliceMemorySize, "Invalid memory size (resulting in memory corruption!)");

    return VolumeAtomic<Voxel>(decompressedBuffer.release(), getSliceDimensions());
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
    return LZ4WriteableSlice<Voxel>(*this, sliceNumber);
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
    return volumeInConstruction_.getWritableSlice(numSlicesPushed_-1);
}

template<typename Voxel>
LZ4SliceVolume<Voxel> LZ4SliceVolumeBuilder<Voxel>::finalize(LZ4SliceVolumeBuilder<Voxel>&& builder) {
    auto tmp = std::move(builder);
    tmp.volumeInConstruction_.metadata_.save();
    return tmp.volumeInConstruction_;
}


}
