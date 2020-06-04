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

#ifndef VRN_SLICEREADER_H
#define VRN_SLICEREADER_H

#include "volumefilter.h"

#include "voreen/core/datastructures/volume/volumeminmax.h"
#include "voreen/core/datastructures/volume/volumebase.h"
#include "voreen/core/datastructures/volume/volumefactory.h"
#include "modules/hdf5/io/hdf5filevolume.h"
#include "voreen/core/io/progressreporter.h"

#include <functional>

namespace voreen {

class SliceReaderMetaData {
public:
    // Does _not_ copy min/max values, since most of the time they are not valid for a filter/reader on top
    static SliceReaderMetaData fromBase(const SliceReaderMetaData& base);
    // Copies min/max values if available. Use only when you are sure that the previous min/max values are accurate!
    static SliceReaderMetaData fromBaseAccurate(const SliceReaderMetaData& base);
    static SliceReaderMetaData fromVolume(const VolumeBase& vol);
    static SliceReaderMetaData fromHDF5Volume(const HDF5FileVolume& volume);


    SliceReaderMetaData(RealWorldMapping rwm, tgt::svec3 dimensions, size_t numChannels, std::string baseType, tgt::vec3 spacing);
    SliceReaderMetaData() = delete;
    SliceReaderMetaData(const SliceReaderMetaData&) = delete;
    SliceReaderMetaData& operator=(const SliceReaderMetaData&) = delete;
    SliceReaderMetaData(SliceReaderMetaData&&) = default;
    SliceReaderMetaData& operator=(SliceReaderMetaData&&) = default;

    void setMinMax(std::vector<tgt::vec2> minmax);
    void setMinMaxNormalized(std::vector<tgt::vec2> minmaxNorm);
    void setDimensions(tgt::svec3 dimensions);
    void setNumChannels(size_t numChannels);
    void setBaseType(std::string baseType);
    void setRealWorldMapping(RealWorldMapping rwm);
    void setSpacing(tgt::vec3 spacing);

    // Have a best guess for valid min/max values. If available, all actual
    // min/max channels are considered. Otherwise the RWM will be used as a
    // fallback.
    tgt::vec2 estimateMinMax() const;

    // Result may be null (no valid min/max data available)
    std::unique_ptr<VolumeMinMax> getVolumeMinMax() const;
    const boost::optional<std::vector<tgt::vec2>>& getMinMax() const;
    const RealWorldMapping& getRealworldMapping() const;
    tgt::svec3 getDimensions() const;
    size_t getNumChannels() const;
    std::string getBaseType() const;
    RealWorldMapping getRealWorldMapping() const;
    tgt::vec3 getSpacing() const;

private:
    RealWorldMapping rwm_;
    tgt::svec3 dimensions_;
    boost::optional<std::vector<tgt::vec2>> minmax_;
    size_t numChannels_;
    std::string baseType_;
    tgt::vec3 spacing_;
};

class SliceReader {
public:
    SliceReader(const tgt::ivec3& signedDim, SliceReaderMetaData&& metadata);
    virtual ~SliceReader() {}

    virtual void advance() = 0;
    virtual void seek(int z) = 0;
    virtual int getCurrentZPos() const = 0;
    virtual const VolumeRAM* getCurrentSlice() const = 0;
    virtual std::string getBaseType() const = 0;
    virtual size_t getNumChannels() const = 0;

    const tgt::ivec3& getSignedDimensions() const;
    tgt::svec3 getDimensions() const;
    const SliceReaderMetaData& getMetaData() const;

    virtual float getVoxelNormalized(const tgt::ivec3& xyz, size_t channel = 0) const = 0;

protected:
    tgt::ivec3 dim_;
    SliceReaderMetaData metadata_;
};


class CachingSliceReader : public SliceReader {
public:
    CachingSliceReader(std::unique_ptr<SliceReader>&& base, int neighborhoodSize);

    virtual ~CachingSliceReader();
    void advance();
    void seek(int z);
    int getCurrentZPos() const;

    float getVoxelNormalized(const tgt::ivec3& xyz, size_t channel = 0) const;
    const VolumeRAM* getCurrentSlice() const;
    std::string getBaseType() const;
    size_t getNumChannels() const;

    int getZExtent() const;

    VolumeRAM* const& getSlice(int dz) const;
protected:
    VolumeRAM*& getSlice(int dz);

    std::unique_ptr<SliceReader> base_;
    std::vector<VolumeRAM*> slices_;
    const int neighborhoodSize_; // has to be >= 0, but is declared (signed) int to avoid lots of casting
};


class VolumeSliceReader : public SliceReader {
public:
    VolumeSliceReader(const VolumeBase& volume);

    virtual ~VolumeSliceReader();
    void advance();
    void seek(int z);
    int getCurrentZPos() const;
    const VolumeRAM* getCurrentSlice() const;
    std::string getBaseType() const;
    float getVoxelNormalized(const tgt::ivec3& xyz, size_t channel = 0) const;
    size_t getNumChannels() const;

protected:

    const VolumeBase& volume_;
    int currentZPos_;
    std::unique_ptr<VolumeRAM> currentSlice_;
    size_t numChannels_; // cache, because getting it from volume can be veeeeeeeery slow
};

class HDF5VolumeSliceReader : public SliceReader {
public:
    HDF5VolumeSliceReader(const HDF5FileVolume& volume);

    virtual ~HDF5VolumeSliceReader();
    void advance();
    void seek(int z);
    int getCurrentZPos() const;
    const VolumeRAM* getCurrentSlice() const;
    std::string getBaseType() const;
    float getVoxelNormalized(const tgt::ivec3& xyz, size_t channel = 0) const;
    size_t getNumChannels() const;

protected:

    const HDF5FileVolume& volume_;
    int currentZPos_;
    std::unique_ptr<VolumeRAM> currentSlice_;
    size_t numChannels_; // cache, because getting it from volume can be veeeeeeeery slow
};


class FilteringSliceReader : public SliceReader {
public:
    FilteringSliceReader(std::unique_ptr<CachingSliceReader>, std::unique_ptr<VolumeFilter>);
    virtual ~FilteringSliceReader() {}

    void advance();
    void seek(int z);
    int getCurrentZPos() const;

    float getVoxelNormalized(const tgt::ivec3& xyz, size_t channel = 0) const;
    const VolumeRAM* getCurrentSlice() const;
    std::string getBaseType() const;
    size_t getNumChannels() const;

protected:
    void updateCurrentSlice();
    int nearestBaseZ(int thisZ);

    std::unique_ptr<CachingSliceReader> base_;
    std::unique_ptr<VolumeRAM> currentSlice_;
    std::unique_ptr<VolumeFilter> filter_;
    int z_;
    float thisToBaseScale_;
    float thisToBaseOffset_;
};

// VolumeFilterStackBuilder -------------------------------------------------------------------------------------

class VolumeFilterStackBuilder {
public:
    VolumeFilterStackBuilder(const VolumeBase& volume);

    VolumeFilterStackBuilder& addLayer(std::unique_ptr<VolumeFilter> conv);
    std::unique_ptr<SliceReader> build(int initZPos);

    std::unique_ptr<CachingSliceReader> buildCaching(int initZPos, int neighborhoodSize);

private:
    std::unique_ptr<SliceReader> top_;
};

// Utility functions -------------------------------------------------------------------------------------

void writeSlicesToHDF5File(SliceReader& reader, HDF5FileVolume& file, ProgressReporter* progress = nullptr);

} //namespace voreen

#endif // VRN_SLICEREADER_H
