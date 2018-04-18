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

#ifndef VRN_SLICEREADER_H
#define VRN_SLICEREADER_H

#include "volumefilter.h"

#include "voreen/core/datastructures/volume/volumebase.h"
#include "voreen/core/datastructures/volume/volumefactory.h"
#include "modules/hdf5/io/hdf5filevolume.h"
#include "voreen/core/io/progressreporter.h"

#include <functional>

namespace voreen {


class SliceReader {
public:
    SliceReader(const tgt::ivec3& signedDim);
    virtual ~SliceReader() {}

    virtual void advance() = 0;
    virtual void seek(int z) = 0;
    virtual int getCurrentZPos() const = 0;
    virtual const tgt::svec3& getDimensions() const = 0;
    virtual const VolumeRAM* getCurrentSlice() const = 0;
    virtual std::string getBaseType() const = 0;
    virtual size_t getNumChannels() const = 0;

    const tgt::ivec3& getSignedDimensions() const;

    virtual float getVoxelNormalized(const tgt::ivec3& xyz, size_t channel = 0) const = 0;

protected:
    tgt::ivec3 dim_; // cache to avoid vtable lookups
};


class CachingSliceReader : public SliceReader {
public:
    CachingSliceReader(std::unique_ptr<SliceReader>&& base, int neighborhoodSize);

    virtual ~CachingSliceReader();
    void advance();
    void seek(int z);
    int getCurrentZPos() const;
    const tgt::svec3& getDimensions() const;

    float getVoxelNormalized(const tgt::ivec3& xyz, size_t channel = 0) const;
    const VolumeRAM* getCurrentSlice() const;
    std::string getBaseType() const;
    size_t getNumChannels() const;

    int getZExtent() const;

protected:
    VolumeRAM*& getSlice(int dz);
    VolumeRAM* const& getSlice(int dz) const;

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
    const tgt::svec3& getDimensions() const;
    const VolumeRAM* getCurrentSlice() const;
    std::string getBaseType() const;
    float getVoxelNormalized(const tgt::ivec3& xyz, size_t channel = 0) const;
    size_t getNumChannels() const;

protected:

    const VolumeBase& volume_;
    tgt::svec3 dimensions_;
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
    const tgt::svec3& getDimensions() const;

    float getVoxelNormalized(const tgt::ivec3& xyz, size_t channel = 0) const;
    const VolumeRAM* getCurrentSlice() const;
    std::string getBaseType() const;
    size_t getNumChannels() const;

protected:
    void updateCurrentSlice();

    std::unique_ptr<CachingSliceReader> base_;
    std::unique_ptr<VolumeRAM> currentSlice_;
    std::unique_ptr<VolumeFilter> filter_;
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
