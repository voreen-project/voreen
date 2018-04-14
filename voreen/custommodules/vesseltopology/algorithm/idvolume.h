/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2015 University of Muenster, Germany.                        *
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

#include <string>
#include <iostream>
#include <fstream>
#include <boost/iostreams/device/mapped_file.hpp>

#include "tgt/vector.h"
#include "surface.h"
#include "../datastructures/protovesselgraph.h"
#include "voreen/core/datastructures/volume/volumebase.h"

namespace voreen {

class IdVolumeStorage;
class IdVolumeStorageInitializer;
struct IdVolumeInitializationReader;

class IdVolume {
public:
    //Usable Values are in the range of [0,UNLABELED_FOREGROUND_VALUE)!
    typedef uint32_t Value;

    const static Value BACKGROUND_VALUE;
    const static Value UNLABELED_FOREGROUND_VALUE;

    IdVolume(IdVolume&&);
    //IdVolume(IdVolumeInitializationReader& initializationReader);
    IdVolume(IdVolumeStorageInitializer&& storage, StoredSurface surface, tgt::svec3 dimensions, size_t numUnlabeledForegroundVoxels);
    ~IdVolume();

    void floodFromLabels(ProgressReporter& progress, size_t maxIt);
    void floodIteration(size_t& numberOfFloodedVoxels, ProgressReporter& progress);

    tgt::svec3 getDimensions() const;

    //static std::unique_ptr<VolumeBase> intoVolume(IdVolume&&);

    static const std::string loggerCat_;

    std::unique_ptr<IdVolumeStorage> data_; // Never null
private:
    StoredSurface surfaceFile_;
    size_t numUnlabeledForegroundVoxels_;
};


class IdVolumeStorage {
public:
    IdVolumeStorage(IdVolumeStorageInitializer&& initializer, tgt::svec3 dimensions);

    ~IdVolumeStorage();

    void set(const tgt::svec3& pos, IdVolume::Value val);

    IdVolume::Value get(const tgt::svec3& pos) const;

    tgt::svec3 dimensions_;
    std::string filename_;
private:
    boost::iostreams::mapped_file file_;
};


struct IdVolumeInitializationReader {
    IdVolumeInitializationReader(const HDF5FileVolume& branchIds, const HDF5FileVolume& holeIds)
        : holeIds_(holeIds)
        , branchIds_(branchIds)
        , holeIdSlices_({nullptr, nullptr, nullptr})
        , branchIdSlices_({nullptr, nullptr, nullptr})
        , z_(-1)
    {
        tgtAssert(holeIds.getDimensions() == branchIds.getDimensions(), "Invalid volume dimensions");
        tgtAssert(holeIds.getBaseType() == "uint32", "Invalid volume format");
        tgtAssert(branchIds.getBaseType() == "uint32", "Invalid volume format");

        seek(-1);
    }

    bool isLabeled(const tgt::ivec3& pos) const {
        return getBranchId(pos) != 0
            && getHoleId(pos) == 0; // holeIdBrick might override regions that are cut off and must be relabeled
    }

    bool isObject(const tgt::ivec3& pos, uint32_t idToFill) const {
        return isLabeled(pos) || getHoleId(pos) == idToFill;
    }

    uint32_t getBranchId(const tgt::ivec3& pos) const {
        int sid = pos.z-z_+ 1;
        tgtAssert(0 <= sid && sid < 3, "invalid z pos");
        return branchIdSlices_[sid]->voxel(tgt::svec3(pos.xy(), 0));
    }

    uint32_t getHoleId(const tgt::ivec3& pos) const {
        int sid = pos.z-z_+ 1;
        tgtAssert(0 <= sid && sid < 3, "invalid z pos");
        return holeIdSlices_[sid]->voxel(tgt::svec3(pos.xy(), 0));
    }

    void seek(int z) {
        z_ = z;

        holeIdSlices_[0].reset(loadSlice(z_-1, holeIds_));
        holeIdSlices_[1].reset(loadSlice(z_  , holeIds_));
        holeIdSlices_[2].reset(loadSlice(z_+1, holeIds_));

        branchIdSlices_[0].reset(loadSlice(z_-1, branchIds_));
        branchIdSlices_[1].reset(loadSlice(z_  , branchIds_));
        branchIdSlices_[2].reset(loadSlice(z_+1, branchIds_));
    }

    void advance() {
        ++z_;
        std::swap(holeIdSlices_[1], holeIdSlices_[0]);
        std::swap(holeIdSlices_[2], holeIdSlices_[1]);

        std::swap(branchIdSlices_[1], branchIdSlices_[0]);
        std::swap(branchIdSlices_[2], branchIdSlices_[1]);

        holeIdSlices_[2].reset(loadSlice(z_+1, holeIds_));
        branchIdSlices_[2].reset(loadSlice(z_+1, branchIds_));
    }

    tgt::svec3 getDimensions() const {
        return holeIds_.getDimensions();
    }

private:
    VolumeRAM_UInt32* loadSlice(int z, const HDF5FileVolume& source) {
        auto dim = getDimensions();
        if(0 <= z && z < dim.z) {
            tgt::svec3 sliceSize(dim.xy(), 1);
            tgt::svec3 sliceOffset(0,0,z);
            auto slice = dynamic_cast<VolumeRAM_UInt32*>(source.loadBrick(sliceOffset, sliceSize));

            tgtAssert(slice, "Invalid holeIdVolume format");
            return slice;
        } else {
            return nullptr;
        }
    }

    const HDF5FileVolume& holeIds_;
    const HDF5FileVolume& branchIds_;
    std::array<std::unique_ptr<VolumeRAM_UInt32>, 3> holeIdSlices_;
    std::array<std::unique_ptr<VolumeRAM_UInt32>, 3> branchIdSlices_;
    int z_;
};

struct IdVolumeReader {
    IdVolumeReader(const IdVolume& vol)
        : idVolume_(vol)
        , z_(0)
    {
    }

    uint64_t getId(const tgt::svec3& pos) const {
        return idVolume_.data_->get(pos);
    }

    void advance() {
        ++z_;
    }

    const IdVolume& idVolume_;
    int z_;
};

class IdVolumeStorageInitializer {
public:
    IdVolumeStorageInitializer(std::string filename);

    IdVolumeStorageInitializer(IdVolumeStorageInitializer&& other);

    ~IdVolumeStorageInitializer();

    void push(IdVolume::Value val);

    const std::string filename_;

private:
    std::ofstream file_;
};

}
