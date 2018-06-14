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
    IdVolume(IdVolumeStorage&& storage, StoredSurface surface, size_t numUnlabeledForegroundVoxels);
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
    IdVolumeStorage(const LZ4SliceVolume<IdVolume::Value>& compressedData, std::string filename);
    IdVolumeStorage(IdVolumeStorage&& other);

    ~IdVolumeStorage();

    void set(const tgt::svec3& pos, IdVolume::Value val);

    IdVolume::Value get(const tgt::svec3& pos) const;
    VolumeAtomic<IdVolume::Value> getSlice(size_t z) const;

    tgt::svec3 dimensions_;
    std::string filename_;
private:
    friend class IdVolume;
    boost::iostreams::mapped_file file_;
};


struct IdVolumeInitializationReader {
    IdVolumeInitializationReader(const LZ4SliceVolume<uint32_t>& branchIds, const LZ4SliceVolume<uint32_t>& holeIds)
        : holeIds_(holeIds)
        , branchIds_(branchIds)
    {
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
        auto id = branchIds_.getVoxel(pos);
        tgtAssert(id, "Invalid voxel pos");
        return *id;
    }

    uint32_t getHoleId(const tgt::ivec3& pos) const {
        auto id = holeIds_.getVoxel(pos);
        tgtAssert(id, "Invalid voxel pos");
        return *id;
    }

    void seek(int z) {
        branchIds_.seek(z);
        holeIds_.seek(z);
    }

    void advance() {
        branchIds_.advance();
        holeIds_.advance();
    }

    tgt::svec3 getDimensions() const {
        return holeIds_.getVolume().getDimensions();
    }

private:

    LZ4SliceVolumeReader<uint32_t, 1> holeIds_;
    LZ4SliceVolumeReader<uint32_t, 1> branchIds_;
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
    void push(IdVolume::Value* vals, size_t number_of_vals);

    const std::string filename_;

private:
    std::ofstream file_;
};

}
