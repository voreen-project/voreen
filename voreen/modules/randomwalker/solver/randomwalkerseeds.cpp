/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#include "randomwalkerseeds.h"

#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "tgt/vector.h"

inline size_t volumeCoordsToIndex(int x, int y, int z, const tgt::ivec3& dim) {
    return z*dim.y*dim.x + y*dim.x + x;
}

inline size_t volumeCoordsToIndex(const tgt::ivec3& coords, const tgt::ivec3& dim) {
    return coords.z*dim.y*dim.x + coords.y*dim.x + coords.x;
}

inline tgt::ivec3 indexToVolumeCoords(size_t index, const tgt::ivec3& dim) {
    tgt::ivec3 result;
    result.z = static_cast<int>(index)/(dim.y*dim.x);
    result.y = (static_cast<int>(index)-(result.z*dim.y*dim.x)) / dim.x;
    result.x = index % dim.x;
    return result;
}


namespace voreen {

void RandomWalkerSeeds::initialize(const VolumeRAM* volume) {
    volume_ = volume;
    volDim_ = volume->getDimensions();
    numVoxels_ = volume->getNumVoxels();
}

size_t RandomWalkerSeeds::getNumSeeds() const {
    return numSeeds_;
}

tgt::vec2 RandomWalkerSeeds::getSeedRange() const {
    return seedRange_;
}

//---------------------------------------------------------------------------------------

RandomWalkerTwoLabelSeeds::RandomWalkerTwoLabelSeeds(
        const PointSegmentListGeometryVec3& foregroundSeedList,
        const PointSegmentListGeometryVec3& backgroundSeedList,
        const VolumeRAM_UInt8* foregroundSeedVolume,
        const VolumeRAM_UInt8* backgroundSeedVolume,
        const tgt::ivec3& clipLLF,
        const tgt::ivec3& clipURB) :
    foregroundSeedList_(foregroundSeedList),
    backgroundSeedList_(backgroundSeedList),
    foregroundSeedVolume_(foregroundSeedVolume),
    backgroundSeedVolume_(backgroundSeedVolume),
    clipLLF_(clipLLF),
    clipURB_(clipURB),
    numForegroundSeeds_(0),
    numBackgroundSeeds_(0),
    seedBuffer_(0)
{ }

RandomWalkerTwoLabelSeeds::~RandomWalkerTwoLabelSeeds() {
    delete[] seedBuffer_;
}

bool RandomWalkerTwoLabelSeeds::isSeedPoint(size_t index) const {
        return seedBuffer_[index];
}

bool RandomWalkerTwoLabelSeeds::isSeedPoint(const tgt::ivec3& voxel) const {
    return isSeedPoint(volumeCoordsToIndex(voxel, volDim_));
}

float RandomWalkerTwoLabelSeeds::getSeedValue(size_t index) const {
    if (seedBuffer_[index] == 2)
        return 1.f;
    else if (seedBuffer_[index] == 1)
        return 0.f;
    else
        return -1.f;
}

float RandomWalkerTwoLabelSeeds::getSeedValue(const tgt::ivec3& voxel) const {
    return getSeedValue(volumeCoordsToIndex(voxel, volDim_));
}

void RandomWalkerTwoLabelSeeds::setClippingPlanes(const tgt::ivec3& clipLLF, const tgt::ivec3& clipURB) {
    clipLLF_ = clipLLF;
    clipURB_ = clipURB;
}

size_t RandomWalkerTwoLabelSeeds::getNumForegroundSeeds() const {
    return numForegroundSeeds_;
}

size_t RandomWalkerTwoLabelSeeds::getNumBackgroundSeeds() const {
    return numBackgroundSeeds_;
}

void RandomWalkerTwoLabelSeeds::initialize(const VolumeRAM* volume) {
    RandomWalkerSeeds::initialize(volume);

    // create seed buffer
    try {
        seedBuffer_ = new char[numVoxels_];
    }
    catch (std::bad_alloc&) {
        throw VoreenException("Bad allocation during creation of seed buffer");
    }

    // initialize seed buffer
    for (size_t i=0; i<numVoxels_; i++) {
        seedBuffer_[i] = 0;
    }

    size_t tNumForegroundSeeds = 0;
    size_t tNumBackgroundSeeds = 0;

    // apply clipping planes: define all voxels outside of roi as background seeds
    if (clipLLF_ != tgt::ivec3(-1) && clipURB_ != tgt::ivec3(-1)) {
        tgt::ivec3 voxel;
        for (voxel.z=0; voxel.z<volDim_.z; voxel.z++) {
            for (voxel.y=0; voxel.y<volDim_.y; voxel.y++) {
                for (voxel.x=0; voxel.x<volDim_.x; voxel.x++) {

                    if (tgt::hor(tgt::lessThan(voxel, clipLLF_)) || tgt::hor(tgt::greaterThan(voxel, clipURB_))) {
                        size_t index = volumeCoordsToIndex(voxel, volDim_);
                        tgtAssert(index < numVoxels_, "Invalid index");
                        if (!seedBuffer_[index]) {
                            seedBuffer_[index] = 1; //< background seed
                            tNumBackgroundSeeds++;
                        }
                    }
                }
            }
        }
    }

    // foreground geometry seeds
    for (int m=0; m<foregroundSeedList_.getNumSegments(); m++) {
        const std::vector<tgt::vec3>& foregroundPoints = foregroundSeedList_.getData()[m];
        if (foregroundPoints.empty())
            continue;
        for (size_t i=0; i<foregroundPoints.size()-1; i++) {
            tgt::vec3 left = foregroundPoints[i];
            tgt::vec3 right = foregroundPoints[i+1];
            tgt::vec3 dir = tgt::normalize(right - left);
            for (float t=0.f; t<tgt::length(right-left); t += 1.f) {
                tgt::vec3 point = tgt::clamp(tgt::iround(left + t*dir), tgt::ivec3(0), volDim_-1);
                size_t index = volumeCoordsToIndex(point, volDim_);
                tgtAssert(index < numVoxels_, "Invalid index");
                if (!seedBuffer_[index]) {
                    seedBuffer_[index] = 2; //< foreground seed
                    tNumForegroundSeeds++;
                }
            }
        }
    }

    // foreground seeds from volume
    if (foregroundSeedVolume_) {
        tgt::ivec3 seedDim = foregroundSeedVolume_->getDimensions();
        tgt::vec3 seedVolScale = tgt::vec3(seedDim)/tgt::vec3(volDim_);
        tgt::ivec3 voxel;
        #ifdef VRN_MODULE_OPENMP
        #pragma omp parallel for reduction(+:tNumForegroundSeeds)
        #endif
        for (int z=0; z<volDim_.z; z++) {
            for (int y=0; y<volDim_.y; y++) {
                for (int x=0; x<volDim_.x; x++) {
                    size_t index = volumeCoordsToIndex(x, y, z, volDim_);
                    if (!seedBuffer_[index]) {
                        tgt::ivec3 fgVoxel = tgt::ifloor(tgt::vec3((float)x,(float)y,(float)z)*seedVolScale + 0.45f);
                        if (foregroundSeedVolume_->voxel(fgVoxel)) {
                            seedBuffer_[index] = 2; //< foreground seed
                            tNumForegroundSeeds++;
                        }
                    }
                }
            }
        }
    }

    // background geometry seeds
    for (int m=0; m<backgroundSeedList_.getNumSegments(); m++) {
        const std::vector<tgt::vec3>& backgroundPoints = backgroundSeedList_.getData()[m];
        if (backgroundPoints.empty())
            continue;
        for (size_t i=0; i<backgroundPoints.size()-1; i++) {
            tgt::vec3 left = backgroundPoints[i];
            tgt::vec3 right = backgroundPoints[i+1];
            tgt::vec3 dir = tgt::normalize(right - left);
            for (float t=0.f; t<tgt::length(right-left); t += 1.f) {
                tgt::vec3 point = tgt::clamp(tgt::iround(left + t*dir), tgt::ivec3(0), volDim_-1);
                size_t index = volumeCoordsToIndex(point, volDim_);
                tgtAssert(index < numVoxels_, "Invalid index");
                if (!seedBuffer_[index]) {
                    seedBuffer_[index] = 1; //< background seed
                    tNumBackgroundSeeds++;
                }
            }
        }
    }

    // background seeds from volume
    if (backgroundSeedVolume_) {
        tgt::ivec3 seedDim = backgroundSeedVolume_->getDimensions();
        tgt::vec3 seedVolScale = tgt::vec3(seedDim)/tgt::vec3(volDim_);
        tgt::ivec3 voxel;
        #ifdef VRN_MODULE_OPENMP
        #pragma omp parallel for reduction(+:tNumBackgroundSeeds)
        #endif
        for (int z=0; z<volDim_.z; z++) {
            for (int y=0; y<volDim_.y; y++) {
                for (int x=0; x<volDim_.x; x++) {
                    size_t index = volumeCoordsToIndex(x, y, z, volDim_);
                    if (!seedBuffer_[index]) {
                        tgt::ivec3 bgVoxel = tgt::ifloor(tgt::vec3((float)x,(float)y,(float)z)*seedVolScale + 0.45f);
                        if (backgroundSeedVolume_->voxel(bgVoxel)) {
                            seedBuffer_[index] = 1; //< background seed
                            tNumBackgroundSeeds++;
                        }
                    }
                }
            }
        }
    }

    numForegroundSeeds_ = tNumForegroundSeeds;
    numBackgroundSeeds_ = tNumBackgroundSeeds;

    numSeeds_ = numForegroundSeeds_ + numBackgroundSeeds_;

    if (numSeeds_ == 0) {
        LWARNINGC("voreen.RandomWalker.RandomWalkerTwoLabelSeeds", "no seeds");
    }
    else if (numForegroundSeeds_ == 0) {
        LWARNINGC("voreen.RandomWalker.RandomWalkerTwoLabelSeeds", "no foreground seeds");
        seedRange_ = tgt::vec2(0.f, 0.f);
    }
    else if (numBackgroundSeeds_ == 0) {
        LWARNINGC("voreen.RandomWalker.RandomWalkerTwoLabelSeeds", "no background seeds");
        seedRange_ = tgt::vec2(1.f, 1.f);
    }
    else {
        LDEBUGC("voreen.RandomWalker.RandomWalkerTwoLabelSeeds", "Foreground/Background seeds: " <<
            numForegroundSeeds_ << "/" << numBackgroundSeeds_);
        seedRange_ = tgt::vec2(0.f, 1.f);
    }
}

}   // namespace
