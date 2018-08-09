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

#ifndef VRN_RANDOMWALKERSEEDS_H
#define VRN_RANDOMWALKERSEEDS_H

#include <string>

#include "randomwalkersolver.h"

#include "voreen/core/datastructures/geometry/pointsegmentlistgeometry.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"

namespace voreen {

/**
 * Abstract base class for Random Walker seed point definition,
 * used by RandomWalkerSolver.
 */
class RandomWalkerSeeds {

public:
    virtual ~RandomWalkerSeeds() {}

    virtual void initialize(const VolumeRAM* volume);

    virtual bool isSeedPoint(size_t index) const = 0;

    virtual bool isSeedPoint(const tgt::ivec3& voxel) const = 0;

    virtual float getSeedValue(size_t index) const = 0;

    virtual float getSeedValue(const tgt::ivec3& voxel) const = 0;

    size_t getNumSeeds() const;

    tgt::vec2 getSeedRange() const;

protected:
    const VolumeRAM* volume_;
    tgt::ivec3 volDim_;
    size_t numVoxels_;

    size_t numSeeds_;
    tgt::vec2 seedRange_;
};

//---------------------------------------------------------------------------------------

/**
 * Basic seed point definer that receives lists of fore- and background voxels
 * and assigns 0.0 to background seeds and 1.0 to foreground seeds.
 *
 * Additionally, seed volumes can be passed, in which all non-zero elements are considered seed voxels.
 *
 * Furthermore, a clipping range can be specified: All voxels outside the range are considered background seeds.
 */
class RandomWalkerTwoLabelSeeds : public RandomWalkerSeeds {

public:
    RandomWalkerTwoLabelSeeds(const PointSegmentListGeometryVec3& foregroundSeedList,
        const PointSegmentListGeometryVec3& backgroundSeedList,
        const VolumeRAM_UInt8* foregroundSeedVolume = 0, const VolumeRAM_UInt8* backgroundSeedVolume = 0,
        const tgt::ivec3& clipLLF = tgt::ivec3(-1), const tgt::ivec3& clipURB = tgt::ivec3(-1));

    ~RandomWalkerTwoLabelSeeds();

    void setClippingPlanes(const tgt::ivec3& clipLLF, const tgt::ivec3& clipURB);

    virtual void initialize(const VolumeRAM* volume);

    virtual bool isSeedPoint(size_t index) const;

    virtual bool isSeedPoint(const tgt::ivec3& voxel) const;

    virtual float getSeedValue(size_t index) const;

    virtual float getSeedValue(const tgt::ivec3& voxel) const;

    size_t getNumForegroundSeeds() const;

    size_t getNumBackgroundSeeds() const;

protected:
    const PointSegmentListGeometryVec3 foregroundSeedList_;
    const PointSegmentListGeometryVec3 backgroundSeedList_;
    const VolumeRAM_UInt8* foregroundSeedVolume_;
    const VolumeRAM_UInt8* backgroundSeedVolume_;
    tgt::ivec3 clipLLF_;
    tgt::ivec3 clipURB_;

    size_t numForegroundSeeds_;
    size_t numBackgroundSeeds_;

    char* seedBuffer_;
};

} //namespace

#endif
