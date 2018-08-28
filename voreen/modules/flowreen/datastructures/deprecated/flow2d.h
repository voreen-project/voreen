/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#ifndef VRN_FLOW2D_H
#define VRN_FLOW2D_H

#include "tgt/vector.h"
#include "tgt/matrix.h"

#include <deque>
#include <vector>

namespace voreen {

class Flow2D {
public:
    enum DataOrientation {XY, YX};

public:
    Flow2D(const tgt::vec2* const flow2D, const tgt::ivec2& dimensions, const DataOrientation& orientation,
        const float minValue = 0.0f, const float maxValue = 1.0f, const float maxMagnitude = 1.0f,
        bool freeOnDestruction = false);

    ~Flow2D();

    void free() { delete [] flow2D_; }

    tgt::ivec2 cellNumberToPos(const size_t n) const;

    tgt::ivec2 flowPosToSlicePos(const tgt::vec2& fp, const tgt::ivec2& sliceSize,
        const tgt::ivec2& offset = tgt::ivec2(0, 0)) const;

    std::vector<tgt::ivec2> flowPosToSlicePos(std::vector<tgt::vec2>& fps,
        const tgt::ivec2& sliceSize, const tgt::ivec2& offset = tgt::ivec2(0, 0)) const;

    std::deque<tgt::ivec2> flowPosToSlicePos(std::deque<tgt::vec2>& fps,
        const tgt::ivec2& sliceSize, const tgt::ivec2& offset = tgt::ivec2(0, 0)) const;

    static tgt::ivec2 getAxisPermutation(const DataOrientation& orientation);

    tgt::ivec2 getNearestCellPos(const tgt::vec2& r) const;

    /**
     * Returns the flow data after "normalization":
     * The data are mapped componentwise to the range of [0.0, 1.0] by calculating
     *
     * (v[i] - minValue_) / (maxValue_ - minValue_)
     *
     * where minValue_(maxValue_) are the smallest(greatest) values among all
     * components of all vectors.
     * This is necessary when using the flow data within OpenGL textures as
     * OpenGL requires normalized float values.
     *
     * NOTE: ensure the deletion of the returned pointer by using delete [] when
     * it is no longer needed.
     */
    tgt::vec3* getNormalizedFlow() const;

    tgt::mat3 getTransformationMatrix(const tgt::vec2& r, const float scaling = 1.0f) const;

    bool isInsideBoundings(const tgt::vec2& r) const;

    inline const tgt::vec2& lookupFlow(const tgt::vec2& r) const {
        return lookupFlowNearest(r);
    }

    //inline tgt::vec2 lookupFlow(const tgt::vec2& r) const {
    //    return lookupFlowBilinear(r);
    //}

    inline const tgt::vec2& lookupFlow(const tgt::ivec2& r) const {
        return flow2D_[posToCellNumber(r)];
    }

    size_t posToCellNumber(const tgt::ivec2& pos) const;

    tgt::ivec2 slicePosToFlowPos(const tgt::ivec2& vp, const tgt::ivec2& sliceSize,
        tgt::ivec2* const error = 0) const;

public: // for fast and easy access...
    const tgt::vec2* const flow2D_;
    const tgt::ivec2 dimensions_;
    const DataOrientation orientation_;
    const tgt::ivec2 axisPermutation_;
    const float minValue_;
    const float maxValue_;
    const float maxMagnitude_;

private:
    inline const tgt::vec2& lookupFlowNearest(const tgt::vec2& r) const {
        return flow2D_[posToCellNumber(getNearestCellPos(r))];
    }

    tgt::vec2 lookupFlowBilinear(const tgt::vec2& r) const;

private:
    const bool freeOnDestruction_;
};

}   // namespace

#endif
