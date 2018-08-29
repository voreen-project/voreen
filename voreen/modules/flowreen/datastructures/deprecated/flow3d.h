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

#ifndef VRN_FLOW3D_H
#define VRN_FLOW3D_H

#include "flow2d.h"

namespace voreen {

class Flow3D {
public:
    enum VoxelOrder { XZY, XYZ, YXZ, YZX, ZXY, ZYX };

public:
    Flow3D(const tgt::vec3* const flow3D, const tgt::ivec3& dimensions, const VoxelOrder& voxelOrder,
        const float minValue = 0.0f, const float maxValue = 1.0f, const float maxMagnitude = 1.0f);

    ~Flow3D() {}

    void free() { delete [] flow3D_; }

    inline const tgt::vec3& lookupFlow(const tgt::vec3& r) const {
        return lookupFlowNearest(r);
    }

    //inline tgt::vec3 lookupFlow(const tgt::vec3& r) const {
    //    return lookupFlowTrilinear(r);
    //}

    inline const tgt::vec3& lookupFlow(const tgt::ivec3& r) const {
        return flow3D_[posToVoxelNumber(r, axisPermutation_, dimensions_)];
    }

    /**
     * Clamps the passed point to the flow's dimensions.
     */
    tgt::ivec3 clampToFlowDimensions(const tgt::ivec3& r) const;

    Flow2D extractSlice(const tgt::ivec3& planeComponents, const size_t sliceNo) const;

    /**
     * Transforms the given position in coordinates within the flow volume to
     * viewport coordinates (projection) according to the plane alignment described
     * by the given permutation of components.
     */
    tgt::ivec2 flowPosToSlicePos(const tgt::vec3& fp, const tgt::ivec2& sliceSize,
        const tgt::ivec3& components, const tgt::ivec2& offset = tgt::ivec2(0, 0)) const;

     /**
     * Transforms the given positions in coordinates within the flow volume to
     * viewport coordinates (projection) according to the plane alignment described
     * by the given permutation of components.
     */
    std::vector<tgt::ivec2> flowPosToSlicePos(std::vector<tgt::vec3>& fps,
        const tgt::ivec2& sliceSize, const tgt::ivec3& components,
        const tgt::ivec2& offset = tgt::ivec2(0, 0)) const;

    /**
     * Transforms the given positions in coordinates within the flow volume to
     * viewport coordinates (projection) according to the plane alignment described
     * by the given permutation of components.
     */
    std::deque<tgt::ivec2> flowPosToSlicePos(std::deque<tgt::vec3>& fps,
        const tgt::ivec2& sliceSize, const tgt::ivec3& components,
        const tgt::ivec2& offset = tgt::ivec2(0, 0)) const;

    /**
     * Returns the dimensions of the slice through the flow's volume, based
     * on the passed permutation of comonents (e.g. [2, 1, 0] for the sagittal
     * plane)
     */
    tgt::ivec2 getFlowSliceDimensions(const tgt::ivec3& planeComponents) const;

    tgt::mat4 getTransformationMatrix(const tgt::vec3& r, const float scaling = 1.0f) const;

    /**
     * Returns the coordinates of the nearest voxel within this volume for the
     * given position.
     */
    tgt::ivec3 getNearestVoxelPos(const tgt::vec3& r) const;

    /**
     * Determines whether the given point lies within the flow's volume boundaries.
     */
    bool isInsideBoundings(const tgt::vec3& r) const;
    bool isInsideBoundings(const tgt::ivec3& r) const;

    bool isNull(const tgt::vec3& r) const;

    /**
     * Transforms the given position from coordinates on the given slice to coordinates
     * within the flow volume according to the given slice number and the plane alignemnt
     * described by the given permutation of components.
     */
    tgt::ivec3 slicePosToFlowPos(const tgt::ivec2& vp, const tgt::ivec2& sliceSize,
        const tgt::ivec3& components, const size_t sliceNo, tgt::ivec2* const error = 0) const;

    tgt::vec3 toFlowPosition(const tgt::ivec3& r, const tgt::ivec3& textureSize,
        tgt::ivec3* const error) const;
    tgt::ivec3 toTexturePosition(const tgt::vec3& r, const tgt::ivec3& textureSize,
        const tgt::ivec3& offset = tgt::ivec3(0, 0, 0)) const;
    std::deque<tgt::ivec3> toTexturePosition(std::deque<tgt::vec3>& fps,
        const tgt::ivec3& textureSize, const tgt::ivec3& offset = tgt::ivec3(0, 0, 0)) const;

    // ------------------------------------------------------------------------

    /**
     * Returns the permutation of the axis which will affect the order
     * of the linearization of the data stored within the volume for the
     * given voxelOrder (e.g. XZY results in [0, 2, 1]).
     */
    static tgt::ivec3 getAxisPermutation(const VoxelOrder& voxelOrder);

    static tgt::ivec3 getAxisPermutation(const unsigned char& voxelOrder) {
        return getAxisPermutation(getVoxelOrder(voxelOrder));
    }

    /**
     * Returns the equivalent VoxelOrder type / enum for the given
     * unsigned char.
     */
    static VoxelOrder getVoxelOrder(const unsigned char& voxelOrder);

    /**
     * Permutes the components of the given input vec3 according to the given
     * permutation.
     */
    static tgt::vec3 permute(const tgt::vec3& input, const tgt::ivec3& permutation);

    /**
     */
    static size_t posToVoxelNumber(const tgt::ivec3& pos, const tgt::ivec3& permutation,
        const tgt::ivec3& dimensions);

    /**
     */
    static tgt::ivec3 voxelNumberToPos(const size_t n, const tgt::ivec3& permutation,
        const tgt::ivec3& dimensions);

public: // for fast and easy access...
    const tgt::vec3* const flow3D_;
    const tgt::ivec3 dimensions_;
    const VoxelOrder voxelOrder_;
    const tgt::ivec3 axisPermutation_;
    const float minValue_;
    const float maxValue_;
    const float maxMagnitude_;

private:
    const tgt::vec3& lookupFlowNearest(const tgt::vec3& r) const;
    tgt::vec3 lookupFlowTrilinear(const tgt::vec3& r) const;
    tgt::vec3 linearInterpolation(const float value, const tgt::ivec3& r, const size_t direction) const;
};

}   // namespace

#endif
