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

#include "flow3d.h"
#include "modules/flowreen/utils/flowmath.h"
#include <limits>

namespace voreen {

Flow3D::Flow3D(const tgt::vec3* const flow3D, const tgt::ivec3& dimensions, const VoxelOrder& voxelOrder,
               const float minValue, const float maxValue, const float maxMagnitude)
    : flow3D_(flow3D),
    dimensions_(dimensions),
    voxelOrder_(voxelOrder),
    axisPermutation_(getAxisPermutation(voxelOrder)),
    minValue_(minValue),
    maxValue_(maxValue),
    maxMagnitude_(maxMagnitude)
{
}

tgt::ivec3 Flow3D::clampToFlowDimensions(const tgt::ivec3& r) const {
    tgt::ivec3 ir(r);
    for (size_t i = 0; i < 3; ++i) {
        if (ir[i] < 0)
            ir[i] = 0;
        else if (ir[i] >= dimensions_[i])
            ir[i] = (dimensions_[i] - 1);
    }
    return ir;
}

Flow2D Flow3D::extractSlice(const tgt::ivec3& planeComponents, const size_t sliceNo) const {
    const int& i = planeComponents.x;
    const int& j = planeComponents.y;
    const int& k = planeComponents.z;

    const tgt::ivec2 sliceDim(dimensions_[i], dimensions_[j]);
    tgt::vec2* slice = new tgt::vec2[sliceDim.x * sliceDim.y];
    memset(slice, 0, sliceDim.x * sliceDim.y * sizeof(tgt::vec2));

    tgt::ivec3 voxelPos(0, 0, 0);
    voxelPos[k] = static_cast<int>(sliceNo);

    for (int y = 0; y < sliceDim.y; ++y) {
        int a = y * sliceDim.x;
        voxelPos[j] = y;
        for (int x = 0; x < sliceDim.x; ++x) {
            voxelPos[i] = x;
            size_t voxelNo = posToVoxelNumber(voxelPos, axisPermutation_, dimensions_);
            const tgt::vec3& v = flow3D_[voxelNo];

            size_t sliceElem = a + x;
            slice[sliceElem].x = v[i];
            slice[sliceElem].y = v[j];
        }
    }

    return Flow2D(slice, sliceDim, Flow2D::XY, minValue_, maxValue_, maxMagnitude_, true);
}

tgt::ivec2 Flow3D::flowPosToSlicePos(const tgt::vec3& fp, const tgt::ivec2& sliceSize,
                                        const tgt::ivec3& components, const tgt::ivec2& offset) const
{
    const int& i = components.x;
    const int& j = components.y;

    tgt::vec3 aux(fp / static_cast<tgt::vec3>(dimensions_));
    tgt::ivec2 vp(0, 0);
    vp.x = static_cast<int>(tgt::round(aux[i] * sliceSize.x) + offset.x) % sliceSize.x;
    vp.y = static_cast<int>(tgt::round(aux[j] * sliceSize.y) + offset.y) % sliceSize.y;
    return vp;
}

std::vector<tgt::ivec2> Flow3D::flowPosToSlicePos(std::vector<tgt::vec3>& fps,
        const tgt::ivec2& sliceSize, const tgt::ivec3& components, const tgt::ivec2& offset) const
{
    const int& i = components.x;
    const int& j = components.y;

    std::vector<tgt::ivec2> sps;    // slice positions
    while (fps.empty() == false) {
        std::vector<tgt::vec3>::iterator it = fps.begin();

        // no call to overloaded flowPosToViewportPos() for performance
        //
        tgt::vec3 aux(*it / static_cast<tgt::vec3>(dimensions_));
        tgt::ivec2 sp(
            (static_cast<int>(tgt::round(aux[i] * sliceSize.x)) + offset.x) % sliceSize.x,
            (static_cast<int>(tgt::round(aux[j] * sliceSize.y)) + offset.y) % sliceSize.y);

        sps.push_back(sp);
        fps.erase(it);
    }
    return sps;
}

std::deque<tgt::ivec2> Flow3D::flowPosToSlicePos(std::deque<tgt::vec3>& fps,
        const tgt::ivec2& sliceSize, const tgt::ivec3& components, const tgt::ivec2& offset) const
{
    const int& i = components.x;
    const int& j = components.y;

    std::deque<tgt::ivec2> sps;    // slice positions
    while (fps.empty() == false) {
        std::deque<tgt::vec3>::iterator it = fps.begin();

        // no call to overloaded flowPosToViewportPos() for performance
        //
        tgt::vec3 aux(*it / static_cast<tgt::vec3>(dimensions_));
        tgt::ivec2 sp(
            (static_cast<int>(tgt::round(aux[i] * sliceSize.x)) + offset.x) % sliceSize.x,
            (static_cast<int>(tgt::round(aux[j] * sliceSize.y)) + offset.y) % sliceSize.y);

        sps.push_back(sp);
        fps.erase(it);
    }
    return sps;
}

tgt::ivec2 Flow3D::getFlowSliceDimensions(const tgt::ivec3& planeComponents) const {
    return tgt::ivec2(dimensions_[planeComponents.x], dimensions_[planeComponents.y]);
}


tgt::mat4 Flow3D::getTransformationMatrix(const tgt::vec3& r, const float scaling) const {
    const tgt::vec3& v = lookupFlow(r);
    if (v == tgt::vec3::zero)
        return tgt::mat4::identity;
    tgt::vec3 front = FlowMath::normalize(v);

    tgt::vec3 temp(0.0f, 1.0f, 0.0f);
    if (fabsf(tgt::dot(temp, front)) >= 0.99f)
        temp = tgt::vec3(0.0f, 0.0f, 1.0f);
    tgt::vec3 up = FlowMath::normalize(tgt::cross(temp, front));
    tgt::vec3 right = FlowMath::normalize(tgt::cross(front, up));

    if (scaling != 1.0f) {
        front *= scaling;
        up *= scaling;
        right *= scaling;
    }

    return tgt::mat4(right.x, up.x, front.x, r.x,
        right.y, up.y, front.y, r.y,
        right.z, up.z, front.z, r.z,
        0.0f, 0.0f, 0.0f, 1.0f);
}

tgt::ivec3 Flow3D::getNearestVoxelPos(const tgt::vec3& r) const {
    tgt::ivec3 ir(0, 0, 0);
    for (size_t i = 0; i < 3; ++i) {
        //ir[i] = (r[i] <= 0.0f) ? 0 : static_cast<int>(r[i]);// cut off fractional part (floor)
        ir[i] = (r[i] <= 0.0f) ? 0 : static_cast<int>((r[i] + 0.5f));
        if (ir[i] >= dimensions_[i])
            ir[i] = (dimensions_[i] - 1);
    }
    return ir;
}

bool Flow3D::isInsideBoundings(const tgt::vec3& r) const {
    for (size_t i = 0; i < 3; ++i) {
        if (r[i] < 0.0f)
            return false;
        else if (static_cast<int>(ceilf(r[i])) >= dimensions_[i])
            return false;
    }
    return true;
}

bool Flow3D::isInsideBoundings(const tgt::ivec3& r) const {
    for (size_t i = 0; i < 3; ++i) {
        if (r[i] < 0)
            return false;
        else if (r[i] >= dimensions_[i])
            return false;
    }
    return true;
}

bool Flow3D::isNull(const tgt::vec3& r) const {
    return (lookupFlow(r) == tgt::vec3::zero);
}

tgt::ivec3 Flow3D::slicePosToFlowPos(const tgt::ivec2& vp, const tgt::ivec2& sliceSize,
                                        const tgt::ivec3& components, const size_t sliceNo,
                                        tgt::ivec2* const error) const
{
    tgt::vec2 c(static_cast<float>(vp.x % sliceSize.x),
        static_cast<float>(vp.y % sliceSize.y));

    const int& i = components.x;
    const int& j = components.y;
    const int& k = components.z;

    tgt::vec3 p(0, 0, 0);
    p[i] = (c.x / sliceSize.x);
    p[j] = (c.y / sliceSize.y);
    p[k] = float(sliceNo) / dimensions_[k];
    tgt::ivec3 r = getNearestVoxelPos(p * static_cast<tgt::vec3>(dimensions_));

    if (error != 0)
        *error = (vp - flowPosToSlicePos(r, sliceSize, components));

    return r;
}

tgt::vec3 Flow3D::toFlowPosition(const tgt::ivec3& r, const tgt::ivec3& textureSize,
                                 tgt::ivec3* const error) const
{
    const tgt::vec3 dim(static_cast<tgt::vec3>(dimensions_));
    const tgt::vec3 textureScaling(static_cast<tgt::vec3>(textureSize) / dim);

    tgt::vec3 rFlow(r.x / textureScaling.x, r.y / textureScaling.y, r.z / textureScaling.z);
    rFlow = tgt::clamp(rFlow, tgt::vec3::zero, dim);

    if (error != 0)
        *error = (r - toTexturePosition(rFlow, textureSize));
    return rFlow;
}

tgt::ivec3 Flow3D::toTexturePosition(const tgt::vec3& r, const tgt::ivec3& textureSize,
        const tgt::ivec3& offset) const
{
    const tgt::vec3 dim(static_cast<tgt::vec3>(dimensions_));
    const tgt::vec3 textureScaling(static_cast<tgt::vec3>(textureSize) / dim);
    tgt::ivec3 temp(
        (static_cast<int>(tgt::round(r.x * textureScaling.x)) + offset.x),
        (static_cast<int>(tgt::round(r.y * textureScaling.y)) + offset.y),
        (static_cast<int>(tgt::round(r.z * textureScaling.z)) + offset.z));
    temp = tgt::clamp(temp, tgt::ivec3::zero, textureSize);
    return temp;
}

std::deque<tgt::ivec3> Flow3D::toTexturePosition(std::deque<tgt::vec3>& fps, const tgt::ivec3& textureSize,
                                         const tgt::ivec3& offset) const
{
    const tgt::vec3 dim(static_cast<tgt::vec3>(dimensions_));
    const tgt::vec3 textureScaling(static_cast<tgt::vec3>(textureSize) / dim);

    std::deque<tgt::ivec3> tps;
    while (fps.empty() == false) {
        std::deque<tgt::vec3>::iterator it = fps.begin();
        const tgt::vec3& r = *it;
        tgt::ivec3 temp(
            (static_cast<int>(tgt::round(r.x * textureScaling.x)) + offset.x),
            (static_cast<int>(tgt::round(r.y * textureScaling.y)) + offset.y),
            (static_cast<int>(tgt::round(r.z * textureScaling.z)) + offset.z));
        temp = tgt::clamp(temp, tgt::ivec3::zero, textureSize);
        tps.push_back(temp);

        fps.erase(it);
    }
    return tps;
}

// public static methods
//

tgt::ivec3 Flow3D::getAxisPermutation(const VoxelOrder& voxelOrder) {
    switch (voxelOrder) {
        case XZY:
            return tgt::ivec3(0, 2, 1);
        case YXZ:
            return tgt::ivec3(1, 0, 2);
        case YZX:
            return tgt::ivec3(2, 0, 1);
        case ZXY:
            return tgt::ivec3(1, 2, 0);
        case ZYX:
            return tgt::ivec3(2, 1, 0);
        case XYZ:
        default:
            break;
    }
    return tgt::ivec3(0, 1, 2);
}

Flow3D::VoxelOrder Flow3D::getVoxelOrder(const unsigned char& voxelOrder) {
    switch (voxelOrder) {
        case XZY:
            return XZY;
        case YXZ:
            return YXZ;
        case YZX:
            return YZX;
        case ZXY:
            return ZXY;
        case ZYX:
            return ZYX;
        case XYZ:
        default:
            break;
    }
    return XYZ;
}

tgt::vec3 Flow3D::permute(const tgt::vec3& input, const tgt::ivec3& permutation) {
    return tgt::vec3(input[permutation[0]], input[permutation[1]], input[permutation[2]]);
}

size_t Flow3D::posToVoxelNumber(const tgt::ivec3& pos, const tgt::ivec3& permutation,
                                const tgt::ivec3& dimensions)
{
    const int& i = permutation[0];
    const int& j = permutation[1];
    const int& k = permutation[2];
    return (pos[i] + pos[j] * dimensions[i] + pos[k] * dimensions[i] * dimensions[j]);
}

tgt::ivec3 Flow3D::voxelNumberToPos(const size_t n, const tgt::ivec3& permutation,
                                    const tgt::ivec3& dimensions)
{
    tgt::ivec3 pos(0, 0, 0);
    const int& i = permutation[0];
    const int& j = permutation[1];
    const int& k = permutation[2];

    pos[i] = (n % dimensions[i]);
    pos[j] = ((n / dimensions[i]) % dimensions[j]);
    pos[k] = (static_cast<int>(n) / (dimensions[i] * dimensions[j]));
    return pos;
}

// private methods
//

const tgt::vec3& Flow3D::lookupFlowNearest(const tgt::vec3& r) const {
    size_t n = posToVoxelNumber(getNearestVoxelPos(r), axisPermutation_, dimensions_);
    return flow3D_[n];
}

tgt::vec3 Flow3D::lookupFlowTrilinear(const tgt::vec3& r) const {
    tgt::ivec3 ir1(static_cast<int>(floorf(r.x)), static_cast<int>(floorf(r.y)),
        static_cast<int>(floorf(r.z)));
    ir1 = clampToFlowDimensions(ir1);
    tgt::ivec3 ir2 = clampToFlowDimensions(ir1 + tgt::ivec3(1));

    // interpolate in x-direction
    //
    tgt::vec3 v1 = linearInterpolation(r.x, ir1, 0);
    tgt::vec3 v2 = linearInterpolation(r.x, tgt::ivec3(0, ir2.y, ir1.z), 0);
    tgt::vec3 v3 = linearInterpolation(r.x, tgt::ivec3(0, ir1.y, ir2.z), 0);
    tgt::vec3 v4 = linearInterpolation(r.x, ir2, 0);

    // interpolates value from x-direction interpolation in y-direction
    //
    float fintegral = 0.0f;
    float b = modff(r.y, &fintegral);
    float a = 1.0f - b;
    tgt::vec3 v5 = ((v1 * a) + (v2 * b));
    tgt::vec3 v6 = ((v3 * a) + (v4 * b));

    // interpolate values from y-direction interpolation in z-drection
    //
    b = modff(r.z, &fintegral);
    a = 1.0f - b;
    return ((v5 * a) + (v6 * b));
}

tgt::vec3 Flow3D::linearInterpolation(const float value, const tgt::ivec3& r, const size_t direction) const
{
    const size_t component = direction;

    float fintegral = 0.0f;
    float b = modff(value, &fintegral);
    float a = 1.0f - b;

    int integral = static_cast<int>(fintegral);
    if (integral < 0)
        integral = 0;
    else if (integral >= dimensions_[component])
        integral = (dimensions_[component] - 1);

    tgt::ivec3 ir(r);
    ir[component] = integral;
    size_t voxel1 = posToVoxelNumber(ir, axisPermutation_, dimensions_);
    ++(ir[component]);
    ir = clampToFlowDimensions(ir);
    size_t voxel2 = posToVoxelNumber(ir, axisPermutation_, dimensions_);

    return ((flow3D_[voxel1] * a) + (flow3D_[voxel2] * b));
}

}   // namespace
