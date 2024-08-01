/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#ifndef VRN_OCTREEUTILS_H
#define VRN_OCTREEUTILS_H

#include "voreen/core/utils/stringutils.h"

#include "tgt/assert.h"
#include "tgt/logmanager.h"
#include "tgt/tgt_math.h"
#include "tgt/stopwatch.h"
#include "tgt/vector.h"
#include "tgt/matrix.h"

#include <sstream>

using tgt::svec3;
using tgt::vec3;

namespace { // helper functions

/// Converts a normalized sample position (within unit cube) to an integer voxel position.
inline svec3 samplePosToVoxel(const vec3& pos, const svec3& volumeDim) {
    return tgt::clamp(svec3(tgt::ifloor(pos*vec3(volumeDim))), svec3::zero, volumeDim-svec3(1));
}

/// Returns whether a is a multiple of b.
inline bool isMultipleOf(uint32_t a, uint32_t b) {
    return (a==0 || b==0 ? false : (a/b)*b == a);
}
inline bool isMultipleOf(uint64_t a, uint64_t b) {
    return (a==0 || b==0 ? false : (a/b)*b == a);
}
#ifdef __APPLE__
inline bool isMultipleOf(long unsigned int a, long unsigned int b) {
    return (a==0 || b==0 ? false : (a/b)*b == a);
}
#endif

/// Returns whether a is a component-wise multiple of b.
inline bool isMultipleOf(const svec3& a, const svec3& b) {
    return isMultipleOf(a.x, b.x) && isMultipleOf(a.y, b.y) && isMultipleOf(a.z, b.z);
}

template<class T>
inline T absDiff(T a, T b) {
    return (a > b ? a-b : b-a);
}

template<class T>
bool inRange(T val, T min, T max) {
    return (min <= val && val <= max);
}

template<class T>
bool inRange(const tgt::Vector2<T>& val, const tgt::Vector2<T>& min, const tgt::Vector2<T>& max) {
    return (tgt::hand(tgt::greaterThanEqual(val, min))) &&
           (tgt::hand(tgt::lessThanEqual(val, max)));
}

template<class T>
bool inRange(const tgt::Vector3<T>& val, const tgt::Vector3<T>& min, const tgt::Vector3<T>& max) {
    return (tgt::hand(tgt::greaterThanEqual(val, min))) &&
           (tgt::hand(tgt::lessThanEqual(val, max)));
}

inline bool isCubicAndPot(const svec3& dim) {
    return (dim.x == dim.y && dim.x == dim.z && tgt::isPowerOfTwo((int)dim.x));
}

inline size_t cubicCoordToLinear(const svec3& coord, const svec3& dim) {
    tgtAssert(tgt::hand(tgt::lessThan(coord, dim)), "coordinates outside dimensions");
    return coord.z*dim.x*dim.y + coord.y*dim.x + coord.x;
}

inline svec3 linearCoordToCubic(size_t coord, const svec3& dim) {
    tgtAssert(coord < dim.x*dim.y*dim.z, "linear coordinate outside dimensions");

    svec3 result;
    result.z = coord / (dim.x*dim.y);
    coord = coord % (dim.x*dim.y);
    result.y = coord / dim.x;
    result.x = coord % dim.x;

    return result;
}

inline uint16_t computeAvgValue(const uint16_t* buffer, const svec3& dim) {
    tgtAssert(buffer, "null pointer passed");

    size_t numElems = tgt::hmul(dim);
    uint64_t result = 0;
    for (size_t i=0; i<numElems; i++)
        result += buffer[i];
    result /= numElems;

    return static_cast<uint16_t>(result);
}


/**
 * Determine the maximum pixel size of a projected voxel on the viewport.
 */
inline float computeVoxelProjectionSize(const tgt::vec3& sample, const tgt::svec3& volumeDim, const tgt::mat4& textureToWorldMatrix,
    const tgt::vec3& camPos, const tgt::mat4& projectionMatrix, const tgt::ivec2& viewportSize)
{
    tgtAssert(inRange(sample, tgt::vec3(0.f), tgt::vec3(1.f)), "sample pos outside range [0.0;1.0]");

    // compute diagonal length of voxel in world space
    tgt::vec3 voxelDim = tgt::vec3(1.f) / static_cast<tgt::vec3>(volumeDim);
    tgt::vec4 voxelLLFWorld = textureToWorldMatrix*tgt::vec4(sample, 1.f);
    voxelLLFWorld /= voxelLLFWorld.w;
    tgt::vec4 voxelURBWorld = textureToWorldMatrix*tgt::vec4(sample+voxelDim, 1.f);
    voxelURBWorld /= voxelURBWorld.w;
    float nodeDiagonalWorld = tgt::length(voxelLLFWorld.xyz() - voxelURBWorld.xyz());

    // project world diagonal to view plane, assuming that diagonal is perpendicular to cam viewing direction
    float camDist = tgt::length(voxelLLFWorld.xyz() - camPos);
    tgt::vec4 projectedLLF = projectionMatrix * tgt::vec4(-nodeDiagonalWorld/2.f, 0.f, camDist, 1.f);
    projectedLLF /= projectedLLF.w;
    tgt::vec4 projectedURB = projectionMatrix * tgt::vec4(+nodeDiagonalWorld/2.f, 0.f, camDist, 1.f);
    projectedURB /= projectedURB.w;
    float nodeDiagonalNDC = tgt::length(projectedLLF.xyz() - projectedURB.xyz());

    // compute projected diagonal length in viewport coordinates
    float nodeDiagonalViewport = (nodeDiagonalNDC/2.f) * viewportSize.x;

    return nodeDiagonalViewport;
}

/**
 * Derives the required node level from the voxel projection size such that one voxel projects to approx. one pixel.
 */
inline size_t computeNodeLevel(const tgt::vec3& sample, const tgt::svec3& volumeDim,
    const tgt::mat4& textureToWorldMatrix, size_t numLevels,
    const tgt::vec3& camPos, const tgt::mat4& projectionMatrix, const tgt::ivec2& viewportSize)
{
    // viewport size (in pixel) of a voxel located at the sampling point
    float voxelProjectionSize = computeVoxelProjectionSize(sample, volumeDim, textureToWorldMatrix,
        camPos, projectionMatrix, viewportSize);

    // choose node level such that a brick voxel projects to approximately one pixel.
    int nodeLevel = ((int)numLevels-1) - tgt::iround((std::log(1.f/voxelProjectionSize) / std::log(2.f)));
    nodeLevel = tgt::clamp(nodeLevel, 0, (int)numLevels-1);

    return nodeLevel;
}

} // namespace anonymous

#endif
