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

#include "voreen/core/datastructures/volume/volumeram.h"

#include "voreen/core/datastructures/volume/volumefactory.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/utils/hashing.h"

#include "tgt/plane.h"
#include <sstream>

using tgt::vec3;
using tgt::ivec3;
using tgt::svec3;

namespace voreen {

const std::string VolumeRAM::loggerCat_("voreen.VolumeRAM");

VolumeRAM::VolumeRAM(const svec3& dimensions)
    : VolumeRepresentation(dimensions)
{}

VolumeRAM::VolumeRAM(const VolumeRAM* vol)
    : VolumeRepresentation(vol->getDimensions())
{}

std::string VolumeRAM::getFormat() const {
    VolumeFactory vf;
    return vf.getFormat(this);
}

std::string VolumeRAM::getBaseType() const {
    VolumeFactory vf;
    return vf.getBaseType(getFormat());
}

/*
 * getters and setters
 */

template<class T>
inline typename T::VoxelType* VolumeRAM::getData(T* v) {
    return (typename T::VoxelType*) v->getData();
}

/*
 * abstract access of the voxels
 */

float VolumeRAM::getVoxelNormalized(size_t x, size_t y, size_t z, size_t channel /*= 0*/) const {
    return getVoxelNormalized(svec3(x, y, z), channel);
}

float VolumeRAM::getVoxelNormalizedLinear(const vec3& pos, size_t channel /*= 0*/) const {
    vec3 posAbs = tgt::max(pos/* - 0.5f*/, vec3(0.0f));
    vec3 p = posAbs - floor(posAbs); // get decimal part
    svec3 llb = svec3(posAbs);
    svec3 urf = svec3(ceil(posAbs));
    urf = min(urf, dimensions_ - svec3(1)); // clamp so the lookups do not exceed the dimensions
    llb = min(llb, dimensions_ - svec3(1)); // dito

    /*
        interpolate linearly
    */
    return  getVoxelNormalized(llb.x, llb.y, llb.z, channel) * (1.f-p.x)*(1.f-p.y)*(1.f-p.z) // llB
          + getVoxelNormalized(urf.x, llb.y, llb.z, channel) * (    p.x)*(1.f-p.y)*(1.f-p.z) // lrB
          + getVoxelNormalized(urf.x, urf.y, llb.z, channel) * (    p.x)*(    p.y)*(1.f-p.z) // urB
          + getVoxelNormalized(llb.x, urf.y, llb.z, channel) * (1.f-p.x)*(    p.y)*(1.f-p.z) // ulB
          + getVoxelNormalized(llb.x, llb.y, urf.z, channel) * (1.f-p.x)*(1.f-p.y)*(    p.z) // llF
          + getVoxelNormalized(urf.x, llb.y, urf.z, channel) * (    p.x)*(1.f-p.y)*(    p.z) // lrF
          + getVoxelNormalized(urf.x, urf.y, urf.z, channel) * (    p.x)*(    p.y)*(    p.z) // urF
          + getVoxelNormalized(llb.x, urf.y, urf.z, channel) * (1.f-p.x)*(    p.y)*(    p.z);// ulF
}

float VolumeRAM::getVoxelNormalizedCubic(const tgt::vec3& pos, size_t channel /*= 0*/) const {
    vec3 posCorrected = tgt::max(pos/* - 0.5f*/, vec3(0.0f));
    vec3 p = posCorrected - floor(posCorrected); // get decimal part of sampling point

    //
    // 0. Compute Catmull-Rom spline base functions for sampling point p
    //
    const tgt::vec4 CINT_X(p.x*((2.f-p.x)*p.x-1.f),
                           p.x*p.x*(3.f*p.x-5.f) + 2.f,
                           p.x*((4.f-3.f*p.x)*p.x+1.f),
                           (p.x-1.f)*p.x*p.x);

    const tgt::vec4 CINT_Y(p.y*((2.f-p.y)*p.y-1.f),
                           p.y*p.y*(3.f*p.y-5.f) + 2.f,
                           p.y*((4.f-3.f*p.y)*p.y+1.f),
                           (p.y-1.f)*p.y*p.y);

    const tgt::vec4 CINT_Z(p.z*((2.f-p.z)*p.z-1.f),
                           p.z*p.z*(3.f*p.z-5.f) + 2.f,
                           p.z*((4.f-3.f*p.z)*p.z+1.f),
                           (p.z-1.f)*p.z*p.z);

    //
    // 1. Retrieve neighbored grid values (4x4x4 neighborhood).
    // Instead of using a three-dimensional array, we wrap the third dimension in vec4 vectors,
    // which can directly be used in the subsequent dot product computations
    //
    tgt::vec4 ENV[4][4];

    // boundaries of the 4x4x4 neighborhood used for interpolation
    ivec3 llb = ivec3(posCorrected) - 1;
    ivec3 urf = llb + 3;

    // only clamp indices, if llb or urf lie outside the grid (saves a lot)
    if (tgt::hor(tgt::lessThan(llb, tgt::ivec3(0))) || tgt::hor(tgt::greaterThanEqual(svec3(urf), dimensions_))) {
        const tgt::ivec3 maxDim = ivec3(dimensions_) - 1;
        tgt::ivec3 offset;
        for (offset.z=0; offset.z<4; ++offset.z) {
            for (offset.y=0; offset.y<4; ++offset.y) {
                for (offset.x=0; offset.x<4; ++offset.x) {
                    svec3 index = svec3(tgt::clamp(llb + offset, tgt::ivec3(0), maxDim));  //< clamping
                    ENV[offset.x][offset.y][offset.z] = getVoxelNormalized(index, channel);
                }
            }
        }
    }
    else { // no clamping necessary
        tgt::ivec3 offset;
        for (offset.z=0; offset.z<4; ++offset.z) {
            for (offset.y=0; offset.y<4; ++offset.y) {
                for (offset.x=0; offset.x<4; ++offset.x) {
                    tgt::svec3 index = svec3(llb + offset);
                    ENV[offset.x][offset.y][offset.z] = getVoxelNormalized(index, channel);
                }
            }
        }
    }

    //
    // 2. Interpolate in z direction, again wrapping the last dimension in vec4s
    //
    tgt::vec4 t[4];
    for (int x=0; x<4; x++) {
        for (int y=0; y<4; y++) {
            t[x][y] = tgt::dot(CINT_Z, ENV[x][y]);
        }
    }

    //
    // 3. Interpolate in y direction, using the z-interpolated values
    //
    tgt::vec4 u;
    for (int i=0; i<4; i++) {
        u[i] = tgt::dot(CINT_Y, t[i]);
    }

    //
    // 4. Interpolate in x direction to obtain final result
    //
    float result = tgt::clamp(tgt::dot(CINT_X, u) / 8.f, 0.f, 1.f);

    return result;
}

std::string VolumeRAM::getVoxelValueAsString(const tgt::svec3& voxelPos, const RealWorldMapping* mapping, const size_t channel) const {

    const size_t sMinusOne = -1; //fixes compiler warnings

    // voxel pos outside volume dimensions
    if (tgt::hor(tgt::greaterThanEqual(voxelPos, getDimensions())))
        return "";

    //return empty string if channel is larger than num channels
    if(channel != sMinusOne /*-1*/ && channel >= getNumChannels())
        return "";

    std::stringstream str;

    if (mapping) {
        // apply real-world mapping to the normalized voxel value
        if(channel == sMinusOne /* -1 */) {
            size_t numChannels = getNumChannels();
            if (numChannels > 1)
                str << "[";
            for (size_t i=0; i<numChannels; i++) {
                str << mapping->normalizedToRealWorld(getVoxelNormalized(voxelPos, i));
                if (i < numChannels-1)
                    str << " ";
            }
            if (numChannels > 1)
                str << "]";
            str << " " << mapping->getUnit();
        } else {
            str << mapping->normalizedToRealWorld(getVoxelNormalized(voxelPos, channel));
        }
    }
    else {
        // no mapping => return native voxel value

        // scalar
        if (dynamic_cast<const VolumeAtomic<uint8_t>*>(this)) {
            str << static_cast<int>(static_cast<const VolumeAtomic<uint8_t>*>(this)->voxel(voxelPos));
        }
        else if (dynamic_cast<const VolumeAtomic<int8_t>*>(this)) {
            str << static_cast<int>(static_cast<const VolumeAtomic<int8_t>*>(this)->voxel(voxelPos));
        }
        else if (dynamic_cast<const VolumeAtomic<uint16_t>*>(this)) {
            str << static_cast<const VolumeAtomic<uint16_t>*>(this)->voxel(voxelPos);
        }
        else if (dynamic_cast<const VolumeAtomic<int16_t>*>(this)) {
            str << static_cast<const VolumeAtomic<int16_t>*>(this)->voxel(voxelPos);
        }
        else if (dynamic_cast<const VolumeAtomic<uint32_t>*>(this)) {
            str << static_cast<const VolumeAtomic<uint32_t>*>(this)->voxel(voxelPos);
        }
        else if (dynamic_cast<const VolumeAtomic<int32_t>*>(this)) {
            str << static_cast<const VolumeAtomic<int32_t>*>(this)->voxel(voxelPos);
        }
        else if (dynamic_cast<const VolumeAtomic<uint64_t>*>(this)) {
            str << static_cast<const VolumeAtomic<uint64_t>*>(this)->voxel(voxelPos);
        }
        else if (dynamic_cast<const VolumeAtomic<int64_t>*>(this)) {
            str << static_cast<const VolumeAtomic<int64_t>*>(this)->voxel(voxelPos);
        }
        else if (dynamic_cast<const VolumeAtomic<float>*>(this)) {
            str << static_cast<const VolumeAtomic<float>*>(this)->voxel(voxelPos);
        }
        else if (dynamic_cast<const VolumeAtomic<double>*>(this)) {
            str << static_cast<const VolumeAtomic<double>*>(this)->voxel(voxelPos);
        }
        // vec2
        else if (dynamic_cast<const VolumeAtomic<tgt::Vector2<uint8_t> >*>(this)) {
            str << static_cast<const VolumeAtomic<tgt::Vector2<uint8_t> >*>(this)->voxel(voxelPos);
        }
        else if (dynamic_cast<const VolumeAtomic<tgt::Vector2<int8_t> >*>(this)) {
            str << static_cast<const VolumeAtomic<tgt::Vector2<int8_t> >*>(this)->voxel(voxelPos);
        }
        else if (dynamic_cast<const VolumeAtomic<tgt::Vector2<uint16_t> >*>(this)) {
            str << static_cast<const VolumeAtomic<tgt::Vector2<uint16_t> >*>(this)->voxel(voxelPos);
        }
        else if (dynamic_cast<const VolumeAtomic<tgt::Vector2<int16_t> >*>(this)) {
            str << static_cast<const VolumeAtomic<tgt::Vector2<int16_t> >*>(this)->voxel(voxelPos);
        }
        else if (dynamic_cast<const VolumeAtomic<tgt::Vector2<uint32_t> >*>(this)) {
            str << static_cast<const VolumeAtomic<tgt::Vector2<uint32_t> >*>(this)->voxel(voxelPos);
        }
        else if (dynamic_cast<const VolumeAtomic<tgt::Vector2<int32_t> >*>(this)) {
            str << static_cast<const VolumeAtomic<tgt::Vector2<int32_t> >*>(this)->voxel(voxelPos);
        }
        else if (dynamic_cast<const VolumeAtomic<tgt::Vector2<uint64_t> >*>(this)) {
            str << static_cast<const VolumeAtomic<tgt::Vector2<uint64_t> >*>(this)->voxel(voxelPos);
        }
        else if (dynamic_cast<const VolumeAtomic<tgt::Vector2<int64_t> >*>(this)) {
            str << static_cast<const VolumeAtomic<tgt::Vector2<int64_t> >*>(this)->voxel(voxelPos);
        }
        else if (dynamic_cast<const VolumeAtomic<tgt::Vector2<float> >*>(this)) {
            str << static_cast<const VolumeAtomic<tgt::Vector2<float> >*>(this)->voxel(voxelPos);
        }
        else if (dynamic_cast<const VolumeAtomic<tgt::Vector2<double> >*>(this)) {
            str << static_cast<const VolumeAtomic<tgt::Vector2<double> >*>(this)->voxel(voxelPos);
        }
        // vec3
        else if (dynamic_cast<const VolumeAtomic<tgt::Vector3<uint8_t> >*>(this)) {
            str << static_cast<const VolumeAtomic<tgt::Vector3<uint8_t> >*>(this)->voxel(voxelPos);
        }
        else if (dynamic_cast<const VolumeAtomic<tgt::Vector3<int8_t> >*>(this)) {
            str << static_cast<const VolumeAtomic<tgt::Vector3<int8_t> >*>(this)->voxel(voxelPos);
        }
        else if (dynamic_cast<const VolumeAtomic<tgt::Vector3<uint16_t> >*>(this)) {
            str << static_cast<const VolumeAtomic<tgt::Vector3<uint16_t> >*>(this)->voxel(voxelPos);
        }
        else if (dynamic_cast<const VolumeAtomic<tgt::Vector3<int16_t> >*>(this)) {
            str << static_cast<const VolumeAtomic<tgt::Vector3<int16_t> >*>(this)->voxel(voxelPos);
        }
        else if (dynamic_cast<const VolumeAtomic<tgt::Vector3<uint32_t> >*>(this)) {
            str << static_cast<const VolumeAtomic<tgt::Vector3<uint32_t> >*>(this)->voxel(voxelPos);
        }
        else if (dynamic_cast<const VolumeAtomic<tgt::Vector3<int32_t> >*>(this)) {
            str << static_cast<const VolumeAtomic<tgt::Vector3<int32_t> >*>(this)->voxel(voxelPos);
        }
        else if (dynamic_cast<const VolumeAtomic<tgt::Vector3<uint64_t> >*>(this)) {
            str << static_cast<const VolumeAtomic<tgt::Vector3<uint64_t> >*>(this)->voxel(voxelPos);
        }
        else if (dynamic_cast<const VolumeAtomic<tgt::Vector3<int64_t> >*>(this)) {
            str << static_cast<const VolumeAtomic<tgt::Vector3<int64_t> >*>(this)->voxel(voxelPos);
        }
        else if (dynamic_cast<const VolumeAtomic<tgt::Vector3<float> >*>(this)) {
            str << static_cast<const VolumeAtomic<tgt::Vector3<float> >*>(this)->voxel(voxelPos);
        }
        else if (dynamic_cast<const VolumeAtomic<tgt::Vector3<double> >*>(this)) {
            str << static_cast<const VolumeAtomic<tgt::Vector3<double> >*>(this)->voxel(voxelPos);
        }
        // vec4
        else if (dynamic_cast<const VolumeAtomic<tgt::Vector4<uint8_t> >*>(this)) {
            str << static_cast<const VolumeAtomic<tgt::Vector4<uint8_t> >*>(this)->voxel(voxelPos);
        }
        else if (dynamic_cast<const VolumeAtomic<tgt::Vector4<int8_t> >*>(this)) {
            str << static_cast<const VolumeAtomic<tgt::Vector4<uint8_t> >*>(this)->voxel(voxelPos);
        }
        else if (dynamic_cast<const VolumeAtomic<tgt::Vector4<uint16_t> >*>(this)) {
            str << static_cast<const VolumeAtomic<tgt::Vector4<uint8_t> >*>(this)->voxel(voxelPos);
        }
        else if (dynamic_cast<const VolumeAtomic<tgt::Vector4<int16_t> >*>(this)) {
            str << static_cast<const VolumeAtomic<tgt::Vector4<uint8_t> >*>(this)->voxel(voxelPos);
        }
        else if (dynamic_cast<const VolumeAtomic<tgt::Vector4<uint32_t> >*>(this)) {
            str << static_cast<const VolumeAtomic<tgt::Vector4<uint8_t> >*>(this)->voxel(voxelPos);
        }
        else if (dynamic_cast<const VolumeAtomic<tgt::Vector4<int32_t> >*>(this)) {
            str << static_cast<const VolumeAtomic<tgt::Vector4<uint8_t> >*>(this)->voxel(voxelPos);
        }
        else if (dynamic_cast<const VolumeAtomic<tgt::Vector4<uint64_t> >*>(this)) {
            str << static_cast<const VolumeAtomic<tgt::Vector4<uint8_t> >*>(this)->voxel(voxelPos);
        }
        else if (dynamic_cast<const VolumeAtomic<tgt::Vector4<int64_t> >*>(this)) {
            str << static_cast<const VolumeAtomic<tgt::Vector4<uint8_t> >*>(this)->voxel(voxelPos);
        }
        else if (dynamic_cast<const VolumeAtomic<tgt::Vector4<float> >*>(this)) {
            str << static_cast<const VolumeAtomic<tgt::Vector4<uint8_t> >*>(this)->voxel(voxelPos);
        }
        else if (dynamic_cast<const VolumeAtomic<tgt::Vector4<double> >*>(this)) {
            str << static_cast<const VolumeAtomic<tgt::Vector4<uint8_t> >*>(this)->voxel(voxelPos);
        }
        else {
            str << getVoxelNormalized(voxelPos);
        }
    }

    return str.str();
}

float VolumeRAM::minNormalizedValue(size_t channel) const {
    float minValue = getVoxelNormalized(0, channel);
    for (size_t i=0; i<numVoxels_; i++) {
        minValue = std::min(minValue, getVoxelNormalized(i, channel));
    }
    return minValue;
}

float VolumeRAM::maxNormalizedValue(size_t channel) const {
    float maxValue = getVoxelNormalized(0, channel);
    for (size_t i=0; i<numVoxels_; i++) {
        maxValue = std::max(maxValue, getVoxelNormalized(i, channel));
    }
    return maxValue;
}

} // namespace voreen
