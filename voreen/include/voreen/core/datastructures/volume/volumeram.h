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

#ifndef VRN_VOLUMERAM_H
#define VRN_VOLUMERAM_H

// Note: please ensure that no OpenGL dependencies are added to this file

#include "voreen/core/datastructures/volume/volumerepresentation.h"
#include "voreen/core/datastructures/geometry/meshgeometry.h"
#include "voreen/core/datastructures/meta/realworldmappingmetadata.h"

#include <stdexcept>
namespace voreen {

/**
 * OpenGL-independent base class for volumetric data sets.
 *
 * This class stores the raw data as well as the
 * required meta information about it. It does,
 * however, neither perform any OpenGL operations nor
 * does it hold any OpenGL-related properties.
 *
 * @see VolumeGL
 */
class VRN_CORE_API VolumeRAM : public VolumeRepresentation {
public:
    enum Filter {
        NEAREST,
        LINEAR,
        CUBIC
    };

    VolumeRAM(const tgt::svec3& dimensions);
    VolumeRAM(const VolumeRAM* vol);

    virtual ~VolumeRAM() {}

    /// Use this as a kind of a virtual constructor.
    virtual VolumeRAM* clone() const = 0;

    /**
     * Use this as a kind of a virtual constructor that does _NOT_ copy over
     * the voxels but uses the given pointer instead as a voxel data. If this
     * pointer is 0, an empty volume without voxel data is created.
     */
    virtual VolumeRAM* clone(void* data) const = 0;
    virtual VolumeRAM* createNew(const tgt::svec3& dimensions, bool allocMem = false) const = 0;

    /// Create new volume which contains part of the data of the current volume.
    virtual VolumeRAM* getSubVolume(tgt::svec3 dimensions, tgt::svec3 offset) const = 0;

    /// Returns the number of bytes held in the \a data_ array.
    virtual size_t getNumBytes() const = 0;

    /// Returns the number of channels of this volume.
    virtual size_t getNumChannels() const = 0;

    /// Returns the number of bits that are allocated by each voxel.
    /// For technical reasons, it may exceed the volume's bit depth.
    virtual size_t getBitsAllocated() const = 0;

    /// Returns the number of bytes that are allocated for each voxel.
    virtual size_t getBytesPerVoxel() const = 0;

    /// Returns whether the volume's data type is a signed type.
    virtual bool isSigned() const = 0;

    /// Returns whether the volume's data type is an integer type.
    virtual bool isInteger() const = 0;

    virtual std::string getFormat() const;
    virtual std::string getBaseType() const;

    /**
     * Invalidates cached values (e.g. min/max), should be called when the volume was modified.
     */
    virtual void invalidate() const = 0;
    //------------------------------------------------------

    /**
     * Returns the minimum value contained by the specified channel converted to float.
     */
    virtual float minNormalizedValue(size_t channel = 0) const;

    /**
     * Returns the maximum value contained by the specified channel converted to float.
     */
    virtual float maxNormalizedValue(size_t channel = 0) const;

    /**
     * Returns the value with the smallest magnitude (absolute value) contained in the volume;
     * in case of multiple channels, this means that the minimum euclidian vector length is returned.
     * The normalized values are considered. For the real world magnitudes use VolumeMinMaxMagnitude.
     */
    virtual float minNormalizedMagnitude() const = 0;

    /**
     * Returns the value with the largest magnitude (absolute value) contained in the volume;
     * in case of multiple channels, this means that the maximum euclidian vector length is returned.
     * The normalized values are considered. For the real world magnitudes use VolumeMinMaxMagnitude.
     */
    virtual float maxNormalizedMagnitude() const = 0;

    //------------------------------------------------------

    /**
     * Returns the data set's minimal and maximal possible element values
     * according to its data type converted to float.
     */
    virtual tgt::vec2 elementRange() const = 0;

    /**
     * Returns the voxel of a given position and channel, normalized to float:
     *  - for \em unsigned integer types, the value range is mapped linearly to [0.0;1.0]
     *  - for \em signed integer types, the value range is mapped linearly to [-1.0;1.0]
     *  - floating point types are not mapped
     *
     * @param pos the position of the voxel
     * @param channel the channel of the voxel
     */
    virtual float getVoxelNormalized(const tgt::svec3& pos, size_t channel = 0) const = 0;

    /**
     * Returns the voxel of a given position and channel, normalized to float:
     *  - for \em unsigned integer types, the value range is mapped linearly to [0.0;1.0]
     *  - for \em signed integer types, the value range is mapped linearly to [-1.0;1.0]
     *  - floating point types are not mapped
     *
     * @param x the x-coordinate of the voxel
     * @param y the y-coordinate of the voxel
     * @param z the z-coordinate of the voxel
     * @param channel the channel of the voxel
     */
    virtual float getVoxelNormalized(size_t x, size_t y, size_t z, size_t channel = 0) const = 0;

    /**
     * Returns the voxel of a given position and channel, normalized to float:
     *  - for \em unsigned integer types, the value range is mapped linearly to [0.0;1.0]
     *  - for \em signed integer types, the value range is mapped linearly to [-1.0;1.0]
     *  - floating point types are not mapped
     *
     * @param index the index of the voxel
     * @param channel the channel of the voxel
     */
    virtual float getVoxelNormalized(size_t index, size_t channel = 0) const = 0;

    /**
     * Sets the voxel of a given postion and channel, de-normalized from a float:
     *  - for \em unsigned integer types, the value range [0.0;1.0] is mapped linearly to the integer range
     *  - for \em signed integer types, the value range [-1.0;1.0] is mapped linearly to the integer range
     *  - for floating point types, no mapping is performed
     *
     * @param value The float value to be set.
     * @param pos the position of the voxel
     * @param channel the channel of the voxel
     */
    virtual void setVoxelNormalized(float value, const tgt::svec3& pos, size_t channel = 0) = 0;

    /**
     * Sets the voxel of a given postion and channel, de-normalized from a float:
     *  - for \em unsigned integer types, the value range [0.0;1.0] is mapped linearly to the integer range
     *  - for \em signed integer types, the value range [-1.0;1.0] is mapped linearly to the integer range
     *  - for floating point types, no mapping is performed
     *
     * @param value The float value to be set.
     * @param x the x-coordinate of the voxel
     * @param y the y-coordinate of the voxel
     * @param z the z-coordinate of the voxel
     * @param channel the channel of the voxel
     */
    virtual void setVoxelNormalized(float value, size_t x, size_t y, size_t z, size_t channel = 0) = 0;

    /**
     * Sets the voxel of a given postion and channel, de-normalized from a float:
     *  - for \em unsigned integer types, the value range [0.0;1.0] is mapped linearly to the integer range
     *  - for \em signed integer types, the value range [-1.0;1.0] is mapped linearly to the integer range
     *  - for floating point types, no mapping is performed
     *
     * @param value The float value to be set.
     * @param index the index of the voxel
     * @param channel the channel of the voxel
     */
    virtual void setVoxelNormalized(float value, size_t index, size_t channel = 0) = 0;

    /**
     * Returns the trilinearly interpolated voxel of a given
     * postion and channel, normalized to a float.
     *
     * @see getVoxelFloat
     *
     * @param pos the position of the voxel
     * @param channel the channel of the voxel
     */
    virtual float getVoxelNormalizedLinear(const tgt::vec3& pos, size_t channel = 0) const;

    /**
     * Returns the tricubic interpolated voxel of a given position and channel,
     * converted to a normalized float. A Catmull-Rom spline is used for the interpolation.
     *
     * @see getVoxelFloat
     *
     * @param pos the position of the voxel
     * @param channel the channel of the voxel
     */
    virtual float getVoxelNormalizedCubic(const tgt::vec3& pos, size_t channel = 0) const;

    /**
     * Helper function that returns the voxel value at the passed voxel position as string (nearest filtering).
     *
     * The default implementation returns for scalar volumes the "stringified" scalar value.
     * For vector types the tgt::Vector stream operator is used, which yields "[x y z]" for Vector3 types, for example.
     * If the passed voxel position lies outside the volume dimensions, an empty string is returned.
     * If a real-world mapping is passed, the mapped normalized voxel value (float) is used instead of the native one.
     *
     * @param voxelPos The position in voxel coordinates whose value is to be returned.
     * @param mapping Optional real-world mapping to apply to the voxel value.
     *
     * @note Non-standard volume types should override this function.
     */
    virtual std::string getVoxelValueAsString(const tgt::svec3& voxelPos, const RealWorldMapping* mapping = 0, const size_t channel = -1) const;

    /// Set all volume data to zero
    virtual void clear() = 0;

    /// Gets a void* to the data stored with this Volume
    virtual const void* getData() const = 0;
    virtual void* getData() = 0;

    /// Gets a void* to the requested data. New memory will be allocated for the return buffer.
    virtual void* getBrickData(const tgt::svec3& offset, const tgt::svec3& dimension) const = 0;
    virtual void* getSliceData(const size_t firstSlice, const size_t lastSlice) const = 0;

    /**
     * Use this as type safe wrapper in order to get a proper typed pointer.
     */
    template<class T>
    inline static typename T::VoxelType* getData(T* v);
protected:
    // protected default constructor
    VolumeRAM() {}

    static const std::string loggerCat_;
};


/**
 * You can use this macro if you want to iterate over all voxels and it is
 * important for you that you do it dimensionwise, for instance: <br>
\code
VRN_FOR_EACH_VOXEL(i, svec3(0, 0, 0), svec(15, 20, 30))
    voxel(i) = i.x * i.y + i.z;
\endcode
 * or like this:
\code
VRN_FOR_EACH_VOXEL(i, svec3(0, 0, 0), svec(15, 20, 30)) {
    voxel(i) = i.x * i.y + i.z;
    foo();
}
\endcode
 * If you want to do sth just with each voxel use this simple for-loop
\code
for (size_t i = 0; i \< numVoxels_; ++i)
     voxel(i) = i;
\endcode
*/
#define VRN_FOR_EACH_VOXEL(INDEX, POS, SIZE) \
    for (auto INDEX = (POS); INDEX.z < (SIZE).z; ++INDEX.z)\
        for (INDEX.y = (POS).y; INDEX.y < (SIZE).y; ++INDEX.y)\
            for (INDEX.x = (POS).x; INDEX.x < (SIZE).x; ++INDEX.x)


/**
 * As an alternative to VRN_FOR_EACH_VOXEL, but don't need the voxel position
 * you can use VRN_FOR_EACH_VOXEL_INDEX, which is typically faster. This index
 * can be used for example with the voxel()-method of VolumeRAM.
 * In constrast tot VRN_FOR_EACH_VOXEL, you will have to provide the size of
 * the whole volume AND the body as a parameters to the macro.
*/
#define VRN_FOR_EACH_VOXEL_INDEX(INDEX, BEGIN, END, SIZE, BODY) \
    size_t vrn_fevi_z_slice_size = (SIZE).x*(SIZE).y;\
    size_t vrn_fevi_z_begin = (BEGIN).z*vrn_fevi_z_slice_size;\
    size_t vrn_fevi_z_end = (END).z*vrn_fevi_z_slice_size;\
    for (size_t vrn_fevi_z = vrn_fevi_z_begin; vrn_fevi_z < vrn_fevi_z_end ; vrn_fevi_z += vrn_fevi_z_slice_size) {\
        size_t vrn_fevi_y_begin = vrn_fevi_z + (BEGIN).y*(SIZE).x;\
        size_t vrn_fevi_y_end = vrn_fevi_z + (END).y*(SIZE).x;\
        for (size_t vrn_fevi_y = vrn_fevi_y_begin; vrn_fevi_y < vrn_fevi_y_end ; vrn_fevi_y += (SIZE).x){\
            size_t vrn_fevi_x_begin = vrn_fevi_y + (BEGIN).x;\
            size_t vrn_fevi_x_end = vrn_fevi_y + (END).x;\
            for (size_t INDEX = vrn_fevi_x_begin; INDEX < vrn_fevi_x_end; ++INDEX){\
                BODY\
            }\
        }\
    }\


} // namespace voreen

#endif // VRN_VOLUMERAM_H
