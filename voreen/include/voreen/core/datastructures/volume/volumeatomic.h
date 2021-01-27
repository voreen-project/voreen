/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#ifndef VRN_VOLUMEATOMIC_H
#define VRN_VOLUMEATOMIC_H

#include <algorithm>
#include <cstring>
#include <fstream>

#include "tgt/assert.h"
#include "tgt/logmanager.h"

#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volumeelement.h"

#include "voreen/core/datastructures/tensor.h"

#include <typeinfo>

namespace voreen {

template<class T>
class VolumeAtomic : public VolumeRAM {

public:
    typedef T VoxelType;

    enum {
        BYTES_PER_VOXEL = sizeof(T),
        BITS_PER_VOXEL  = BYTES_PER_VOXEL * 8
    };

    /**
     * While using this constructor the class will automatically allocate
     * an appropriate chunk of memory. This memory will be deleted by this class.
     * If allocMem is false, no memory will be allocated.
     */
    VolumeAtomic(const tgt::svec3& dimensions, bool allocMem=true);

    /**
     * While using this constructor the class will use an preallocated chunk
     * of memory given in \p data. If takeOwnership is true, this memory will be deleted by this class.
     */
    VolumeAtomic(T* data, const tgt::svec3& dimensions, bool takeOwnership=true);

    /**
     * Move constructor.
     */
    VolumeAtomic(VolumeAtomic<T>&& other);

    /**
     * Move assignment.
     */
    VolumeAtomic<T>& operator=(VolumeAtomic<T>&& other);

    /// Deletes the \a data_ array
    virtual ~VolumeAtomic();

    virtual VolumeAtomic<T>* clone() const;
    virtual VolumeAtomic<T>* clone(void* data) const;
    virtual VolumeAtomic<T>* createNew(const tgt::svec3& dimensions, bool allocMem = false) const;
    virtual VolumeAtomic<T>* getSubVolume(tgt::svec3 dimensions, tgt::svec3 offset = tgt::svec3(0,0,0)) const;

    VolumeAtomic<T> copy() const;

    virtual size_t getBitsAllocated() const;

    virtual size_t getNumChannels() const;

    virtual size_t getBytesPerVoxel() const;

    virtual bool isSigned() const;

    virtual bool isInteger() const;

    /**
     * Returns the minimal and maximal element data values that can be contained
     * by the volume according to its data type.
     *
     * @note for unbounded types, (0.0, 1.0) is returned
     */
    virtual tgt::vec2 elementRange() const;

    /**
     * Returns the minimum data value that is actually contained by the volume.
     *
     * @note Calculated min/max values are cached. Call invalidate() in order to
     *      enforce a recalculation.
     */
    T min() const;

    /**
     * Returns the minimum value contained by the specified channel converted to float.
     * This function chooses either a fast scalar implementation (if possible) or the
     * slower parent method otherwise.
     */
    virtual float minNormalizedValue(size_t channel = 0) const;

    /**
     * Returns the maximal data value that is actually contained by the volume.
     *
     * @note Calculated min/max values are cached. Call invalidate() in order to
     *      enforce a recalculation.
     */
    T max() const;

    /**
     * Returns the maximum value contained by the specified channel converted to float.
     * This function chooses either a fast scalar implementation (if possible) or the
     * slower parent method otherwise.
     */
    virtual float maxNormalizedValue(size_t channel = 0) const;

    /**
     * Returns the value with minimum magnitude contained in the volume, converted to float.
     * This function chooses between a scalar implementation or a euclidian vector length method,
     * in case the volume has multiple channels.
     */
    virtual float minNormalizedMagnitude() const;

    /**
     * Returns the value with maximum magnitude contained in the volume, converted to float.
     * This function chooses between a scalar implementation or a euclidian vector length method,
     * in case the volume has multiple channels.
     */
    virtual float maxNormalizedMagnitude() const;

    virtual size_t getNumBytes() const;

    /*
     * Helpers for calculating the position in 3d
     */
    inline static size_t calcPos(const tgt::svec3& dimensions, size_t x, size_t y, size_t z);

    inline static size_t calcPos(const tgt::svec3& dimensions, const tgt::svec3& pos);

    /*
      Methods for accessing the voxels
    */

    /// just get a proper pointer for read and write access
    inline T* voxel();

    /// just get a proper pointer for read access only
    inline const T* voxel() const;

    /// get or set voxel
    inline T& voxel(size_t x, size_t y, size_t z);

    /// get voxel
    inline const T& voxel(size_t x, size_t y, size_t z) const;

    /// get or set voxel
    inline T& voxel(const tgt::svec3& pos);

    /// get voxel
    inline const T& voxel(const tgt::svec3& pos) const;

    /// get or set voxel
    inline T& voxel(size_t i);

    /// get voxel
    inline const T& voxel(size_t i) const;


    /*
     * getVoxelNormalized and setVoxelNormalized
     */
    virtual float getVoxelNormalized(const tgt::svec3& pos, size_t channel = 0) const;
    virtual float getVoxelNormalized(size_t x, size_t y, size_t z, size_t channel = 0) const;
    virtual float getVoxelNormalized(size_t index, size_t channel = 0) const;
    virtual void setVoxelNormalized(float value, const tgt::svec3& pos, size_t channel = 0);
    virtual void setVoxelNormalized(float value, size_t x, size_t y, size_t z, size_t channel = 0);
    virtual void setVoxelNormalized(float value, size_t index, size_t channel = 0);

    /**
     * TODO: implement this function for nearest/cubic in all volume representations
     * TODO: use channel right
     * @param allowOutsidePosition true => clamp to edge, false => tgtAssert
     */
    const T getVoxelLinear(const tgt::vec3& pos, size_t channel = 0, bool allowOutsidePosition = false) const;
    /**
     * TODO: implement this function for nearest/cubic in all volume representations
     * TODO: use channel right
     * @param allowOutsidePosition true => clamp to edge, false => tgtAssert
     */
    const T getVoxelNearest(const tgt::vec3& pos, size_t channel = 0, bool allowOutsidePosition = false) const;


    virtual void clear();
    void fill(T value);
    virtual const void* getData() const;
    virtual void* getData();

    virtual void* getBrickData(const tgt::svec3& offset, const tgt::svec3& dimensions) const;
    virtual void* getSliceData(const size_t firstSlice, const size_t lastSlice) const;

    /**
     * Invalidates cached values (e.g. min/max), should be called when the volume was modified.
     */
    virtual void invalidate() const;

    //-------------------------------------------------------------------
protected:
    // small utility
    template<bool> struct IsScalar{};

    // Helper method for scalar minimum float value
    float minNormalizedImpl(size_t channel, IsScalar<true>) const;

    // Helper method for non-scalar minimum float value
    float minNormalizedImpl(size_t channel, IsScalar<false>) const;

    // Helper method for scalar maximum float value
    float maxNormalizedImpl(size_t channel, IsScalar<true>) const;

    // Helper method for non-scalar maximum float value
    float maxNormalizedImpl(size_t channel, IsScalar<false>) const;

    // Helper method for scalar maximum magnitude
    float maxNormalizedMagnitudeImpl(IsScalar<true>) const;

    // Helper method for non-scalar maximum magnitude (euclidian vector length)
    float maxNormalizedMagnitudeImpl(IsScalar<false>) const;

    // Helper method for scalar minimum magnitude
    float minNormalizedMagnitudeImpl(IsScalar<true>) const;

    // Helper method for non-scalar minimum magnitude (euclidian vector length)
    float minNormalizedMagnitudeImpl(IsScalar<false>) const;

    T* data_;
    bool ownsData_;

    tgt::vec2 elementRange_;

    mutable T minValue_;
    mutable T maxValue_;
    mutable bool minMaxValid_;
private:
    /*
     *  Automatically generated copy constructor copies data_,
     * we therefore disable it. Use clone()/copy() instead.
     *
     * The implementation is intentionally omitted to provoke a linker error
     * when accidentally using this constructor.
     */
    VolumeAtomic(const VolumeAtomic<T>&);
};

/*
 * typedefs for easy usage
 */

typedef VolumeAtomic<uint8_t>   VolumeRAM_UInt8;
typedef VolumeAtomic<uint16_t>  VolumeRAM_UInt16;
typedef VolumeAtomic<uint32_t>  VolumeRAM_UInt32;
typedef VolumeAtomic<uint64_t>  VolumeRAM_UInt64;

typedef VolumeAtomic<int8_t>    VolumeRAM_Int8;
typedef VolumeAtomic<int16_t>   VolumeRAM_Int16;
typedef VolumeAtomic<int32_t>   VolumeRAM_Int32;
typedef VolumeAtomic<int64_t>   VolumeRAM_Int64;

typedef VolumeAtomic<float>     VolumeRAM_Float;
typedef VolumeAtomic<double>    VolumeRAM_Double;

typedef VolumeAtomic<tgt::Vector2<uint8_t> >    VolumeRAM_2xUInt8;
typedef VolumeAtomic<tgt::Vector2< int8_t> >    VolumeRAM_2xInt8;
typedef VolumeAtomic<tgt::Vector2<uint16_t> >   VolumeRAM_2xUInt16;
typedef VolumeAtomic<tgt::Vector2< int16_t> >   VolumeRAM_2xInt16;
typedef VolumeAtomic<tgt::Vector2<uint32_t> >   VolumeRAM_2xUInt32;
typedef VolumeAtomic<tgt::Vector2< int32_t> >   VolumeRAM_2xInt32;
typedef VolumeAtomic<tgt::Vector2<uint64_t> >   VolumeRAM_2xUInt64;
typedef VolumeAtomic<tgt::Vector2< int64_t> >   VolumeRAM_2xInt64;
typedef VolumeAtomic<tgt::vec2>                 VolumeRAM_2xFloat;
typedef VolumeAtomic<tgt::dvec2>                VolumeRAM_2xDouble;

typedef VolumeAtomic<tgt::Vector3<uint8_t> >    VolumeRAM_3xUInt8;
typedef VolumeAtomic<tgt::Vector3<int8_t> >     VolumeRAM_3xInt8;
typedef VolumeAtomic<tgt::Vector3<uint16_t> >   VolumeRAM_3xUInt16;
typedef VolumeAtomic<tgt::Vector3<int16_t> >    VolumeRAM_3xInt16;
typedef VolumeAtomic<tgt::Vector3<uint32_t> >   VolumeRAM_3xUInt32;
typedef VolumeAtomic<tgt::Vector3<int32_t> >    VolumeRAM_3xInt32;
typedef VolumeAtomic<tgt::Vector3<uint64_t> >   VolumeRAM_3xUInt64;
typedef VolumeAtomic<tgt::Vector3<int64_t> >    VolumeRAM_3xInt64;
typedef VolumeAtomic<tgt::vec3>                 VolumeRAM_3xFloat;
typedef VolumeAtomic<tgt::dvec3>                VolumeRAM_3xDouble;

typedef VolumeAtomic<tgt::Vector4<uint8_t> >    VolumeRAM_4xUInt8;
typedef VolumeAtomic<tgt::Vector4<int8_t> >     VolumeRAM_4xInt8;
typedef VolumeAtomic<tgt::Vector4<uint16_t> >   VolumeRAM_4xUInt16;
typedef VolumeAtomic<tgt::Vector4<int16_t> >    VolumeRAM_4xInt16;
typedef VolumeAtomic<tgt::Vector4<uint32_t> >   VolumeRAM_4xUInt32;
typedef VolumeAtomic<tgt::Vector4<int32_t> >    VolumeRAM_4xInt32;
typedef VolumeAtomic<tgt::Vector4<uint64_t> >   VolumeRAM_4xUInt64;
typedef VolumeAtomic<tgt::Vector4<int64_t> >    VolumeRAM_4xInt64;
typedef VolumeAtomic<tgt::vec4>                 VolumeRAM_4xFloat;
typedef VolumeAtomic<tgt::dvec4>                VolumeRAM_4xDouble;

typedef VolumeAtomic<tgt::mat3>  VolumeRAM_Mat3Float;

typedef VolumeAtomic<Tensor2<float> >  VolumeRAM_Tensor2Float;

//------------------------------------------------------------------------------
//  implementation
//------------------------------------------------------------------------------

/*
 * constructors and destructor
 */

template<class T>
VolumeAtomic<T>::VolumeAtomic(const tgt::svec3& dimensions, bool allocMem)
    : VolumeRAM(dimensions)
    , data_(0)
    , ownsData_(true)
    , elementRange_(static_cast<float>(VolumeElement<T>::rangeMinElement()),
        static_cast<float>(VolumeElement<T>::rangeMaxElement()))
    , minMaxValid_(false)
{
    if (allocMem) {
        try {
            data_ = new T[numVoxels_];
        }
        catch (std::bad_alloc&) {
            LERROR("Failed to allocate memory: bad allocation");
            throw; // throw it to the caller
        }
    }
}

template<class T>
VolumeAtomic<T>::VolumeAtomic(T* data,
                              const tgt::svec3& dimensions,
                              bool takeOwnership)
    : VolumeRAM(dimensions)
    , data_(data)
    , ownsData_(takeOwnership)
    , elementRange_(static_cast<float>(VolumeElement<T>::rangeMinElement()),
         static_cast<float>(VolumeElement<T>::rangeMaxElement()))
    , minMaxValid_(false)
{
}

template<class T>
VolumeAtomic<T>::VolumeAtomic(VolumeAtomic<T>&& other)
    : VolumeRAM(other.dimensions_)
    , data_(other.data_)
    , ownsData_(other.ownsData_)
    , elementRange_(other.elementRange_)
    , minMaxValid_(other.minMaxValid_)
{
    other.data_ = nullptr;
}

template<typename T>
VolumeAtomic<T>& VolumeAtomic<T>::operator=(VolumeAtomic<T>&& other) {
    // This looks quite funky, but is actually safe:

    // Destruct the current object, but keep the memory.
    this->~VolumeAtomic();

    // Call the move constructor on the memory region of the current object.
    new(this) VolumeAtomic(std::move(other));

    return *this;
}

template<class T>
VolumeAtomic<T>* VolumeAtomic<T>::clone() const {
    // create clone
    VolumeAtomic<T>* newVolume = 0;
    try {
        newVolume = new VolumeAtomic<T>(dimensions_); // allocate a chunk of data
    }
    catch (const std::bad_alloc&) {
        LERROR("Failed to clone volume: bad allocation");
        throw; // throw it to the caller
    }

    // copy over the voxel data
    T* dst = newVolume->data_;
    const T* src = data_;
    size_t len = getNumVoxels();
    std::copy(src, src+len, dst);

    return newVolume;
}

template<class T>
VolumeAtomic<T>* VolumeAtomic<T>::clone(void* data) const {
    // create clone
    VolumeAtomic<T>* newVolume = 0;
    if (data) {
        // use preallocated data
        newVolume = new VolumeAtomic<T>((T*) data, dimensions_);
    }
    else {
        // create volume without allocating memory
        newVolume = new VolumeAtomic<T>(dimensions_, false);
    }

    return newVolume;
}

template<class T>
VolumeAtomic<T> VolumeAtomic<T>::copy() const {
    // create clone
    VolumeAtomic<T> newVolume(dimensions_);

    // copy over the voxel data
    T* dst = newVolume.data_;
    const T* src = data_;
    size_t len = getNumVoxels();
    std::copy(src, src+len, dst);

    return newVolume;
}

template<class T>
VolumeAtomic<T>* VolumeAtomic<T>::createNew(const tgt::svec3& dimensions, bool allocMem) const {
    // create new volume
    VolumeAtomic<T>* newVolume = new VolumeAtomic<T>(dimensions, allocMem);

    return newVolume;
}

template<class T>
VolumeAtomic<T>* VolumeAtomic<T>::getSubVolume(tgt::svec3 dimensions, tgt::svec3 offset) const {
    //check wrong parameters
    if(tgt::hmul(dimensions) == 0)
        throw std::invalid_argument("requested subvolume dimensions are zero!");
    if(!tgt::hand(tgt::lessThanEqual(offset+dimensions,getDimensions())))
        throw std::invalid_argument("requested subvolume outside volume date!");

    // create new volume
    VolumeAtomic<T>* newVolume = new VolumeAtomic<T>(dimensions, true);
    T* data = reinterpret_cast<T*>(newVolume->getData());

    // determine parameters
    tgt::svec3 dataDims = getDimensions();
    size_t initialStartPos = (offset.z * dataDims.x * dataDims.y)+(offset.y * dataDims.x) + offset.x;


    // copy each row for every slice to form sub volume
    size_t volumePos;
    size_t subVolumePos;
    for (size_t i=0; i < dimensions.z; i++) {
        for (size_t j=0; j < dimensions.y; j++) {
            volumePos = (j*dataDims.x) + (i*dataDims.x*dataDims.y);
            subVolumePos = (j*dimensions.x) + (i*dimensions.x*dimensions.y);
            T* dst = data + subVolumePos;
            const T* src = data_ + volumePos + initialStartPos;
            size_t len = dimensions.x;
            std::copy(src, src+len, dst);
        }
    }

    return newVolume;
}

template<class T>
VolumeAtomic<T>::~VolumeAtomic() {
    if(ownsData_) {
        delete[] data_;
    }
}

template<class T>
size_t VolumeAtomic<T>::getNumChannels() const {
    return VolumeElement<T>::getNumChannels();
}

template<class T>
size_t VolumeAtomic<T>::getBitsAllocated() const {
    return BITS_PER_VOXEL;
}

template<class T>
size_t VolumeAtomic<T>::getBytesPerVoxel() const {
    return BYTES_PER_VOXEL;
}

template<class T>
size_t VolumeAtomic<T>::getNumBytes() const {
    return sizeof(T) * numVoxels_;
}

template<class T>
bool VolumeAtomic<T>::isSigned() const {
    return VolumeElement<T>::isSigned();
}

template<class T>
bool VolumeAtomic<T>::isInteger() const {
    return VolumeElement<T>::isInteger();
}

/*
 * getVoxelNormalized and setVoxelNormalized
 */

template<class T>
float VolumeAtomic<T>::getVoxelNormalized(const tgt::svec3& pos, size_t channel) const {
    return getVoxelNormalized(calcPos(dimensions_, pos), channel);
}

template<class T>
float VolumeAtomic<T>::getVoxelNormalized(size_t x, size_t y, size_t z, size_t channel) const {
    return getVoxelNormalized(tgt::svec3(x, y, z), channel);
}

template<class T>
float voreen::VolumeAtomic<T>::getVoxelNormalized(size_t index, size_t channel) const {
    typedef typename VolumeElement<T>::BaseType Base;

    Base value = VolumeElement<T>::getChannel(voxel(index), channel);
    return getTypeAsFloat(value);
}

template<class T>
void VolumeAtomic<T>::setVoxelNormalized(float value, const tgt::svec3& pos, size_t channel) {
    setVoxelNormalized(value, calcPos(dimensions_, pos), channel);
}

template<class T>
void VolumeAtomic<T>::setVoxelNormalized(float value, size_t x, size_t y, size_t z, size_t channel) {
    setVoxelNormalized(value, tgt::svec3(x, y, z), channel);
}

template<class T>
void voreen::VolumeAtomic<T>::setVoxelNormalized(float value, size_t index, size_t channel) {
    typedef typename VolumeElement<T>::BaseType Base;

    VolumeElement<T>::setChannel(getFloatAsType<Base>(value), voxel(index), channel);
}

template<class T>
const T VolumeAtomic<T>::getVoxelLinear(const tgt::vec3& pos, size_t channel,  bool allowOutsidePosition) const {
    tgtAssert(channel == 0, "not supported yet");

    //clamp to dimension
    tgt::vec3 clampedPos = tgt::clamp(pos,tgt::vec3::zero,tgt::vec3(getDimensions()-tgt::svec3::one));
    tgtAssert(allowOutsidePosition || clampedPos == pos, "position outside volume");

    //get values
    tgt::vec3 decimal = clampedPos - tgt::floor(clampedPos); // get decimal part
    tgt::svec3 llf = tgt::svec3(tgt::floor(clampedPos));
    tgt::svec3 urb = tgt::svec3(tgt::ceil(clampedPos));


    //interpolation
    return  voxel(llf.x, llf.y, llf.z) * (1.f-decimal.x)*(1.f-decimal.y)*(1.f-decimal.z)
          + voxel(urb.x, llf.y, llf.z) * (    decimal.x)*(1.f-decimal.y)*(1.f-decimal.z)
          + voxel(llf.x, urb.y, llf.z) * (1.f-decimal.x)*(    decimal.y)*(1.f-decimal.z)
          + voxel(urb.x, urb.y, llf.z) * (    decimal.x)*(    decimal.y)*(1.f-decimal.z)
          + voxel(llf.x, llf.y, urb.z) * (1.f-decimal.x)*(1.f-decimal.y)*(    decimal.z)
          + voxel(urb.x, llf.y, urb.z) * (    decimal.x)*(1.f-decimal.y)*(    decimal.z)
          + voxel(llf.x, urb.y, urb.z) * (1.f-decimal.x)*(    decimal.y)*(    decimal.z)
          + voxel(urb.x, urb.y, urb.z) * (    decimal.x)*(    decimal.y)*(    decimal.z);
}

template<class T>
const T VolumeAtomic<T>::getVoxelNearest(const tgt::vec3& pos, size_t channel,  bool allowOutsidePosition) const {
    tgtAssert(channel == 0, "not supported yet");

    //clamp to dimension
    tgt::vec3 clampedPos = tgt::clamp(pos,tgt::vec3::zero,tgt::vec3(getDimensions()-tgt::svec3::one));
    tgtAssert(allowOutsidePosition || clampedPos == pos, "position outside volume");

    //get values
    tgt::ivec3 intPos = tgt::iround(clampedPos);

    return voxel(tgt::svec3(intPos));
}

/*
 * Helpers for calculating the position in 3d
 */

template<class T>
inline size_t VolumeAtomic<T>::calcPos(const tgt::svec3& dimensions, size_t x, size_t y, size_t z) {
    //we need to work with size_t to avoid problems when working with very large datasets
    return z*dimensions.x*dimensions.y + y*dimensions.x + x;
}

template<class T>
inline size_t VolumeAtomic<T>::calcPos(const tgt::svec3& dimensions, const tgt::svec3& pos) {
    //we need to work with size_t to avoid problems when working with very large datasets
    return pos.z*dimensions.x*dimensions.y + pos.y*dimensions.x + pos.x;
}

/*
 * further methods
 */

template<class T>
tgt::vec2 VolumeAtomic<T>::elementRange() const {
    return elementRange_;
}

template<class T>
T VolumeAtomic<T>::min() const {
    if (!minMaxValid_) {
        auto minMax = std::minmax_element(data_, data_ + getNumVoxels());
        minValue_ = *minMax.first;
        maxValue_ = *minMax.second;
        minMaxValid_ = true;
    }

    return minValue_;

}

template<class T>
T VolumeAtomic<T>::max() const {
    if (!minMaxValid_) {
        auto minMax = std::minmax_element(data_, data_ + getNumVoxels());
        minValue_ = *minMax.first;
        maxValue_ = *minMax.second;
        minMaxValid_ = true;
    }

    return maxValue_;
}

// scalar (fast)
template<typename T>
float VolumeAtomic<T>::minNormalizedImpl(size_t /*channel*/, IsScalar<true>) const {
    return getTypeAsFloat(min());
}

// non-scalar (we have to use the slow way via getVoxelNormalized())
template<typename T>
float VolumeAtomic<T>::minNormalizedImpl(size_t channel, IsScalar<false>) const {
    return VolumeRAM::minNormalizedValue(channel);
}

template<class T>
float VolumeAtomic<T>::minNormalizedValue(size_t channel) const {
    tgtAssert(channel < (size_t)getNumChannels(), "channel number too large");
    return minNormalizedImpl(channel, IsScalar<std::numeric_limits<T>::is_specialized>());
}

// scalar (fast)
template<typename T>
float VolumeAtomic<T>::maxNormalizedImpl(size_t /*channel*/, IsScalar<true>) const {
    return getTypeAsFloat(max());
}

// non-scalar (we have to use the slow way via getVoxelNormalized())
template<typename T>
float VolumeAtomic<T>::maxNormalizedImpl(size_t channel, IsScalar<false>) const {
    return VolumeRAM::maxNormalizedValue(channel);
}

template<class T>
float VolumeAtomic<T>::maxNormalizedValue(size_t channel) const {
    return maxNormalizedImpl(channel, IsScalar<std::numeric_limits<T>::is_specialized>());
}

// scalar
template<typename T>
float VolumeAtomic<T>::minNormalizedMagnitudeImpl(IsScalar<true>) const {
    float minMagnitude = std::numeric_limits<float>::max();
    for (size_t i=0; i<numVoxels_; i++) {
        float val = getVoxelNormalized(i);
        minMagnitude = std::min(minMagnitude, std::abs(val));
    }
    return minMagnitude;
}

// scalar
template<typename T>
float VolumeAtomic<T>::maxNormalizedMagnitudeImpl(IsScalar<true>) const {
    float maxMagnitude = 0.f;
    for (size_t i=0; i<numVoxels_; i++) {
        float val = getVoxelNormalized(i);
        maxMagnitude = std::max(maxMagnitude, std::abs(val));
    }
    return maxMagnitude;
}

// non-scalar
template<typename T>
float VolumeAtomic<T>::minNormalizedMagnitudeImpl(IsScalar<false>) const {
    float minMagnitude = std::numeric_limits<float>::max();
    for (size_t i=0; i<numVoxels_; i++) {
        float magnitude = 0.0f;
        for(size_t c=0; c<getNumChannels(); c++) {
            float value = getVoxelNormalized(i, c);
            magnitude += value * value;
        }
        minMagnitude = std::min(minMagnitude, magnitude);
    }
    return std::sqrt(minMagnitude);
}

// non-scalar
template<typename T>
float VolumeAtomic<T>::maxNormalizedMagnitudeImpl(IsScalar<false>) const {
    float maxMagnitude = 0.f;
    for (size_t i=0; i<numVoxels_; i++) {
        float magnitude = 0.0f;
        for(size_t c=0; c<getNumChannels(); c++) {
            float value = getVoxelNormalized(i, c);
            magnitude += value * value;
        }
        maxMagnitude = std::max(maxMagnitude, magnitude);
    }
    return std::sqrt(maxMagnitude);
}

template<class T>
float VolumeAtomic<T>::minNormalizedMagnitude() const {
    return minNormalizedMagnitudeImpl(IsScalar<std::numeric_limits<T>::is_specialized>());
}

template<class T>
float VolumeAtomic<T>::maxNormalizedMagnitude() const {
    return maxNormalizedMagnitudeImpl(IsScalar<std::numeric_limits<T>::is_specialized>());
}

template<class T>
void VolumeAtomic<T>::clear() {
    fill(T(0.0));
}

template<class T>
void VolumeAtomic<T>::fill(T value) {
    std::fill(data_, data_+getNumVoxels(), value);
    invalidate();
}

template<class T>
void VolumeAtomic<T>::invalidate() const {
    minMaxValid_ = false;
}

template<class T>
const void* VolumeAtomic<T>::getData() const {
    return reinterpret_cast<const void*>(data_);
}

template<class T>
void* VolumeAtomic<T>::getData() {
    return reinterpret_cast<void*>(data_);
}

template<class T>
void* VolumeAtomic<T>::getBrickData(const tgt::svec3& offset, const tgt::svec3& dimensions) const {
    //check wrong parameters
    if(tgt::hmul(dimensions) == 0)
        throw std::invalid_argument("requested subvolume dimensions are zero!");
    if(!tgt::hand(tgt::lessThanEqual(offset+dimensions,getDimensions())))
        throw std::invalid_argument("requested subvolume outside volume date!");

    // create buffer
    T* data = new T[tgt::hmul(dimensions)];

    // determine parameters
    tgt::svec3 dataDims = getDimensions();
    size_t initialStartPos = (offset.z * dataDims.x * dataDims.y)+(offset.y * dataDims.x) + offset.x;

    // copy each row for every slice to form sub volume
    size_t volumePos;
    size_t bufferPos;
    for (size_t z=0; z < dimensions.z; z++) {
        for (size_t y=0; y < dimensions.y; y++) {
            volumePos = (y*dataDims.x) + (z*dataDims.x*dataDims.y);
            bufferPos = (y*dimensions.x) + (z*dimensions.x*dimensions.y);
            T* dst = data + bufferPos;
            const T* src = data_ + volumePos + initialStartPos;
            size_t len = dimensions.x;
            std::copy(src, src+len, dst);
        }
    }

    return reinterpret_cast<void*>(data);
}

template<class T>
void* VolumeAtomic<T>::getSliceData(const size_t firstSlice, const size_t lastSlice) const {
    //check for wrong parameter
    if(getDimensions().z <= lastSlice)
        throw std::invalid_argument("lastSlice is out of volume dimension!!!");
    if(firstSlice > lastSlice)
        throw std::invalid_argument("firstSlice has to be less or equal lastSlice!!!");

    //create buffer
    tgt::svec3 dataDims = getDimensions();
    size_t bufferSize = dataDims.x*dataDims.y*(lastSlice-firstSlice+1);
    T* data = new T[bufferSize];
    //determine parameters
    size_t initialStartPos = dataDims.x * dataDims.y * firstSlice;
    //copy data
    T* dst = data;
    const T* src = data_ + initialStartPos;
    size_t len = bufferSize;
    std::copy(src, src+len, dst);

    return reinterpret_cast<void*>(data);
}

/*
  Methods for accessing the voxels
*/

/// just get a proper pointer for read and write access
template<class T>
inline T* VolumeAtomic<T>::voxel() {
    return data_;
}

/// just get a proper pointer for read access only
template<class T>
inline const T* VolumeAtomic<T>::voxel() const {
    return data_;
}

/// get or set voxel
template<class T>
inline T& VolumeAtomic<T>::voxel(size_t x, size_t y, size_t z) {
    tgtAssert(x < size_t(dimensions_.x), "x index out of bounds");
    tgtAssert(y < size_t(dimensions_.y), "y index out of bounds");
    tgtAssert(z < size_t(dimensions_.z), "z index out of bounds");

    return data_[calcPos(dimensions_, x, y, z)];
}

/// get voxel
template<class T>
inline const T& VolumeAtomic<T>::voxel(size_t x, size_t y, size_t z) const {
    tgtAssert(x < size_t(dimensions_.x), "x index out of bounds");
    tgtAssert(y < size_t(dimensions_.y), "y index out of bounds");
    tgtAssert(z < size_t(dimensions_.z), "z index out of bounds");

    return data_[calcPos(dimensions_, x, y, z)];
}

/// get or set voxel
template<class T>
inline T& VolumeAtomic<T>::voxel(const tgt::svec3& pos) {
    tgtAssert(pos.x < dimensions_.x, "x index out of bounds");
    tgtAssert(pos.y < dimensions_.y, "y index out of bounds");
    tgtAssert(pos.z < dimensions_.z, "z index out of bounds");

    return data_[calcPos(dimensions_, pos)];
}

/// get voxel
template<class T>
inline const T& VolumeAtomic<T>::voxel(const tgt::svec3& pos) const {
    tgtAssert(pos.x < dimensions_.x, "x index out of bounds");
    tgtAssert(pos.y < dimensions_.y, "y index out of bounds");
    tgtAssert(pos.z < dimensions_.z, "z index out of bounds");

    return data_[calcPos(dimensions_, pos)];
}

/// get or set voxel
template<class T>
inline T& VolumeAtomic<T>::voxel(size_t i) {
    tgtAssert(i < numVoxels_, "index out of bounds");

    return data_[i];
}

/// get voxel
template<class T>
inline const T& VolumeAtomic<T>::voxel(size_t i) const {
    tgtAssert(i < numVoxels_, "index out of bounds");

    return data_[i];
}


} // namespace voreen

#endif // VRN_VOLUMEATOMIC_H
