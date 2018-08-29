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

#ifndef VRN_VOLUMEFUSION_H
#define VRN_VOLUMEFUSION_H

#include <algorithm>
#include <new>
#include <cstring>

#include "tgt/assert.h"

#include "voreen/core/datastructures/volume/volumeatomic.h"

namespace voreen {

/**
 * This class handles several VolumeAtomic objects of the same type as one
 * merged multichannel Volume.
 * All \a voxel methods have an additonal parameter which specifies the channel.
 */
template<class T, size_t N>
class VRN_CORE_API VolumeFusion : public VolumeRAM {
public:
    typedef T* VoxelType; ///< So that VolumeRAM::getData\< VolumeFusion\<T, N> will return T**.

    //
    // constructors
    //

    /// Creates a VolumeFusion with N new VolumeAtomic objects.
    VolumeFusion();

    /// Creates a VolumeFusion with N existing VolumeAtomic objects
    VolumeFusion(VolumeAtomic<T>* volumes[N]);

    virtual VolumeFusion<T, N>* clone() const;

    /// The \p data parameter must be a T**.
    virtual VolumeFusion<T, N>* clone(void* data) const;

    virtual VolumeFusion<T, N>* createNew(const tgt::svec3& dimensions, bool allocMem = false) const;

    virtual VolumeFusion<T, N>* getSubVolume(tgt::svec3 dimensions, tgt::svec3 offset = tgt::svec3::zero) const;
    virtual void setSubVolume(const VolumeRAM* vol, tgt::svec3 offset = tgt::svec3::zero);

    //
    // getters and setters.
    // inline here since the definition must be known to the compiler anyway
    //

    /// get a volume channel
    inline VolumeAtomic<T>* getVolumeAtomic(size_t channel) {
        return volumes_[channel];
    }
    /// set a volume channel
    inline void setVolumeAtomic(VolumeAtomic<T>* volume, size_t channel) {
        volumes_[channel] = volume;
    }

    //
    // methods for accessing the voxels
    //

    /// just get a proper pointer for read and write access
    inline T* voxel(size_t channel) {
        tgtAssert(channel < N, "channel must be less than N");
        return volumes_[channel]->voxel();
    }
    /// just get a proper pointer for read access only
    inline const T* voxel(size_t channel) const {
        tgtAssert(channel < N, "channel must be less than N");
        return volumes_[channel]->voxel();
    }
    /// get or set voxel
    inline T& voxel(size_t x, size_t y, size_t z, size_t channel) {
        tgtAssert(channel < N, "channel must be less than N");
        return volumes_[channel]->voxel(x, y, z);
    }
    /// get voxel
    inline const T& voxel(size_t x, size_t y, size_t z, size_t channel) const {
        tgtAssert(channel < N, "channel must be less than N");
        return volumes_[channel]->voxel(x, y, z);
    }
    /// get or set voxel
    inline T& voxel(const tgt::ivec3& pos, size_t channel) {
        tgtAssert(channel < N, "channel must be less than N");
        return volumes_[channel]->voxel(pos);
    }
    /// get voxel
    inline const T& voxel(const tgt::ivec3& pos, size_t channel) const {
        tgtAssert(channel < N, "channel must be less than N");
        return volumes_[channel]->voxel(pos);
    }
    /// get or set voxel
    inline T& voxel(size_t i, size_t channel) {
        tgtAssert(channel < N, "channel must be less than N");
        return volumes_[channel]->voxel(i);
    }
    /// get voxel
    inline const T& voxel(size_t i, size_t channel) const {
        tgtAssert(channel < N, "channel must be less than N");
        return volumes_[channel]->voxel(i);
    }

    // getVoxelNormalized and setVoxelNormalized
    virtual float getVoxelNormalized(const tgt::ivec3& pos, size_t channel = 0) const;
    virtual float getVoxelNormalized(size_t x, size_t y, size_t z, size_t channel = 0) const;
    virtual void setVoxelNormalized(float value, const tgt::ivec3& pos, size_t channel = 0);
    virtual void setVoxelNormalized(float value, size_t x, size_t y, size_t z, size_t channel = 0);

    //
    // getters and setters
    //

    virtual size_t getBitsAllocated() const;
    virtual size_t getNumChannels() const;
    virtual size_t getBytesPerVoxel() const;
    /// Returns the number of bytes held in one channel.
    virtual size_t getNumBytes() const;

    //
    // further methods
    //

    virtual void clear();

    /// Returns the minimum value of the given channel.
    T min(size_t channel) const;
    /// Returns the maximum value of the given channel.
    T max(size_t channel) const;
    /// Returns the minimum value of all channels.
    T min() const;
    /// Returns the maximum value of all channels.
    T max() const;

    /**
     * Returns an array of size N of type T*. So cast the returning void* to T**.
     * Because this array is allocated on the Heap the caller must destroy it.
     *
     * @return A T** which holds pointers to the data (an T* array of size N).
     *      In the very unlikely case, that not enough memory is available for
     *      the result array (this array isn't huge) 0 is returned.
     */
    virtual void* getData();
    virtual VolumeFusion<T, N>* createSubset(const tgt::ivec3 &pos, const tgt::ivec3 &size) const;
    virtual VolumeFusion<T, N>* resample(const tgt::ivec3& newDims, Filter filter) const;
    virtual VolumeFusion<T, N>* mirrorZ() const;

protected:
    VolumeAtomic<T>* volumes_[N]; ///< N VolumeAtomic instances of type T in an array.
};

//------------------------------------------------------------------------------
//  non inline implementation
//------------------------------------------------------------------------------

//
// constructors
//

template<class T, size_t N>
VolumeFusion<T, N>::VolumeFusion() {
    try {
        for (size_t i = 0; i < N; ++i)
            volumes_[i] = new VolumeAtomic<T>();
    }
    catch (std::bad_alloc) {
        throw; // throw it to the caller
    }
}

template<class T, size_t N>
VolumeFusion<T, N>::VolumeFusion(VolumeAtomic<T>* volumes[N]) {
    memcpy(volumes, volumes_, sizeof(VolumeAtomic<T>*) * N);
}

template<class T, size_t N>
VolumeFusion<T, N>* VolumeFusion<T, N>::clone() const {
    VolumeAtomic<T>* volumes[N];

    try {
        for (size_t i = 0; i < N; ++i)
            volumes[i] = volumes_[i]->clone();
    }
    catch (std::bad_alloc) {
        throw; // throw it to the caller
    }

    return new VolumeFusion<T, N>(volumes);
}

template<class T, size_t N>
VolumeFusion<T, N>* VolumeFusion<T, N>::clone(void* data) const {
    VolumeAtomic<T>* volumes[N];
    T** t = (T**) data;

    try {
        for (size_t i = 0; i < N; ++i) {
            if (t == 0)
                volumes[i] = volumes_[i]->clone(0);
            else
                volumes[i] = volumes_[i]->clone(t[i]);
        }
    }
    catch (std::bad_alloc) {
        throw; // throw it to the caller
    }

    return new VolumeFusion<T, N>(volumes);
}

template<class T, size_t N>
VolumeFusion<T, N>* VolumeFusion<T, N>::createNew(const tgt::svec3& dimensions, bool allocMem) const {
    VolumeAtomic<T>* volumes[N];

    // create new volumes
    for (size_t i = 0; i < N; ++i) {
        volumes[i] = volumes_[i]->createNew(dimensions, allocMem);
    }

    return new VolumeFusion<T, N>(volumes);
}

template<class T, size_t N>
VolumeFusion<T, N>* VolumeFusion<T, N>::getSubVolume(tgt::svec3 dimensions, tgt::svec3 offset) const {
    VolumeAtomic<T>* volumes[N];

    try {
        for (size_t i = 0; i < N; ++i) {
                volumes[i] = volumes_[i]->getSubVolume(dimensions, offset);
        }
    }
    catch (std::bad_alloc) {
        throw; // throw it to the caller
    }

    return new VolumeFusion<T, N>(volumes);
}

template<class T, size_t N>
void VolumeFusion<T, N>::setSubVolume(const VolumeRAM* vol, tgt::svec3 offset) {
    const VolumeFusion<T, N>* vols = static_cast<const VolumeFusion<T, N>*>(vol);

    for (size_t i = 0; i < N; ++i) {
        volumes_[i]->setSubVolume(vols->getVolumeAtomic(i), offset);
    }
}

//
// getVoxelNormalized and setVoxelNormalized
//

template<class T, size_t N>
float VolumeFusion<T, N>::getVoxelNormalized(const tgt::ivec3& pos, size_t channel /*= 0*/) const {
    return volumes_[channel]->getVoxelNormalized(pos);
}

template<class T, size_t N>
float VolumeFusion<T, N>::getVoxelNormalized(size_t x, size_t y, size_t z, size_t channel /*= 0*/) const {
    return volumes_[channel]->getVoxelNormalized(x, y, z);
}

template<class T, size_t N>
void VolumeFusion<T, N>::setVoxelNormalized(float value, const tgt::ivec3& pos, size_t channel /*= 0*/) {
    volumes_[channel]->setVoxelNormalized(value, pos);
}

template<class T, size_t N>
void VolumeFusion<T, N>::setVoxelNormalized(float value, size_t x, size_t y, size_t z, size_t channel /*= 0*/) {
    volumes_[channel]->setVoxelNormalized(value, x, y, z);
}

//
// getters and setters
//

template<class T, size_t N>
size_t VolumeFusion<T, N>::getBitsAllocated() const {
    return volumes_[0]->getBitsAllocated();
}

template<class T, size_t N>
size_t VolumeFusion<T, N>::getNumChannels() const {
    return N;
}

template<class T, size_t N>
size_t  VolumeFusion<T, N>::getBytesPerVoxel() const {
    return volumes_[0]->getBytesPerVoxel();
}

template<class T, size_t N>
size_t VolumeFusion<T, N>::getNumBytes() const {
    return volumes_[0]->getNumBytes();
}

template<class T, size_t N>
void* VolumeFusion<T, N>::getData() {
    // allocate the array which holds the result, get 0 if memory isn't available.
    T** result = new(std::nothrow) T*[N];

    if (result) {
        // fill the array with data pointers
        for (size_t i = 0; i < N; ++i)
            result[i] = VolumeRAM::getData< VolumeAtomic<T> >(volumes_[i]);
    }

    return reinterpret_cast<void*>(result);
}

//
// further methods
//

template<class T, size_t N>
T VolumeFusion<T, N>::min(size_t channel) const {
    VolumeAtomic<T>* v = volumes_[channel];
    return std::min_element(v->voxel(), v->voxel() + v->getNumVoxels());
}

template<class T, size_t N>
T VolumeFusion<T, N>::max(size_t channel) const {
    VolumeAtomic<T>* v = volumes_[channel];
    return std::max_element(v->voxel(), v->voxel() + v->getNumVoxels());
}

template<class T, size_t N>
T VolumeFusion<T, N>::min() const {
    T mins[N];
    for (size_t i = 0; i < N; ++i)
        mins[i] = min(i);

    return std::min_element(mins, mins + N);
}

template<class T, size_t N>
T VolumeFusion<T, N>::max() const {
    T maxis[N];
    for (size_t i = 0; i < N; ++i)
        maxis[i] = max(i);

    return std::max_element(maxis, maxis + N);
}

template<class T, size_t N>
void VolumeFusion<T, N>::clear() {
    for (size_t i = 0; i < N; ++i)
        volumes_[i]->clear();
}

template<class T, size_t N>
VolumeFusion<T, N>* VolumeFusion<T, N>::createSubset(const tgt::ivec3 &pos, const tgt::ivec3 &size) const {
    VolumeAtomic<T>* volumes[N];

    try {
        for (size_t i = 0; i < N; ++i)
            volumes[i] = volumes_[i]->createSubset(pos, size);
    }
    catch (std::bad_alloc) {
        throw; // throw it to the caller
    }

    return new VolumeFusion<T, N>(volumes);
}

template<class T, size_t N>
VolumeFusion<T, N>* VolumeFusion<T, N>::resample(const tgt::ivec3& newDims, Filter filter) const {
    VolumeAtomic<T>* volumes[N];

    try {
        for (size_t i = 0; i < N; ++i)
            volumes[i] = volumes_[i]->resample(newDims, filter);
    }
    catch (std::bad_alloc) {
        throw; // throw it to the caller
    }

    return new VolumeFusion<T, N>(volumes);
}

template<class T, size_t N>
VolumeFusion<T, N>* VolumeFusion<T, N>::mirrorZ() const {
    VolumeAtomic<T>* volumes[N];

    for (size_t i = 0; i < N; ++i)
        volumes[i] = volumes_[i]->mirrorZ();

    return new VolumeFusion<T, N>(volumes);
}

} // namespace voreen

#endif // VRN_VOLUMEFUSION_H
