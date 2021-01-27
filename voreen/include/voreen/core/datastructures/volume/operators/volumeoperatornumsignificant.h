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

#ifndef VRN_VOLUMEOPERATORNUMSIGNIFICANT_H
#define VRN_VOLUMEOPERATORNUMSIGNIFICANT_H

#include "voreen/core/datastructures/volume/volumeoperator.h"

namespace voreen {

/**
 * Returns the number of significant voxels, i.e., the number of voxels with
 * a value greater than the smallest possible value that can be stored by the volume.
 */
class VRN_CORE_API VolumeOperatorNumSignificantBase : public UnaryVolumeOperatorBase {
public:
    virtual size_t apply(const VolumeBase* volume) const = 0;
};

// ============================================================================

// Generic implementation:
template<typename T>
class VolumeOperatorNumSignificantGeneric : public VolumeOperatorNumSignificantBase {
public:
    virtual size_t apply(const VolumeBase* volume) const;
    //Implement isCompatible using a handy macro:
    IS_COMPATIBLE
};

template<typename T>
size_t VolumeOperatorNumSignificantGeneric<T>::apply(const VolumeBase* vh) const {
    const VolumeRAM* v = vh->getRepresentation<VolumeRAM>();
    if(!v)
        return 0;

    const VolumeAtomic<T>* volume = dynamic_cast<const VolumeAtomic<T>*>(v);
    if(!volume)
        return 0;

    size_t result = 0;
    VRN_FOR_EACH_VOXEL(pos, tgt::svec3::zero, volume->getDimensions()) {
        if (volume->voxel(pos) > VolumeElement<T>::rangeMin())
            result++;
    }
    return result;
}

// ============================================================================

template<typename S>
class VolumeOperatorNumSignificantGeneric<tgt::Vector2<S> >: public VolumeOperatorNumSignificantBase {
public:
    virtual size_t apply(const VolumeBase* volume) const;
    bool isCompatible(const VolumeBase* volume) const {
        const VolumeRAM* v = volume->getRepresentation<VolumeRAM>();
        if(!v)
            return false;
        const VolumeAtomic<tgt::Vector2<S> >* va = dynamic_cast<const VolumeAtomic<tgt::Vector2<S> >*>(v);
        if(!va)
            return false;
        return true;
    }
};

template<typename S>
size_t VolumeOperatorNumSignificantGeneric<tgt::Vector2<S> >::apply(const VolumeBase* vh) const {
    const VolumeRAM* v = vh->getRepresentation<VolumeRAM>();
    if(!v)
        return 0;

    const VolumeAtomic<tgt::Vector2<S> >* volume = dynamic_cast<const VolumeAtomic<tgt::Vector2<S> >*>(v);
    if(!volume)
        return 0;

    size_t result = 0;
    VRN_FOR_EACH_VOXEL(pos, tgt::svec3::zero, volume->getDimensions()) {
        if (tgt::hadd(volume->voxel(pos)) > VolumeElement<S>::rangeMin())
            result++;
    }
    return result;
}

// ============================================================================

template<typename S>
class VolumeOperatorNumSignificantGeneric<tgt::Vector3<S> >: public VolumeOperatorNumSignificantBase {
public:
    virtual size_t apply(const VolumeBase* volume) const;
    bool isCompatible(const VolumeBase* volume) const {
        const VolumeRAM* v = volume->getRepresentation<VolumeRAM>();
        if(!v)
            return false;
        const VolumeAtomic<tgt::Vector3<S> >* va = dynamic_cast<const VolumeAtomic<tgt::Vector3<S> >*>(v);
        if(!va)
            return false;
        return true;
    }
};

template<typename S>
size_t VolumeOperatorNumSignificantGeneric<tgt::Vector3<S> >::apply(const VolumeBase* vh) const {
    const VolumeRAM* v = vh->getRepresentation<VolumeRAM>();
    if(!v)
        return 0;

    const VolumeAtomic<tgt::Vector3<S> >* volume = dynamic_cast<const VolumeAtomic<tgt::Vector3<S> >*>(v);
    if(!volume)
        return 0;

    size_t result = 0;
    VRN_FOR_EACH_VOXEL(pos, tgt::svec3::zero, volume->getDimensions()) {
        if (tgt::hadd(volume->voxel(pos)) > VolumeElement<S>::rangeMin())
            result++;
    }
    return result;
}


// ============================================================================

template<typename S>
class VolumeOperatorNumSignificantGeneric<tgt::Vector4<S> >: public VolumeOperatorNumSignificantBase {
public:
    virtual size_t apply(const VolumeBase* volume) const;
    bool isCompatible(const VolumeBase* volume) const {
        const VolumeRAM* v = volume->getRepresentation<VolumeRAM>();
        if(!v)
            return false;
        const VolumeAtomic<tgt::Vector4<S> >* va = dynamic_cast<const VolumeAtomic<tgt::Vector4<S> >*>(v);
        if(!va)
            return false;
        return true;
    }
};

template<typename S>
size_t VolumeOperatorNumSignificantGeneric<tgt::Vector4<S> >::apply(const VolumeBase* vh) const {
    const VolumeRAM* v = vh->getRepresentation<VolumeRAM>();
    if(!v)
        return 0;

    const VolumeAtomic<tgt::Vector4<S> >* volume = dynamic_cast<const VolumeAtomic<tgt::Vector4<S> >*>(v);
    if(!volume)
        return 0;

    size_t result = 0;
    VRN_FOR_EACH_VOXEL(pos, tgt::svec3::zero, volume->getDimensions()) {
        if (tgt::hadd(volume->voxel(pos)) > VolumeElement<S>::rangeMin())
            result++;
    }
    return result;
}


// ============================================================================

template<typename S>
class VolumeOperatorNumSignificantGeneric<Tensor2<S> >: public VolumeOperatorNumSignificantBase {
public:
    virtual size_t apply(const VolumeBase* volume) const;
    bool isCompatible(const VolumeBase* volume) const {
        const VolumeRAM* v = volume->getRepresentation<VolumeRAM>();
        if(!v)
            return false;
        const VolumeAtomic<Tensor2<S> >* va = dynamic_cast<const VolumeAtomic<Tensor2<S> >*>(v);
        if(!va)
            return false;
        return true;
    }
};

template<typename S>
size_t VolumeOperatorNumSignificantGeneric<Tensor2<S> >::apply(const VolumeBase* vh) const {
    const VolumeRAM* v = vh->getRepresentation<VolumeRAM>();
    if(!v)
        return 0;

    const VolumeAtomic<Tensor2<S> >* volume = dynamic_cast<const VolumeAtomic<Tensor2<S> >*>(v);
    if(!volume)
        return 0;

    size_t result = 0;
    VRN_FOR_EACH_VOXEL(pos, tgt::svec3::zero, volume->getDimensions()) {
        const Tensor2<S>& t = volume->voxel(pos);
        if ((t.Dxx + t.Dxy + t.Dxz + t.Dyy + t.Dyz + t.Dzz) > VolumeElement<S>::rangeMin())
            result++;
    }
    return result;
}

// ============================================================================

typedef UniversalUnaryVolumeOperatorGeneric<VolumeOperatorNumSignificantBase> VolumeOperatorNumSignificant;

} // namespace

#endif // VRN_VOLUMEOPERATOR_H
