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

#ifndef VRN_VOLUMEOPERATORMIRROR_H
#define VRN_VOLUMEOPERATORMIRROR_H

#include "voreen/core/datastructures/volume/volumeoperator.h"

namespace voreen {

///Mirrors the volume on the X axis.
class VRN_CORE_API VolumeOperatorMirrorXBase : public UnaryVolumeOperatorBase {
public:
    virtual Volume* apply(const VolumeBase* vh) const = 0;
};

template<typename T>
class VolumeOperatorMirrorXGeneric : public VolumeOperatorMirrorXBase {
public:
    virtual Volume* apply(const VolumeBase* vh) const;
    IS_COMPATIBLE
};

template<typename T>
Volume* VolumeOperatorMirrorXGeneric<T>::apply(const VolumeBase* vh) const {
    const VolumeRAM* v = vh->getRepresentation<VolumeRAM>();
    if(!v)
        return 0;

    const VolumeAtomic<T>* va = dynamic_cast<const VolumeAtomic<T>*>(v);
    if(!va)
        return 0;

    VolumeAtomic<T>* mirror = va->clone();

    //VRN_FOR_EACH_VOXEL_WITH_PROGRESS(i, tgt::ivec3(0), mirror->getDimensions(), progressBar_) //FIXME
    VRN_FOR_EACH_VOXEL(i, tgt::svec3::zero, mirror->getDimensions())
        mirror->voxel(i) = va->voxel(mirror->getDimensions().x - i.x - 1, i.y, i.z);

    //if (progressBar_) //FIXME
        //progressBar_->setProgress(1.f);

    return new Volume(mirror, vh);
}

typedef UniversalUnaryVolumeOperatorGeneric<VolumeOperatorMirrorXBase> VolumeOperatorMirrorX;

// ============================================================================
///Mirrors the volume on the Y axis.
class VolumeOperatorMirrorYBase : public UnaryVolumeOperatorBase {
public:
    virtual Volume* apply(const VolumeBase* volume) const = 0;
};

template<typename T>
class VolumeOperatorMirrorYGeneric : public VolumeOperatorMirrorYBase {
public:
    virtual Volume* apply(const VolumeBase* vh) const;
    IS_COMPATIBLE
};

template<typename T>
Volume* VolumeOperatorMirrorYGeneric<T>::apply(const VolumeBase* vh) const {
    const VolumeRAM* v = vh->getRepresentation<VolumeRAM>();
    if(!v)
        return 0;

    const VolumeAtomic<T>* va = dynamic_cast<const VolumeAtomic<T>*>(v);
    if(!va)
        return 0;

    VolumeAtomic<T>* mirror = va->clone();

    //VRN_FOR_EACH_VOXEL_WITH_PROGRESS(i, tgt::ivec3(0), mirror->getDimensions(), progressBar_) //FIXME
    VRN_FOR_EACH_VOXEL(i, tgt::svec3::zero, mirror->getDimensions())
        mirror->voxel(i) = va->voxel(i.x, va->getDimensions().y - i.y - 1, i.z);

    //if (progressBar_) //FIXME
        //progressBar_->setProgress(1.f);

    return new Volume(mirror, vh);
}

typedef UniversalUnaryVolumeOperatorGeneric<VolumeOperatorMirrorYBase> VolumeOperatorMirrorY;

// ============================================================================
///Mirrors the volume on the Z axis.
class VolumeOperatorMirrorZBase : public UnaryVolumeOperatorBase {
public:
    virtual Volume* apply(const VolumeBase* volume) const = 0;
};

template<typename T>
class VolumeOperatorMirrorZGeneric : public VolumeOperatorMirrorZBase {
public:
    virtual Volume* apply(const VolumeBase* vh) const;
    IS_COMPATIBLE
};

template<typename T>
Volume* VolumeOperatorMirrorZGeneric<T>::apply(const VolumeBase* vh) const {
    const VolumeRAM* v = vh->getRepresentation<VolumeRAM>();
    if(!v)
        return 0;

    const VolumeAtomic<T>* va = dynamic_cast<const VolumeAtomic<T>*>(v);
    if(!va)
        return 0;

    VolumeAtomic<T>* mirror = va->clone();

    //VRN_FOR_EACH_VOXEL_WITH_PROGRESS(i, tgt::ivec3(0), mirror->getDimensions(), progressBar_)
    VRN_FOR_EACH_VOXEL(i, tgt::svec3::zero, mirror->getDimensions())
        mirror->voxel(i) = va->voxel(i.x, i.y, mirror->getDimensions().z - i.z - 1);

    //if (progressBar_) //FIXME
        //progressBar_->setProgress(1.f);

    return new Volume(mirror, vh);
}

typedef UniversalUnaryVolumeOperatorGeneric<VolumeOperatorMirrorZBase> VolumeOperatorMirrorZ;

// ============================================================================

} // namespace

#endif // VRN_VOLUMEOPERATOR_H
