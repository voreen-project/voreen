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

#ifndef VRN_VOLUMEOPERATORTRANSPOSE_H
#define VRN_VOLUMEOPERATORTRANSPOSE_H

#include "voreen/core/datastructures/volume/volumeoperator.h"

namespace voreen {

class VRN_CORE_API VolumeOperatorTransposeBase : public UnaryVolumeOperatorBase {
public:
    virtual Volume* apply(const VolumeBase* vh, int a, int b) const = 0;
};

template<typename T>
class VolumeOperatorTransposeGeneric : public VolumeOperatorTransposeBase {
public:
    virtual Volume* apply(const VolumeBase* vh, int a, int b) const;
    IS_COMPATIBLE
};

template<typename T>
Volume* VolumeOperatorTransposeGeneric<T>::apply(const VolumeBase* vh, int a, int b) const {
    const VolumeRAM* v = vh->getRepresentation<VolumeRAM>();
    if(!v)
        return 0;

    const VolumeAtomic<T>* va = dynamic_cast<const VolumeAtomic<T>*>(v);
    if(!va)
        return 0;

    tgt::ivec3 tp(0,1,2);
    tp[a] = b;
    tp[b] = a;

    tgt::svec3 dims = va->getDimensions();
    tgt::svec3 newDims(dims[tp[0]], dims[tp[1]], dims[tp[2]]);

    VolumeAtomic<T>* transposed = new VolumeAtomic<T>(newDims);

    //VRN_FOR_EACH_VOXEL_WITH_PROGRESS(i, tgt::ivec3(0), mirror->getDimensions(), progressBar_) //FIXME
    VRN_FOR_EACH_VOXEL(i, tgt::svec3(0,0,0), dims)
        transposed->voxel(tgt::svec3(i[tp[0]], i[tp[1]], i[tp[2]])) = va->voxel(i);

    //if (progressBar_) //FIXME
        //progressBar_->setProgress(1.f);

    Volume* ret = new Volume(transposed, vh);
    tgt::vec3 sp = ret->getSpacing();
    ret->setSpacing(tgt::vec3(sp[tp[0]], sp[tp[1]], sp[tp[2]]));
    return ret;
}

typedef UniversalUnaryVolumeOperatorGeneric<VolumeOperatorTransposeBase> VolumeOperatorTranspose;

} // namespace

#endif
