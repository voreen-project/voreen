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

#ifndef VRN_VOLUMEOPERATORGRADIENTDESCENT_H
#define VRN_VOLUMEOPERATORGRADIENTDESCENT_H

#include "voreen/core/datastructures/volume/volumeoperator.h"

#define VRN_FOR_EACH_VOXEL_INT(INDEX, START, END)\
    for (tgt::ivec3 (INDEX) = (START); (INDEX).z < (END).z; ++(INDEX).z)\
        for ((INDEX).y = (START).y; (INDEX).y < (END).y; ++(INDEX).y)\
            for ((INDEX).x = (START).x; (INDEX).x < (END).x; ++(INDEX).x)

namespace voreen {

// Generic implementation Squared Euclidean distance transform:
class VRN_CORE_API VolumeOperatorGradientDescentBase : public UnaryVolumeOperatorBase {
public:
    virtual Volume* apply(const VolumeBase* volume, ProgressReporter* progressReporter = 0) const = 0;
};

//specific implementation squared euclidean distance transform
template<typename T>
class VolumeOperatorGradientDescentGeneric : public VolumeOperatorGradientDescentBase {
public:
    virtual Volume* apply(const VolumeBase* v1, ProgressReporter* progressReporter = 0) const;

    IS_COMPATIBLE
};

template<typename T>
Volume* VolumeOperatorGradientDescentGeneric<T>::apply(const VolumeBase* vb, ProgressReporter* pR) const{
    const VolumeRAM* v = vb->getRepresentation<VolumeRAM>();
    if (!v)
        return 0;

    const VolumeAtomic<T>* va = dynamic_cast<const VolumeAtomic<T>*>(v);
    if (!va)
        return 0;

    tgt::svec3 dim = va->getDimensions();
    VolumeRAM_UInt16* descendants = new VolumeRAM_UInt16(dim);
    memset(descendants->getData(), 0, descendants->getNumBytes());

    VRN_FOR_EACH_VOXEL(idx, tgt::svec3::zero, dim){
        if(va->voxel(idx) > T(0)){
            //descent Voxel in vol
            tgt::svec3 oldPos = idx;
            tgt::svec3 newPos = idx;
            do{
                oldPos = newPos;
                T max = va->voxel(newPos);
                VRN_FOR_EACH_VOXEL_INT(sub_idx, tgt::ivec3(-1), tgt::ivec3(2)){
                    if(sub_idx != tgt::ivec3::zero){
                        long long x = newPos.x + sub_idx.x;
                        long long y = newPos.y + sub_idx.y;
                        long long z = newPos.z + sub_idx.z;

                        tgt::svec3 index = tgt::clamp(tgt::svec3((size_t)x,(size_t)y,(size_t)z), tgt::svec3::zero, dim - tgt::svec3::one);

                        T value = va->voxel(index);
                        if(value > max){
                            max = value;
                            newPos = index;
                        }
                    }
                }
            }while(oldPos != newPos);
            descendants->voxel(newPos)++;
        }
    }

    return new Volume(descendants, vb);
}


typedef UniversalUnaryVolumeOperatorGeneric<VolumeOperatorGradientDescentBase> VolumeOperatorGradientDescent;

} // namespace


#endif //VRN_VOLUMEOPERATORGRADIENTDESCENT_H
