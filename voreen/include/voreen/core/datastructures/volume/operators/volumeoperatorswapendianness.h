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

#ifndef VRN_VOLUMEOPERATORSWAPENDIANNESS_H
#define VRN_VOLUMEOPERATORSWAPENDIANNESS_H

#include "voreen/core/datastructures/volume/volumeoperator.h"

namespace voreen {

// Base class, defines interface for the operator (-> apply):
class VRN_CORE_API VolumeOperatorSwapEndiannessBase : public UnaryVolumeOperatorBase {
public:
    virtual void apply(Volume* volume) const = 0;
};

// Generic implementation:
template<typename T>
class VolumeOperatorSwapEndiannessGeneric : public VolumeOperatorSwapEndiannessBase {
public:
    virtual void apply(Volume* volume) const;
    //Implement isCompatible using a handy macro:
    IS_COMPATIBLE
};

template<typename T>
void VolumeOperatorSwapEndiannessGeneric<T>::apply(Volume* vh) const {
    VolumeRAM* v = vh->getWritableRepresentation<VolumeRAM>();
    if(!v)
        return;

    if(v->getBytesPerVoxel() == 1)
        return; //no need to swap

    VolumeAtomic<T>* va = dynamic_cast<VolumeAtomic<T>*>(v);
    if(!va)
        return;

    T* voxel = va->voxel();
    for (size_t i = 0; i < va->getNumVoxels(); ++i)
        voxel[i] = VolumeElement<T>::swapEndianness(voxel[i]);
}

typedef UniversalUnaryVolumeOperatorGeneric<VolumeOperatorSwapEndiannessBase> VolumeOperatorSwapEndianness;

} // namespace

#endif // VRN_VOLUMEOPERATOR_H
