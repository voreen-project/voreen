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

#ifndef VRN_VOLUMEOPERATORISUNIFORM_H
#define VRN_VOLUMEOPERATORISUNIFORM_H

#include "voreen/core/datastructures/volume/volumeoperator.h"

namespace voreen {

// Base class, defines interface for the operator (-> apply):
class VRN_CORE_API VolumeOperatorIsUniformBase : public UnaryVolumeOperatorBase {
public:
    virtual bool apply(const VolumeBase* volume) const = 0;
};

// Generic implementation:
template<typename T>
class VolumeOperatorIsUniformGeneric : public VolumeOperatorIsUniformBase {
public:
    virtual bool apply(const VolumeBase* volume) const;
    //Implement isCompatible using a handy macro:
    IS_COMPATIBLE
};

template<typename T>
bool VolumeOperatorIsUniformGeneric<T>::apply(const VolumeBase* vh) const {
    const VolumeRAM* v = vh->getRepresentation<VolumeRAM>();
    if(!v)
        return 0;

    const VolumeAtomic<T>* volume = dynamic_cast<const VolumeAtomic<T>*>(v);
    if(!volume)
        return 0;

    T firstVoxel = volume->voxel(0);
    bool allVoxelsEqual=true;
    for (size_t i=1; i < volume->getNumVoxels(); i++) {
        T currentVoxel = volume->voxel(i);
        if (firstVoxel != currentVoxel) {
            allVoxelsEqual = false;
            break;
        }
    }
    return allVoxelsEqual;
}

typedef UniversalUnaryVolumeOperatorGeneric<VolumeOperatorIsUniformBase> VolumeOperatorIsUniform;

} // namespace

#endif // VRN_VOLUMEOPERATOR_H
