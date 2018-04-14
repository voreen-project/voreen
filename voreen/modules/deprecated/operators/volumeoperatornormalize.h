/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#ifndef VRN_VOLUMEOPERATORNORMALIZE_H
#define VRN_VOLUMEOPERATORNORMALIZE_H

#include "voreen/core/datastructures/volume/volumeoperator.h"

namespace voreen {

class VolumeOperatorNormalizeBase : public UnaryVolumeOperatorBase {
public:
    virtual Volume* apply(const VolumeBase* vh) const = 0;
};

template<typename T>
class VolumeOperatorNormalizeGeneric : public VolumeOperatorNormalizeBase {
public:
    Volume* apply(const VolumeBase* vh) const;
    IS_COMPATIBLE
};

template<typename T>
Volume* VolumeOperatorNormalizeGeneric<T>::apply(const VolumeBase* vh) const {
    const VolumeRAM* v = vh->getRepresentation<VolumeRAM>();
    if (!v)
        return 0;

    const VolumeAtomic<T>* va = dynamic_cast<const VolumeAtomic<T>*>(v);
    if (!va)
        return 0;

    VolumeAtomic<T>* normalized = va->clone();

    T minLocalValue = va->min();
    T maxLocalValue = va->max();
    T maxGlobalValue = static_cast<T>(1 << va->getBitsAllocated());

    VRN_FOR_EACH_VOXEL(i, tgt::ivec3(0), va->getDimensions()) {
        T value = va->voxel(i);
        normalized->voxel(i) = T((float(value - minLocalValue) / float(maxLocalValue - minLocalValue)) * maxGlobalValue);
    }

    return new Volume(normalized, vh);
}

typedef UniversalUnaryVolumeOperatorGeneric<VolumeOperatorNormalizeBase> VolumeOperatorNormalize;

} // namespace

#endif // VRN_VOLUMEOPERATORNORMALIZE_H
