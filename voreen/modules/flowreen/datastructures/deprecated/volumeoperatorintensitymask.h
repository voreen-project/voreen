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

#ifndef VRN_VOLUMEOPERATORINTENSITYMASK_H
#define VRN_VOLUMEOPERATORINTENSITYMASK_H

#include "tgt/vector.h"
#include "voreen/core/datastructures/volume/volumeoperator.h"

namespace voreen {

/**
 * A VolumeOperator which returns all voxel positions having values within the
 * thresholds being passed to the ctor.
 * Note that this currently works only for volumes containing intensities, i.e.
 * the template parameter T for the underlying VolumeAtomic<T> may not be of any type
 * which cannot be compared to or casted to float.
 */
class VolumeOperatorIntensityMaskBase : public UnaryVolumeOperatorBase {
public:
    virtual Volume* apply(const VolumeBase* volume, tgt::vec2 thresholds) const = 0;
};

// Generic implementation:
template<typename T>
class VolumeOperatorIntensityMaskGeneric : public VolumeOperatorIntensityMaskBase {
public:
    virtual Volume* apply(const VolumeBase* volume, tgt::vec2 thresholds) const;
    //Implement isCompatible using a handy macro:
    IS_COMPATIBLE
};

template<typename T>
Volume* VolumeOperatorIntensityMaskGeneric<T>::apply(const VolumeBase* vh, tgt::vec2 thresholds) const {
    const VolumeRAM* v = vh->getRepresentation<VolumeRAM>();
    if(!v)
        return 0;

    const VolumeAtomic<T>* volume = dynamic_cast<const VolumeAtomic<T>*>(v);
    if(!volume)
        return 0;


    const tgt::ivec3 dim = volume->getDimensions();
    const T* const data = volume->voxel();

    VolumeAtomic<bool>* subVolume = new VolumeAtomic<bool>(dim);
    bool* const binVoxels = subVolume->voxel();

    for (size_t i = 0; i < volume->getNumVoxels(); ++i) {
        float value = static_cast<float>(data[i]);
        if ((value >= thresholds.x) && (value <= thresholds.y))
            binVoxels[i] = true;
        else
            binVoxels[i] = false;
    }   // for (i

    return new Volume(subVolume, vh);
}

typedef UniversalUnaryVolumeOperatorGeneric<VolumeOperatorIntensityMaskBase> VolumeOperatorIntensityMask;
}   // namespace

#endif
