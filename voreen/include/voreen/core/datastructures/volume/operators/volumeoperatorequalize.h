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

#ifndef VRN_VOLUMEOPERATOREQUALIZE_H
#define VRN_VOLUMEOPERATOREQUALIZE_H

#include "voreen/core/datastructures/volume/operators/volumeoperatorminmax.h"
#include "voreen/core/datastructures/volume/volumeoperator.h"

namespace voreen {

/**
 * Performs a linear histogram equalization of the input volume.
 */
class VRN_CORE_API VolumeOperatorEqualizeBase : public UnaryVolumeOperatorBase {
public:
    virtual Volume* apply(const VolumeBase* volume, ProgressBar* progressBar = 0) const = 0;
};

template<typename T>
class VolumeOperatorEqualizeGeneric : public VolumeOperatorEqualizeBase {
public:
    virtual Volume* apply(const VolumeBase* volume, ProgressBar* progressBar = 0) const;
    //Implement isCompatible using a handy macro:
    IS_COMPATIBLE
};

template<typename T>
Volume* VolumeOperatorEqualizeGeneric<T>::apply(const VolumeBase* vh, ProgressBar* progressBar) const {
    const VolumeRAM* v = vh->getRepresentation<VolumeRAM>();
    if(!v)
        return 0;

    const VolumeAtomic<T>* va = dynamic_cast<const VolumeAtomic<T>*>(v);
    if (!va)
        return 0;

    VolumeAtomic<T>* out = va->clone();
    RealWorldMapping mapping = vh->getRealWorldMapping();

    if (VolumeElement<T>::isInteger()) {
        T min = VolumeOperatorMinValue::apply(va);
        T max = VolumeOperatorMaxValue::apply(va);
        T rangeMin = VolumeElement<T>::rangeMin();
        T rangeMax = VolumeElement<T>::rangeMax();
        double scale = static_cast<double>(VolumeElement<T>::rangeMax() - VolumeElement<T>::rangeMin()) / static_cast<double>(max - min);
        float invScale = static_cast<float>(1.0/scale);

        double offset;
        float relOffset;
        if (VolumeElement<T>::isSigned()) {
            offset = min/2.0 + max/2.0;
            relOffset = static_cast<float>(offset*2.0 / (rangeMax-rangeMin));
        }
        else { // unsigned
            offset = static_cast<double>(min);
            relOffset = static_cast<float>(offset / (rangeMax-rangeMin));
        }

        VRN_FOR_EACH_VOXEL_WITH_PROGRESS(index, tgt::svec3(0, 0, 0), out->getDimensions(), progressBar) {
            out->voxel(index) = static_cast<T>( ((double)out->voxel(index) - offset) * scale );
        }

        mapping = RealWorldMapping::combine(RealWorldMapping(invScale, relOffset, ""), mapping);
    }
    else {
        LWARNINGC("voreen.core.VolumeOperatorEqualize", "Histogram equalization only supported for integer volume types");
    }

    if (progressBar)
        progressBar->setProgress(1.f);

    Volume* result = new Volume(out, vh);
    result->setRealWorldMapping(mapping);
    return result;
}

typedef UniversalUnaryVolumeOperatorGeneric<VolumeOperatorEqualizeBase> VolumeOperatorEqualize;

} // namespace

#endif // VRN_VOLUMEOPERATOREQUALIZE_H
