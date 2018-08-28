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

#ifndef VRN_VOLUMEOPERATORCALCERROR_H
#define VRN_VOLUMEOPERATORCALCERROR_H

#include "voreen/core/datastructures/volume/volumeoperator.h"

namespace voreen {

/**
 * Calculates the root mean square error between the two volumes.
 */
class VRN_CORE_API VolumeOperatorCalcErrorBase : public BinaryVolumeOperatorBase {
public:
    virtual float apply(const VolumeBase* vh1, const VolumeBase* vh2) const = 0;
};

// Generic implementation:
template<typename T>
class VolumeOperatorCalcErrorGeneric : public VolumeOperatorCalcErrorBase {
public:
    virtual float apply(const VolumeBase* vh1, const VolumeBase* vh2) const;
    virtual bool isCompatible(const VolumeBase* vh1, const VolumeBase* vh2) const {
        const VolumeRAM* v1 = vh1->getRepresentation<VolumeRAM>();
        if(!v1)
            return false;
        const VolumeAtomic<T>* va1 = dynamic_cast<const VolumeAtomic<T>*>(v1);
        if(!va1)
            return false;

        const VolumeRAM* v2 = vh2->getRepresentation<VolumeRAM>();
        if(!v2)
            return false;
        const VolumeAtomic<T>* va2 = dynamic_cast<const VolumeAtomic<T>*>(v2);
        if(!va2)
            return false;

        return true;
    }
};

template<typename T>
float VolumeOperatorCalcErrorGeneric<T>::apply(const VolumeBase* vh1, const VolumeBase* vh2) const {
    const VolumeRAM* vol1 = vh1->getRepresentation<VolumeRAM>();
    if(!vol1)
        return 0;

    const VolumeRAM* vol2 = vh2->getRepresentation<VolumeRAM>();
    if(!vol2)
        return 0;

    const VolumeAtomic<T>* v1 = dynamic_cast<const VolumeAtomic<T>*>(vol1);
    if(!v1)
        return 0;

    const VolumeAtomic<T>* v2 = dynamic_cast<const VolumeAtomic<T>*>(vol2);
    if(!v2)
        return 0;

    double errorSum = 0.0;
    tgt::svec3 factor = v1->getDimensions() / v2->getDimensions();

    for (size_t i=0; i<v1->getDimensions().z; i++) {
        for (size_t j=0; j<v1->getDimensions().y; j++) {
            for (size_t k=0; k<v1->getDimensions().x; k++) {
                tgt::svec3 currentPos = tgt::svec3(k,j,i);
                tgt::svec3 smallVolumePos;
                smallVolumePos.x = static_cast<int>( floor(currentPos.x / (float)factor.x));
                smallVolumePos.y = static_cast<int>( floor(currentPos.y / (float)factor.y));
                smallVolumePos.z = static_cast<int>( floor(currentPos.z / (float)factor.z));
                T origVoxel = v1->voxel(currentPos);
                T errVoxel = v2->voxel(smallVolumePos);

                errorSum = errorSum + VolumeElement<T>::calcSquaredDifference(origVoxel, errVoxel);
            }
        }
    }

    errorSum = errorSum / (float)(v1->getNumVoxels());
    errorSum = sqrt(errorSum);
    errorSum = errorSum / (double)VolumeElement<T>::rangeMaxElement();

    return (float)errorSum;
}

typedef UniversalBinaryVolumeOperatorGeneric<VolumeOperatorCalcErrorBase> VolumeOperatorCalcError;

} // namespace

#endif // VRN_VOLUMEOPERATOR_H
