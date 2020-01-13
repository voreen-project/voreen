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

#ifndef VRN_VOLUMEOPERATORRESIZE_H
#define VRN_VOLUMEOPERATORRESIZE_H

#include "voreen/core/datastructures/volume/volumeoperator.h"

namespace voreen {

/**
 * Returns a resized copy of the passed input volume by keeping its
 * remaining properties.
 *
 * @note The volume data is not copied by this operation.
 *
 * @see VolumeOperatorResample
 *
 * @return the resized volume
 */
class VRN_CORE_API VolumeOperatorResizeBase : public UnaryVolumeOperatorBase {
public:
    /**
     * @param newDims The target dimensions
     * @param allocMem If true, a new data buffer is allocated
     */
    virtual Volume* apply(const VolumeBase* volume, tgt::ivec3 newDims, bool allocMem = true) const = 0;
};

// Generic implementation:
template<typename T>
class VolumeOperatorResizeGeneric : public VolumeOperatorResizeBase {
public:
    virtual Volume* apply(const VolumeBase* volume, tgt::ivec3 newDims, bool allocMem = true) const;
    //Implement isCompatible using a handy macro:
    IS_COMPATIBLE
};

template<typename T>
Volume* VolumeOperatorResizeGeneric<T>::apply(const VolumeBase* vh, tgt::ivec3 newDims, bool allocMem) const {
    const VolumeRAM* vol = vh->getRepresentation<VolumeRAM>();
    if(!vol)
        return 0;

    const VolumeAtomic<T>* volume = dynamic_cast<const VolumeAtomic<T>*>(vol);
    if(!volume)
        return 0;

    LINFOC("voreen.VolumeOperatorResize", "Resizing from dimensions " << volume->getDimensions() << " to " << newDims);

    // build target volume
    VolumeAtomic<T>* result;
    try {
        result = new VolumeAtomic<T>(newDims, allocMem);
    }
    catch (std::bad_alloc&) {
        throw; // throw it to the caller
    }

    return new Volume(result, vh);
}

typedef UniversalUnaryVolumeOperatorGeneric<VolumeOperatorResizeBase> VolumeOperatorResize;

} // namespace

#endif // VRN_VOLUMEOPERATOR_H
