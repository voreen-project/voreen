/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#ifndef VRN_VOLUMEOPERATORRESIZEPOWEROFTWO_H
#define VRN_VOLUMEOPERATORRESIZEPOWEROFTWO_H

#include "voreen/core/datastructures/volume/operators/volumeoperatorresize.h"

#include "tgt/tgt_math.h"

namespace voreen {

/**
 * Returns a copy of the passed input volume that has been
 * resized to the next largest power-of-two.
 */
class VRN_CORE_API VolumeOperatorResizePowerOfTwoBase : public UnaryVolumeOperatorBase {
public:
    /**
     * @param volume the input volume
     * @param isotropic if true, the volume is resized to a cube. Otherwise, the next power-of-two dimensions
     *  are determined component-wise.
     * @param if true, the volume data of the source volume is copied over to the returned volume
     */
    virtual Volume* apply(const VolumeBase* volume, bool isotropic = false, bool copyVolumeData = true) const = 0;
};

// Generic implementation:
template<typename T>
class VolumeOperatorResizePowerOfTwoGeneric : public VolumeOperatorResizePowerOfTwoBase {
public:
    virtual Volume* apply(const VolumeBase* volume, bool isotropic = false, bool copyVolumeData = true) const;
    //Implement isCompatible using a handy macro:
    IS_COMPATIBLE
};

template<typename T>
Volume* VolumeOperatorResizePowerOfTwoGeneric<T>::apply(const VolumeBase* vh, bool isotropic, bool copyVolumeData) const {
    // extract VolumeAtomic from passed input volume
    const VolumeRAM* vol = vh->getRepresentation<VolumeRAM>();
    if (!vol) {
        LWARNINGC("voreen.VolumeOperatorResizePowerOfTwo", "No RAM volume available");
        return 0;
    }
    const VolumeAtomic<T>* volume = dynamic_cast<const VolumeAtomic<T>*>(vol);
    if (!volume) {
        LWARNINGC("voreen.VolumeOperatorResizePowerOfTwo", "Failed to cast VolumeRAM into VolumeAtomic.");
        return 0;
    }
    tgt::svec3 srcDim = volume->getDimensions();

    // determine power-of-two dimensions
    tgt::svec3 nextPotDim;
    if (isotropic) {
        // next power-of-two of largest dimension component
        size_t maxDim = tgt::max(srcDim);
        size_t nextPot = tgt::nextLargerPowerOfTwo((int)maxDim);
        nextPotDim = tgt::svec3(nextPot);
    }
    else {
        // component-wise next power-of-two
        nextPotDim = tgt::svec3(tgt::nextLargerPowerOfTwo((int)srcDim.x), tgt::nextLargerPowerOfTwo((int)srcDim.y), tgt::nextLargerPowerOfTwo((int)srcDim.z));
    }
    tgtAssert(tgt::hand(tgt::greaterThanEqual(nextPotDim, srcDim)), "invalid next-power-of-two");
    LINFOC("voreen.VolumeOperatorResizePowerOfTwo", "Resizing from " << srcDim << " to " << nextPotDim);

    // created resized volume
    Volume* result = VolumeOperatorResizeGeneric<T>().apply(vh, tgt::svec3(nextPotDim), true);
    VolumeAtomic<T>* resultAtomic = dynamic_cast<VolumeAtomic<T>* >(result->getWritableRepresentation<VolumeRAM>());
    if (!resultAtomic) {
        LERRORC("voreen.VolumeOperatorResizePowerOfTwo", "VolumeOperatorResize returned unexpected volume type");
        delete result;
        return 0;
    }
    resultAtomic->clear();

    // copy volume data (centered)
    if (copyVolumeData) {
        tgt::svec3 offset = (nextPotDim - srcDim) / tgt::svec3(2, 2, 2);
        VRN_FOR_EACH_VOXEL(i, tgt::svec3(0, 0, 0), volume->getDimensions()) {
            resultAtomic->voxel(offset + i) = volume->voxel(i);
        }
    }

    return result;
}

typedef UniversalUnaryVolumeOperatorGeneric<VolumeOperatorResizePowerOfTwoBase> VolumeOperatorResizePowerOfTwo;

} // namespace

#endif // VRN_VOLUMEOPERATORRESIZEPOWEROFTWO_H
