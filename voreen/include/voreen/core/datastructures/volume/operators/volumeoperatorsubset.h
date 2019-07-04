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

#ifndef VRN_VOLUMEOPERATORSUBSET_H
#define VRN_VOLUMEOPERATORSUBSET_H

#include "voreen/core/datastructures/volume/volumeoperator.h"
#include "voreen/core/datastructures/volume/volumedisk.h"

namespace voreen {

/**
 * Returns a volume containing the subset [pos, pos+size[ of the passed input volume.
 */
class VRN_CORE_API VolumeOperatorSubsetBase : public UnaryVolumeOperatorBase {
public:
    virtual Volume* apply(const VolumeBase* volume, tgt::ivec3 pos, tgt::ivec3 size, ProgressReporter* progressReporter = 0) const = 0;
};

// Generic implementation:
template<typename T>
class VolumeOperatorSubsetGeneric : public VolumeOperatorSubsetBase {
public:
    virtual Volume* apply(const VolumeBase* volume, tgt::ivec3 pos, tgt::ivec3 size, ProgressReporter* progressReporter = 0) const;
    //Implement isCompatible using a handy macro:
    IS_COMPATIBLE
};

template<typename T>
Volume* VolumeOperatorSubsetGeneric<T>::apply(const VolumeBase* vh, tgt::ivec3 pos, tgt::ivec3 size, ProgressReporter* progressReporter) const {
    const VolumeRAM* volRam = 0;
    const VolumeDisk* volDisk = 0;

    if (vh->hasRepresentation<VolumeRAM>())
        volRam = vh->getRepresentation<VolumeRAM>();
    else if (vh->hasRepresentation<VolumeDisk>())
        volDisk = vh->getRepresentation<VolumeDisk>();

    if(!volRam && !volDisk) {
        LERRORC("voreen.VolumeOperatorCreateSubset", "No VolumeRam or VolumeDisk present!");
        return 0;
    }

    VolumeRAM* subsetRAM = 0;
    tgt::svec3 start = tgt::svec3(tgt::max(pos, tgt::ivec3::zero));// clamp values
    tgt::svec3 end   = tgt::svec3(tgt::min(pos + size, tgt::ivec3(vh->getDimensions())));    // clamp values
    tgt::svec3 diff  = end - start;

    if(volRam) {
        const VolumeAtomic<T>* volume = dynamic_cast<const VolumeAtomic<T>*>(volRam);
        if(!volume)
            return 0;

        VolumeAtomic<T>* subset;
        try {
            subset = new VolumeAtomic<T>(size);
        }
        catch (std::bad_alloc&) {
            LERRORC("voreen.VolumeOperatorCreateSubset", "Failed to create subset: bad allocation");
            throw; // throw it to the caller
        }

        LINFOC("voreen.VolumeOperatorCreateSubset", "Creating subset " << size << " from position " << pos);

        // create values for ranges less than zero and greater equal dimensions_
        subset->clear(); // TODO: This can be optimized by avoiding to clear the values in range

        VRN_FOR_EACH_VOXEL_WITH_PROGRESS(index, tgt::svec3(0, 0, 0), diff, progressReporter)
            subset->voxel(index) = volume->voxel(index+start);

        subsetRAM = subset;
    } else { //case volDisk
        subsetRAM = volDisk->loadBrick(start,diff);
    }

    Volume* newvh = new Volume(subsetRAM, vh);
    newvh->setOffset(vh->getOffset() + (tgt::vec3(start) * vh->getSpacing()));
    return newvh;
}

typedef UniversalUnaryVolumeOperatorGeneric<VolumeOperatorSubsetBase> VolumeOperatorSubset;

} // namespace

#endif // VRN_VOLUMEOPERATOR_H
