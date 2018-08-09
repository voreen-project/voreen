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

#ifndef VRN_VOLUMEOPERATORHALFSAMPLE_H
#define VRN_VOLUMEOPERATORHALFSAMPLE_H

#include "voreen/core/datastructures/volume/volumeoperator.h"

namespace voreen {

/**
 * Reduces the Volumes resolution by half, by linearly downsampling 8 voxels
 * to 1 voxel. This does not necessarily happen when using the resample(..) function.
 *
 * @return the resampled volume
 */
class VRN_CORE_API VolumeOperatorHalfsampleBase : public UnaryVolumeOperatorBase {
public:
    virtual Volume* apply(const VolumeBase* volume, ProgressReporter* progressReporter = 0) const = 0;
};

// Generic implementation:
template<typename T>
class VolumeOperatorHalfsampleGeneric : public VolumeOperatorHalfsampleBase {
public:
    virtual Volume* apply(const VolumeBase* volume, ProgressReporter* progressReporter = 0) const;
    //Implement isCompatible using a handy macro:
    IS_COMPATIBLE
};

template<typename T>
Volume* VolumeOperatorHalfsampleGeneric<T>::apply(const VolumeBase* vh, ProgressReporter* progressReporter) const {
    const VolumeRAM* v = vh->getRepresentation<VolumeRAM>();
    if(!v)
        return 0;

    const VolumeAtomic<T>* volume = dynamic_cast<const VolumeAtomic<T>*>(v);
    if(!volume)
        return 0;

    tgt::svec3 halfDims = volume->getDimensions() / tgt::svec3(2);

    VolumeAtomic<T>* newVolume = new VolumeAtomic<T>(halfDims);

    typedef typename VolumeElement<T>::DoubleType Double;
    VRN_FOR_EACH_VOXEL_WITH_PROGRESS(index, tgt::svec3(0, 0, 0), halfDims, progressReporter) {
        tgt::svec3 pos = index*tgt::svec3(2); // tgt::ivec3(2*x,2*y,2*z);
        newVolume->voxel(index) =
            T(  Double(volume->voxel(pos.x, pos.y, pos.z))          * (1.0/8.0) //LLF
              + Double(volume->voxel(pos.x, pos.y, pos.z+1))        * (1.0/8.0) //LLB
              + Double(volume->voxel(pos.x, pos.y+1, pos.z))        * (1.0/8.0) //ULF
              + Double(volume->voxel(pos.x, pos.y+1, pos.z+1))      * (1.0/8.0) //ULB
              + Double(volume->voxel(pos.x+1, pos.y, pos.z))        * (1.0/8.0) //LRF
              + Double(volume->voxel(pos.x+1, pos.y, pos.z+1))      * (1.0/8.0) //LRB
              + Double(volume->voxel(pos.x+1, pos.y+1, pos.z))      * (1.0/8.0) //URF
              + Double(volume->voxel(pos.x+1, pos.y+1, pos.z+1))    * (1.0/8.0)); //URB
    }
    if (progressReporter)
        progressReporter->setProgress(1.f);

    Volume* ret = new Volume(newVolume, vh);
    ret->setSpacing(vh->getSpacing()*2.f);
    return ret;
}

typedef UniversalUnaryVolumeOperatorGeneric<VolumeOperatorHalfsampleBase> VolumeOperatorHalfsample;

} // namespace

#endif // VRN_VOLUMEOPERATOR_H
