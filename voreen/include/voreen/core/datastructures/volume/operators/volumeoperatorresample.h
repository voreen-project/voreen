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

#ifndef VRN_VOLUMEOPERATORRESAMPLE_H
#define VRN_VOLUMEOPERATORRESAMPLE_H

#include "voreen/core/datastructures/volume/volumeoperator.h"

namespace voreen {

/**
 * Returns a copy of the input volume that has been resampled to the specified dimensions
 * by using the given filtering mode.
 *
 * @return the resampled volume
 */
class VRN_CORE_API VolumeOperatorResampleBase : public UnaryVolumeOperatorBase {
public:
    /**
     * @param newDims the target dimensions
     * @param filter The filtering mode to use for calculating the resampled values.
     */
    virtual Volume* apply(const VolumeBase* volume, tgt::ivec3 newDims, VolumeRAM::Filter filter, ProgressReporter* progressReporter = 0) const = 0;
};

// Generic implementation:
template<typename T>
class VolumeOperatorResampleGeneric : public VolumeOperatorResampleBase {
public:
    virtual Volume* apply(const VolumeBase* volume, tgt::ivec3 newDims, VolumeRAM::Filter filter, ProgressReporter* progressReporter = 0) const;
    //Implement isCompatible using a handy macro:
    IS_COMPATIBLE
};

template<typename T>
Volume* VolumeOperatorResampleGeneric<T>::apply(const VolumeBase* vh, tgt::ivec3 newDims, VolumeRAM::Filter filter, ProgressReporter* progressReporter) const {
    const VolumeRAM* vol = vh->getRepresentation<VolumeRAM>();
    if(!vol)
        return 0;

    const VolumeAtomic<T>* volume = dynamic_cast<const VolumeAtomic<T>*>(vol);
    if(!volume)
        return 0;

    using tgt::vec3;
    using tgt::ivec3;
    using tgt::svec3;

    LDEBUGC("voreen.VolumeOperatorResample", "Resampling from dimensions " << volume->getDimensions() << " to " << newDims);

    vec3 ratio = vec3(volume->getDimensions()) / vec3(newDims);
    // here is the actually correct rescaling formula for geting the position in the old volume:
    // (pos - (dim_new - 1) / 2) * ratio + (dim_old - 1) / 2
    // define d_a := (dim_new - 1) / 2, d_b := (dim_old - 1) / 2
    vec3 d_a = vec3(newDims - tgt::ivec3::one) / vec3(2.f);
    vec3 d_b = vec3(volume->getDimensions() - tgt::svec3::one) / vec3(2.f);

    ivec3 pos = ivec3::zero; // iteration variable
    vec3 nearest; // knows the new position of the target volume

    // build target volume
    VolumeAtomic<T>* v;
    try {
         v = new VolumeAtomic<T>(newDims);
    }
    catch (std::bad_alloc&) {
        throw; // throw it to the caller
    }

    if (progressReporter)
        progressReporter->setProgress(0.f);

    /*
        Filter from the source volume to the target volume.
    */
    for (pos.z = 0; pos.z < newDims.z; ++pos.z) {

        if (progressReporter)
            progressReporter->setProgress(static_cast<float>(pos.z) / static_cast<float>(newDims.z));

        nearest.z = (static_cast<float>(pos.z) - d_a.z) * ratio.z + d_b.z;

        for (pos.y = 0; pos.y < newDims.y; ++pos.y) {
            nearest.y = (static_cast<float>(pos.y) - d_a.y) * ratio.y + d_b.y;

            for (pos.x = 0; pos.x < newDims.x; ++pos.x) {
                nearest.x = (static_cast<float>(pos.x) - d_a.x) * ratio.x + d_b.x;

                // switch between filtering options
                switch (filter) {
                    case VolumeRAM::NEAREST: {
                        svec3 index = tgt::clamp(svec3(tgt::iround(nearest)), svec3(0, 0, 0), volume->getDimensions() - svec3(1, 1, 1));
                        v->voxel(pos) = volume->voxel(index); // round and do the lookup
                        break; //switch-case
                    }
                    case VolumeRAM::LINEAR: {
                        // clamp to volume dimensions
                        nearest = tgt::clamp(nearest, tgt::vec3::zero, tgt::vec3(volume->getDimensions() - tgt::svec3::one));

                        // get decimal part and lower / upper voxel
                        vec3 p = nearest - floor(nearest);
                        ivec3 llb = ivec3(nearest);
                        ivec3 urf = ivec3(ceil(nearest));

                        // clamp again for safety so the lookups do not exceed the dimensions
                        llb = tgt::max(llb, tgt::ivec3::zero);
                        urf = tgt::min(urf, ivec3(volume->getDimensions()) - 1);

                        /*
                          interpolate linearly
                        */
                        typedef typename VolumeElement<T>::DoubleType Double;
                        v->voxel(pos) =
                            T(  Double(volume->voxel(llb.x, llb.y, llb.z)) * static_cast<double>((1.f-p.x)*(1.f-p.y)*(1.f-p.z))  // llB
                              + Double(volume->voxel(urf.x, llb.y, llb.z)) * static_cast<double>((    p.x)*(1.f-p.y)*(1.f-p.z))  // lrB
                              + Double(volume->voxel(urf.x, urf.y, llb.z)) * static_cast<double>((    p.x)*(    p.y)*(1.f-p.z))  // urB
                              + Double(volume->voxel(llb.x, urf.y, llb.z)) * static_cast<double>((1.f-p.x)*(    p.y)*(1.f-p.z))  // ulB
                              + Double(volume->voxel(llb.x, llb.y, urf.z)) * static_cast<double>((1.f-p.x)*(1.f-p.y)*(    p.z))  // llF
                              + Double(volume->voxel(urf.x, llb.y, urf.z)) * static_cast<double>((    p.x)*(1.f-p.y)*(    p.z))  // lrF
                              + Double(volume->voxel(urf.x, urf.y, urf.z)) * static_cast<double>((    p.x)*(    p.y)*(    p.z))  // urF
                              + Double(volume->voxel(llb.x, urf.y, urf.z)) * static_cast<double>((1.f-p.x)*(    p.y)*(    p.z)));// ulF

                        break; // switch-case
                    }
                    case VolumeRAM::CUBIC: {
                        // clamp to volume dimensions
                        nearest = tgt::clamp(nearest, tgt::vec3::zero, tgt::vec3(volume->getDimensions() - tgt::svec3::one));
                        v->setVoxelNormalized(volume->getVoxelNormalizedCubic(nearest), pos);
                        break; // switch-case
                    }
                }
            }
        }
    }

    if (progressReporter)
        progressReporter->setProgress(1.f);

    Volume* h = new Volume(v, vh);
    h->setSpacing(vh->getSpacing() * ratio);
    return h;
}

typedef UniversalUnaryVolumeOperatorGeneric<VolumeOperatorResampleBase> VolumeOperatorResample;

} // namespace

#endif // VRN_VOLUMEOPERATOR_H
