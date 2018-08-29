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

#ifndef VRN_VOLUMEOPERATORGAUSSIAN_H
#define VRN_VOLUMEOPERATORGAUSSIAN_H

#include "voreen/core/datastructures/volume/volumeoperator.h"

namespace voreen {

// Base class, defines interface for the operator (-> apply):
class VRN_CORE_API VolumeOperatorGaussianBase : public UnaryVolumeOperatorBase {
public:
    virtual Volume* apply(const VolumeBase* volume, size_t kernelSize = 3, float sigma = 1.f, ProgressReporter* progressReporter = 0) const = 0;
};

// Generic implementation:
template<typename T>
class VolumeOperatorGaussianGeneric : public VolumeOperatorGaussianBase {
public:
    virtual Volume* apply(const VolumeBase* volume, size_t kernelSize = 3, float sigma = 1.f, ProgressReporter* progressReporter = 0) const;
    //Implement isCompatible using a handy macro:
    IS_COMPATIBLE
};

template<typename T>
Volume* VolumeOperatorGaussianGeneric<T>::apply(const VolumeBase* vh, size_t kernelSize, float sigma, ProgressReporter* progressReporter) const {
    const VolumeRAM* v = vh->getRepresentation<VolumeRAM>();
    if (!v)
        return 0;

    const VolumeAtomic<T>* va = dynamic_cast<const VolumeAtomic<T>*>(v);
    if (!va)
        return 0;

    // TODO: this needs two volumes for ping-pong, optionally use just one using a 3D filter kernel instead of a separated computation in each direction
    VolumeAtomic<T>* tmp = va->clone();
    VolumeAtomic<T>* output = va->clone();

    // compute sigma from kernel size
    size_t kernelRadius = kernelSize / 2;
    tgtAssert(kernelSize >= 3, "invalid kernel size");
    tgtAssert(kernelRadius > 0, std::string("invalid kernel radius: " + ftos(kernelRadius)).c_str());
    tgtAssert(sigma > 0, std::string("Invalid sigma: " + ftos(sigma)).c_str());

    float* gaussKernel = new float[kernelRadius + 1];

    // compute (half) 1D gauss kernel
    for (size_t i=0; i<=kernelRadius; i++)
        gaussKernel[i] = exp(-static_cast<float>(i*i)/(2.f*sigma*sigma));

    // compute normalization factor
    float norm = 0.0;
    for (size_t i=1; i<=kernelRadius; i++)
        norm += gaussKernel[i];
    // so far we have just computed norm for one half without x=0
    norm = 2.f * norm + gaussKernel[0];

    // use 1D filter kernel in each direction
    tgt::svec3 volDim = vh->getDimensions();
    const VolumeAtomic<T>* currentSource = va;
    VolumeAtomic<T>* currentDestination = output;
    for (size_t dim = 0; dim < 3; ++dim) {
        // volume ping-pong
        if (dim == 1) {
            currentSource = output;
            currentDestination = tmp;
        }
        else if (dim == 2) {
            currentSource = tmp;
            currentDestination = output;
        }

        VRN_FOR_EACH_VOXEL_WITH_PROGRESS(pos, tgt::ivec3(0), volDim, progressReporter) {
            // for (u = -kernelRadius, u <= kernelRadius) getVoxel, multiply by kernel[u], add to sum and write to a volume (switch input and output volume due to ping-pong)
            float sum = 0;
            for (int u = -static_cast<int>(kernelRadius); u <= static_cast<int>(kernelRadius); ++u) {
                // calculate position
                tgt::ivec3 currentVoxel = tgt::ivec3(pos);
                currentVoxel[dim] += u;
                currentVoxel = tgt::clamp(currentVoxel, tgt::ivec3::zero, tgt::ivec3(volDim - tgt::svec3::one));
                sum += gaussKernel[std::abs(u)] * static_cast<float>(currentSource->voxel(tgt::svec3(currentVoxel)));
            }
            // FIXME: should round for integer values
            currentDestination->voxel(pos) = static_cast<T>(sum / norm);
        }
    }

    delete[] gaussKernel;
    delete tmp;

    if (progressReporter)
        progressReporter->setProgress(1.f);

    return new Volume(output, vh);
}

typedef UniversalUnaryVolumeOperatorGeneric<VolumeOperatorGaussianBase> VolumeOperatorGaussian;

} // namespace

#endif
