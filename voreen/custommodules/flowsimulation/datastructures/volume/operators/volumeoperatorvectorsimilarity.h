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

#ifndef VRN_VOLUMEOPERATORVECTORSIMILARITY_H
#define VRN_VOLUMEOPERATORVECTORSIMILARITY_H

#include "voreen/core/datastructures/volume/volumeoperator.h"
#include "voreen/core/datastructures/volume/volumefactory.h"
#include "voreen/core/datastructures/volume/volumeminmaxmagnitude.h"
#include "voreen/core/utils/stringutils.h"

namespace voreen {

// Base class, defines interface for the operator (-> apply):
class VRN_CORE_API VolumeOperatorVectorSimilarityBase : public UnaryVolumeOperatorBase {
public:
    virtual Volume* apply(const VolumeBase* volume, float p, ProgressReporter* progressReporter = 0) const = 0;
};

// Generic implementation:
template<typename T>
class VolumeOperatorVectorSimilarityGeneric : public VolumeOperatorVectorSimilarityBase {
public:
    virtual Volume* apply(const VolumeBase* volume, float p, ProgressReporter* progressReporter = 0) const;
    //Implement isCompatible using a handy macro:
    IS_COMPATIBLE
};

template<typename T>
Volume* VolumeOperatorVectorSimilarityGeneric<T>::apply(const VolumeBase* vh, float p, ProgressReporter* progressReporter) const {

    VolumeRAM* out = VolumeFactory().create(getBaseTypeFromType<T>(), vh->getDimensions());
    if (!out)
        return nullptr;

    const VolumeRAM* volume = vh->getRepresentation<VolumeRAM>();
    if(!volume)
        return nullptr;

    RealWorldMapping rwm = vh->getRealWorldMapping();
    float maxMagnitude = vh->getDerivedData<VolumeMinMaxMagnitude>()->getMaxMagnitude();

    VRN_FOR_EACH_VOXEL_WITH_PROGRESS(index, tgt::svec3::zero, out->getDimensions(), progressReporter) {

        T voxel;
        for(size_t channel=0; channel<vh->getNumChannels(); channel++) {
            voxel[channel] = rwm.normalizedToRealWorld(volume->getVoxelNormalized(index, channel));
        }

        //float angle = std::acos(tgt::dot(xAxis, voxel) / (tgt::length(xAxis) * tgt::length(voxel) )
        float angle = std::acos(voxel.x / tgt::length(voxel)) / (2*tgt::PIf);
        float magnitude = tgt::length(voxel) / maxMagnitude;

        float value = (1 - p) * angle  + p * magnitude;
        out->setVoxelNormalized(value, index);
    }
    if (progressReporter)
        progressReporter->setProgress(1.f);

    return new Volume(out, vh);
}

typedef UniversalUnaryVolumeOperatorGeneric<VolumeOperatorVectorSimilarityBase> VolumeOperatorVectorSimilarity;

} // namespace

#endif // VRN_VOLUMEOPERATORVECTORSIMILARITY_H
