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

#ifndef VRN_VOLUMEOPERATORFASTVOLUMECOMBINE_H
#define VRN_VOLUMEOPERATORFASTVOLUMECOMBINE_H

#include "voreen/core/datastructures/volume/volumeoperator.h"
#include <queue>

namespace voreen {

/// Voxel-wise combine operation.
enum CombineOperation {
    OP_ADD,
    OP_A_MINUS_B,
    OP_B_MINUS_A,
    OP_MULT,
    OP_AVG,
    OP_MAX,
    OP_MIN,
    OP_WEIGHTED_SUM,
    OP_WEIGHTED_SUM_2P,
    OP_BLEND,
    OP_MASK_A_BY_B,
    OP_MASK_B_BY_A,
    OP_PRIORITY_FIRST,
    OP_PRIORITY_SECOND,
    OP_TAKE_FIRST,
    OP_TAKE_SECOND
};

// Generic implementation Squared Euclidean distance transform:
class VRN_CORE_API VolumeOperatorFastVolumeCombineBase : public BinaryVolumeOperatorBase {
public:
    virtual Volume* apply(const VolumeBase* volume1, const VolumeBase* volume2, CombineOperation co, float FactorC, float FactorD, bool nearest, bool linear, bool cubic, ProgressReporter* pR = 0) const = 0;
};

//specific implementation squared euclidean distance transform
template<typename T>
class VolumeOperatorFastVolumeCombineGeneric : public VolumeOperatorFastVolumeCombineBase {
public:
    virtual Volume* apply(const VolumeBase* volume1, const VolumeBase* volume2, CombineOperation co, float FactorC, float FactorD, bool nearest, bool linear, bool cubic,  ProgressReporter* pR = 0) const;
    //Implement isCompatible using a handy macro:

    virtual bool isCompatible(const VolumeBase* volume1, const VolumeBase* volume2) const{
        const VolumeRAM* v1 = volume1->getRepresentation<VolumeRAM>();
        if(!v1)
            return false;
        const VolumeAtomic<T>* va1 = dynamic_cast<const VolumeAtomic<T>*>(v1);
        if(!va1)
            return false;

        const VolumeRAM* v2 = volume2->getRepresentation<VolumeRAM>();
        if(!v2)
            return false;
        const VolumeAtomic<T>* va2 = dynamic_cast<const VolumeAtomic<T>*>(v2);
        if(!va2)
            return false;

        return true;
    }
};

template<typename T>
Volume* VolumeOperatorFastVolumeCombineGeneric<T>::apply(const VolumeBase* vol1, const VolumeBase* vol2, CombineOperation co, float factorC, float factorD, bool nearest, bool linear, bool cubic, ProgressReporter* pR) const{
    const VolumeRAM* v = vol1->getRepresentation<VolumeRAM>();
    if (!v)
        return 0;

    const VolumeAtomic<T>* va = dynamic_cast<const VolumeAtomic<T>*>(v);
    if (!va)
        return 0;

    VolumeAtomic<T>* combinedVolume = va->clone();

    // determine which volume is the other (non reference) one
    const VolumeBase* otherVolume = vol2;

    //
    // voxel-wise combination
    //
    tgt::vec3 dimFirst(vol1->getDimensions() - tgt::svec3(1));
    tgt::vec3 dimSecond(vol2->getDimensions() - tgt::svec3(1));
    const float c = factorC;
    const float d = factorD;
    const bool nearestFiltering = nearest;
    const bool linearFiltering = linear;
    const bool cubicFiltering = cubic;

    const VolumeRAM* v1 = vol1->getRepresentation<VolumeRAM>();
    const VolumeRAM* v2 = vol2->getRepresentation<VolumeRAM>();

    VRN_FOR_EACH_VOXEL_WITH_PROGRESS(pos, tgt::ivec3(0), combinedVolume->getDimensions(), pR) {
        for (size_t referenceChannel = 0; referenceChannel < combinedVolume->getNumChannels(); ++referenceChannel) {
            size_t otherChannel = (combinedVolume->getNumChannels() == otherVolume->getNumChannels()) ? referenceChannel : 0;

            // sample input volumes, if transformed voxel position lies inside respective volume
            float valFirst = 0.f;
            float valSecond = 0.f;
            if (nearestFiltering) {
                valFirst = v1->getVoxelNormalized(pos, referenceChannel);
                valSecond = v2->getVoxelNormalized(pos, otherChannel);
            }
            else if (linearFiltering) {
                valFirst = v1->getVoxelNormalizedLinear(pos, referenceChannel);
                valSecond = v2->getVoxelNormalizedLinear(pos, otherChannel);

            }
            else if (cubicFiltering) {
                valFirst = v1->getVoxelNormalizedCubic(pos, referenceChannel);
                valSecond = v2->getVoxelNormalizedCubic(pos, otherChannel);
            }
            else {
                return nullptr;
            }

            // apply operation to sampled values (note: a switch-block within a volume traversal loop
            // should normally be avoided, however in this case the main operations are way more expensive)
            float result = 0.f;
            switch (co) {
                case OP_MAX:
                    result = std::max(valFirst, valSecond);
                    break;
                case OP_MIN:
                    result = std::min(valFirst, valSecond);
                    break;
                case OP_ADD:
                    result = valFirst + valSecond;
                    break;
                case OP_A_MINUS_B:
                    result = valFirst - valSecond;
                    break;
                case OP_B_MINUS_A:
                    result = valSecond - valFirst;
                    break;
                case OP_MULT:
                    result = valFirst * valSecond;
                    break;
                case OP_AVG:
                    result = (valFirst + valSecond) / 2.f;
                    break;
                case OP_WEIGHTED_SUM:
                    result = c*valFirst + (1.f-c)*valSecond;
                    break;
                case OP_WEIGHTED_SUM_2P:
                    result = c*valFirst + d*valSecond;
                    break;
                case OP_BLEND:
                    result = valFirst + valSecond*(1.f - valFirst);
                    break;
                case OP_MASK_A_BY_B:
                    result = (valSecond > 0.f) ? valFirst : 0.f;
                    break;
                case OP_MASK_B_BY_A:
                    result = (valFirst > 0.f) ? valSecond : 0.f;
                    break;
                case OP_PRIORITY_FIRST:
                    result = (valFirst > 0.f) ? valFirst : valSecond;
                    break;
                case OP_PRIORITY_SECOND:
                    result = (valSecond > 0.f) ? valSecond : valFirst;
                    break;
                case OP_TAKE_FIRST:
                    result = valFirst;
                    break;
                case OP_TAKE_SECOND:
                    result = valSecond;
                    break;
                default:
                    return nullptr;
            }

            // assign clamped result to combined volume
            combinedVolume->setVoxelNormalized(tgt::clamp(result, 0.f, 1.f), pos, referenceChannel);
        }
    } // VRN_FOR_EACH_VOXEL_WITH_PROGRESS

    if(pR){
        pR->setProgress(1.0f);
    }
    return new Volume(combinedVolume, vol1);
}


typedef UniversalBinaryVolumeOperatorGeneric<VolumeOperatorFastVolumeCombineBase> VolumeOperatorFastVolumeCombine;

} // namespace


#endif
