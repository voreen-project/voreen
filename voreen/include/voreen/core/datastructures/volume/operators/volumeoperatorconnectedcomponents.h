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

#ifndef VRN_VOLUMEOPERATORCONNECTEDCOMPONENTANALYSIS_H
#define VRN_VOLUMEOPERATORCONNECTEDCOMPONENTANALYSIS_H

#include "voreen/core/datastructures/volume/volumeoperator.h"

namespace voreen {

/// Computes connected component analysis and returns a volume (of type uint16) containing the component IDs at each voxel
class VRN_CORE_API VolumeOperatorConnectedComponentAnalysisBase : public UnaryVolumeOperatorBase {
public:
    virtual Volume* apply(const VolumeBase* volume, int connectivity, bool stretchLabels, ProgressReporter* progressReporter = 0) const = 0;
};

template<typename T>
class VolumeOperatorConnectedComponentAnalysisGeneric : public VolumeOperatorConnectedComponentAnalysisBase {
public:
    virtual Volume* apply(const VolumeBase* volume, int connectivity, bool stretchLabels, ProgressReporter* progressReporter = 0) const;
    //Implement isCompatible using a handy macro:
    IS_COMPATIBLE
};

template<typename T>
Volume* VolumeOperatorConnectedComponentAnalysisGeneric<T>::apply(const VolumeBase* vb, int connectivity, bool stretchLabels, ProgressReporter* pR) const{
    const VolumeRAM* v = vb->getRepresentation<VolumeRAM>();
    if (!v)
        return 0;

    const VolumeAtomic<T>* va = dynamic_cast<const VolumeAtomic<T>*>(v);
    if (!va)
        return 0;

    tgt::svec3 dim = va->getDimensions();

    VolumeAtomic<uint16_t>* marked = new VolumeAtomic<uint16_t>(dim, true);
    memset(marked->getData(), 0, marked->getNumBytes());

    std::vector<tgt::svec3> stack;
    stack.reserve(va->getNumVoxels());

    uint16_t label = 1;
    VRN_FOR_EACH_VOXEL(idx, tgt::svec3::zero, dim){
        if(marked->voxel(idx) == 0 && va->voxel(idx) != T(0)){
            stack.push_back(idx);
            marked->voxel(idx) = label;
            while(!stack.empty()){
                tgt::svec3 back = stack.back();
                stack.pop_back();


                for (tgt::ivec3 sub_idx = tgt::ivec3(-1,-1,-1); sub_idx.z < 2; ++sub_idx.z){
                    for (sub_idx.y = -1; sub_idx.y < 2; ++sub_idx.y){
                        for (sub_idx.x = -1; sub_idx.x < 2; ++sub_idx.x){
                            long long x = back.x + sub_idx.x;
                            long long y = back.y + sub_idx.y;
                            long long z = back.z + sub_idx.z;

                            tgt::svec3 index = tgt::clamp(tgt::svec3((size_t)x,(size_t)y,(size_t)z), tgt::svec3::zero, dim - tgt::svec3::one);
                            if(sub_idx != tgt::ivec3::zero && marked->voxel(index) == 0 && va->voxel(index) != T(0)){
                                stack.push_back(index);
                                marked->voxel(index) = label;
                            }
                        }
                    }
                }
            }
            label++;
        }
    }

    label--;
    if(stretchLabels){
        VRN_FOR_EACH_VOXEL(idx, tgt::svec3::zero, dim){
            marked->voxel(idx) = marked->voxel(idx)*static_cast<uint16_t>(65535.f/static_cast<float>(label));
        }
    }

    if (pR)
        pR->setProgress(1.f);

    return new Volume(marked, vb);
}


typedef UniversalUnaryVolumeOperatorGeneric<VolumeOperatorConnectedComponentAnalysisBase> VolumeOperatorConnectedComponentAnalysis;

} // namespace


#endif
