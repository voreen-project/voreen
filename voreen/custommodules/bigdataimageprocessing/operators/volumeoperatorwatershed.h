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

#ifndef VRN_VOLUMEOPERATORWATERSHEDTRANSFORM_H
#define VRN_VOLUMEOPERATORWATERSHEDTRANSFORM_H

#include "ternaryvolumeoperator.h"
#include <queue>

namespace Watershed {

template<typename T>
struct Voxel{
    tgt::svec3 position;
    T label;
    T height;
};

template<typename T>
bool VoxelCompare(const Voxel<T>& v1, const Voxel<T>& v2){
    return (v1.height < v2.height);
}

//To resolve VS bug with decltype()
template <typename T>
T identity(T);

};

#define VRN_FOR_EACH_VOXEL_INT(INDEX, START, END)\
    for (tgt::ivec3 (INDEX) = (START); (INDEX).z < (END).z; ++(INDEX).z)\
        for ((INDEX).y = (START).y; (INDEX).y < (END).y; ++(INDEX).y)\
            for ((INDEX).x = (START).x; (INDEX).x < (END).x; ++(INDEX).x)

namespace voreen {

// Generic implementation Squared Euclidean distance transform:
class VRN_CORE_API VolumeOperatorWatershedTransformBase : public TernaryVolumeOperatorBase {
public:
    virtual Volume* apply(const VolumeBase* image, const VolumeBase* marker, const VolumeBase* mask, ProgressReporter* progressReporter = 0) const = 0;
};

//specific implementation squared euclidean distance transform
template<typename T>
class VolumeOperatorWatershedTransformGeneric : public VolumeOperatorWatershedTransformBase {
public:
    virtual Volume* apply(const VolumeBase* image, const VolumeBase* marker, const VolumeBase* mask, ProgressReporter* progressReporter = 0) const;
    //Implement isCompatible using a handy macro:

    virtual bool isCompatible(const VolumeBase* volume1, const VolumeBase* volume2, const VolumeBase* volume3) const{
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

        const VolumeRAM* v3 = volume3->getRepresentation<VolumeRAM>();
        if(!v2)
            return false;
        const VolumeAtomic<T>* va3 = dynamic_cast<const VolumeAtomic<T>*>(v3);
        if(!va3)
            return false;

        return true;
    }
};

template<typename T>
Volume* VolumeOperatorWatershedTransformGeneric<T>::apply(const VolumeBase* image, const VolumeBase* marker, const VolumeBase* mask, ProgressReporter* pR) const{
    const VolumeRAM* v = image->getRepresentation<VolumeRAM>();
    if (!v){
        return 0;
    }

    const VolumeAtomic<T>* vaImage = dynamic_cast<const VolumeAtomic<T>*>(v);
    if (!vaImage){
        return 0;
    }

    v = marker->getRepresentation<VolumeRAM>();
    if(!v){
        return 0;
    }

    const VolumeAtomic<T>* vaMarker = dynamic_cast<const VolumeAtomic<T>*>(v);
    if(!vaMarker){
        return 0;
    }

    v = mask->getRepresentation<VolumeRAM>();
    if(!v){
        return 0;
    }

    const VolumeAtomic<T>* vaMask = dynamic_cast<const VolumeAtomic<T>*>(v);
    if(!vaMask){
        return 0;
    }

    if(vaImage->getDimensions() != vaMask->getDimensions() || vaMask->getDimensions() != vaMarker->getDimensions()){
        return 0;
    }

    tgt::svec3 dim = vaImage->getDimensions();

    VolumeAtomic<T>* output = new VolumeAtomic<T>(dim, true);
    memset(output->getData(), 0, output->getNumBytes());
    std::priority_queue<Watershed::Voxel<T>, std::vector<Watershed::Voxel<T> >, decltype(Watershed::identity(&Watershed::VoxelCompare<T>))> WatershedHeap(&Watershed::VoxelCompare<T>);

    VRN_FOR_EACH_VOXEL(idx, tgt::svec3(0,0,0), dim){
        T label = vaMarker->voxel(idx);
        if(vaMask->voxel(idx) && label){
            Watershed::Voxel<T> v = {};
            v.position = idx;
            v.label = label;
            v.height = vaImage->voxel(idx);
            WatershedHeap.push(v);

            output->voxel(idx) = label;
        }
    }

    while(!WatershedHeap.empty()){
        Watershed::Voxel<T> v = WatershedHeap.top();
        WatershedHeap.pop();

        output->voxel(v.position) = v.label;

        VRN_FOR_EACH_VOXEL_INT(sub_idx, tgt::ivec3(-1), tgt::ivec3(2)){
            long long x = v.position.x + sub_idx.x;
            long long y = v.position.y + sub_idx.y;
            long long z = v.position.z + sub_idx.z;

            tgt::svec3 index = tgt::clamp(tgt::svec3((size_t)x,(size_t)y,(size_t)z), tgt::svec3::zero, dim - tgt::svec3::one);
            if(output->voxel(index) == 0 && vaMask->voxel(index)){
                output->voxel(index) = v.label;
                Watershed::Voxel<T> u = {};
                u.position = index;
                u.label = v.label;
                u.height = vaImage->voxel(index);
                WatershedHeap.push(u);
            }
        }
    }


    return new Volume(output, image);
}


typedef UniversalTernaryVolumeOperatorGeneric<VolumeOperatorWatershedTransformBase> VolumeOperatorWatershedTransform;

} // namespace


#endif
