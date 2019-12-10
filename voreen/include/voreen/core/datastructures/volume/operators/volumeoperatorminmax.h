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

#ifndef VRN_VOLUMEOPERATORMINMAX_H
#define VRN_VOLUMEOPERATORMINMAX_H

#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/volume.h"

namespace voreen {


///Returns the minimum voxel value in the volume.
class VRN_CORE_API VolumeOperatorMinValue {
public:
    template<typename T>
    static T apply(const VolumeAtomic<T>* volume);

    template<typename S>
    static S apply(const VolumeAtomic<tgt::Vector2<S> >* volume);

    template<typename S>
    static S apply(const VolumeAtomic<tgt::Vector3<S> >* volume);

    template<typename S>
    static S apply(const VolumeAtomic<tgt::Vector4<S> >* volume);

    template<typename S>
    static S apply(const VolumeAtomic<Tensor2<S> >* volume);
};

///Returns the maximum voxel value in the volume.
class VRN_CORE_API VolumeOperatorMaxValue {
public:
    template<typename T>
    static T apply(const VolumeAtomic<T>* volume);

    template<typename S>
    static tgt::Vector2<S> apply(const VolumeAtomic<tgt::Vector2<S> >* volume);

    template<typename S>
    static tgt::Vector3<S> apply(const VolumeAtomic<tgt::Vector3<S> >* volume);

    template<typename S>
    static tgt::Vector4<S> apply(const VolumeAtomic<tgt::Vector4<S> >* volume);

    template<typename S>
    static Tensor2<S> apply(const VolumeAtomic<Tensor2<S> >* volume);
};

// ============================================================================

template<typename T>
T VolumeOperatorMinValue::apply(const VolumeAtomic<T>* volume) {
    const T* voxel = volume->voxel();
    T min = std::numeric_limits<T>::max();
    for (size_t i = 0; i < volume->getNumVoxels(); ++i) {
        if (voxel[i] < min)
            min = voxel[i];
    }
    return min;
}

// specialized template version for tgt::VectorX clases which do no implement operator>().
//

template<typename S>
S VolumeOperatorMinValue::apply(const VolumeAtomic<tgt::Vector2<S> >* volume) {
    const tgt::Vector2<S>* voxel = volume->voxel();
    S min = std::numeric_limits<S>::max();
    S voxelMin = min;

    for (size_t i = 0; i < volume->getNumVoxels(); ++i) {
        voxelMin = tgt::min(voxel[i]);
        if (voxelMin < min)
            min = voxelMin;
    }
    return min;
}

template<typename S>
S VolumeOperatorMinValue::apply(const VolumeAtomic<tgt::Vector3<S> >* volume) {
    const tgt::Vector3<S>* voxel = volume->voxel();
    S min = std::numeric_limits<S>::max();
    S voxelMin = min;

    for (size_t i = 0; i < volume->getNumVoxels(); ++i) {
        voxelMin = tgt::min(voxel[i]);
        if (voxelMin < min)
            min = voxelMin;
    }
    return min;
}

template<typename S>
S VolumeOperatorMinValue::apply(const VolumeAtomic<tgt::Vector4<S> >* volume) {
    const tgt::Vector4<S>* voxel = volume->voxel();
    S min = std::numeric_limits<S>::max();
    S voxelMin = min;

    for (size_t i = 0; i < volume->getNumVoxels(); ++i) {
        voxelMin = tgt::min(voxel[i]);
        if (voxelMin < min)
            min = voxelMin;
    }
    return min;
}

template<typename S>
S VolumeOperatorMinValue::apply(const VolumeAtomic<Tensor2<S> >* volume) {
    const Tensor2<S>* voxel = volume->voxel();
    S min = std::numeric_limits<S>::max();

    for (size_t i = 0; i < volume->getNumVoxels(); ++i) {
        for (int j = 0; j < 6; ++j) {
            if (voxel[i].elem[j] < min)
                min = voxel[i].elem[j];
        }
    }
    return min;
}

// ============================================================================

template<typename T>
T VolumeOperatorMaxValue::apply(const VolumeAtomic<T>* volume) {
    const T* voxel = volume->voxel();
    T max = std::numeric_limits<T>::min();
    for (size_t i = 0; i < volume->getNumVoxels(); ++i) {
        if (voxel[i] > max)
            max = voxel[i];
    }
    return max;
}

// specialized template version for tgt::VectorX classes which do no implement operator>().

template<typename S>
tgt::Vector2<S> VolumeOperatorMaxValue::apply(const VolumeAtomic<tgt::Vector2<S> >* volume) {
    const tgt::Vector2<S>* voxel = volume->voxel();
    S max = std::numeric_limits<S>::min();
    S voxelMax = max;

    for (size_t i = 0; i < volume->getNumVoxels(); ++i) {
        voxelMax = tgt::max(voxel[i]);
        if (voxelMax > max)
            max = voxelMax;
    }
    return tgt::Vector2<S>(max);
}

template<typename S>
tgt::Vector3<S> VolumeOperatorMaxValue::apply(const VolumeAtomic<tgt::Vector3<S> >* volume) {
    const tgt::Vector3<S>* voxel = volume->voxel();
    S max = std::numeric_limits<S>::min();
    S voxelMax = max;

    for (size_t i = 0; i < volume->getNumVoxels(); ++i) {
        voxelMax = tgt::max(voxel[i]);
        if (voxelMax > max)
            max = voxelMax;
    }
    return tgt::Vector3<S>(max);
}

template<typename S>
tgt::Vector4<S> VolumeOperatorMaxValue::apply(const VolumeAtomic<tgt::Vector4<S> >* volume) {
    const tgt::Vector4<S>* voxel = volume->voxel();
    S max = std::numeric_limits<S>::min();
    S voxelMax = max;

    for (size_t i = 0; i < volume->getNumVoxels(); ++i) {
        voxelMax = tgt::max(voxel[i]);
        if (voxelMax > max)
            max = voxelMax;
    }
    return tgt::Vector4<S>(max);
}

template<typename S>
Tensor2<S> VolumeOperatorMaxValue::apply(const VolumeAtomic<Tensor2<S> >* volume) {
    const Tensor2<S>* voxel = volume->voxel();
    S max = std::numeric_limits<S>::min();

    for (size_t i = 0; i < volume->getNumVoxels(); ++i) {
        for (int j = 0; j < 6; ++j) {
            if (voxel[i].elem[j] > max)
                max = voxel[i].elem[j];
        }
    }
    return Tensor2<S>(max);
}

} // namespace

#endif // VRN_VOLUMEOPERATOR_H
