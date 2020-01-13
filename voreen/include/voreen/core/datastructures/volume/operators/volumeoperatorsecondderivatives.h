/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#ifndef VRN_VOLUMEOPERATORSECONDDERIVATIVES_H
#define VRN_VOLUMEOPERATORSECONDDERIVATIVES_H

//#include "voreen/core/datastructures/volume/volumeoperator.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "tgt/vector.h"

namespace voreen {

class VRN_CORE_API VolumeOperatorSecondDerivatives {
public:
    VolumeOperatorSecondDerivatives(){};

    /**
     * Computes an simple approximation of the second directional derivative along the gradient direction at each voxel.
     * The calculation is done by applying the Laplacian operator.
     *
     * Use uint8_t or uint16_t as U template argument in order to generate 8 or 16 bit datasets.
     */
    template<typename U>
    Volume* apply(const VolumeBase* srcVolume) ;
private:
    /** T = Input    U = Output */
    template<typename T, typename U>
    Volume* calc2ndDerivatives(const VolumeBase* handle);
};

//---------------------------------------------------------------------------------------------
//      apply function
//---------------------------------------------------------------------------------------------
    template<typename U>
    Volume* VolumeOperatorSecondDerivatives::apply(const VolumeBase* srcVolume) {
        //case uint8_t
        if(dynamic_cast<const VolumeAtomic<uint8_t>*>(srcVolume->getRepresentation<VolumeRAM>())){
            return calc2ndDerivatives<uint8_t,U>(srcVolume);
        }  //case uint16_t
        else if (dynamic_cast<const VolumeAtomic<uint16_t>*>(srcVolume->getRepresentation<VolumeRAM>())){
            return calc2ndDerivatives<uint16_t,U>(srcVolume);
        } //case uint32_t
        else if (dynamic_cast<const VolumeAtomic<uint32_t>*>(srcVolume->getRepresentation<VolumeRAM>())){
            return calc2ndDerivatives<uint32_t,U>(srcVolume);
        } else {
            LERRORC("calc2ndDerivatives", "Unknown or unsupported input volume type");
            return 0;
        }
    }
//---------------------------------------------------------------------------------------------
//      second derivatives
//---------------------------------------------------------------------------------------------
    template<typename T, typename U>
    Volume* VolumeOperatorSecondDerivatives::calc2ndDerivatives(const VolumeBase* handle) {

        VolumeAtomic<U>* result = new VolumeAtomic<U>(handle->getDimensions());

        float maxValueT;
        if ( typeid(*(handle->getRepresentation<VolumeRAM>())) == typeid(VolumeRAM_UInt8)  ||
             typeid(*(handle->getRepresentation<VolumeRAM>())) == typeid(VolumeRAM_UInt16) ||
             typeid(*(handle->getRepresentation<VolumeRAM>())) == typeid(VolumeRAM_UInt32) )
            maxValueT = static_cast<float>( (1 << handle->getRepresentation<VolumeRAM>()->getBitsAllocated()) - 1);
        else if ( typeid(*(handle->getRepresentation<VolumeRAM>())) == typeid(VolumeRAM_Float) || typeid(*(handle->getRepresentation<VolumeRAM>())) == typeid(VolumeRAM_Double))
            maxValueT = 1.f;
        else {
            LERRORC("calc2ndDerivatives", "Unknown or unsupported input volume type");
            tgtAssert(false, "Unknown or unsupported input volume type");
            return new Volume(result,handle);
        }

        float maxValueU;
        if ( typeid(*result) == typeid(VolumeRAM_UInt8)  ||
             typeid(*result) == typeid(VolumeRAM_UInt16) ||
             typeid(*result) == typeid(VolumeRAM_UInt32) )
            maxValueU = static_cast<float>( (1 << result->getBitsAllocated()) - 1);
        else if ( typeid(*result) == typeid(VolumeRAM_Float) || typeid(*result) == typeid(VolumeRAM_Double))
            maxValueU = 1.f;
        else {
            LERRORC("calc2ndDerivatives", "Unknown or unsupported output volume type");
            tgtAssert(false, "Unknown or unsupported output volume type");
            return new Volume(result,handle);
        }


        float derivative;
        tgt::ivec3 pos;
        tgt::ivec3 dim = handle->getRepresentation<VolumeRAM>()->getDimensions();
        for (pos.z = 0; pos.z < dim.z; pos.z++) {
            for (pos.y = 0; pos.y < dim.y; pos.y++) {
                for (pos.x = 0; pos.x < dim.x; pos.x++) {

                    if (pos.x >= 1 && pos.x < dim.x-1 &&
                        pos.y >= 1 && pos.y < dim.y-1 &&
                        pos.z >= 1 && pos.z < dim.z-1)
                    {
                        //left plane
                        //T v000 = input->voxel(pos + ivec3(-1,-1,-1));
                        //T v001 = input->voxel(pos + ivec3(-1, -1, 0));
                        //T v002 = input->voxel(pos + ivec3(-1, -1, 1));
                        //T v010 = input->voxel(pos + ivec3(-1, 0, -1));
                        T v011 = handle->getRepresentation<VolumeAtomic<T> >()->voxel(pos + tgt::ivec3(-1, 0, 0));
                        //T v012 = input->voxel(pos + ivec3(-1, 0, 1));
                        //T v020 = input->voxel(pos + ivec3(-1, 1, -1));
                        //T v021 = input->voxel(pos + ivec3(-1, 1, 0));
                        //T v022 = input->voxel(pos + ivec3(-1, 1, 1));

                        //mid plane
                        //T v100 = input->voxel(pos + ivec3(0, -1, -1));
                        T v101 = handle->getRepresentation<VolumeAtomic<T> >()->voxel(pos + tgt::ivec3(0, -1, 0));
                        //T v102 = input->voxel(pos + ivec3(0, -1, 1));
                        T v110 = handle->getRepresentation<VolumeAtomic<T> >()->voxel(pos + tgt::ivec3(0, 0, -1));
                        T v111 = handle->getRepresentation<VolumeAtomic<T> >()->voxel(pos + tgt::ivec3(0, 0, 0));
                        T v112 = handle->getRepresentation<VolumeAtomic<T> >()->voxel(pos + tgt::ivec3(0, 0, 1));
                        //T v120 = input->voxel(pos + ivec3(0, 1, -1));
                        T v121 = handle->getRepresentation<VolumeAtomic<T> >()->voxel(pos + tgt::ivec3(0, 1, 0));
                        //T v122 = input->voxel(pos + ivec3(0, 1, 1));

                        //right plane
                        //T v200 = input->voxel(pos + ivec3(1, -1, -1));
                        //T v201 = input->voxel(pos + ivec3(1, -1, 0));
                        //T v202 = input->voxel(pos + ivec3(1, -1, 1));
                        //T v210 = input->voxel(pos + ivec3(1, 0, -1));
                        T v211 = handle->getRepresentation<VolumeAtomic<T> >()->voxel(pos + tgt::ivec3(1, 0, 0));
                        //T v212 = input->voxel(pos + ivec3(1, 0, 1));
                        //T v220 = input->voxel(pos + ivec3(1, 1, -1));
                        //T v221 = input->voxel(pos + ivec3(1, 1, 0));
                        //T v222 = input->voxel(pos + ivec3(1, 1, 1));

                        // Original Laplace operator
                        derivative = static_cast<float>( -6.0*v111                  +
                                                          v101 + v110 + v112 + v121  +
                                                          v011 + v211     );

                        // Simple Laplacian of Gaussian
                        /*derivative = static_cast<double>( -36.0*v111 +
                                                          4.0*v101 + 4.0*v110 + 4.0*v112 + 4.0*v121  +
                                                          4.0*v011 + 4.0*v211 +
                                                          v100 + v102 + v120 + v122 +
                                                          v001 + v010 + v012 + v021 +
                                                          v201 + v210 + v212 + v221      ) / 8.0;  */

                    }
                    else {
                        derivative = 0.f;
                    }

                    // map value from [-maxT:maxT] to [0:maxU] since we expect an unsigned volume as input/output type
                    result->voxel(pos) = static_cast<U>( ( (derivative / maxValueT) / 2.f + 0.5f ) * maxValueU );
                }
            }
        }
        return new Volume(result,handle);
    }

} // namespace

#endif // VRN_VOLUMEOPERATORSECONDDERIVATIVES_H
