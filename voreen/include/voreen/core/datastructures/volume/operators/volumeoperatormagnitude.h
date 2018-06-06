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

#ifndef VRN_VOLUMEOPERATORMAGNITUDE_H
#define VRN_VOLUMEOPERATORMAGNITUDE_H

//#include "voreen/core/datastructures/volume/volumeoperator.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "tgt/vector.h"

namespace voreen {

class VRN_CORE_API VolumeOperatorMagnitude {
public:
    VolumeOperatorMagnitude(){};

    /**
     * Calculates gradient magnitudes from a gradient volume.
     *
     * Use uint8_t or uint16_t as U template argument in order to generate 8 or 16 bit datasets.
     */
    template<typename U>
    Volume* apply(const VolumeBase* srcVolume) ;
private:
    /** T = Input   U = Output */
    template<typename T, typename U>
    Volume* calcGradientMagnitudesGeneric(const VolumeBase* handle);
};

//---------------------------------------------------------------------------------------------
//      apply function
//---------------------------------------------------------------------------------------------
    template<typename U>
    Volume* VolumeOperatorMagnitude::apply(const VolumeBase* srcVolume) {
        if (typeid(*(srcVolume->getRepresentation<VolumeRAM>())) == typeid(VolumeRAM_3xUInt8))
            return calcGradientMagnitudesGeneric<tgt::Vector3<uint8_t>, U>(srcVolume);
        else
        if (typeid(*(srcVolume->getRepresentation<VolumeRAM>())) == typeid(VolumeRAM_3xUInt16))
            return calcGradientMagnitudesGeneric<tgt::Vector3<uint16_t>, U >(srcVolume);
        else
        if (typeid(*(srcVolume->getRepresentation<VolumeRAM>())) == typeid(VolumeRAM_4xUInt8))
            return calcGradientMagnitudesGeneric<tgt::Vector4<uint8_t>, U>(srcVolume);
        else
        if (typeid(*(srcVolume->getRepresentation<VolumeRAM>())) == typeid(VolumeRAM_4xUInt16))
            return calcGradientMagnitudesGeneric<tgt::Vector4<uint16_t>, U>(srcVolume);
        else
        if (typeid(*(srcVolume->getRepresentation<VolumeRAM>())) == typeid(VolumeRAM_3xFloat))
            return calcGradientMagnitudesGeneric<tgt::Vector3<float>, U>(srcVolume);
        else
        if (typeid(*(srcVolume->getRepresentation<VolumeRAM>())) == typeid(VolumeRAM_3xDouble))
            return calcGradientMagnitudesGeneric<tgt::Vector3<double>, U>(srcVolume);
        else
        if (typeid(*(srcVolume->getRepresentation<VolumeRAM>())) == typeid(VolumeRAM_4xFloat))
            return calcGradientMagnitudesGeneric<tgt::Vector4<float>, U>(srcVolume);
        else
        if (typeid(*(srcVolume->getRepresentation<VolumeRAM>())) == typeid(VolumeRAM_4xDouble))
            return calcGradientMagnitudesGeneric<tgt::Vector4<double>, U>(srcVolume);
        else {
            LERRORC("calcGradientMagnitudes", "Unhandled type!");
            return 0;
        }
    }
//---------------------------------------------------------------------------------------------
//      magnitude function
//---------------------------------------------------------------------------------------------
    template<typename T, typename U>
    Volume* VolumeOperatorMagnitude::calcGradientMagnitudesGeneric(const VolumeBase* handle) {
        const VolumeAtomic<T>* input = dynamic_cast<const VolumeAtomic<T>*>(handle->getRepresentation<VolumeRAM>());

        VolumeAtomic<U>* result = new VolumeAtomic<U>(input->getDimensions());

        //float maxValueT;
        //if ( typeid(*input) == typeid(VolumeRAM_3xUInt8)  ||
             //typeid(*input) == typeid(VolumeRAM_3xUInt16) ||
             //typeid(*input) == typeid(VolumeRAM_4xUInt8)  ||
             //typeid(*input) == typeid(VolumeRAM_4xUInt16)    ) {
             //int bitsT = input->getBitsStored() / input->getNumChannels();
            //maxValueT = static_cast<float>( (1 << bitsT) - 1);
        //}
        //else if ( typeid(*input) == typeid(VolumeRAM_3xFloat) || typeid(*input) == typeid(VolumeRAM_3xDouble) ||
            //typeid(*input) == typeid(VolumeRAM_4xFloat) || typeid(*input) == typeid(VolumeRAM_4xDouble) ) {
            //maxValueT = 1.f;
        //}
        //else {
            //LERRORC("calcGradientMagnitudes", "Unknown or unsupported input volume type");
            //tgtAssert(false, "Unknown or unsupported input volume type");
            //return result;
        //}

        float maxValueU;
        if ( typeid(*result) == typeid(VolumeRAM_UInt8)  ||
             typeid(*result) == typeid(VolumeRAM_UInt16) ||
             typeid(*result) == typeid(VolumeRAM_UInt32) )
            maxValueU = static_cast<float>( (1 << result->getBitsAllocated()) - 1);
        else if ( typeid(*result) == typeid(VolumeRAM_Float) || typeid(*result) == typeid(VolumeRAM_Double))
            maxValueU = 1.f;
        else {
            LERRORC("calcGradientMagnitudes", "Unknown or unsupported output volume type");
            tgtAssert(false, "Unknown or unsupported output volume type");
            return new Volume(result,handle);
        }


        tgt::ivec3 pos;
        tgt::ivec3 dim = input->getDimensions();
        for (pos.z = 0; pos.z < dim.z; pos.z++) {
            for (pos.y = 0; pos.y < dim.y; pos.y++) {
                for (pos.x = 0; pos.x < dim.x; pos.x++) {

                    tgt::vec3 gradient;
                    gradient.x = input->getVoxelNormalized(pos, 0);
                    gradient.y = input->getVoxelNormalized(pos, 1);
                    gradient.z = input->getVoxelNormalized(pos, 2);

                    // input value range is [0:maxValue] with (maxValue/2.f) corresponding to zero
                    gradient = (gradient*2.f)-1.f;

                    float gradientMagnitude = tgt::length(gradient);

                    //result->voxel(pos) = static_cast<U>( ( (derivative / maxValueT) / 2.f + 0.5f ) * maxValueU );
                    result->voxel(pos) = static_cast<U>( gradientMagnitude * maxValueU );

                }
            }
        }
        return new Volume(result,handle);
    }

} // namespace

#endif // VRN_VOLUMEOPERATORMAGNITUDE_H
