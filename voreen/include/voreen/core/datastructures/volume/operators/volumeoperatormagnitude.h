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
        if (typeid(*(srcVolume->getRepresentation<VolumeRAM>())) == typeid(VolumeRAM_3xInt8))
            return calcGradientMagnitudesGeneric<tgt::Vector3<int8_t>, U>(srcVolume);
        else
        if (typeid(*(srcVolume->getRepresentation<VolumeRAM>())) == typeid(VolumeRAM_3xInt16))
            return calcGradientMagnitudesGeneric<tgt::Vector3<int16_t>, U >(srcVolume);
        else
        if (typeid(*(srcVolume->getRepresentation<VolumeRAM>())) == typeid(VolumeRAM_4xUInt8))
            return calcGradientMagnitudesGeneric<tgt::Vector4<uint8_t>, U>(srcVolume);
        else
        if (typeid(*(srcVolume->getRepresentation<VolumeRAM>())) == typeid(VolumeRAM_4xUInt16))
            return calcGradientMagnitudesGeneric<tgt::Vector4<uint16_t>, U>(srcVolume);
        else
        if (typeid(*(srcVolume->getRepresentation<VolumeRAM>())) == typeid(VolumeRAM_4xInt8))
            return calcGradientMagnitudesGeneric<tgt::Vector4<int8_t>, U>(srcVolume);
        else
        if (typeid(*(srcVolume->getRepresentation<VolumeRAM>())) == typeid(VolumeRAM_4xInt16))
            return calcGradientMagnitudesGeneric<tgt::Vector4<int16_t>, U >(srcVolume);
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
        // get RAM representation of the input volume and relevant meta data
        const VolumeAtomic<T>* input = dynamic_cast<const VolumeAtomic<T>*>(handle->getRepresentation<VolumeRAM>());
        RealWorldMapping rwm = handle->getRealWorldMapping();  
        tgt::svec3 dim = input->getDimensions();

        // create the output volume
        VolumeAtomic<U>* result = new VolumeAtomic<U>(dim);

        // we need to find out the max vector magnitude of the input first
        float maxMagnitude = 0.f;
        tgt::svec3 pos;
        for (pos.z = 0; pos.z < dim.z; pos.z++) {
            for (pos.y = 0; pos.y < dim.y; pos.y++) {
                for (pos.x = 0; pos.x < dim.x; pos.x++) {
                    // get value of all three channels for x, y, and z vector values and transform using real-world mapping
                    tgt::vec3 gradient;
                    gradient.x = rwm.normalizedToRealWorld(input->getVoxelNormalized(pos, 0));
                    gradient.y = rwm.normalizedToRealWorld(input->getVoxelNormalized(pos, 1));
                    gradient.z = rwm.normalizedToRealWorld(input->getVoxelNormalized(pos, 2));

                    float gradientMagnitude = tgt::length(gradient);
                    
                    maxMagnitude = std::max(maxMagnitude, gradientMagnitude);
                }
            }
        }

        LDEBUGC("calcGradientMagnitudes", "Found max. magnitude of " << maxMagnitude);

        // now we can normalize magnitudes to the range [0, 1] for storing floats
        for (pos.z = 0; pos.z < dim.z; pos.z++) {
            for (pos.y = 0; pos.y < dim.y; pos.y++) {
                for (pos.x = 0; pos.x < dim.x; pos.x++) {
                    tgt::vec3 gradient;
                    gradient.x = rwm.normalizedToRealWorld(input->getVoxelNormalized(pos, 0));
                    gradient.y = rwm.normalizedToRealWorld(input->getVoxelNormalized(pos, 1));
                    gradient.z = rwm.normalizedToRealWorld(input->getVoxelNormalized(pos, 2));

                    float gradientMagnitude = tgt::length(gradient);

                    // magnitude is always positive, so we can just normalize using the maximum
                    result->setVoxelNormalized(gradientMagnitude / maxMagnitude, pos);
                }
            }
        }
        
        Volume* magnitudeVolume = new Volume(result, handle);

        // Overwrite real-world mapping to rescale the normalized values
        RealWorldMapping rescaleMapping(maxMagnitude, 0.f, "");
        magnitudeVolume->setRealWorldMapping(rescaleMapping);

        return  magnitudeVolume;
    }

} // namespace

#endif // VRN_VOLUMEOPERATORMAGNITUDE_H
