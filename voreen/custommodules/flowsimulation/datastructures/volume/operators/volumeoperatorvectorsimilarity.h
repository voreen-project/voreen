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

//#include "voreen/core/datastructures/volume/volumeoperator.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "tgt/vector.h"

namespace voreen {

class VRN_CORE_API VolumeOperatorVectorSimilarity {
public:
    VolumeOperatorVectorSimilarity(){};

    /**
     * Calculates voxel-wise similarity between
     *
     * Use uint8_t or uint16_t as U template argument in order to generate 8 or 16 bit datasets.
     */
    template<typename U>
    Volume* apply(const VolumeBase* srcVolume, float p);
private:
    /** T = Input   U = Output */
    template<typename T, typename U>
    Volume* calcVectorSimilarityGeneric(const VolumeBase* handle, float p);
};

//---------------------------------------------------------------------------------------------
//      apply function
//---------------------------------------------------------------------------------------------
template<typename U>
Volume* VolumeOperatorVectorSimilarity::apply(const VolumeBase* srcVolume, float p) {
    if (typeid(*(srcVolume->getRepresentation<VolumeRAM>())) == typeid(VolumeRAM_3xUInt8))
        return calcVectorSimilarityGeneric<tgt::Vector3<uint8_t>, U>(srcVolume, p);
    else
    if (typeid(*(srcVolume->getRepresentation<VolumeRAM>())) == typeid(VolumeRAM_3xUInt16))
        return calcVectorSimilarityGeneric<tgt::Vector3<uint16_t>, U >(srcVolume, p);
    else
    if (typeid(*(srcVolume->getRepresentation<VolumeRAM>())) == typeid(VolumeRAM_3xInt8))
        return calcVectorSimilarityGeneric<tgt::Vector3<int8_t>, U>(srcVolume, p);
    else
    if (typeid(*(srcVolume->getRepresentation<VolumeRAM>())) == typeid(VolumeRAM_3xInt16))
        return calcVectorSimilarityGeneric<tgt::Vector3<int16_t>, U >(srcVolume, p);
    else
    if (typeid(*(srcVolume->getRepresentation<VolumeRAM>())) == typeid(VolumeRAM_4xUInt8))
        return calcVectorSimilarityGeneric<tgt::Vector4<uint8_t>, U>(srcVolume, p);
    else
    if (typeid(*(srcVolume->getRepresentation<VolumeRAM>())) == typeid(VolumeRAM_4xUInt16))
        return calcVectorSimilarityGeneric<tgt::Vector4<uint16_t>, U>(srcVolume, p);
    else
    if (typeid(*(srcVolume->getRepresentation<VolumeRAM>())) == typeid(VolumeRAM_4xInt8))
        return calcVectorSimilarityGeneric<tgt::Vector4<int8_t>, U>(srcVolume, p);
    else
    if (typeid(*(srcVolume->getRepresentation<VolumeRAM>())) == typeid(VolumeRAM_4xInt16))
        return calcVectorSimilarityGeneric<tgt::Vector4<int16_t>, U >(srcVolume, p);
    else
    if (typeid(*(srcVolume->getRepresentation<VolumeRAM>())) == typeid(VolumeRAM_3xFloat))
        return calcVectorSimilarityGeneric<tgt::Vector3<float>, U>(srcVolume, p);
    else
    if (typeid(*(srcVolume->getRepresentation<VolumeRAM>())) == typeid(VolumeRAM_3xDouble))
        return calcVectorSimilarityGeneric<tgt::Vector3<double>, U>(srcVolume, p);
    else
    if (typeid(*(srcVolume->getRepresentation<VolumeRAM>())) == typeid(VolumeRAM_4xFloat))
        return calcVectorSimilarityGeneric<tgt::Vector4<float>, U>(srcVolume, p);
    else
    if (typeid(*(srcVolume->getRepresentation<VolumeRAM>())) == typeid(VolumeRAM_4xDouble))
        return calcVectorSimilarityGeneric<tgt::Vector4<double>, U>(srcVolume, p);
    else {
        LERRORC("calcGradientMagnitudes", "Unhandled type!");
        return 0;
    }
}
//---------------------------------------------------------------------------------------------
//      similarity function
//---------------------------------------------------------------------------------------------
template<typename T, typename U>
Volume* VolumeOperatorVectorSimilarity::calcVectorSimilarityGeneric(const VolumeBase* handle, float p) {
    // get RAM representation of the input volume and relevant meta data
    const VolumeAtomic<T>* input = dynamic_cast<const VolumeAtomic<T>*>(handle->getRepresentation<VolumeRAM>());
    tgt::svec3 dim = input->getDimensions();

    // create the output volume
    VolumeAtomic<U>* result = new VolumeAtomic<U>(dim);

    tgt::svec3 pos;
    for (pos.z = 0; pos.z < dim.z; pos.z++) {
        for (pos.y = 0; pos.y < dim.y; pos.y++) {
            for (pos.x = 0; pos.x < dim.x; pos.x++) {
                // TODO:
            }
        }
    }

    Volume* volume = new Volume(result, handle);
    return  volume;
}

} // namespace

#endif // VRN_VOLUMEOPERATORVECTORSIMILARITY_H
