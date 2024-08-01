/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#ifndef VRN_VOLUMEOPERATORGRADIENT_H
#define VRN_VOLUMEOPERATORGRADIENT_H

#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "tgt/vector.h"

namespace voreen {

class VRN_CORE_API VolumeOperatorGradient {
public:
    enum GradientType{
        VOG_CENTRAL_DIFFERENCE = 0,
        VOG_LINEAR_REGRESSION = 1,
        VOG_SOBEL = 2,
    };

    VolumeOperatorGradient(){};

    template<typename U>
    Volume* apply(const VolumeBase* srcVolume, GradientType gt) ;
private:
/** applys    T = Input    U = Output */
    template<typename T, typename U>
    Volume* calcGradientsCentralDifferences(const VolumeBase* handle);
    template<typename T, typename U>
    Volume* calcGradientsLinearRegression(const VolumeBase* handle);
    template<typename T, typename U>
    Volume* calcGradientsSobel(const VolumeBase* handle);

//utils
    //Max gradient = max change over min distance:
    template<typename T>
    float getMaxGradientLength(tgt::vec3 spacing);
    //Store gradient in volume.
    template<typename T>
    void storeGradient(tgt::vec3 gradient, const tgt::svec3& pos, VolumeAtomic<tgt::Vector3<T> >* result);

public:
    /**
     * Calculates gradients with neighborhood of 26, using the central differences.
     * Note: Boundary layers are set to zero.
     */
    template<typename T>
    static tgt::vec3 calcGradientCentralDifferences(const VolumeAtomic<T>* input, const tgt::vec3& spacing, const tgt::svec3& pos);

    /**
     * Calculates gradients with neighborhood of 26, using the central differences.
     * Note: Boundary layers are set to zero.
     * Note: Uses normalized values, hence this implementation is slower than the templated version.
     */
    static tgt::vec3 calcGradientCentralDifferences(const VolumeRAM* input, const tgt::vec3& spacing, const tgt::svec3& pos, size_t channel = 0);

    /**
     * Calculates gradients using linear regression according to Neumann et al.
     *
     \verbatim
        L. Neumann, B. Csbfalvi, A. Koenig, and M. E. Groeller.
        Gradient estimation in volume data using 4D linear regression.
        In Proceedings of Eurographics 2000, pages 351-358, 2000.
     \endverbatim
     *
     * The neighboring voxels are weighted by their reciprocal Euclidean distance.
     *
     * Note: Boundary layers are set to zero.
     */
    template<typename T>
    tgt::vec3 calcGradientLinearRegression(const VolumeAtomic<T>* input, const tgt::vec3& spacing, const tgt::svec3& pos);

   /**
    * Calculates gradients with neighborhood of 26, using the Sobel filter.
    * Note: Boundary layers are set to zero.
    */
    template<class T>
    static tgt::vec3 calcGradientSobel(const VolumeAtomic<T>* input, const tgt::vec3& spacing, const tgt::svec3& pos);

    /**
     * Calculates gradients with neighborhood of 26, using the Sobel filter.
     * Note: Boundary layers are set to zero.
     * Note: Uses normalized values, hence this implementation is slower than the templated version.
     */
    static tgt::vec3 calcGradientSobel(const VolumeRAM* input, const tgt::vec3& spacing, const tgt::svec3& pos, const size_t channel);
};

//---------------------------------------------------------------------------------------------
//      apply function
//---------------------------------------------------------------------------------------------
template<typename U>
Volume* VolumeOperatorGradient::apply(const VolumeBase* srcVolume, GradientType gt) {
    switch(gt){
    case VOG_CENTRAL_DIFFERENCE:
        //case uint8_t
        if(dynamic_cast<const VolumeAtomic<uint8_t>*>(srcVolume->getRepresentation<VolumeRAM>())){
            return calcGradientsCentralDifferences<uint8_t,U>(srcVolume);
        }  //case uint16_t
        else if (dynamic_cast<const VolumeAtomic<uint16_t>*>(srcVolume->getRepresentation<VolumeRAM>())){
            return calcGradientsCentralDifferences<uint16_t,U>(srcVolume);
        } //case float
        else if (dynamic_cast<const VolumeAtomic<float>*>(srcVolume->getRepresentation<VolumeRAM>())){
            return calcGradientsCentralDifferences<float,U>(srcVolume);
        } //case double
        else if (dynamic_cast<const VolumeAtomic<double>*>(srcVolume->getRepresentation<VolumeRAM>())) {
            return calcGradientsCentralDifferences<double, U>(srcVolume);
        } //unsupported input
        else {
            LERRORC("calcGradientsCentralDifferences", "Unsupported input");
            return 0;
        }
        break;
    case VOG_LINEAR_REGRESSION:
        //case uint8_t
        if(dynamic_cast<const VolumeAtomic<uint8_t>*>(srcVolume->getRepresentation<VolumeRAM>())){
            return calcGradientsLinearRegression<uint8_t,U>(srcVolume);
        }  //case uint16_t
        else if (dynamic_cast<const VolumeAtomic<uint16_t>*>(srcVolume->getRepresentation<VolumeRAM>())){
            return calcGradientsLinearRegression<uint16_t,U>(srcVolume);
        } //case float
        else if (dynamic_cast<const VolumeAtomic<float>*>(srcVolume->getRepresentation<VolumeRAM>())) {
            return calcGradientsLinearRegression<float, U>(srcVolume);
        } //case double
        else if (dynamic_cast<const VolumeAtomic<double>*>(srcVolume->getRepresentation<VolumeRAM>())) {
            return calcGradientsLinearRegression<double, U>(srcVolume);
        } //unsupported input
        else {
            LERRORC("calcGradientsLinearRegression", "calcGradientsLinearRegression needs a 8-, 12- or 16-bit dataset as input");
            return 0;
        }
        break;
    case VOG_SOBEL:
        //case uint8_t
        if(dynamic_cast<const VolumeAtomic<uint8_t>*>(srcVolume->getRepresentation<VolumeRAM>())){
            return calcGradientsSobel<uint8_t,U>(srcVolume);
        }  //case uint16_t
        else if (dynamic_cast<const VolumeAtomic<uint16_t>*>(srcVolume->getRepresentation<VolumeRAM>())){
            return calcGradientsSobel<uint16_t,U>(srcVolume);
        } //case float
        else if (dynamic_cast<const VolumeAtomic<float>*>(srcVolume->getRepresentation<VolumeRAM>())) {
            return calcGradientsSobel<float, U>(srcVolume);
        } //case double
        else if (dynamic_cast<const VolumeAtomic<double>*>(srcVolume->getRepresentation<VolumeRAM>())) {
            return calcGradientsSobel<double, U>(srcVolume);
        } //unsupported input
        else {
            LERRORC("calcGradientsSobel", "calcGradientsSobel needs a 8-, 12- or 16-bit dataset as input");
            return 0;
        }
        break;
    default:
        return 0;
    }
    //should not get here
    return 0;
}

//---------------------------------------------------------------------------------------------
//      gradient functions
//---------------------------------------------------------------------------------------------
    template<typename T, typename U>
    Volume* VolumeOperatorGradient::calcGradientsCentralDifferences(const VolumeBase* handle) {
        const VolumeAtomic<T>* input = dynamic_cast<const VolumeAtomic<T>*>(handle->getRepresentation<VolumeRAM>());
        VolumeAtomic<tgt::Vector3<U> >* result = new VolumeAtomic<tgt::Vector3<U> >(input->getDimensions());

        tgt::vec3 gradient;
        tgt::svec3 pos;
        tgt::svec3 dim = input->getDimensions();
        tgt::vec3 spacing = handle->getSpacing();

        //We normalize gradients for integer datasets:
        bool normalizeGradient = VolumeElement<T>::isInteger();
        float maxGradientLength = 1.0f;
        if(normalizeGradient) {
            maxGradientLength = getMaxGradientLength<T>(handle->getSpacing());
        }

        for (pos.z = 0; pos.z < dim.z; pos.z++) {
            for (pos.y = 0; pos.y < dim.y; pos.y++) {
                for (pos.x = 0; pos.x < dim.x; pos.x++) {
                    gradient = calcGradientCentralDifferences(input, spacing, pos);

                    if(normalizeGradient)
                        gradient /= maxGradientLength;

                    storeGradient(gradient, pos, result);
                }
            }
        }
        return new Volume(result, handle);
    }

    template<typename T, typename U>
    Volume* VolumeOperatorGradient::calcGradientsLinearRegression(const VolumeBase* handle) {
        const VolumeAtomic<T>* input = dynamic_cast<const VolumeAtomic<T>*>(handle->getRepresentation<VolumeRAM>());
        VolumeAtomic<tgt::Vector3<U> >* result = new VolumeAtomic<tgt::Vector3<U> >(input->getDimensions());

        //We normalize gradients for integer datasets:
        bool normalizeGradient = VolumeElement<T>::isInteger();
        float maxGradientLength = 1.0f;
        if(normalizeGradient) {
            maxGradientLength = getMaxGradientLength<T>(handle->getSpacing());
        }

        tgt::vec3 gradient;
        tgt::svec3 pos;
        tgt::svec3 dim = input->getDimensions();
        tgt::vec3 spacing = handle->getSpacing();

        for (pos.z = 0; pos.z < dim.z; pos.z++) {
            for (pos.y = 0; pos.y < dim.y; pos.y++) {
                for (pos.x = 0; pos.x < dim.x; pos.x++) {
                    gradient = calcGradientLinearRegression(input, spacing, pos);

                    if(normalizeGradient)
                        gradient /= maxGradientLength;

                    storeGradient(gradient, pos, result);
                }
            }
        }

        return new Volume(result, handle);
    }

    template<typename T, typename U>
    Volume* VolumeOperatorGradient::calcGradientsSobel(const VolumeBase* handle) {
        const VolumeAtomic<T>* input = dynamic_cast<const VolumeAtomic<T>*>(handle->getRepresentation<VolumeRAM>());
        VolumeAtomic<tgt::Vector3<U> >* result = new VolumeAtomic<tgt::Vector3<U> >(input->getDimensions());

        //We normalize gradients for integer datasets:
        bool normalizeGradient = VolumeElement<T>::isInteger();
        float maxGradientLength = 1.0f;
        if(normalizeGradient) {
            maxGradientLength = getMaxGradientLength<T>(handle->getSpacing());
        }

        tgt::vec3 gradient;
        tgt::vec3 spacing = handle->getSpacing();
        tgt::svec3 pos;
        tgt::svec3 dim = input->getDimensions();

        for (pos.z = 0; pos.z < dim.z; pos.z++) {
            for (pos.y = 0; pos.y < dim.y; pos.y++) {
                for (pos.x = 0; pos.x < dim.x; pos.x++) {
                    gradient = calcGradientSobel(input, spacing, pos);

                    if(normalizeGradient)
                        gradient /= maxGradientLength;

                    storeGradient(gradient, pos, result);
                }
            }
        }
        return new Volume(result, handle);
    }

//---------------------------------------------------------------------------------------------
//      private functions
//---------------------------------------------------------------------------------------------
    //Max gradient = max change over min distance:
    template<typename T>
    float VolumeOperatorGradient::getMaxGradientLength(tgt::vec3 spacing) {
        return (VolumeElement<T>::rangeMaxElement() - VolumeElement<T>::rangeMinElement()) / min(spacing);
    }

    //Store gradient in volume.
    template<typename T>
    void VolumeOperatorGradient::storeGradient(tgt::vec3 gradient, const tgt::svec3& pos, VolumeAtomic<tgt::Vector3<T> >* result) {
        if(VolumeElement<tgt::Vector3<T> >::isInteger()) {
            //map to [0,1]:
            gradient += 1.0f;
            gradient /= 2.0f;

            //...and to [minElement,maxElement]:
            gradient *= VolumeElement<T>::rangeMaxElement() - VolumeElement<T>::rangeMinElement();
            gradient += VolumeElement<T>::rangeMinElement();

            result->voxel(pos) = tgt::Vector3<T>(static_cast<T>(gradient.x), static_cast<T>(gradient.y), static_cast<T>(gradient.z));
        }
        else {
            //floating point output, no remapping:
            result->voxel(pos) = tgt::Vector3<T>(static_cast<T>(gradient.x), static_cast<T>(gradient.y), static_cast<T>(gradient.z));
        }
    }

    template<typename T>
    tgt::vec3 VolumeOperatorGradient::calcGradientCentralDifferences(const VolumeAtomic<T>* input, const tgt::vec3& spacing, const tgt::svec3& pos) {
        if (pos != tgt::clamp(pos, tgt::svec3::one, input->getDimensions() - tgt::svec3::two)) {
            return tgt::vec3::zero;
        }

        T v0 = input->voxel(pos + tgt::svec3(1, 0, 0));
        T v1 = input->voxel(pos + tgt::svec3(0, 1, 0));
        T v2 = input->voxel(pos + tgt::svec3(0, 0, 1));
        T v3 = input->voxel(pos - tgt::svec3(1, 0, 0));
        T v4 = input->voxel(pos - tgt::svec3(0, 1, 0));
        T v5 = input->voxel(pos - tgt::svec3(0, 0, 1));

        tgt::vec3 gradient(static_cast<float>(v3 - v0), static_cast<float>(v4 - v1), static_cast<float>(v5 - v2));
        gradient /= (spacing * 2.0f);

        return gradient;
    }

    template<typename T>
    tgt::vec3 VolumeOperatorGradient::calcGradientLinearRegression(const VolumeAtomic<T>* input, const tgt::vec3& spacing, const tgt::svec3& pos) {
        if (pos != tgt::clamp(pos, tgt::svec3::one, input->getDimensions() - tgt::svec3::two)) {
            return tgt::vec3::zero;
        }

        // Euclidean weights for voxels with Manhattan distances of 1/2/3
        const float w_1 = 1.f;
        const float w_2 = 0.5f;
        const float w_3 = 1.f/3.f;

        const float w_A = 1.f / (8.f + 2.f/3.f);
        const float w_B = w_A;
        const float w_C = w_A;

        //left plane
        T v000 = input->voxel(pos + tgt::svec3(-1,-1,-1));
        T v001 = input->voxel(pos + tgt::svec3(-1, -1, 0));
        T v002 = input->voxel(pos + tgt::svec3(-1, -1, 1));
        T v010 = input->voxel(pos + tgt::svec3(-1, 0, -1));
        T v011 = input->voxel(pos + tgt::svec3(-1, 0, 0));
        T v012 = input->voxel(pos + tgt::svec3(-1, 0, 1));
        T v020 = input->voxel(pos + tgt::svec3(-1, 1, -1));
        T v021 = input->voxel(pos + tgt::svec3(-1, 1, 0));
        T v022 = input->voxel(pos + tgt::svec3(-1, 1, 1));

        //mid plane
        T v100 = input->voxel(pos + tgt::svec3(0, -1, -1));
        T v101 = input->voxel(pos + tgt::svec3(0, -1, 0));
        T v102 = input->voxel(pos + tgt::svec3(0, -1, 1));
        T v110 = input->voxel(pos + tgt::svec3(0, 0, -1));
        //T v111 = input->voxel(pos + svec3(0, 0, 0));
        T v112 = input->voxel(pos + tgt::svec3(0, 0, 1));
        T v120 = input->voxel(pos + tgt::svec3(0, 1, -1));
        T v121 = input->voxel(pos + tgt::svec3(0, 1, 0));
        T v122 = input->voxel(pos + tgt::svec3(0, 1, 1));

        //right plane
        T v200 = input->voxel(pos + tgt::svec3(1, -1, -1));
        T v201 = input->voxel(pos + tgt::svec3(1, -1, 0));
        T v202 = input->voxel(pos + tgt::svec3(1, -1, 1));
        T v210 = input->voxel(pos + tgt::svec3(1, 0, -1));
        T v211 = input->voxel(pos + tgt::svec3(1, 0, 0));
        T v212 = input->voxel(pos + tgt::svec3(1, 0, 1));
        T v220 = input->voxel(pos + tgt::svec3(1, 1, -1));
        T v221 = input->voxel(pos + tgt::svec3(1, 1, 0));
        T v222 = input->voxel(pos + tgt::svec3(1, 1, 1));

        tgt::vec3 gradient = tgt::vec3::zero;

        gradient.x = static_cast<float>( w_1 * ( v211 - v011 )               +
                w_2 * ( v201 + v210 + v212 + v221
                    -v001 - v010 - v012 - v021 ) +
                w_3 * ( v200 + v202 + v220 + v222
                    -v000 - v002 - v020 - v022 )   );

        gradient.y = static_cast<float>( w_1 * ( v121 - v101 )               +
                w_2 * ( v021 + v120 + v122 + v221
                    -v001 - v100 - v102 - v201 ) +
                w_3 * ( v020 + v022 + v220 + v222
                    -v000 - v002 - v200 - v202 )   );

        gradient.z = static_cast<float>( w_1 * ( v112 - v110 )               +
                w_2 * ( v012 + v102 + v122 + v212
                    -v010 - v100 - v120 - v210 ) +
                w_3 * ( v002 + v022 + v202 + v222
                    -v000 - v020 - v200 - v220 )   );

        gradient.x *= w_A;
        gradient.y *= w_B;
        gradient.z *= w_C;

        gradient /= spacing;
        gradient *= -1.f;

        return gradient;
    }

    template<class T>
    tgt::vec3 VolumeOperatorGradient::calcGradientSobel(const VolumeAtomic<T>* input, const tgt::vec3& spacing, const tgt::svec3& pos) {
        if (pos != tgt::clamp(pos, tgt::svec3::one, input->getDimensions() - tgt::svec3::two)) {
            return tgt::vec3::zero;
        }

        //left plane
        T v000 = input->voxel(pos + tgt::svec3(-1, -1, -1));
        T v001 = input->voxel(pos + tgt::svec3(-1, -1, 0));
        T v002 = input->voxel(pos + tgt::svec3(-1, -1, 1));
        T v010 = input->voxel(pos + tgt::svec3(-1, 0, -1));
        T v011 = input->voxel(pos + tgt::svec3(-1, 0, 0));
        T v012 = input->voxel(pos + tgt::svec3(-1, 0, 1));
        T v020 = input->voxel(pos + tgt::svec3(-1, 1, -1));
        T v021 = input->voxel(pos + tgt::svec3(-1, 1, 0));
        T v022 = input->voxel(pos + tgt::svec3(-1, 1, 1));
        //mid plane
        T v100 = input->voxel(pos + tgt::svec3(0, -1, -1));
        T v101 = input->voxel(pos + tgt::svec3(0, -1, 0));
        T v102 = input->voxel(pos + tgt::svec3(0, -1, 1));
        T v110 = input->voxel(pos + tgt::svec3(0, 0, -1));
        //T v111 = input->voxel(pos + svec3(0, 0, 0)); //not needed for calculation
        T v112 = input->voxel(pos + tgt::svec3(0, 0, 1));
        T v120 = input->voxel(pos + tgt::svec3(0, -1, -1));
        T v121 = input->voxel(pos + tgt::svec3(0, -1, 0));
        T v122 = input->voxel(pos + tgt::svec3(0, -1, 1));
        //right plane
        T v200 = input->voxel(pos + tgt::svec3(1, -1, -1));
        T v201 = input->voxel(pos + tgt::svec3(1, -1, 0));
        T v202 = input->voxel(pos + tgt::svec3(1, -1, 1));
        T v210 = input->voxel(pos + tgt::svec3(1, 0, -1));
        T v211 = input->voxel(pos + tgt::svec3(1, 0, 0));
        T v212 = input->voxel(pos + tgt::svec3(1, 0, 1));
        T v220 = input->voxel(pos + tgt::svec3(1, 1, -1));
        T v221 = input->voxel(pos + tgt::svec3(1, 1, 0));
        T v222 = input->voxel(pos + tgt::svec3(1, 1, 1));

        tgt::vec3 gradient = tgt::vec3::zero;

        //filter x-direction
        gradient.x += -1 * v000;
        gradient.x += -3 * v010;
        gradient.x += -1 * v020;
        gradient.x += 1 * v200;
        gradient.x += 3 * v210;
        gradient.x += 1 * v220;
        gradient.x += -3 * v001;
        gradient.x += -6 * v011;
        gradient.x += -3 * v021;
        gradient.x += +3 * v201;
        gradient.x += +6 * v211;
        gradient.x += +3 * v221;
        gradient.x += -1 * v002;
        gradient.x += -3 * v012;
        gradient.x += -1 * v022;
        gradient.x += +1 * v202;
        gradient.x += +3 * v212;
        gradient.x += +1 * v222;

        //filter y-direction
        gradient.y += -1 * v000;
        gradient.y += -3 * v100;
        gradient.y += -1 * v200;
        gradient.y += +1 * v020;
        gradient.y += +3 * v120;
        gradient.y += +1 * v220;
        gradient.y += -3 * v001;
        gradient.y += -6 * v101;
        gradient.y += -3 * v201;
        gradient.y += +3 * v021;
        gradient.y += +6 * v121;
        gradient.y += +3 * v221;
        gradient.y += -1 * v002;
        gradient.y += -3 * v102;
        gradient.y += -1 * v202;
        gradient.y += +1 * v022;
        gradient.y += +3 * v122;
        gradient.y += +1 * v222;

        //filter z-direction
        gradient.z += -1 * v000;
        gradient.z += -3 * v100;
        gradient.z += -1 * v200;
        gradient.z += +1 * v002;
        gradient.z += +3 * v102;
        gradient.z += +1 * v202;
        gradient.z += -3 * v010;
        gradient.z += -6 * v110;
        gradient.z += -3 * v210;
        gradient.z += +3 * v012;
        gradient.z += +6 * v112;
        gradient.z += +3 * v212;
        gradient.z += -1 * v020;
        gradient.z += -3 * v120;
        gradient.z += -1 * v220;
        gradient.z += +1 * v022;
        gradient.z += +3 * v122;
        gradient.z += +1 * v222;

        gradient /= 22.f;   // sum of all positive weights
        gradient /= 2.f;    // this mask has a step length of 2 voxels
        gradient /= spacing;
        gradient *= -1.f;

        return gradient;
    }

} // namespace

#endif // VRN_VOLUMEOPERATORGRADIENT_H
