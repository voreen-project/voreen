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

#ifndef VRN_VOLUMEOPERATORVORTICITY_H
#define VRN_VOLUMEOPERATORVORTICITY_H

#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "volumeoperatorgradient.h"
#include "tgt/vector.h"

namespace voreen {

class VRN_CORE_API VolumeOperatorVorticity {
public:

    VolumeOperatorVorticity(){};

    template<typename U>
    Volume* apply(const VolumeBase* srcVolume, VolumeOperatorGradient::GradientType gt);

private:
/** applys    T = Input    U = Output */
    template<typename T, typename U>
    Volume* calcVorticities(const VolumeBase* handle, VolumeOperatorGradient::GradientType gt);

//utils
    //Max vorticity = max change over min distance:
    template<typename T>
    float getMaxVorticityLength(tgt::vec3 spacing);
    //Store vorticity in volume.
    template<typename T>
    void storeVorticity(tgt::vec3 vorticity, const tgt::ivec3& pos, VolumeAtomic<tgt::Vector3<T> >* result);
};

//---------------------------------------------------------------------------------------------
//      apply function
//---------------------------------------------------------------------------------------------
template<typename U>
Volume* VolumeOperatorVorticity::apply(const VolumeBase* srcVolume, VolumeOperatorGradient::GradientType gt) {
    switch(gt){
        case VolumeOperatorGradient::VOG_CENTRAL_DIFFERENCE:
            //case uint8_t
            if(dynamic_cast<const VolumeAtomic<tgt::Vector3<uint8_t> >*>(srcVolume->getRepresentation<VolumeRAM>())){
                return calcVorticities<uint8_t,U>(srcVolume, VolumeOperatorGradient::VOG_CENTRAL_DIFFERENCE);
            }  //case uint16_t
            else if (dynamic_cast<const VolumeAtomic<tgt::Vector3<uint16_t> >*>(srcVolume->getRepresentation<VolumeRAM>())){
                return calcVorticities<uint16_t,U>(srcVolume, VolumeOperatorGradient::VOG_CENTRAL_DIFFERENCE);
            } //case float
            else if (dynamic_cast<const VolumeAtomic<tgt::vec3>*>(srcVolume->getRepresentation<VolumeRAM>())){
                return calcVorticities<float,U>(srcVolume, VolumeOperatorGradient::VOG_CENTRAL_DIFFERENCE);
            } //case double
            else if (dynamic_cast<const VolumeAtomic<tgt::dvec3>*>(srcVolume->getRepresentation<VolumeRAM>())) {
                return calcVorticities<double, U>(srcVolume, VolumeOperatorGradient::VOG_CENTRAL_DIFFERENCE);
            } //unsupported input
            else {
                LERRORC("calcVorticitiesCentralDifferences", "Unsupported input");
                return 0;
            }
            break;
        case VolumeOperatorGradient::VOG_LINEAR_REGRESSION:
            //case uint8_t
            if(dynamic_cast<const VolumeAtomic<tgt::Vector3<uint8_t> >*>(srcVolume->getRepresentation<VolumeRAM>())){
                return calcVorticities<uint8_t,U>(srcVolume, VolumeOperatorGradient::VOG_LINEAR_REGRESSION);
            }  //case uint16_t
            else if (dynamic_cast<const VolumeAtomic<tgt::Vector3<uint16_t> >*>(srcVolume->getRepresentation<VolumeRAM>())){
                return calcVorticities<uint16_t,U>(srcVolume, VolumeOperatorGradient::VOG_LINEAR_REGRESSION);
            } //case float
            else if (dynamic_cast<const VolumeAtomic<tgt::vec3>*>(srcVolume->getRepresentation<VolumeRAM>())) {
                return calcVorticities<float, U>(srcVolume, VolumeOperatorGradient::VOG_LINEAR_REGRESSION);
            } //case double
            else if (dynamic_cast<const VolumeAtomic<tgt::dvec3>*>(srcVolume->getRepresentation<VolumeRAM>())) {
                return calcVorticities<double, U>(srcVolume, VolumeOperatorGradient::VOG_LINEAR_REGRESSION);
            } //unsupported input
            else {
                LERRORC("calcVorticitiesLinearRegression", "calcVorticitiesLinearRegression needs a 8-, 12- or 16-bit dataset as input");
                return 0;
            }
            break;
        case VolumeOperatorGradient::VOG_SOBEL:
            //case uint8_t
            if(dynamic_cast<const VolumeAtomic<tgt::Vector3<uint8_t> >*>(srcVolume->getRepresentation<VolumeRAM>())){
                return calcVorticities<uint8_t,U>(srcVolume, VolumeOperatorGradient::VOG_SOBEL);
            }  //case uint16_t
            else if (dynamic_cast<const VolumeAtomic<tgt::Vector3<uint16_t> >*>(srcVolume->getRepresentation<VolumeRAM>())){
                return calcVorticities<uint16_t,U>(srcVolume, VolumeOperatorGradient::VOG_SOBEL);
            } //case float
            else if (dynamic_cast<const VolumeAtomic<tgt::vec3>*>(srcVolume->getRepresentation<VolumeRAM>())) {
                return calcVorticities<float, U>(srcVolume, VolumeOperatorGradient::VOG_SOBEL);
            } //case double
            else if (dynamic_cast<const VolumeAtomic<tgt::dvec3>*>(srcVolume->getRepresentation<VolumeRAM>())) {
                return calcVorticities<double, U>(srcVolume, VolumeOperatorGradient::VOG_SOBEL);
            } //unsupported input
            else {
                LERRORC("calcVorticitiesSobel", "calcVorticitiesSobel needs a 8-, 12- or 16-bit dataset as input");
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
//      vorticity function
//---------------------------------------------------------------------------------------------

    template<typename T, typename U>
    Volume* VolumeOperatorVorticity::calcVorticities(const VolumeBase* handle, VolumeOperatorGradient::GradientType gt) {
        const VolumeAtomic<tgt::Vector3<T> >* input = dynamic_cast<const VolumeAtomic<tgt::Vector3<T> >*>(handle->getRepresentation<VolumeRAM>());
        VolumeAtomic<tgt::Vector3<U> >* result = new VolumeAtomic<tgt::Vector3<U> >(input->getDimensions());

        tgt::vec3 vorticity;
        tgt::svec3 pos;
        tgt::svec3 dim = input->getDimensions();
        //tgt::vec3 spacing = handle->getSpacing();

        //We normalize vorticities for integer datasets:
        bool normalizeVorticity = VolumeElement<T>::isInteger();
        float maxVorticityLength = 1.0f;
        if(normalizeVorticity) {
            maxVorticityLength = getMaxVorticityLength<T>(handle->getSpacing());
        }

        VolumeAtomic<T>* xVals = new VolumeAtomic<T>(input->getDimensions());
        VolumeAtomic<T>* yVals = new VolumeAtomic<T>(input->getDimensions());
        VolumeAtomic<T>* zVals = new VolumeAtomic<T>(input->getDimensions());

        for (pos.z = 0; pos.z < dim.z; pos.z++) {
            for (pos.y = 0; pos.y < dim.y; pos.y++) {
                for (pos.x = 0; pos.x < dim.x; pos.x++) {
                    tgt::Vector3<T> inputVec = input->voxel(pos);
                    xVals->voxel(pos) = inputVec.x;
                    yVals->voxel(pos) = inputVec.y;
                    zVals->voxel(pos) = inputVec.z;
                }
            }
        }

        Volume* xValsVol = new Volume(xVals, handle);
        Volume* yValsVol = new Volume(yVals, handle);
        Volume* zValsVol = new Volume(zVals, handle);

        VolumeOperatorGradient voOpGr;

        Volume* xvh = voOpGr.apply<U>(xValsVol, gt);
        const VolumeAtomic<tgt::Vector3<U> >* xGrad = xvh->getRepresentation<VolumeAtomic<tgt::Vector3<U> > >();
        delete xValsVol;

        Volume* yvh = voOpGr.apply<U>(yValsVol, gt);
        const VolumeAtomic<tgt::Vector3<U> >* yGrad = yvh->getRepresentation<VolumeAtomic<tgt::Vector3<U> > >();
        delete yValsVol;

        Volume* zvh = voOpGr.apply<U>(zValsVol, gt);
        const VolumeAtomic<tgt::Vector3<U> >* zGrad = zvh->getRepresentation<VolumeAtomic<tgt::Vector3<U> > >();
        delete zValsVol;

        for (pos.z = 0; pos.z < dim.z; pos.z++) {
            for (pos.y = 0; pos.y < dim.y; pos.y++) {
                for (pos.x = 0; pos.x < dim.x; pos.x++) {
                    tgt::Vector3<U> gradientX = xGrad->voxel(pos);
                    if (gradientX != tgt::Vector3<U>::zero) {
                        gradientX = tgt::normalize(gradientX);
                    }
                    tgt::Vector3<U> gradientY = yGrad->voxel(pos);
                    if (gradientY != tgt::Vector3<U>::zero) {
                        gradientY = tgt::normalize(gradientY);
                    }
                    tgt::Vector3<U> gradientZ = zGrad->voxel(pos);
                    if (gradientZ != tgt::Vector3<U>::zero) {
                        gradientZ = tgt::normalize(gradientZ);
                    }

                    vorticity.x = static_cast<float>(gradientZ.y - gradientY.z);
                    vorticity.y = static_cast<float>(gradientX.z - gradientZ.x);
                    vorticity.z = static_cast<float>(gradientY.x - gradientX.y);

                    if(normalizeVorticity)
                        vorticity /= maxVorticityLength;

                    storeVorticity(vorticity, pos, result);
                }
            }
        } 
        delete xvh;
        delete yvh;
        delete zvh;

        return new Volume(result, handle);
    }

//---------------------------------------------------------------------------------------------
//      private functions
//---------------------------------------------------------------------------------------------
    //Max vorticity = max change over min distance:
    template<typename T>
    float VolumeOperatorVorticity::getMaxVorticityLength(tgt::vec3 spacing) {
        return (VolumeElement<T>::rangeMaxElement() - VolumeElement<T>::rangeMinElement()) / min(spacing);
    }

    //Store vorticity in volume.
    template<typename T>
    void VolumeOperatorVorticity::storeVorticity(tgt::vec3 vorticity, const tgt::ivec3& pos, VolumeAtomic<tgt::Vector3<T> >* result) {
        if(VolumeElement<tgt::Vector3<T> >::isInteger()) {
            //map to [0,1]:
            vorticity += 1.0f;
            vorticity /= 2.0f;

            //...and to [minElement,maxElement]:
            vorticity *= VolumeElement<T>::rangeMaxElement() - VolumeElement<T>::rangeMinElement();
            vorticity += VolumeElement<T>::rangeMinElement();

            result->voxel(pos) = tgt::Vector3<T>(static_cast<T>(vorticity.x), static_cast<T>(vorticity.y), static_cast<T>(vorticity.z));
        }
        else {
            //floating point output, no remapping:
            result->voxel(pos) = tgt::Vector3<T>(static_cast<T>(vorticity.x), static_cast<T>(vorticity.y), static_cast<T>(vorticity.z));
        }
    }


} // namespace

#endif // VRN_VOLUMEOPERATORVORTICITY_H
