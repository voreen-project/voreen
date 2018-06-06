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

#ifndef VRN_VOLUMEOPERATORCURVATURE_H
#define VRN_VOLUMEOPERATORCURVATURE_H

//#include "voreen/core/datastructures/volume/volumeoperator.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "tgt/vector.h"

namespace voreen {

class VRN_CORE_API VolumeOperatorCurvature {
public:
    VolumeOperatorCurvature(){};

    /**
     * Calculates curvature.
     *
     * Use uint8_t or col4 as template argument in order to generate 8 or 4x8 bit datasets.
     * Use uint16_t or Vector4<uint16_t> as template argument in order
     * to generate 16 or 4x16 bit datasets.
     */
    template<typename U>
    Volume* apply(const VolumeBase* srcVolume, unsigned int curvatureType = 0) ;
private:
    /** T = Input   U = Output */
    template<typename T, typename U>
    Volume* calcCurvature(const VolumeBase* handle, unsigned int curvatureType);
};

//---------------------------------------------------------------------------------------------
//      apply function
//---------------------------------------------------------------------------------------------
    template<typename U>
    Volume* VolumeOperatorCurvature::apply(const VolumeBase* srcVolume, unsigned int curvatureType) {
        if (srcVolume->getRepresentation<VolumeRAM>()->getNumChannels() != 1) {
            LERRORC("calcCurvature", "calcCurvature needs an input volume with 1 intensity channel");
            return 0;
        }
        //case uint8_t
        if(dynamic_cast<const VolumeAtomic<uint8_t>*>(srcVolume->getRepresentation<VolumeRAM>())){
            return calcCurvature<uint8_t,U>(srcVolume,curvatureType);
        }  //case uint16_t
        else if (dynamic_cast<const VolumeAtomic<uint16_t>*>(srcVolume->getRepresentation<VolumeRAM>())){
            return calcCurvature<uint16_t,U>(srcVolume,curvatureType);
        } else {
            LERRORC("calcCuvature", "calcCuvature needs a 8-, 12- or 16-bit dataset as input");
            return 0;
        }
    }
//---------------------------------------------------------------------------------------------
//      curvature function
//---------------------------------------------------------------------------------------------
    template<typename T, typename U>
    Volume* VolumeOperatorCurvature::calcCurvature(const VolumeBase* handle, unsigned int curvatureType) {
        const VolumeAtomic<T>* input = dynamic_cast<const VolumeAtomic<T>*>(handle->getRepresentation<VolumeRAM>());
        VolumeAtomic<U>* result = new VolumeAtomic<U>(input->getDimensions());

        tgt::ivec3 pos;
        tgt::ivec3 dim = input->getDimensions();
        float minCurvature = std::numeric_limits<float>::max();
        float maxCurvature = std::numeric_limits<float>::min();
        for (pos.z = 0; pos.z < dim.z; pos.z++) {
            for (pos.y = 0; pos.y < dim.y; pos.y++) {
                for (pos.x = 0; pos.x < dim.x; pos.x++) {
                    if (pos.x >= 2 && pos.x < dim.x-2 &&
                        pos.y >= 2 && pos.y < dim.y-2 &&
                        pos.z >= 2 && pos.z < dim.z-2)
                    {

                        // fetch necessary data
                        float c = input->getVoxelNormalized(pos);

                        float r0 = input->getVoxelNormalized(pos+tgt::ivec3(1,0,0));
                        float r1 = input->getVoxelNormalized(pos+tgt::ivec3(2,0,0));
                        float l0 = input->getVoxelNormalized(pos+tgt::ivec3(-1,0,0));
                        float l1 = input->getVoxelNormalized(pos+tgt::ivec3(-2,0,0));

                        float u0 = input->getVoxelNormalized(pos+tgt::ivec3(0,1,0));
                        float u1 = input->getVoxelNormalized(pos+tgt::ivec3(0,2,0));
                        float d0 = input->getVoxelNormalized(pos+tgt::ivec3(0,-1,0));
                        float d1 = input->getVoxelNormalized(pos+tgt::ivec3(0,-2,0));

                        float f0 = input->getVoxelNormalized(pos+tgt::ivec3(0,0,1));
                        float f1 = input->getVoxelNormalized(pos+tgt::ivec3(0,0,2));
                        float b0 = input->getVoxelNormalized(pos+tgt::ivec3(0,0,-1));
                        float b1 = input->getVoxelNormalized(pos+tgt::ivec3(0,0,-2));

                        float ur0 = input->getVoxelNormalized(pos+tgt::ivec3(1,1,0));
                        float dr0 = input->getVoxelNormalized(pos+tgt::ivec3(1,-1,0));
                        float ul0 = input->getVoxelNormalized(pos+tgt::ivec3(-1,1,0));
                        float dl0 = input->getVoxelNormalized(pos+tgt::ivec3(-1,-1,0));

                        float fr0 = input->getVoxelNormalized(pos+tgt::ivec3(1,0,1));
                        float br0 = input->getVoxelNormalized(pos+tgt::ivec3(1,0,-1));
                        float fl0 = input->getVoxelNormalized(pos+tgt::ivec3(-1,0,1));
                        float bl0 = input->getVoxelNormalized(pos+tgt::ivec3(-1,0,-1));

                        float uf0 = input->getVoxelNormalized(pos+tgt::ivec3(0,1,1));
                        float ub0 = input->getVoxelNormalized(pos+tgt::ivec3(0,1,-1));
                        float df0 = input->getVoxelNormalized(pos+tgt::ivec3(0,-1,1));
                        float db0 = input->getVoxelNormalized(pos+tgt::ivec3(0,-1,-1));

                        tgt::vec3 gradient = tgt::vec3(l0-r0,d0-u0,b0-f0);

                        float gradientLength = length(gradient);
                        if (gradientLength == 0.0f) gradientLength = 1.0f;

                        tgt::vec3 n = -gradient / gradientLength;

                        tgt::mat3 nxn; // matrix to hold the outer product of n and n^T
                        for (int i=0; i<3; ++i)
                            for (int j=0; j<3; ++j)
                                nxn[i][j] = n[i]*n[j];

                        tgt::mat3 P = tgt::mat3::identity - nxn;

                        // generate Hessian matrix
                        float fxx = (((r1-c)/2.0f)-((c-l1)/2.0f))/2.0f;
                        float fyy = (((u1-c)/2.0f)-((c-d1)/2.0f))/2.0f;
                        float fzz = (((f1-c)/2.0f)-((c-b1)/2.0f))/2.0f;
                        float fxy = (((ur0-ul0)/2.0f)-((dr0-dl0)/2.0f))/2.0f;
                        float fxz = (((fr0-fl0)/2.0f)-((br0-bl0)/2.0f))/2.0f;
                        float fyz = (((uf0-ub0)/2.0f)-((df0-db0)/2.0f))/2.0f;
                        tgt::mat3 H;
                        H[0][0] = fxx;
                        H[0][1] = fxy;
                        H[0][2] = fxz;
                        H[1][0] = fxy;
                        H[1][1] = fyy;
                        H[1][2] = fyz;
                        H[2][0] = fxz;
                        H[2][1] = fyz;
                        H[2][2] = fzz;

                        tgt::mat3 G = -P*H*P / gradientLength;

                        // compute trace of G
                        float trace = G.t00 + G.t11 + G.t22;

                        // compute Frobenius norm of G
                        float F = 0.0f;
                        for (int i=0; i<3; ++i)
                            for (int j=0; j<3; ++j)
                                F += powf(std::abs(G[i][j]), 2.0f);
                        F = sqrt(F);

                        float kappa1 = (trace + sqrtf(2.0f * powf(F,2.0f) - powf(trace,2.0f))) / 2.0f;
                        float kappa2 = (trace - sqrtf(2.0f * powf(F,2.0f) - powf(trace, 2.0f))) / 2.0f;

                        float curvature = 0.f;
                        if (curvatureType == 0) // first principle
                            curvature = kappa1;
                        else if (curvatureType == 1) // second principle
                            curvature = kappa2;
                        else if (curvatureType == 2) // mean
                            curvature = (kappa1+kappa2)/2.0f;
                        else if (curvatureType == 3) // Gaussian
                            curvature = kappa1*kappa2;
                        result->voxel(pos) = static_cast<U>(curvature);

                        if (curvature < minCurvature) minCurvature = curvature;
                        else if (curvature > maxCurvature) maxCurvature = curvature;
                    } else
                        result->voxel(pos) = static_cast<U>(0.0f);
                }
            }
        }

        // scale curvature to lie in interval [0.0,1.0], where 0.5 equals zero curvature
        for (pos.z = 0; pos.z < dim.z; pos.z++) {
            for (pos.y = 0; pos.y < dim.y; pos.y++) {
                for (pos.x = 0; pos.x < dim.x; pos.x++) {
                    float c = result->getVoxelNormalized(pos, 0);
                    if (c < 0.0f) c /= -minCurvature;
                    else if (c >= 0.0f) c /= maxCurvature;
                    c /= 2.0f;
                    c += 0.5f;
                    result->voxel(pos) = static_cast<U>(c);
                }
            }
        }
        return new Volume(result, handle);
    }

} // namespace

#endif // VRN_VOLUMEOPERATORCURVATURE_H
