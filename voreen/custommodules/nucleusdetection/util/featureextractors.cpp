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

#include <omp.h>
#include "cosinetransform.h"
#include "featureextractors.h"
#include <cstring>

void samplePlanefrom3D(size_t N, size_t plane, float* data, float* result){
    // A sample at position [x,y,z] has index x + y*N + z*N*N in the data set.
    // To extract samples from the YZ-plane, one has to loop over y and z, keeping x constant.
    // This function loops over a and b, while each of them is multiplied with one of the coefficients 1, N or N*N
    // for the YZ-plane (plane=1), coeffA will be N, coeffB will be N*N and x will have a constant value of N/2
    // plane=2 is the XZ-plane, plane=0 is the XY-plane
    size_t layerOffset;
    size_t coeffA;
    size_t coeffB;
    switch (plane) {
        case 1: coeffA=N; coeffB=N*N; layerOffset=N/2; break;
        case 2: coeffA=1; coeffB=N*N; layerOffset=N/2 *N; break;
        default: coeffA=1; coeffB=N; layerOffset=N/2 * N*N;
    }
    for (size_t b=0; b<N; ++b){
        for (size_t a=0; a<N; ++a){
            result[b*N+a]=data[layerOffset + b*coeffB + a*coeffA];
        }
    }
}

size_t sample2_5raw(size_t N, float* data, float* feature){
    size_t layerSize = N*N;
    size_t featureLength = layerSize*3;
//#pragma omp parallel for shared(N, data, result, layerSize) num_threads(3)
    for (size_t plane = 0; plane < 3; ++plane) {
        samplePlanefrom3D(N, plane, data, feature + plane * layerSize);
    }
    return featureLength;
}

size_t sample3raw(size_t N, float* data, float* feature){
    size_t featureLength = N*N*N;
    std::memcpy(feature, data, featureLength * sizeof(float)); 
    return featureLength;
}

size_t sample3filtered(size_t N, float* data, size_t K, float* filters, float* feature){
    size_t featureLength = K; //number of filters
    size_t filterSize = N*N*N; //elements per filter
    
#pragma omp parallel for shared(data, filters, feature)
    for (size_t kernel = 0; kernel < K; ++kernel) {
        size_t offset = filterSize * kernel;
        float current = 0.f;
        for (size_t kernelEntry = 0; kernelEntry < filterSize; ++kernelEntry) {
            current += data[kernelEntry] * filters[ offset + kernelEntry];
        }
        // normalize feature by number of voxels in patch
        feature[kernel] = current / static_cast<float>(filterSize);
    }
    //end omp parallel
    return featureLength;
}

size_t sample2_5dct(size_t N, float* data, float* dctCoefficients, float* feature){
    size_t featureLength = 3*N*N;
    performDCTmulti(N, 2, 3, data, dctCoefficients, feature);
    return featureLength;
}

size_t sample1dct(size_t N, float* data, float* dctCoefficients, float* feature){
    size_t featureLength = 3*N;
//    float patches1D[3 * N];
    float* patches1D = new float[3*N];
    int centerIndex = N*N*N/2;
    for (size_t direction=0; direction <3; ++direction){
        int offset = direction * N;
        
        int step = 1;        
        if (direction==1) step=N;
        else if (direction==2) step=N*N;
        
        for (int index=-static_cast<int>(N)/2; index<=static_cast<int>(N)/2; ++index){
            patches1D[offset+index+(N/2)]=data[centerIndex+index*step];
        }
    }
    performDCTmulti(N, 1, 3, patches1D, dctCoefficients, feature);

    delete[] patches1D;

    return featureLength;
}
