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

#include <cmath>
#include <omp.h>

extern "C"
{
   #include <cblas.h>
}

#include "cosinetransform.h"

void generateDCTCoefficients(size_t N, float* coefficients){
    float N_ = M_PI/static_cast<float>(N);
    
    for (size_t k = 0; k < N; ++k) {
        float k_ = static_cast<float>(k);
        for (size_t n = 0; n < N; ++n) {
            float n_ = static_cast<float>(n);
            coefficients[k*N+n] = std::cos(N_*n_*(k_+0.5f));
        }
    }
}


void performDCT(size_t N, size_t D, float* data, float* coefficients, float* result){
    if (D==1) {
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, 1, N, 1.f, coefficients, N, data, 1, 0.f, result, 1);
    } else { //D==2
        float* t = new float[N*N];
        cblas_sgemm(CblasRowMajor, CblasTrans, CblasNoTrans, N, N, N, 1.f, coefficients, N, data, N, 0.f, t, N);
        cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, 1.f, t, N, coefficients, N, 0.f, result, N);
        delete[] t;
    }
}

void performDCTmulti(size_t N, size_t D, size_t L, float* data, float* coefficients, float* result) {
    size_t planeSize = N;
    if (D>1) planeSize *= N;
    
//#pragma omp parallel for shared(N, D, data, coefficients, result, planeSize) num_threads(3)
    for (size_t plane = 0; plane < L; ++plane) {
        performDCT(N, D, &data[plane*planeSize], coefficients, &result[plane*planeSize]);
    }
}

/*
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
*/

void sampleXYfrom3D(size_t N, float* data, float* result){
    //this function is deprecated and can be replaced by samplePlanefrom3D with plane=0
    size_t layerNum = N/2;
    size_t z = layerNum*N*N;
    for (size_t y=0; y<N; ++y){
        for (size_t x=0; x<N; ++x){
            result[y*N+x]=data[z+y*N+x];
        }
    }
}
void sampleYZfrom3D(size_t N, float* data, float* result){
    //this function is deprecated and can be replaced by samplePlanefrom3D with plane=1
    size_t layerNum = N/2;
    size_t x = layerNum;
    for (size_t z=0; z<N; ++z){
        for (size_t y=0; y<N; ++y){
            result[z*N+y]=data[z*N*N+y*N+x];
        }
    }
}
void sampleXZfrom3D(size_t N, float* data, float* result){
    //this function is deprecated and can be replaced by samplePlanefrom3D with plane=2
    size_t layerNum = N/2;
    size_t y = layerNum*N;
    for (size_t z=0; z<N; ++z){
        for (size_t x=0; x<N; ++x){
            result[z*N+x]=data[z*N*N+y+x];
        }
    }
}
/*
void sample2_5raw(size_t N, float* data, float* result) {
    size_t layerSize = N*N;

//#pragma omp parallel for shared(N, data, result, layerSize) num_threads(3)
    for (size_t plane = 0; plane < 3; ++plane) {
        samplePlanefrom3D(N, plane, data, result + plane * layerSize);
    }

}
*/

