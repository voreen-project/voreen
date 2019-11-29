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

#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable

/**
 * Single-precision vector sum: z = alpha*X + Y;
 */
__kernel void sAXPY(
  vrn_size_t n,
  __global const float* X,
  __global const float* Y,
  float alpha,
  __global float* z
)
{
    // get work item id
    size_t gid = get_global_id(0);
    size_t globalSize = get_global_size(0);

    size_t block = 0;
    while ((block*globalSize+gid) < n) {
        z[block*globalSize+gid] = alpha*X[block*globalSize+gid] + Y[block*globalSize+gid];
        block++;
    }
}

/**
 * Single-precision vector sum: z = alpha*(nom/den)*X + Y;
 */
__kernel void sAXPYDiv(
  vrn_size_t n,
  __global const float* X,
  __global const float* Y,
  float alpha,
  __global const float* nom,
  __global const float* den,
  __global float* z
)
{
    // get work item id
    size_t gid = get_global_id(0);
    size_t globalSize = get_global_size(0);

    float m = nom[0];
    float d = den[0];
    float scalar = alpha * (m == d ? 1.0 : m / d); // Special case for nom==den==0.0

    size_t block = 0;
    while ((block*globalSize+gid) < n) {
        z[block*globalSize+gid] = scalar*X[block*globalSize+gid] + Y[block*globalSize+gid];
        block++;
    }
}

/**
 * Single-precision dot product: z = X*Y
 */
__kernel void sDOT(
  vrn_size_t n,
  __global const float* X,
  __global const float* Y,
  __global float* z,
  __global int* mutex,
  __local float* sdata
)
{
    // get local (within workgroup) and global indices of workitem
    size_t tid = get_local_id(0);
    size_t gid = get_global_id(0);
    size_t globalSize = get_global_size(0);

    // Each work-item accumulates as many products as necessary
    // and writes the result to shared workgroup mem
    float partialProduct = 0.0;
    size_t block = 0;
    while ((block*globalSize+gid) < n) {
        partialProduct += X[block*globalSize+gid]*Y[block*globalSize+gid];
        block++;
    }
    sdata[tid] = partialProduct;

    // synchronize all work items before reduction
    barrier(CLK_LOCAL_MEM_FENCE);

    // sum up component products within workgroup by reduction in shared mem
    for (size_t s=1; s < get_local_size(0); s *= 2) {
        if (2*s*tid < get_local_size(0)) {
            sdata[2*s*tid] += sdata[2*s*tid + s];
        }

        // synchronize workgroup after each reduction step
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    // first global item initializes result buffer and releases mutex
    if (gid==0) {
        z[0] = 0.0;
        atom_cmpxchg(mutex,0,1);
    }

    // no sdata[0] of each workgroup contains the sum of the workgroup's component products
    // => first item of each workgroup adds this to the final result
    if (tid == 0) {
        while (atom_cmpxchg(mutex,1,0)==0);  // acquire mutex
        z[0] += sdata[0];
        atom_cmpxchg(mutex,0,1);             // release mutex
    }
}

/**
 * Single-precision vector 2-norm (eclidean distance): y = ||X||_2
 */
__kernel void sNRM2(
  vrn_size_t n,
  __global const float* X,
  __global float* y,
  __global int* mutex,
  __global int* semaphor,
  __local float* sdata
)
{
    // get local (within workgroup) and global indices of workitem
    size_t tid = get_local_id(0);
    size_t gid = get_global_id(0);
    size_t globalSize = get_global_size(0);
    size_t numGroups = get_num_groups(0);

    // Each work-item accumulates as many products as necessary
    // and writes the result to shared workgroup mem
    float partialProduct = 0.0;
    size_t block = 0;
    while ((block*globalSize+gid) < n) {
        partialProduct += X[block*globalSize+gid]*X[block*globalSize+gid];
        block++;
    }
    sdata[tid] = partialProduct;

    // synchronize all work items before reduction
    barrier(CLK_LOCAL_MEM_FENCE);

    // sum up component products within workgroup by reduction in shared mem
    for (size_t s=1; s < get_local_size(0); s *= 2) {
        if (2*s*tid < get_local_size(0)) {
            sdata[2*s*tid] += sdata[2*s*tid + s];
        }

        // synchronize workgroup after each reduction step
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    // first global item initializes result buffer/semaphor and releases mutex
    if (gid==0) {
        y[0] = 0.0;
        semaphor[0] = 0;
        atom_cmpxchg(mutex,0,1);
    }

    // no sdata[0] of each workgroup contains the sum of the workgroup's component products
    // => first item of each workgroup adds this to the final result
    if (tid == 0) {
        while (atom_cmpxchg(mutex,1,0)==0);  // acquire mutex
        y[0] += sdata[0];                    // add workgroup result to final result
        semaphor[0]++;                       // count number of groups
        if (semaphor[0] == numGroups)        // all groups have passed => compute final result
            y[0] = sqrt(y[0]);
        atom_cmpxchg(mutex,0,1);             // release mutex
    }
}

/**
 * Single-precision sparse matrix-vector multiplication
 * for a matrix stored in ELLPACK format: y = M*x;
 */
__kernel void sSpMV_Ell(
  vrn_size_t num_rows,
  vrn_size_t num_cols,
  vrn_size_t num_cols_per_row,
  __global const vrn_size_t* indices,
  __global const float* M,
  __global const float* x,
  __global float* y
)
{

    // get work item id
    size_t gid = get_global_id(0);
    size_t globalSize = get_global_size(0);

    // each work item iterates over a range of rows and computes the dot products
    size_t block = 0;
    while ((block*globalSize+gid) < num_rows) {
        size_t row = block*globalSize+gid;
        float dot = 0.0;
        // iterate over cols of current matrix row
        const size_t start = num_cols_per_row * row;
        const size_t end = start + num_cols_per_row;
        for (size_t index = start; index != end; index++) {
            size_t col = indices[index];
            if(col != -1) {
                dot += M[index] * x[col];
            }
        }
        y[row] = dot;
        block++;
    }
}

/**
 * Sparse matrix-vector multiplication for a
 * half-precision matrix (16 bit integer)
 * stored in ELLPACK format: y = M*x;
 */
__kernel void hSpMV_Ell(
  vrn_size_t num_rows,
  vrn_size_t num_cols,
  vrn_size_t num_cols_per_row,
  __global const vrn_size_t* indices,
  __global const short* M,
  __global const float* x,
  __global float* y
)
{

    // get work item id
    size_t gid = get_global_id(0);
    size_t globalSize = get_global_size(0);

    float valueScale = 1.0 / (float)((1<<15) - 1);

    // each work item iterates over a range of rows and computes the dot products
    size_t block = 0;
    while ((block*globalSize+gid) < num_rows) {
        size_t row = block*globalSize+gid;
        float dot = 0.0;
        // iterate over cols of current matrix row
        const size_t start = num_cols_per_row * row;
        const size_t end = start + num_cols_per_row;
        for (size_t index = start; index != end; index++) {
            size_t col = indices[index];
            if(col != -1) {
                dot += (M[index]*valueScale) * x[col];
            }
        }
        y[row] = dot;
        block++;
    }
}

/**
 * Single-precision general inner product for a sparse matrix
 * stored in ELLPACK format: z = x^t*A*y;
 */
__kernel void sInnerProduct_Ell(
  vrn_size_t num_rows,
  vrn_size_t num_cols,
  vrn_size_t num_cols_per_row,
  __global const vrn_size_t* indices,
  __global const float* A,
  __global const float* x,
  __global const float* y,
  __global float* z,
  __global int* mutex,
  __local float* sdata
)
{

    // get work item id
    size_t tid = get_local_id(0);
    size_t gid = get_global_id(0);
    size_t globalSize = get_global_size(0);

    // each work item iterates over a range of rows and computes the dot products
    size_t block = 0;
    float partialProduct = 0.0;
    while ((block*globalSize+gid) < num_rows) {
        size_t row = block*globalSize+gid;
        float dot = 0.0;
        // iterate over cols of current matrix row
        const size_t start = num_cols_per_row * row;
        const size_t end = start + num_cols_per_row;
        for (size_t index = start; index != end; index++) {
            size_t col = indices[index];
            if(col != -1) {
                dot += A[index] * y[col];
            }
        }
        partialProduct += dot * x[row];
        block++;
    }
    sdata[tid] = partialProduct;

    // synchronize all work items before reduction
    barrier(CLK_LOCAL_MEM_FENCE);

    // sum up component products within workgroup by reduction in shared mem
    for (size_t s=1; s < get_local_size(0); s *= 2) {
        if (2*s*tid < get_local_size(0)) {
            sdata[2*s*tid] += sdata[2*s*tid + s];
        }

        // synchronize workgroup after each reduction step
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    // first global item initializes result buffer and releases mutex
    if (gid==0) {
        z[0] = 0.0;
        atom_cmpxchg(mutex,0,1);
    }

    // no sdata[0] of each workgroup contains the sum of the workgroup's component products
    // => first item of each workgroup adds this to the final result
    if (tid == 0) {
        while (atom_cmpxchg(mutex,1,0)==0);  // acquire mutex
        z[0] += sdata[0];
        atom_cmpxchg(mutex,0,1);             // release mutex
    }
}
