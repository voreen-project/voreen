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

#include "modules/openmp/include/voreenblasmp.h"

#include <cmath>
#include <cstring>

#ifdef WIN32
#include <omp.h>
#endif

namespace voreen {

const std::string VoreenBlasMP::loggerCat_("voreen.VoreenBlasMP");

void VoreenBlasMP::sAXPY(size_t vecSize, const float* vecx, const float* vecy, float alpha, float* result) const {

    #pragma omp parallel for
    for (int i=0; i<static_cast<int>(vecSize); ++i) {
        result[i] = alpha*vecx[i] + vecy[i];
    }

}

float VoreenBlasMP::sDOT(size_t vecSize, const float* vecx, const float* vecy) const {
    float result = 0.f;

    #pragma omp parallel for reduction(+:result)
    for (int i=0; i < static_cast<int>(vecSize); i++)
        result += (vecx[i] * vecy[i]);

    return result;
}

float VoreenBlasMP::sNRM2(size_t vecSize, const float* vecx) const {
    float result = 0.f;

    #pragma omp parallel for reduction(+:result)
    for (int i=0; i < static_cast<int>(vecSize); i++)
        result += (vecx[i] * vecx[i]);

    return sqrt(result);
}

void VoreenBlasMP::sSpMVEll(const EllpackMatrix<float>& mat, const float* vec, float* result) const {

    const int numRows = static_cast<int>(mat.getNumRows());
    const size_t numColsPerRow = mat.getNumColsPerRow();

    #pragma omp parallel for
    for (int row=0; row < numRows; ++row) {
        result[row] = 0.f;
        for (size_t colIndex=0; colIndex < numColsPerRow; ++colIndex)
            result[row] += mat.getValueByIndex(row, colIndex) * vec[mat.getColumn(row, colIndex)];
    }

}

void VoreenBlasMP::hSpMVEll(const EllpackMatrix<int16_t>& mat, const float* vec, float* result) const {

    const float maxValue = static_cast<float>((1<<15) - 1);
    const int numRows = static_cast<int>(mat.getNumRows());
    const size_t numColsPerRow = mat.getNumColsPerRow();

    #pragma omp parallel for
    for (int row=0; row < numRows; ++row) {
        result[row] = 0.f;
        for (size_t colIndex=0; colIndex < numColsPerRow; ++colIndex)
            result[row] += (mat.getValueByIndex(row, colIndex) * vec[mat.getColumn(row, colIndex)]) / maxValue;
    }
}

float VoreenBlasMP::sSpInnerProductEll(const EllpackMatrix<float>& mat, const float* vecx, const float* vecy) const {

    const int numRows = static_cast<int>(mat.getNumRows());
    const size_t numColsPerRow = mat.getNumColsPerRow();

    float result = 0.f;

    #pragma omp parallel for reduction(+:result)
    for (int row=0; row < numRows; ++row) {
        float dot = 0.f;
        for (size_t colIndex=0; colIndex < numColsPerRow; ++colIndex)
            dot += mat.getValueByIndex(row, colIndex) * vecy[mat.getColumn(row, colIndex)];
        result += dot*vecx[row];
    }
    return result;
}

int VoreenBlasMP::sSpConjGradEll(const EllpackMatrix<float>& mat, const float* vec, float* result,
                                 float* initial, ConjGradPreconditioner precond, float threshold, int maxIterations, ProgressReporter& progress) const {

    if (!mat.isSymmetric()) {
        LERROR("Symmetric matrix expected.");
        return -1;
    }

    size_t vecSize = mat.getNumRows();

    bool initialAllocated = false;
    if (initial == 0) {
        initial = new float[vecSize];
        initialAllocated = true;
        for (size_t i=0; i<vecSize; ++i)
            initial[i] = 0.f;
    }

    float* xBuf = result;
    memcpy(xBuf, initial, sizeof(float)*vecSize);

    float* tmpBuf = initial;
    memcpy(tmpBuf, vec, sizeof(float)*vecSize);

    float* rBuf = new float[vecSize];
    float* pBuf = new float[vecSize];

    float* zBuf = 0;

    EllpackMatrix<float>* preconditioner = 0;
    if (precond == Jacobi) {
        preconditioner = new EllpackMatrix<float>(vecSize, vecSize, 1);
        preconditioner->initializeBuffers();

        #pragma omp parallel for
        for (int i=0; i<static_cast<int>(vecSize); i++)
            preconditioner->setValueByIndex(i, i, 0, 1.f / std::max(mat.getValue(i,i), 1e-6f));

        zBuf = new float[vecSize];
    }

    float nominator;
    float denominator;

    int iteration = 0;

    // r <= A*x_0
    sSpMVEll(mat, xBuf, rBuf);

    // p <= -r + b
    sAXPY(vecSize, rBuf, tmpBuf, -1.f, pBuf);

    // r <= -r + b
    sAXPY(vecSize, rBuf, tmpBuf, -1.f, rBuf);

    // preconditioning
    if (precond == Jacobi) {
        sSpMVEll(*preconditioner, rBuf, zBuf);
        memcpy(pBuf, zBuf, vecSize * sizeof(float));
    }


    try {
        while (iteration < maxIterations) {
            progress.setProgress(static_cast<float>(iteration)/maxIterations);

            iteration++;

            // norm(r_k)
            nominator = sDOT(vecSize, rBuf, rBuf);

            if (precond == Jacobi)
                nominator = sDOT(vecSize, rBuf, zBuf);
            else
                nominator = sDOT(vecSize, rBuf, rBuf);

            // tmp <= A * p_k
            sSpMVEll(mat, pBuf, tmpBuf);

            // dot(p_k^T, tmp)
            denominator = sDOT(vecSize, pBuf, tmpBuf);

            float alpha = nominator / denominator;

            // x <= alpha*p + x
            sAXPY(vecSize, pBuf, xBuf, alpha, xBuf);

            // r <= -alpha*tmp + r
            sAXPY(vecSize, tmpBuf, rBuf, -alpha, rBuf);

            float beta;

            // norm(r_k+1)
            if (precond == Jacobi) {
                sSpMVEll(*preconditioner, rBuf, zBuf);
                beta = sDOT(vecSize, rBuf, zBuf);
            }
            else {
                beta = sDOT(vecSize, rBuf, rBuf);
            }

            if (sqrt(beta) < threshold)
                break;

            beta /= nominator;

            // p <= beta*p + r
            if (precond == Jacobi) {
                sAXPY(vecSize, pBuf, zBuf, beta, pBuf);
            }
            else {
                sAXPY(vecSize, pBuf, rBuf, beta, pBuf);
            }
        }
        progress.setProgress(1.0);
    } catch(boost::thread_interrupted& e) {
        delete[] rBuf;
        delete[] pBuf;
        delete[] zBuf;
        delete preconditioner;

        if (initialAllocated)
            delete[] initial;

        throw e;
    }

    delete[] rBuf;
    delete[] pBuf;
    delete[] zBuf;
    delete preconditioner;

    if (initialAllocated)
        delete[] initial;

    return iteration;
}

int VoreenBlasMP::hSpConjGradEll(const EllpackMatrix<int16_t>& mat, const float* vec, float* result,
                                    float* initial, float threshold, int maxIterations) const {

    if (!mat.isSymmetric()) {
        LERROR("Symmetric matrix expected.");
        return -1;
    }

    size_t vecSize = mat.getNumRows();

    bool initialAllocated = false;
    if (initial == 0) {
        initial = new float[vecSize];
        initialAllocated = true;
        for (size_t i=0; i<vecSize; ++i)
            initial[i] = 0.f;
    }

    float* xBuf = result;
    memcpy(xBuf, initial, sizeof(float)*vecSize);

    float* tmpBuf = initial;
    memcpy(tmpBuf, vec, sizeof(float)*vecSize);

    float* rBuf = new float[vecSize];
    float* pBuf = new float[vecSize];

    float nominator;
    float denominator;

    int iteration = 0;

    // r <= A*x_0
    hSpMVEll(mat, xBuf, rBuf);

    // p <= -r + b
    sAXPY(vecSize, rBuf, tmpBuf, -1.f, pBuf);

    // r <= -r + b
    sAXPY(vecSize, rBuf, tmpBuf, -1.f, rBuf);

    while (iteration < maxIterations) {

        iteration++;

        // norm(r_k)
        nominator = sDOT(vecSize, rBuf, rBuf);

        // tmp <= A * p_k
        hSpMVEll(mat, pBuf, tmpBuf);

        // dot(p_k^T, tmp)
        denominator = sDOT(vecSize, pBuf, tmpBuf);

        float alpha = nominator / denominator;

        // x <= alpha*p + x
        sAXPY(vecSize, pBuf, xBuf, alpha, xBuf);

        // r <= -alpha*tmp + r
        sAXPY(vecSize, tmpBuf, rBuf, -alpha, rBuf);

        // norm(r_k+1)
        float beta = sDOT(vecSize, rBuf, rBuf);

        if (sqrt(beta) < threshold)
            break;

        beta /= nominator;

        // p <= beta*p + r
        sAXPY(vecSize, pBuf, rBuf, beta, pBuf);
    }

    delete[] rBuf;
    delete[] pBuf;

    if (initialAllocated)
        delete[] initial;

    return iteration;
}

}   // namespace
