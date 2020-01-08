/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#include "voreen/core/utils/voreenblas/voreenblascpu.h"
#include "voreen/core/voreenapplication.h"

namespace voreen {

const std::string VoreenBlasCPU::loggerCat_("voreen.VoreenBlasCPU");

void VoreenBlasCPU::sAXPY( size_t vecSize, const float* vecx, const float* vecy, float alpha, float* result ) const {
    for (size_t i=0; i<vecSize; ++i) {
        result[i] = alpha*vecx[i] + vecy[i];
    }
}

float VoreenBlasCPU::sDOT(size_t vecSize, const float* vecx, const float* vecy) const {
    float result = 0.f;
    for (size_t i=0; i<vecSize; i++) {
        result += vecx[i]*vecy[i];
    }
    return result;
}

float VoreenBlasCPU::sNRM2(size_t vecSize, const float* vecx) const {
    float result = 0.f;
    for (size_t i=0; i<vecSize; i++) {
        result += vecx[i]*vecx[i];
    }
    return sqrt(result);
}

void VoreenBlasCPU::sSpMVEll(const EllpackMatrix<float>& mat, const float* vec, float* result) const {

    const size_t numRows = mat.getNumRows();
    const size_t numColsPerRow = mat.getNumColsPerRow();

    for (size_t row=0; row < numRows; ++row) {
        result[row] = 0.f;
        for (size_t colIndex=0; colIndex < numColsPerRow; ++colIndex) {
            size_t col = mat.getColumn(row, colIndex);
            if(col != -1) {
                result[row] += mat.getValueByIndex(row, colIndex) * vec[col];
            }
        }
    }
}
void VoreenBlasCPU::hSpMVEll(const EllpackMatrix<int16_t>& mat, const float* vec, float* result) const {

    const float maxValue = static_cast<float>((1<<15) - 1);
    const size_t numRows = mat.getNumRows();
    const size_t numColsPerRow = mat.getNumColsPerRow();

    for (size_t row=0; row < numRows; ++row) {
        result[row] = 0.f;
        for (size_t colIndex=0; colIndex < numColsPerRow; ++colIndex) {
            size_t col = mat.getColumn(row, colIndex);
            if(col != -1) {
                result[row] += (mat.getValueByIndex(row, colIndex) * vec[col]) / maxValue;
            }
        }
    }
}

float VoreenBlasCPU::sSpInnerProductEll(const EllpackMatrix<float>& mat, const float* vecx, const float* vecy) const {

    const size_t numRows = mat.getNumRows();
    const size_t numColsPerRow = mat.getNumColsPerRow();

    float result = 0.f;
    for (size_t row=0; row < numRows; ++row) {
        float dot = 0.f;
        for (size_t colIndex=0; colIndex < numColsPerRow; ++colIndex) {
            size_t col = mat.getColumn(row, colIndex);
            if(col != -1) {
                dot += mat.getValueByIndex(row, colIndex) * vecy[col];
            }
        }
        result += dot*vecx[row];
    }
    return result;
}

int VoreenBlasCPU::sSpConjGradEll(const EllpackMatrix<float>& mat, const float* vec, float* result,
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
        preconditioner = new EllpackMatrix<float>(mat.getNumRows(), mat.getNumRows(), 1);
        preconditioner->initializeBuffers();
        for (size_t i=0; i<mat.getNumRows(); i++)
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

int VoreenBlasCPU::hSpConjGradEll(const EllpackMatrix<int16_t>& mat, const float* vec, float* result,
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
