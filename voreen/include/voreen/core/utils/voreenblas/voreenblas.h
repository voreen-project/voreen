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

#ifndef VRN_VOREENBLAS_H
#define VRN_VOREENBLAS_H

#include <string>
#include <sstream>

#include "voreen/core/utils/voreenblas/ellpackmatrix.h"

namespace voreen {

/**
 * Interface for implementations of the Basic Linear Algebra Subprograms (BLAS).
 *
 * Level 1 and 2 operations in single precision as well as
 * a conjugate gradient solver for sparse matrices are provided.
 *
 * @see VoreenBlasCPU, a basic CPU implementation.
 */
class VRN_CORE_API VoreenBlas {

public:

    /**
     * Preconditioner to be used by
     * the conjugate gradient matrix solver.
     */
    enum ConjGradPreconditioner {
        NoPreconditioner,
        Jacobi
    };

    virtual void sAXPY(size_t vecSize, const float* vecx, const float* vecy, float alpha, float* result) const = 0;

    virtual float sDOT(size_t vecSize, const float* vecx, const float* vecy) const = 0;

    virtual float sNRM2(size_t vecSize, const float* vecx) const = 0;

    virtual void sSpMVEll(const EllpackMatrix<float>& mat, const float* vec, float* result) const = 0;

    virtual void hSpMVEll(const EllpackMatrix<int16_t>& mat, const float* vec, float* result) const = 0;

    virtual float sSpInnerProductEll(const EllpackMatrix<float>& mat, const float* vecx, const float* vecy) const = 0;

    virtual int sSpConjGradEll(const EllpackMatrix<float>& mat, const float* vec, float* result,
        float* initial = 0, ConjGradPreconditioner precond = NoPreconditioner, float threshold = 1e-4f, int maxIterations = 1000) const = 0;

    virtual int hSpConjGradEll(const EllpackMatrix<int16_t>& mat, const float* vec, float* result,
        float* initial = 0, float threshold = 1e-4f, int maxIterations = 1000) const = 0;

};

} // namespace

#endif // VRN_VOREENBLAS_H
