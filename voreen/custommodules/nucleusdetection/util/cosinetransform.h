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

#ifndef COSINETRANSFORM_H
#define COSINETRANSFORM_H

#include <cstddef>
/**
 * @brief Generates a matrix of coefficients for performing the DCT.
 * @param N Number of samples per dimension
 * @param coefficients N*N matrix containing the coefficients for calculating the DCT as M*d (1D) or 
 * transpose(M) * D * M (2D)
 */
void generateDCTCoefficients(size_t N, float* coefficients);

/**
 * @brief Performs DCT on a single 1 or 2 dimensional Dataset.
 * @param N Size of dataset per dimension
 * @param D Dimensions of the dataset (only 1 or 2 are valid values)
 * @param data 1- or 2-D data to calculate DCT on. Must have size pow(N,D)
 * @param coefficients Coefficients matrix
 * @param result DCT coefficients. Will be the same size as data.
 */
void performDCT(size_t N, size_t D, float* data, float* coefficients, float* result);

/**
 * @brief Performs DCT on multiple 1 or 2 dimensional Dataset.
 * @param N Size of dataset per dimension
 * @param D Dimensions of the dataset (only 1 or 2 are valid values)
 * @param L Number of datasets concatenated.
 * @param data 1- or 2-D data to calculate DCT on. Must have size L*pow(N,D).
 * @param coefficients Coefficients matrix.
 * @param result DCT coefficients. Will be the same size as data.
 */
void performDCTmulti(size_t N, size_t D, size_t L, float* data, float* coefficients, float* result);

/**
 * @brief Samples a single axis-orthogonal plane from a 3D dataset.
 * @param N Size of dataset per dimension.
 * @param plane Index of the plane, to sample from (0=XY, 1=YZ, 2=XZ).
 * @param data 3D dataset, must have size N*N*N.
 * @param result 2D dataset, must have size N*N.
 */
//void samplePlanefrom3D(size_t N, size_t plane, float* data, float* result);

/**
 * @brief Samples XY plane from a 3D dataset.
 * @param N Size of dataset per dimension.
 * @param data 3D dataset, must have size N*N*N.
 * @param result 2D dataset, must have size N*N.
 */
void sampleXYfrom3D(size_t N, float* data, float* result);

/**
 * @brief Samples YZ plane from a 3D dataset.
 * @param N Size of dataset per dimension.
 * @param data 3D dataset, must have size N*N*N.
 * @param result 2D dataset, must have size N*N.
 */
void sampleYZfrom3D(size_t N, float* data, float* result);

/**
 * @brief Samples XZ plane from a 3D dataset.
 * @param N Size of dataset per dimension.
 * @param data 3D dataset, must have size N*N*N.
 * @param result 2D dataset, must have size N*N.
 */
void sampleXZfrom3D(size_t N, float* data, float* result);

/**
 * @brief Samples a 2.5D Sample from a 3-D dataset by concatenating axis-orthogonal planes.
 * @param N Size of dataset per dimension.
 * @param data 3D dataset, must have size N*N*N.
 * @param result 2.5D dataset with XY-, YZ- and XZ-plane concatenated in this order. Size is 3*N*N.
 */
//void sample2_5raw(size_t N, float* data, float* result);

#endif	// COSINETRANSFORM_H

