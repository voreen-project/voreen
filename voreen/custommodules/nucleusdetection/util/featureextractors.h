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

#ifndef FEATUREEXTRACTOR_H
#define FEATUREEXTRACTOR_H
/**
* @brief Generates 2.5D patches from a 3D patch
* @param N Number of samples per dimension
* @param data Pointer to the 3D patch (N*N*N elements)
* @param feature Pointer to the destination buffer
* @return the number of floats written to the buffer (i.e. the length of the feature)
*/
size_t sample2_5raw(size_t N, float* data, float* feature);

/**
* @brief Generates a serialised 3D patch 
* @param N Number of samples per dimension
* @param data Pointer to the 3D patch (N*N*N elements)
* @param feature Pointer to the destination buffer
* @return the number of floats written to the buffer (i.e. the length of the feature)
*
* This is just a memcpy
*/
size_t sample3raw(size_t N, float* data, float* feature);

/**
* @brief Calculates filter responses for a whole filter bank
* @param N Number of samples per dimension
* @param data Pointer to the 3D patch (N*N*N elements)
* @param K Number of filters in the filter bank
* @param filters Pointer to the filter bank (K*N*N*N elements)
* @param feature Pointer to the destination buffer
* @return the number of floats written to the buffer (i.e. the length of the feature)
*/
size_t sample3filtered(size_t N, float* data, size_t K, float* filters, float* feature);

/**
* @brief Calculates the 2D-DCTs for a 2.5D Patch
* @param N Number of samples per dimension
* @param data Pointer to the 2.5D patch (N*N*3 elements)
* @param dctCoefficients Pointer to the pre-calculated DCT coefficients (N*N elements)
* @param feature Pointer to the destination buffer
* @return the number of floats written to the buffer (i.e. the length of the feature)
*/
size_t sample2_5dct(size_t N, float* data, float* dctCoefficients, float* feature);

/**
* @brief Calculates the 1D-DCTs for a 3D Patch
* @param N Number of samples per dimension
* @param data Pointer to the 3D patch (N*N*N elements)
* @param dctCoefficients Pointer to the pre-calculated DCT coefficients (N*N elements)
* @param feature Pointer to the destination buffer
* @return the number of floats written to the buffer (i.e. the length of the feature)
*/
size_t sample1dct(size_t N, float* data, float* dctCoefficients, float* feature);

#endif // FEATUREEXTRACTOR_H
