/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Jonas Latt
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
*/

/** \file
 * Descriptor for all types of 2D and 3D lattices. In principle, thanks
 * to the fact that the OpenLB code is generic, it is sufficient to
 * write a new descriptor when a new type of lattice is to be used.
 *  -- header file
 */
#ifndef MRT_LATTICE_DESCRIPTORS_H
#define MRT_LATTICE_DESCRIPTORS_H

#include "latticeDescriptors.h"
#include <vector>

namespace olb {

/// Descriptors for the 2D and 3D lattices.
/** \warning Attention: The lattice directions must always be ordered in
 * such a way that c[i] = -c[i+(q-1)/2] for i=1..(q-1)/2, and c[0] = 0 must
 * be the rest velocity. Furthermore, the velocities c[i] for i=1..(q-1)/2
 * must verify
 *  - in 2D: (c[i][0]<0) || (c[i][0]==0 && c[i][1]<0)
 *  - in 3D: (c[i][0]<0) || (c[i][0]==0 && c[i][1]<0)
 *                       || (c[i][0]==0 && c[i][1]==0 && c[i][2]<0)
 * Otherwise some of the code will work erroneously, because the
 * aformentioned relations are taken as given to enable a few
 * optimizations.
*/
namespace descriptors {

/// MRT D2Q9 lattice. The numbering follows the one in "Viscous flow computations
/// with the method of lattice Boltzmann equation", D. Yu, L.-S. Luo, W. Shi,
/// Progress in Aerospace Sciences 39, (2003), p. 329-367
template <typename T>
struct MRTD2Q9DescriptorBase : public D2Q9DescriptorBase<T> {
  enum { d_ = 2, q_ = 9 };     ///< number of dimensions/distr. functions
  static const T M[q_][q_];    // Matrix of base change between f and moments : moments=M.f
  static const T invM[q_][q_]; // inverse of base change matrix : f=invM.moments
  static const T S[q_];       // relaxation times
  enum {shearIndexes = 2};
  static const int shearViscIndexes[shearIndexes]; // relevant indexes of r. t. for shear viscosity
  static const int bulkViscIndex  = 2; // relevant index of r. t. for bulk viscosity
};

/// MRT D3Q19 lattice. The numbering follows the one in "Multiple-relaxation-
/// time lattice Boltzmann models in three dimensions", D. D'Humières,
/// I. Ginzburg, M. Krafzcyk, P. Lallemand, L.-S. Luo,
/// Phil. Trans. R. Soc. Lond. A (2002) 660, p. 437-451
template <typename T>
struct MRTD3Q19DescriptorBase : public D3Q19DescriptorBase<T> {
  enum { d_ = 3, q_ = 19 };     ///< number of dimensions/distr. functions
  static const T M[q_][q_];    // Matrix of base change between f and moments : moments=M.f
  static const T invM[q_][q_]; // inverse of base change matrix : f=invM.moments
  static const T S[q_];       // relaxation times
  static const T S_2[q_];       // relaxation times
  enum {shearIndexes = 5};
  static const int shearViscIndexes[shearIndexes]; // relevant indexes of r. t. for shear viscosity
  static const int bulkViscIndex  = 1; // relevant index of r. t. for bulk viscosity
};

template <typename T>
struct MRTD2Q9Descriptor
  : public MRTD2Q9DescriptorBase<T>, public NoExternalFieldBase {
};

template <typename T>
struct MRTD3Q19Descriptor
  : public MRTD3Q19DescriptorBase<T>, public NoExternalFieldBase {
};

template <typename T>
struct ForcedMRTD2Q9Descriptor
  : public MRTD2Q9DescriptorBase<T>, public Force2dDescriptorBase {
};

template <typename T>
struct ForcedMRTD3Q19Descriptor
  : public MRTD3Q19DescriptorBase<T>, public Force3dDescriptorBase {
};

}  // namespace descriptors

}  // namespace olb

#endif
