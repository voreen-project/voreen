/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2017 Albert Mink
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
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA 02110-1301, USA.
*/

/** \file
 *  -- header file
 */
#ifndef RTLBM_DESCRIPTORS_H
#define RTLBM_DESCRIPTORS_H

#include "latticeDescriptors.h"

namespace olb {

namespace descriptors {

//TODO AM: get control of BaseDescriptor, since template specialization in lbHelpers is mostly wrong for this stencils

/// D3Q7 lattice for radiative transport problems @2016 A. Mink et al.
template <typename T>
struct D3Q7DescriptorBaseRTLBM {
  typedef D3Q7DescriptorBaseRTLBM<T> BaseDescriptor;
  enum { d = 3, q = 7 };     ///< number of dimensions/distr. functions
  static const int vicinity;  ///< size of neighborhood
  static const int c[q][d];   ///< lattice directions
  static const int opposite[q]; ///< opposite entry
  static const T t[q];        ///< lattice weights
  static const T invCs2;      ///< inverse square of speed of sound
  static const double henyeyPhaseFunction[q][q]; ///<anisotropic discrete scattering coefficient
};

template <typename T>
struct D3Q7DescriptorRTLBM : public D3Q7DescriptorBaseRTLBM<T>, public NoExternalFieldBase { };


/// D3Q19 lattice TODO: AM
template <typename T>
struct D3Q19DescriptorBaseRTLBM {
  typedef D3Q19DescriptorBaseRTLBM<T> BaseDescriptor;
  enum { d = 3, q = 19 };     ///< number of dimensions/distr. functions
  static const int vicinity;  ///< size of neighborhood
  static const int c[q][d];   ///< lattice directions
  static const int opposite[q]; ///< opposite entry
  static const T t[q];        ///< lattice weights
  static const T invCs2;      ///< inverse square of speed of sound
  static const double henyeyPhaseFunction[q][q]; ///<anisotropic discrete scattering coefficient
};

template <typename T>
struct D3Q19DescriptorRTLBM : public D3Q19DescriptorBaseRTLBM<T>, public NoExternalFieldBase {};



/** D3Q27 lattice.
 *  zero direction only need for correct stream process. Contains dummy values.
 */
template <typename T>
struct D3Q27DescriptorBaseRTLBM {
  typedef D3Q27DescriptorBaseRTLBM<T> BaseDescriptor;
  enum { d = 3, q = 27 };     ///< number of dimensions/distr. functions
  static const int vicinity;  ///< size of neighborhood
  static const int c[q][d];   ///< lattice directions
  static const int opposite[q]; ///< opposite entry
  static const T t[q];        ///< lattice weights
  static const T invCs2;      ///< inverse square of speed of sound
  static const double henyeyPhaseFunction[q][q]; ///<anisotropic discrete scattering coefficient
};

template <typename T>
struct D3Q27DescriptorRTLBM : public D3Q27DescriptorBaseRTLBM<T>, public NoExternalFieldBase {};


}  // namespace descriptors

}  // namespace olb

#endif
