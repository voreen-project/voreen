/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014 Mathias J. Krause
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
 * Descriptor for 2D and 3D lattices with dynamic omega. In principle,
 * thanks to the fact that the OpenLB code is generic, it is sufficient
 * to write a new descriptor when a new type of lattice is to be used.
 *  -- header file
 */

#ifndef DYN_OMEGA_DESCRIPTOR_H
#define DYN_OMEGA_DESCRIPTOR_H

#include "latticeDescriptors.h"

namespace olb {

/// Descriptors for 2D and 3D lattices with variable omega.
/** The implementation is to be extended by combination with other
 * lattice descriptors.
 */

namespace descriptors {

struct DynOmegaDescriptor {
  static const int numScalars = 1;
  static const int numSpecies = 1;
  static const int omegaBeginsAt = 0;
  static const int sizeOfOmega = 1;
};

struct DynOmegaDescriptorBase {
  typedef DynOmegaDescriptor ExternalField;
};


/// 2D Descriptors for modells with variable omega

struct ForcedDynOmega2dDescriptor {
  static const int numScalars = 3;
  static const int numSpecies = 2;
  static const int omegaBeginsAt = 0;
  static const int sizeOfOmega = 1;
  static const int forceBeginsAt = 1;
  static const int sizeOfForce   = 2;
};

struct ForcedDynOmega2dDescriptorBase {
  typedef ForcedDynOmega2dDescriptor ExternalField;
};


template <typename T>
struct DynOmegaD2Q9Descriptor
  : public D2Q9DescriptorBase<T>, public DynOmegaDescriptorBase {
};

template <typename T>
struct ForcedDynOmegaD2Q9Descriptor
  : public D2Q9DescriptorBase<T>, public ForcedDynOmega2dDescriptorBase {
};


/// 3D Descriptors for modells with variable omega

struct ForcedDynOmega3dDescriptor {
  static const int numScalars = 4;
  static const int numSpecies = 2;
  static const int omegaBeginsAt = 0;
  static const int sizeOfOmega = 1;
  static const int forceBeginsAt = 1;
  static const int sizeOfForce   = 3;
};

struct ForcedDynOmega3dDescriptorBase {
  typedef ForcedDynOmega3dDescriptor ExternalField;
};

template <typename T>
struct DynOmegaD3Q19Descriptor
  : public D3Q19DescriptorBase<T>, public DynOmegaDescriptorBase {
};

template <typename T>
struct ForcedDynOmegaD3Q19Descriptor
  : public D3Q19DescriptorBase<T>, public ForcedDynOmega3dDescriptorBase {
};


} // namespace descriptors

} // namespace olb

#endif
