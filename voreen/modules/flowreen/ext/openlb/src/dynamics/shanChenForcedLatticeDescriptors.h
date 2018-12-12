/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2008 Andrea Parmigiani, Orestis Malaspinas,
 *  Jonas Latt
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
#ifndef SHAN_CHEN_FORCED_LATTICE_DESCRIPTORS_H
#define SHAN_CHEN_FORCED_LATTICE_DESCRIPTORS_H

#include <vector>

namespace olb {

namespace descriptors {

struct ShanChenForcedExternal2Ddescriptor {
  static const int numScalars = 4;
  static const int numSpecies = 2;
  static const int velocityBeginsAt = 0;
  static const int sizeOfVelocity   = 2;
  static const int forceBeginsAt    = 2;
  static const int sizeOfForce      = 2;
};

/// D2Q9 lattice
template <typename T>
struct ShanChenForcedD2Q9Descriptor {
  typedef ShanChenForcedD2Q9Descriptor<T> BaseDescriptor;
  typedef ShanChenForcedExternal2Ddescriptor ExternalField;
  enum { d = 2, q = 9 };      ///< number of dimensions/distr. functions
  static const int c[q][d];   ///< lattice directions
  static const T t[q];        ///< lattice weights
  static const T invCs2;      ///< inverse square of speed of sound
  static const int vicinity=1;
};

/// 2D Descriptors for modells with variable omega and Forced Shan Chen

struct ShanChenDynOmega2Ddescriptor {
  static const int numScalars = 5;
  static const int numSpecies = 3;
  static const int velocityBeginsAt = 0;
  static const int sizeOfVelocity   = 2;
  static const int forceBeginsAt    = 2;
  static const int sizeOfForce      = 2;
  static const int omegaBeginsAt = 4;
  static const int sizeOfOmega = 1;
};

struct ShanChenDynOmega2DdescriptorBase {
  typedef ShanChenDynOmega2Ddescriptor ExternalField;
};

template <typename T>
struct ShanChenDynOmegaD2Q9Descriptor
  : public D2Q9DescriptorBase<T>, public ShanChenDynOmega2DdescriptorBase {
};

struct ShanChenDynOmegaForced2Ddescriptor {
  static const int numScalars = 7;
  static const int numSpecies = 4;
  static const int velocityBeginsAt = 0;
  static const int sizeOfVelocity   = 2;
  static const int forceBeginsAt    = 2;
  static const int sizeOfForce      = 2;
  static const int externalForceBeginsAt    = 4;
  static const int sizeOfExternalForce      = 2;
  static const int omegaBeginsAt = 6;
  static const int sizeOfOmega = 1;
};

struct ShanChenDynOmegaForced2DdescriptorBase {
  typedef ShanChenDynOmegaForced2Ddescriptor ExternalField;
};

template <typename T>
struct ShanChenDynOmegaForcedD2Q9Descriptor
  : public D2Q9DescriptorBase<T>, public ShanChenDynOmegaForced2DdescriptorBase {
};

/// 2D Descriptors for modells with variable G and Forced Shan Chen

struct ShanChenDynG2Ddescriptor {
  static const int numScalars = 6;
  static const int numSpecies = 3;
  static const int velocityBeginsAt = 0;
  static const int sizeOfVelocity   = 2;
  static const int forceBeginsAt    = 2;
  static const int sizeOfForce      = 2;
  static const int gBeginsAt = 4;
  static const int sizeOfG = 2;
};

struct ShanChenDynG2DdescriptorBase {
  typedef ShanChenDynG2Ddescriptor ExternalField;
};

template <typename T>
struct ShanChenDynGD2Q9Descriptor
  : public D2Q9DescriptorBase<T>, public ShanChenDynG2DdescriptorBase {
};

struct ShanChenDynGForced2Ddescriptor {
  static const int numScalars = 8;
  static const int numSpecies = 4;
  static const int velocityBeginsAt = 0;
  static const int sizeOfVelocity   = 2;
  static const int forceBeginsAt    = 2;
  static const int sizeOfForce      = 2;
  static const int externalForceBeginsAt    = 4;
  static const int sizeOfExternalForce      = 2;
  static const int gBeginsAt = 6;
  static const int sizeOfG = 2;
};

struct ShanChenDynGForced2DdescriptorBase {
  typedef ShanChenDynGForced2Ddescriptor ExternalField;
};

template <typename T>
struct ShanChenDynGForcedD2Q9Descriptor
  : public D2Q9DescriptorBase<T>, public ShanChenDynGForced2DdescriptorBase {
};


struct ShanChenForcedExternal3Ddescriptor {
  static const int numScalars = 6;
  static const int numSpecies = 2;
  static const int velocityBeginsAt = 0;
  static const int sizeOfVelocity   = 3;
  static const int forceBeginsAt    = 3;
  static const int sizeOfForce      = 3;
};

/// D3Q19 lattice
template <typename T> struct ShanChenForcedD3Q19Descriptor {
  typedef ShanChenForcedD3Q19Descriptor<T> BaseDescriptor;
  typedef ShanChenForcedExternal3Ddescriptor ExternalField;
  enum { d = 3, q = 19 };     ///< number of dimensions/distr. functions
  static const int c[q][d];   ///< lattice directions
  static const T t[q];        ///< lattice weights
  static const T invCs2;      ///< inverse square of speed of sound
  static const int vicinity=1;
};

/// 3D Descriptors for modells with variable omega and Forced Shan Chen

struct ShanChenDynOmega3Ddescriptor {
  static const int numScalars = 7;
  static const int numSpecies = 3;
  static const int velocityBeginsAt = 0;
  static const int sizeOfVelocity   = 3;
  static const int forceBeginsAt    = 3;
  static const int sizeOfForce      = 3;
  static const int omegaBeginsAt = 6;
  static const int sizeOfOmega = 1;
};

struct ShanChenDynOmega3DdescriptorBase {
  typedef ShanChenDynOmega3Ddescriptor ExternalField;
};

template <typename T>
struct ShanChenDynOmegaD3Q19Descriptor
  : public D3Q19DescriptorBase<T>, public ShanChenDynOmega3DdescriptorBase {
};

struct ShanChenDynOmegaForced3Ddescriptor {
  static const int numScalars = 10;
  static const int numSpecies = 4;
  static const int velocityBeginsAt = 0;
  static const int sizeOfVelocity   = 3;
  static const int forceBeginsAt    = 3;
  static const int sizeOfForce      = 3;
  static const int externalForceBeginsAt    = 6;
  static const int sizeOfExternalForce      = 3;
  static const int omegaBeginsAt = 9;
  static const int sizeOfOmega = 1;
};

struct ShanChenDynOmegaForced3DdescriptorBase {
  typedef ShanChenDynOmegaForced3Ddescriptor ExternalField;
};

template <typename T>
struct ShanChenDynOmegaForcedD3Q19Descriptor
  : public D3Q19DescriptorBase<T>, public ShanChenDynOmegaForced3DdescriptorBase {
};

}  // namespace descriptors

}  // namespace olb

#endif
