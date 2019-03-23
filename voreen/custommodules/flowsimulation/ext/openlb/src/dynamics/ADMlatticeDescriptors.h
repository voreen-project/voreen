/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Mathias J. Krause, Jonas Latt
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

#ifndef ADM_BGK_DYNAMICS_DESCRIPTOR_H
#define ADM_BGK_DYNAMICS_DESCRIPTOR_H

#include "dynamics/latticeDescriptors.h"
//#include <cmath>


namespace olb {

namespace descriptors {

// 2D Descriptors for ADM


struct ADM2dDescriptor {
  static const int numScalars = 3;
  static const int numSpecies = 2;
  static const int filRhoIsAt      = 0;
  static const int localFilVelXBeginsAt = 1;
  static const int sizeOfFilVelX   = 1;
  static const int localFilVelYBeginsAt = 2;
  static const int sizeOfFilVelY   = 1;
};

struct ADM2dDescriptorBase {
  typedef ADM2dDescriptor ExternalField;
};

template <typename T> struct ADMD2Q9Descriptor
  : public D2Q9DescriptorBase<T>, public ADM2dDescriptorBase {
};




////////////////////////////////////////////////////////////////////////////////
// extended descriptor for ADM

struct ADM3dDescriptor {
  static const int numScalars = 4;
  static const int numSpecies = 2;
  static const int filRhoIsAt      = 0;
  static const int localFilVelXBeginsAt = 1;
  static const int sizeOfFilVelX   = 1;
  static const int localFilVelYBeginsAt = 2;
  static const int sizeOfFilVelY   = 1;
  static const int localFilVelZBeginsAt = 3;
  static const int sizeOfFilVelZ   = 1;
};

struct ADM3dDescriptorBase {
  typedef ADM3dDescriptor ExternalField;
};

template <typename T> struct ADMD3Q19Descriptor
  : public D3Q19DescriptorBase<T>, public ADM3dDescriptorBase {
};


////////////////////////////////////////
/// ADM Descriptors for forced fields

//// Forced 2D ADM scheme
struct ForcedADM2dDescriptor {
  static const int numScalars = 5;
  static const int numSpecies = 3;
  static const int forceBeginsAt = 0;
  static const int sizeOfForce   = 2;
  static const int filRhoIsAt      = 2;
  static const int localFilVelXBeginsAt = 3;
  static const int sizeOfFilVelX   = 1;
  static const int localFilVelYBeginsAt = 4;
  static const int sizeOfFilVelY   = 1;
};

struct ForcedADM2dDescriptorBase {
  typedef ForcedADM2dDescriptor ExternalField;
};

template <typename T> struct ForcedADMD2Q9Descriptor
  : public D2Q9DescriptorBase<T>, public ForcedADM2dDescriptorBase {
};


//// Forced 3D ADM scheme

struct ForcedADM3dDescriptor {
  static const int numScalars = 7;
  static const int numSpecies = 3;
  static const int forceBeginsAt = 0;
  static const int sizeOfForce   = 3;
  static const int filRhoIsAt      = 3;
  static const int localFilVelXBeginsAt = 4;
  static const int sizeOfFilVelX   = 1;
  static const int localFilVelYBeginsAt = 5;
  static const int sizeOfFilVelY   = 1;
  static const int localFilVelZBeginsAt = 6;
  static const int sizeOfFilVelZ   = 1;
};

struct ForcedADM3dDescriptorBase {
  typedef ForcedADM3dDescriptor ExternalField;
};

template <typename T> struct ForcedADMD3Q19Descriptor
  : public D3Q19DescriptorBase<T>, public ForcedADM3dDescriptorBase {
};


//// Forced adapted 3D ADM scheme

struct ForcedAdaptiveADM3dDescriptor {
  static const int numScalars = 10;
  static const int numSpecies = 5;
  static const int forceBeginsAt = 0;
  static const int sizeOfForce   = 3;
  static const int filRhoIsAt      = 3;
  static const int localFilVelXBeginsAt = 4;
  static const int sizeOfFilVelX   = 1;
  static const int localFilVelYBeginsAt = 5;
  static const int sizeOfFilVelY   = 1;
  static const int localFilVelZBeginsAt = 6;
  static const int sizeOfFilVelZ   = 1;
  static const int localAvDissBeginsAt   = 7;
  static const int sizeOfAvDiss   = 1;
  static const int localAvTKEBeginsAt   = 8;
  static const int sizeOfAvTKE   = 1;
  static const int localSigmaADMBeginsAt   = 9;
  static const int sizeOfSigmaADM   = 1;
  static const int localNuEddyBeginsAt   = 10;
  static const int sizeOfNuEddy   = 1;
  static const int tauWIsAt   = 11;
  static const int sizeOfTauW   = 1;
};

struct ForcedAdaptiveADM3dDescriptorBase {
  typedef ForcedAdaptiveADM3dDescriptor ExternalField;
};

template <typename T> struct ForcedAdaptiveADMD3Q19Descriptor
  : public D3Q19DescriptorBase<T>, public ForcedAdaptiveADM3dDescriptorBase {
};



} // namespace descriptors

} // namespace olb

#endif
