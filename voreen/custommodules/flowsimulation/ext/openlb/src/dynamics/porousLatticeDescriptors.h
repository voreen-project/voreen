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

#ifndef POROUS_LATTICE_DESCRIPTOR_H
#define POROUS_LATTICE_DESCRIPTOR_H

#include "latticeDescriptors.h"


namespace olb {

namespace descriptors {

// 2D Descriptors for flow through porous media

struct Porous2dDescriptor {
  static const int numScalars = 1;
  static const int numSpecies = 1;
  static const int porosityIsAt = 0;
};

struct Porous2dDescriptorBase {
  typedef Porous2dDescriptor ExternalField;
};

template <typename T> struct PorousD2Q9Descriptor
  : public D2Q9DescriptorBase<T>, public Porous2dDescriptorBase {
};

////////////////////////////////////////////////////////////////////////////////
// extended descriptor for drag computation - 2D

struct ExtendedPorous2dDescriptor {
  static const int numScalars = 3;
  static const int numSpecies = 2;
  static const int porosityIsAt      = 0;
  static const int localDragBeginsAt = 1;
  static const int sizeOfLocalDrag   = 2;
};

struct ExtendedPorous2dDescriptorBase {
  typedef ExtendedPorous2dDescriptor ExternalField;
};

template <typename T> struct ExtendedPorousD2Q9Descriptor
  : public D2Q9DescriptorBase<T>, public ExtendedPorous2dDescriptorBase {
};

////////////////////////////////////////////////////////////////////////////////
// extended descriptor for porous particles - 2D

struct PorousParticle2dDescriptor {
  static const int numScalars = 4;
  static const int numSpecies = 3;

  static const int porosityIsAt      = 0;
  static const int velNumerator      = 1;
  static const int velDenominator    = 3;
  //  static const int deltaMomentum     = 4;

  static const int sizeOfPorosity    = 1;
  static const int sizeOfVelNum      = 2;
  static const int sizeOfVelDenom    = 1;
  //  static const int sizeOfDeltaMomentum    = 2;
};

struct PorousParticle2dDescriptorBase {
  typedef PorousParticle2dDescriptor ExternalField;
};

template <typename T> struct PorousParticleD2Q9Descriptor
  : public D2Q9DescriptorBase<T>, public PorousParticle2dDescriptorBase {
};


/////////////////////////////////////////////////////////////////////////////////
// 3D Descriptors for flow through porous media

struct Porous3dDescriptor {
  static const int numScalars = 1;
  static const int numSpecies = 1;
  static const int porosityIsAt = 0;
};

struct Porous3dDescriptorBase {
  typedef Porous3dDescriptor ExternalField;
};

template <typename T> struct PorousD3Q19Descriptor
  : public D3Q19DescriptorBase<T>, public Porous3dDescriptorBase {
};



////////////////////////////////////////////////////////////////////////////////
// extended descriptor for drag computation - 3D

struct ExtendedPorous3dDescriptor {
  static const int numScalars = 4;
  static const int numSpecies = 2;
  static const int porosityIsAt      = 0;
  static const int localDragBeginsAt = 1;
  static const int sizeOfLocalDrag   = 3;
};

struct ExtendedPorous3dDescriptorBase {
  typedef ExtendedPorous3dDescriptor ExternalField;
};

template <typename T> struct ExtendedPorousD3Q19Descriptor
  : public D3Q19DescriptorBase<T>, public ExtendedPorous3dDescriptorBase {
};

////////////////////////////////////////////////////////////////////////////////
// extended descriptor for porous particles - 3D

struct PorousParticle3dDescriptor {
  static const int numScalars = 5;
  static const int numSpecies = 3;

  static const int porosityIsAt      = 0;
  static const int velNumerator      = 1;
  static const int velDenominator    = 4;

  static const int sizeOfPorosity    = 1;
  static const int sizeOfVelNum      = 3;
  static const int sizeOfVelDenom    = 1;
};

struct PorousParticle3dDescriptorBase {
  typedef PorousParticle3dDescriptor ExternalField;
};

template <typename T> struct PorousParticleD3Q19Descriptor
  : public D3Q19DescriptorBase<T>, public PorousParticle3dDescriptorBase {
};


} // namespace descriptors

} // namespace olb

#endif
