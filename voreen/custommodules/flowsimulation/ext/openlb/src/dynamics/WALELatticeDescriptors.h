/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2015 Mathias J. Krause, Patrick Nathen
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

#ifndef WALE_LATTICE_DESCRIPTOR_H
#define WALE_LATTICE_DESCRIPTOR_H

#include "dynamics/latticeDescriptors.h"


namespace olb {

namespace descriptors {



/////////////////////////////////////////////////////////////////////////////////
// 3D Descriptors for flow with Wall Adaptive Local Eddy Viscosity (WALE)

struct WALE3dDescriptor {
  static const int numScalars = 10;
  static const int numSpecies = 2;
  static const int EffectiveOmegaIsAt = 0;
  static const int sizeOfEffectiveOmega = 1;
  static const int veloGradIsAt = 1;
  static const int sizeOfVeloGrad = 9;
};

struct WALE3dDescriptorBase {
  typedef WALE3dDescriptor ExternalField;
};

template <typename T> struct WALED3Q19Descriptor
  : public D3Q19DescriptorBase<T>, public WALE3dDescriptorBase {
};




} // namespace descriptors

} // namespace olb

#endif
