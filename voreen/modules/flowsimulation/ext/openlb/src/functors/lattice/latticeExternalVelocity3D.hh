/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Lukas Baron, Tim Dornieden, Mathias J. Krause,
 *  Albert Mink
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

#ifndef LATTICE_EXTERNAL_VELOCITY_3D_HH
#define LATTICE_EXTERNAL_VELOCITY_3D_HH

#include<vector>    // for generic i/o
#include<cmath>     // for lpnorm
#include<math.h>

#include "latticeExternalVelocity3D.h"
#include "superBaseF3D.h"
#include "functors/analytical/indicator/indicatorBaseF3D.h"
#include "indicator/superIndicatorF3D.h"
#include "dynamics/lbm.h"  // for computation of lattice rho and velocity
#include "geometry/superGeometry.h"
#include "blockBaseF3D.h"
#include "communication/mpiManager.h"
#include "utilities/vectorHelpers.h"

namespace olb {

template<typename T, typename DESCRIPTOR>
SuperLatticeExternalVelocity3D<T, DESCRIPTOR>::SuperLatticeExternalVelocity3D(
  SuperLattice<T, DESCRIPTOR>& sLattice)
  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, 3)
{
  this->getName() = "externalVelocity";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticeExternalVelocity3D<T, DESCRIPTOR>(this->_sLattice.getBlock(iC)));
  }
}

template<typename T, typename DESCRIPTOR>
BlockLatticeExternalVelocity3D<T, DESCRIPTOR>::BlockLatticeExternalVelocity3D(
  BlockLattice<T, DESCRIPTOR>& blockLattice)
  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, 3)
{
  this->getName() = "externalVelocity";
}

template<typename T, typename DESCRIPTOR>
bool BlockLatticeExternalVelocity3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  auto ExtVel = this->_blockLattice.get(input[0], input[1], input[2]).template getField<descriptors::VELOCITY>();
  for (int iVel=0; iVel<DESCRIPTOR::d; ++iVel) {
    output[iVel] = ExtVel[iVel];
  }
  return true;
}

}
#endif
