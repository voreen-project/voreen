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

#ifndef LATTICE_PHYS_EXTERNAL_3D_HH
#define LATTICE_PHYS_EXTERNAL_3D_HH

#include<vector>    // for generic i/o
#include<cmath>     // for lpnorm
#include<math.h>

#include "latticePhysExternal3D.h"
#include "superBaseF3D.h"
#include "functors/analytical/indicator/indicatorBaseF3D.h"
#include "indicator/superIndicatorF3D.h"
#include "dynamics/lbHelpers.h"  // for computation of lattice rho and velocity
#include "geometry/superGeometry3D.h"
#include "blockBaseF3D.h"
#include "core/blockLatticeStructure3D.h"
#include "communication/mpiManager.h"
#include "utilities/vectorHelpers.h"

namespace olb {

template<typename T, typename DESCRIPTOR, typename FIELD>
SuperLatticePhysExternal3D<T,DESCRIPTOR,FIELD>::SuperLatticePhysExternal3D(
  SuperLattice3D<T,DESCRIPTOR>& sLattice, T convFactorToPhysUnits)
  : SuperLatticeF3D<T,DESCRIPTOR>(sLattice, 3)
{
  this->getName() = "physExtField";
  const int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(
      new BlockLatticePhysExternal3D<T, DESCRIPTOR, FIELD>(
        this->_sLattice.getBlockLattice(iC), convFactorToPhysUnits)
    );
  }
}

template <typename T, typename DESCRIPTOR, typename FIELD>
BlockLatticePhysExternal3D<T,DESCRIPTOR,FIELD>::BlockLatticePhysExternal3D(
  BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
  T convFactorToPhysUnits)
  : BlockLatticeF3D<T,DESCRIPTOR>(blockLattice, 3),
    _convFactorToPhysUnits(convFactorToPhysUnits)
{
  this->getName() = "physExtField";
}

template <typename T, typename DESCRIPTOR, typename FIELD>
bool BlockLatticePhysExternal3D<T,DESCRIPTOR,FIELD>::operator()(
  T output[], const int input[])
{
  ConstCell<T,DESCRIPTOR> cell = this->_blockLattice.get(input);
  const auto physField = _convFactorToPhysUnits * cell.template getField<FIELD>();
  for (int i = 0; i < DESCRIPTOR::d; ++i) {
    output[i] = physField[i];
  }
  return true;
}

}
#endif
