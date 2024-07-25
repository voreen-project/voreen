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

#ifndef LATTICE_PHYS_TEMPERATURE_3D_HH
#define LATTICE_PHYS_TEMPERATURE_3D_HH

#include<vector>    // for generic i/o
#include<cmath>     // for lpnorm
#include<math.h>

#include "latticePhysTemperature3D.h"
#include "superBaseF3D.h"
#include "functors/analytical/indicator/indicatorBaseF3D.h"
#include "indicator/superIndicatorF3D.h"
#include "dynamics/lbm.h"  // for computation of lattice rho and velocity
#include "geometry/superGeometry.h"
#include "blockBaseF3D.h"
#include "communication/mpiManager.h"
#include "utilities/vectorHelpers.h"

namespace olb {

template <typename T, typename DESCRIPTOR, typename TDESCRIPTOR>
SuperLatticePhysTemperature3D<T,DESCRIPTOR,TDESCRIPTOR>::SuperLatticePhysTemperature3D(
  SuperLattice<T,TDESCRIPTOR>& sLattice, ThermalUnitConverter<T,DESCRIPTOR,TDESCRIPTOR> const& converter)
  : SuperLatticeThermalPhysF3D<T,DESCRIPTOR,TDESCRIPTOR>(sLattice, converter, 1)
{
  this->getName() = "physTemperature";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticePhysTemperature3D<T,DESCRIPTOR,TDESCRIPTOR>(this->_sLattice.getBlock(iC), this->_converter));
  }
}

template <typename T, typename DESCRIPTOR, typename TDESCRIPTOR>
BlockLatticePhysTemperature3D<T,DESCRIPTOR,TDESCRIPTOR>::BlockLatticePhysTemperature3D
(BlockLattice<T,TDESCRIPTOR>& blockLattice, ThermalUnitConverter<T,DESCRIPTOR,TDESCRIPTOR> const& converter)
  : BlockLatticeThermalPhysF3D<T,DESCRIPTOR,TDESCRIPTOR>(blockLattice,converter,1)
{
  this->getName() = "physTemperature";
}


template <typename T, typename DESCRIPTOR, typename TDESCRIPTOR>
bool BlockLatticePhysTemperature3D<T,DESCRIPTOR,TDESCRIPTOR>::operator() (T output[], const int input[])
{
  T latticeTemperature = this->_blockLattice.get( input[0], input[1], input[2]).computeRho();
  output[0] = this->_converter.getPhysTemperature(latticeTemperature);

  return true;
}

}
#endif
