/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2019 Albert Mink, Mathias J. Krause, Lukas Baron
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

#ifndef LATTICE_PHYS_TEMPERATURE_2D_HH
#define LATTICE_PHYS_TEMPERATURE_2D_HH

#include <vector>
#include "utilities/omath.h"
#include <limits>

#include "latticePhysTemperature2D.h"
#include "dynamics/lbm.h"  // for computation of lattice rho and velocity
#include "geometry/superGeometry.h"
#include "indicator/superIndicatorF2D.h"
#include "blockBaseF2D.h"
#include "functors/genericF.h"
#include "functors/analytical/analyticalF.h"
#include "functors/analytical/indicator/indicatorF2D.h"
#include "communication/mpiManager.h"


namespace olb {

template <typename T, typename DESCRIPTOR, typename TDESCRIPTOR>
SuperLatticePhysTemperature2D<T,DESCRIPTOR,TDESCRIPTOR>::SuperLatticePhysTemperature2D(
  SuperLattice<T,TDESCRIPTOR>& sLattice, ThermalUnitConverter<T,DESCRIPTOR,TDESCRIPTOR> const& converter)
  : SuperLatticeThermalPhysF2D<T,DESCRIPTOR,TDESCRIPTOR>(sLattice, converter, 1)
{
  this->getName() = "physTemperature";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticePhysTemperature2D<T,DESCRIPTOR,TDESCRIPTOR>(this->_sLattice.getBlock(iC), this->_converter));
  }
}

template <typename T, typename DESCRIPTOR, typename TDESCRIPTOR>
BlockLatticePhysTemperature2D<T,DESCRIPTOR,TDESCRIPTOR>::BlockLatticePhysTemperature2D
(BlockLattice<T,TDESCRIPTOR>& blockLattice, ThermalUnitConverter<T,DESCRIPTOR,TDESCRIPTOR> const& converter)
  : BlockLatticeThermalPhysF2D<T,DESCRIPTOR,TDESCRIPTOR>(blockLattice,converter,1)
{
  this->getName() = "physTemperature";
}


template <typename T, typename DESCRIPTOR, typename TDESCRIPTOR>
bool BlockLatticePhysTemperature2D<T,DESCRIPTOR,TDESCRIPTOR>::operator() (T output[], const int input[])
{
  T latticeTemperature = this->_blockLattice.get( input[0], input[1] ).computeRho();
  output[0] = this->_converter.getPhysTemperature(latticeTemperature);

  return true;
}

}
#endif
