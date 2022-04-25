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

#ifndef LATTICE_PHYS_HEAT_FLUX_2D_HH
#define LATTICE_PHYS_HEAT_FLUX_2D_HH

#include <vector>
#include "utilities/omath.h"
#include <limits>

#include "latticePhysHeatFlux2D.h"
#include "dynamics/lbm.h"  // for computation of lattice rho and velocity
#include "geometry/superGeometry.h"
#include "indicator/superIndicatorF2D.h"
#include "blockBaseF2D.h"
#include "functors/genericF.h"
#include "functors/analytical/analyticalF.h"
#include "functors/analytical/indicator/indicatorF2D.h"
#include "communication/mpiManager.h"


namespace olb {

template<typename T,typename DESCRIPTOR, typename TDESCRIPTOR>
SuperLatticePhysHeatFlux2D<T,DESCRIPTOR,TDESCRIPTOR>::SuperLatticePhysHeatFlux2D(
  SuperLattice<T,TDESCRIPTOR>& sLattice, const ThermalUnitConverter<T,DESCRIPTOR,TDESCRIPTOR>& converter)
  : SuperLatticeThermalPhysF2D<T,DESCRIPTOR,TDESCRIPTOR>(sLattice, converter, 2)
{
  this->getName() = "physHeatFlux";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticePhysHeatFlux2D<T,DESCRIPTOR,TDESCRIPTOR>(this->_sLattice.getBlock(iC),
                               this->_converter));
  }
}

template <typename T, typename DESCRIPTOR, typename TDESCRIPTOR>
BlockLatticePhysHeatFlux2D<T,DESCRIPTOR,TDESCRIPTOR>::BlockLatticePhysHeatFlux2D
(BlockLattice<T,TDESCRIPTOR>& blockLattice, ThermalUnitConverter<T,DESCRIPTOR,TDESCRIPTOR> const& converter)
  : BlockLatticeThermalPhysF2D<T,DESCRIPTOR,TDESCRIPTOR>(blockLattice,converter,2),
    _temp(converter.getLatticeSpecificHeatCapacity(converter.getPhysSpecificHeatCapacity())*(converter.getLatticeThermalRelaxationTime() - 0.5) / converter.getLatticeThermalRelaxationTime())
{
  this->getName() = "physHeatFlux";

  // std::cout << _temp << " " << converter.getConversionFactorHeatFlux() << " " << 0.426 / _temp / converter.getConversionFactorHeatFlux() << std::endl;
}

template <typename T, typename DESCRIPTOR, typename TDESCRIPTOR>
bool BlockLatticePhysHeatFlux2D<T,DESCRIPTOR,TDESCRIPTOR>::operator() (T output[], const int input[])
{
  T j[2];
  const T temperature = this->_blockLattice.get( input[0], input[1] ).computeRho();
  const auto extVel = this->_blockLattice.get( input[0], input[1] ).template getField<descriptors::VELOCITY>();
  this->_blockLattice.get( input[0], input[1] ).computeJ(j);

  output[0] = this->_converter.getPhysHeatFlux((j[0] - temperature * extVel[0])*_temp);
  output[1] = this->_converter.getPhysHeatFlux((j[1] - temperature * extVel[1])*_temp);

  return true;
}

}
#endif
