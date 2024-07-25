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

#ifndef LATTICE_PHYS_STRAIN_RATE_2D_HH
#define LATTICE_PHYS_STRAIN_RATE_2D_HH

#include <vector>
#include "utilities/omath.h"
#include <limits>

#include "latticePhysStrainRate2D.h"
#include "dynamics/lbm.h"  // for computation of lattice rho and velocity
#include "geometry/superGeometry.h"
#include "indicator/superIndicatorF2D.h"
#include "blockBaseF2D.h"
#include "functors/genericF.h"
#include "functors/analytical/analyticalF.h"
#include "functors/analytical/indicator/indicatorF2D.h"
#include "communication/mpiManager.h"


namespace olb {

template<typename T,typename DESCRIPTOR>
SuperLatticePhysStrainRate2D<T,DESCRIPTOR>::SuperLatticePhysStrainRate2D(
  SuperLattice<T,DESCRIPTOR>& sLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF2D<T,DESCRIPTOR>(sLattice, converter, 4)
{
  this->getName() = "physStrainRate";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticePhysStrainRate2D<T,DESCRIPTOR>(this->_sLattice.getBlock(iC),this->_converter));
  }
}

template <typename T, typename DESCRIPTOR>
BlockLatticePhysStrainRate2D<T,DESCRIPTOR>::BlockLatticePhysStrainRate2D
(BlockLattice<T,DESCRIPTOR>& blockLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF2D<T,DESCRIPTOR>(blockLattice,converter,4)
{
  this->getName() = "strainRate";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticePhysStrainRate2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  T rho, uTemp[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR >::n];
  this->_blockLattice.get( input[0], input[1] ).computeAllMomenta(rho, uTemp, pi);

  T omega = 1. / this->_converter.getLatticeRelaxationTime();
  T dt = this->_converter.getConversionFactorTime();

  output[0] = -pi[0]*omega*descriptors::invCs2<T,DESCRIPTOR>()/rho/2./dt;
  output[1] = -pi[1]*omega*descriptors::invCs2<T,DESCRIPTOR>()/rho/2./dt;
  output[2] = -pi[1]*omega*descriptors::invCs2<T,DESCRIPTOR>()/rho/2./dt;
  output[3] = -pi[2]*omega*descriptors::invCs2<T,DESCRIPTOR>()/rho/2./dt;

  return true;
}

}
#endif
