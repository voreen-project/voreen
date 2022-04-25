/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2019 Albert Mink, Mathias J. Krause, Lukas Baron, Davide Dapelo
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

#ifndef LATTICE_EXTERNAL_SCALAR_FIELD_2D_HH
#define LATTICE_EXTERNAL_SCALAR_FIELD_2D_HH

#include <vector>
#include "utilities/omath.h"
#include <limits>

#include "latticeExternalScalarField2D.h"
#include "dynamics/lbm.h"  // for computation of lattice rho and velocity
#include "geometry/superGeometry.h"
#include "indicator/superIndicatorF2D.h"
#include "blockBaseF2D.h"
#include "functors/genericF.h"
#include "functors/analytical/analyticalF.h"
#include "functors/analytical/indicator/indicatorF2D.h"
#include "communication/mpiManager.h"


namespace olb {

template<typename T,typename DESCRIPTOR, typename FIELD>
SuperLatticeExternalScalarField2D<T,DESCRIPTOR,FIELD>::SuperLatticeExternalScalarField2D(
  SuperLattice<T,DESCRIPTOR>& sLattice)
  : SuperLatticeF2D<T,DESCRIPTOR>(sLattice, 1)
{
  this->getName() = "externalScalarField";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back( new BlockLatticeExternalScalarField2D<T,DESCRIPTOR,FIELD>(this->_sLattice.getBlock(iC)) );
  }
}

template <typename T, typename DESCRIPTOR, typename FIELD>
BlockLatticeExternalScalarField2D<T,DESCRIPTOR,FIELD>::BlockLatticeExternalScalarField2D
(BlockLattice<T,DESCRIPTOR>& blockLattice)
  : BlockLatticeF2D<T,DESCRIPTOR>(blockLattice,1)
{
  this->getName() = "externalScalarField";
}


template <typename T, typename DESCRIPTOR, typename FIELD>
bool BlockLatticeExternalScalarField2D<T,DESCRIPTOR,FIELD>::operator() (T output[], const int input[])
{
  output[0] = static_cast<T>(this->_blockLattice.get( input[0], input[1] ).template getField<FIELD>());
  return true;
}

}
#endif
