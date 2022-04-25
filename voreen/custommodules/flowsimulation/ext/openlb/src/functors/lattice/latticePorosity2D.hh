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

#ifndef LATTICE_POROSITY_2D_HH
#define LATTICE_POROSITY_2D_HH

#include <vector>
#include "utilities/omath.h"
#include <limits>

#include "latticePorosity2D.h"
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
SuperLatticePorosity2D<T,DESCRIPTOR>::SuperLatticePorosity2D(
  SuperLattice<T,DESCRIPTOR>& sLattice, SuperGeometry<T,2>& superGeometry,
  const int material, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF2D<T,DESCRIPTOR>(sLattice, converter, 1),
    _superGeometry(superGeometry), _material(material)
{
  this->getName() = "porosity";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new BlockLatticePorosity2D<T,DESCRIPTOR>(
                                 this->_sLattice.getBlock(iC),
                                 this->_superGeometry.getBlockGeometry(iC),
                                 _material,
                                 converter));
  }
}

template <typename T, typename DESCRIPTOR>
BlockLatticePorosity2D<T,DESCRIPTOR>::BlockLatticePorosity2D
(BlockLattice<T,DESCRIPTOR>& blockLattice, BlockGeometry<T,2>& blockGeometry,
 int material, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF2D<T,DESCRIPTOR>(blockLattice,converter,1), _blockGeometry(blockGeometry), _material(material)
{
  this->getName() = "porosity";
}

// under construction
template <typename T, typename DESCRIPTOR>
bool BlockLatticePorosity2D<T,DESCRIPTOR>::operator()(T output[], const int input[])
{
  return false;
}

}
#endif
