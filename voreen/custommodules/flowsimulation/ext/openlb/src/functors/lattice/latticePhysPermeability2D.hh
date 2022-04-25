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

#ifndef LATTICE_PHYS_PERMEABILITY_2D_HH
#define LATTICE_PHYS_PERMEABILITY_2D_HH

#include <vector>
#include "utilities/omath.h"
#include <limits>

#include "latticePhysPermeability2D.h"
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
SuperLatticePhysPermeability2D<T,DESCRIPTOR>::SuperLatticePhysPermeability2D(
  SuperLattice<T,DESCRIPTOR>& sLattice, SuperGeometry<T,2>& superGeometry,
  const int material, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF2D<T,DESCRIPTOR>(sLattice, converter, 1),
    _superGeometry(superGeometry), _material(material)
{
  this->getName() = "permeability";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back( new BlockLatticePhysPermeability2D<T,DESCRIPTOR>(
                                  this->_sLattice.getBlock(iC),
                                  this->_superGeometry.getBlockGeometry(iC),
                                  _material,
                                  this->getConverter() ) );
  }
}

template<typename T,typename DESCRIPTOR>
bool SuperLatticePhysPermeability2D<T,DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  //  int globIC = input[0];
  //  int locix= input[1];
  //  int lociy= input[2];
  ////  int lociz= input[3];

  //  T* value = new T[1];
  //  this->sLattice.getBlock(this->sLattice.getLoadBalancer().loc(globIC)).get(locix, lociy/*, lociz*/).computeField(0,1,value);
  //  std::vector<T> result(1,this->converter.physPermeability(value[0]));
  //  delete value;
  //  if (!(result[0]<42)&&!(result[0]>42)&&!(result[0]==42)) result[0] = 999999;
  //  if (isinf(result[0])) result[0] = 1e6;
  //  return result;
  return false;
}

template <typename T, typename DESCRIPTOR>
BlockLatticePhysPermeability2D<T,DESCRIPTOR>::BlockLatticePhysPermeability2D
(BlockLattice<T,DESCRIPTOR>& blockLattice, BlockGeometry<T,2>& blockGeometry,
 int material, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF2D<T,DESCRIPTOR>(blockLattice,converter,1), _blockGeometry(blockGeometry), _material(material)
{
  this->getName() = "permeability";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticePhysPermeability2D<T,DESCRIPTOR>::operator()(T output[], const int input[])
{
  return false;
}

}
#endif
