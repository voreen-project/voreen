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

#ifndef LATTICE_PHYS_DARCY_FORCE_2D_HH
#define LATTICE_PHYS_DARCY_FORCE_2D_HH

#include <vector>
#include "utilities/omath.h"
#include <limits>

#include "latticePhysDarcyForce2D.h"
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
SuperLatticePhysDarcyForce2D<T,DESCRIPTOR>::SuperLatticePhysDarcyForce2D(
  SuperLattice<T,DESCRIPTOR>& sLattice, SuperGeometry<T,2>& superGeometry,
  const int material, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF2D<T,DESCRIPTOR>(sLattice, converter, 2),
    _superGeometry(superGeometry), _material(material)
{
  this->getName() = "alphaU";
}

template<typename T,typename DESCRIPTOR>
bool SuperLatticePhysDarcyForce2D<T,DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  //  SuperLatticePhysPermeability2D<T,DESCRIPTOR> permeability(this->sLattice,this->superGeometry,this->material,this->converter);
  //  SuperLatticeVelocity2D<T,DESCRIPTOR> velocity(this->sLattice);

  //  T nu = this->converter.getCharNu();
  //  T K = permeability(input)[0];
  //  std::vector<T> u = velocity(input);

  //  std::vector<T> result(2,T());
  //  result[0] = -nu/K*u[0];
  //  result[1] = -nu/K*u[1];
  ////  result[2] = -nu/K*u[2];

  //  return result;
  return false;
}

template <typename T, typename DESCRIPTOR>
BlockLatticePhysDarcyForce2D<T,DESCRIPTOR>::BlockLatticePhysDarcyForce2D
(BlockLattice<T,DESCRIPTOR>& blockLattice, BlockGeometry<T,2>& blockGeometry,
 int material, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF2D<T,DESCRIPTOR>(blockLattice,converter,2),
    _blockGeometry(blockGeometry), _material(material)
{
  this->getName() = "alphaU";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticePhysDarcyForce2D<T,DESCRIPTOR>::operator()(T output[], const int input[])
{
  BlockLatticePhysPermeability2D<T,DESCRIPTOR> permeability(this->_blockLattice,this->_blockGeometry,this->_material,this->_converter);
  BlockLatticeVelocity2D<T,DESCRIPTOR> velocity(this->_blockLattice);

  T nu = this->_converter.getPhysViscosity();
  permeability(output,input);
  T K = output[0];
  velocity(output,input);

  output[0] *= -nu/K;
  output[1] *= -nu/K;

  return true;
}

}
#endif
