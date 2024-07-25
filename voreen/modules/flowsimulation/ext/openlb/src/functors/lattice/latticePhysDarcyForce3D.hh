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

#ifndef LATTICE_PHYS_DARCY_FORCE_3D_HH
#define LATTICE_PHYS_DARCY_FORCE_3D_HH

#include<vector>    // for generic i/o
#include<cmath>     // for lpnorm
#include<math.h>

#include "latticePhysDarcyForce3D.h"
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
SuperLatticePhysDarcyForce3D<T, DESCRIPTOR>::SuperLatticePhysDarcyForce3D(
  SuperLattice<T, DESCRIPTOR>& sLattice, SuperGeometry<T,3>& superGeometry,
  const int material, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF3D<T, DESCRIPTOR>(sLattice, converter, 3),
    _superGeometry(superGeometry), _material(material)
{
  this->getName() = "alphaU";
}

template<typename T, typename DESCRIPTOR>
bool SuperLatticePhysDarcyForce3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  SuperLatticePhysPermeability3D<T, DESCRIPTOR> permeability(this->_sLattice, this->_converter);
  SuperLatticeVelocity3D<T, DESCRIPTOR> velocity(this->_sLattice);

  T nu = this->_converter.getPhysViscosity();
  T K;
  T u[velocity.getTargetDim()];
  permeability(&K,input);
  velocity(u,input);

  output[0] = -nu / K * u[0];
  output[1] = -nu / K * u[1];
  output[2] = -nu / K * u[2];

  return true;
}

template<typename T, typename DESCRIPTOR>
BlockLatticePhysDarcyForce3D<T, DESCRIPTOR>::BlockLatticePhysDarcyForce3D(
  BlockLattice<T, DESCRIPTOR>& blockLattice, BlockGeometry<T,3>& blockGeometry,
  int material, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF3D<T, DESCRIPTOR>(blockLattice, converter, 3),
    _blockGeometry(blockGeometry),
    _material(material)
{
  this->getName() = "alphaU";
}

template<typename T, typename DESCRIPTOR>
bool BlockLatticePhysDarcyForce3D<T, DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  BlockLatticePhysPermeability3D<T, DESCRIPTOR> permeability(this->_blockLattice, this->_converter);
  BlockLatticeVelocity3D<T, DESCRIPTOR> velocity(this->_blockLattice);

  T nu = this->_converter.getPhysViscosity();
  permeability(output,input);
  T K = output[0];
  velocity(output,input);

  output[0] *= -nu / K;
  output[1] *= -nu / K;
  output[2] *= -nu / K;

  return true;
}

}
#endif
