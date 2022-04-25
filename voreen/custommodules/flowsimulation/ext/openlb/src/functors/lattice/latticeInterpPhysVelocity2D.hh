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

#ifndef LATTICE_INTERP_PHYS_VELOCITY_2D_HH
#define LATTICE_INTERP_PHYS_VELOCITY_2D_HH

#include<vector>    // for generic i/o
#include<cmath>     // for lpnorm
#include<math.h>

#include "latticeInterpPhysVelocity2D.h"
#include "superBaseF2D.h"
#include "functors/analytical/indicator/indicatorBaseF2D.h"
#include "indicator/superIndicatorF2D.h"
#include "dynamics/lbm.h"  // for computation of lattice rho and velocity
#include "geometry/superGeometry.h"
#include "blockBaseF2D.h"
#include "communication/mpiManager.h"
#include "utilities/vectorHelpers.h"

namespace olb {

template<typename T, typename DESCRIPTOR>
SuperLatticeInterpPhysVelocity2D<T, DESCRIPTOR>::SuperLatticeInterpPhysVelocity2D(
  SuperLattice<T, DESCRIPTOR>& sLattice, UnitConverter<T,DESCRIPTOR> const& converter)
  : SuperLatticePhysF2D<T, DESCRIPTOR>(sLattice, converter, 2)
{
  this->getName() = "InterpVelocity";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int lociC = 0; lociC < maxC; lociC++) {
    int globiC = this->_sLattice.getLoadBalancer().glob(lociC);

    this->_blockF.emplace_back(
      new BlockLatticeInterpPhysVelocity2D<T, DESCRIPTOR>(
        sLattice.getBlock(lociC),
        converter,
        &sLattice.getCuboidGeometry().get(globiC))
    );
  }
}

template<typename T, typename DESCRIPTOR>
bool SuperLatticeInterpPhysVelocity2D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  return false;
}

template<typename T, typename DESCRIPTOR>
void SuperLatticeInterpPhysVelocity2D<T, DESCRIPTOR>::operator()(T output[],
    const T input[], const int globiC)
{
  if (this->_sLattice.getLoadBalancer().isLocal(globiC)) {
    static_cast<BlockLatticeInterpPhysVelocity2D<T, DESCRIPTOR>*>(
      this->_blockF[this->_sLattice.getLoadBalancer().loc(globiC)].get()
    )->operator()(output, input);
  }
}

template<typename T, typename DESCRIPTOR>
BlockLatticeInterpPhysVelocity2D<T, DESCRIPTOR>::BlockLatticeInterpPhysVelocity2D(
  BlockLattice<T, DESCRIPTOR>& blockLattice, UnitConverter<T,DESCRIPTOR> const& converter, Cuboid2D<T>* c)
  : BlockLatticePhysF2D<T, DESCRIPTOR>(blockLattice, converter, 2),
    _cuboid(c)
{
  this->getName() = "BlockLatticeInterpVelocity2D";
}

template<typename T, typename DESCRIPTOR>
BlockLatticeInterpPhysVelocity2D<T, DESCRIPTOR>::BlockLatticeInterpPhysVelocity2D(
  const BlockLatticeInterpPhysVelocity2D<T, DESCRIPTOR>& rhs) :
  BlockLatticePhysF2D<T, DESCRIPTOR>(rhs._blockLattice, rhs._converter, 2),
  _cuboid(rhs._cuboid)
{
}

template<typename T, typename DESCRIPTOR>
void BlockLatticeInterpPhysVelocity2D<T, DESCRIPTOR>::operator()(T output[2], const T input[2])
{
  T u[2], rho, volume;
  T d[2], e[2];
  int latIntPos[2] = {0};
  T latPhysPos[2] = {T()};
  _cuboid->getFloorLatticeR(latIntPos, &input[0]);
  _cuboid->getPhysR(latPhysPos, latIntPos);

  T deltaRinv = 1. / _cuboid->getDeltaR();
  d[0] = (input[0] - latPhysPos[0]) * deltaRinv;
  d[1] = (input[1] - latPhysPos[1]) * deltaRinv;

  e[0] = 1. - d[0];
  e[1] = 1. - d[1];

  this->_blockLattice.get(latIntPos[0], latIntPos[1]).computeRhoU(rho, u);
  volume = e[0] * e[1];
  output[0] = u[0] * volume;
  output[1] = u[1] * volume;

  this->_blockLattice.get(latIntPos[0], latIntPos[1] + 1).computeRhoU(rho, u);
  volume = e[0] * d[1];
  output[0] += u[0] * volume;
  output[1] += u[1] * volume;

  this->_blockLattice.get(latIntPos[0] + 1, latIntPos[1]).computeRhoU(rho, u);
  volume = d[0] * e[1];
  output[0] += u[0] * volume;
  output[1] += u[1] * volume;

  this->_blockLattice.get(latIntPos[0] + 1, latIntPos[1] + 1).computeRhoU(rho, u);
  volume = d[0] * d[1];
  output[0] += u[0] * volume;
  output[1] += u[1] * volume;

  output[0] = this->_converter.getPhysVelocity(output[0]);
  output[1] = this->_converter.getPhysVelocity(output[1]);
}

}
#endif
