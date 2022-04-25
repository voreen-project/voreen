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

#ifndef LATTICE_PHYS_EXTERNAL_PARTICLE_VELOCITY_3D_HH
#define LATTICE_PHYS_EXTERNAL_PARTICLE_VELOCITY_3D_HH

#include<vector>    // for generic i/o
#include<cmath>     // for lpnorm
#include<math.h>

#include "latticePhysExternalParticleVelocity3D.h"
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
SuperLatticePhysExternalParticleVelocity3D<T, DESCRIPTOR>::SuperLatticePhysExternalParticleVelocity3D(
  SuperLattice<T, DESCRIPTOR>& sLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : SuperLatticePhysF3D<T, DESCRIPTOR>(sLattice, converter, 2)
{
  this->getName() = "ExtPartVelField";

  for (int iC = 0; iC < sLattice.getLoadBalancer().size(); ++iC) {
    this->_blockF.emplace_back(
      new BlockLatticePhysExternalParticleVelocity3D<T, DESCRIPTOR>(
        sLattice.getBlock(iC),
        converter)
    );
  }
}

template<typename T, typename DESCRIPTOR>
bool SuperLatticePhysExternalParticleVelocity3D<T, DESCRIPTOR>::operator()(
  T output[], const int input[])
{
  auto& load = this->_sLattice.getLoadBalancer();
  const int& globIC = input[0];

  if (load.rank(globIC) == singleton::mpi().getRank()) {
    return this->getBlockF(load.loc(globIC))(output, input+1);
  }
  else {
    return false;
  }
}

template <typename T, typename DESCRIPTOR>
BlockLatticePhysExternalParticleVelocity3D<T,DESCRIPTOR>::BlockLatticePhysExternalParticleVelocity3D
(BlockLattice<T,DESCRIPTOR>& blockLattice, const UnitConverter<T,DESCRIPTOR>& converter)
  : BlockLatticePhysF3D<T,DESCRIPTOR>(blockLattice,converter,2)
{
  this->getName() = "ExtParticleVelocityField";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticePhysExternalParticleVelocity3D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  auto velocity_numerator   = this->blockLattice.get(input).template getFieldPointer<descriptors::VELOCITY_NUMERATOR>();
  auto velocity_denominator = this->blockLattice.get(input).template getFieldPointer<descriptors::VELOCITY_DENOMINATOR>();

  if (velocity_denominator[0] > std::numeric_limits<T>::epsilon()) {
    output[0]=this->_converter.getPhysVelocity(velocity_numerator[0]/velocity_denominator[0]);
    output[1]=this->_converter.getPhysVelocity(velocity_numerator[1]/velocity_denominator[0]);
    output[2]=this->_converter.getPhysVelocity(velocity_numerator[2]/velocity_denominator[0]);
    return true;
  }
  output[0]=this->_converter.getPhysVelocity(velocity_numerator[0]);
  output[1]=this->_converter.getPhysVelocity(velocity_numerator[1]);
  output[2]=this->_converter.getPhysVelocity(velocity_numerator[2]);
  return true;
}

}
#endif
