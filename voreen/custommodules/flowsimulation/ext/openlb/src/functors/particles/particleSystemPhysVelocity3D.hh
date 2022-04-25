/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Nicolas Hafen, Mathias J. Krause
 *  
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

#ifndef PARTICLE_SYSTEM_PHYS_VELOCITY_3D_HH
#define PARTICLE_SYSTEM_PHYS_VELOCITY_3D_HH

#include<vector>    // for generic i/o
#include<cmath>     // for lpnorm
#include<math.h>

//#include "latticePhysVelocity3D.h"
#include "superParticleBaseF3D.h"
//#include "functors/analytical/indicator/indicatorBaseF3D.h"
//#include "indicator/superIndicatorF3D.h"
//#include "dynamics/lbHelpers.h"  // for computation of lattice rho and velocity
//#include "geometry/superGeometry.h"
#include "blockParticleBaseF3D.h"
//#include "core/blockLatticeStructure3D.h"
//#include "communication/mpiManager.h"
//#include "utilities/vectorHelpers.h"

namespace olb {

template<typename T, typename PARTICLETYPE>
SuperParticleSystemPhysVelocity3D<T, PARTICLETYPE>::SuperParticleSystemPhysVelocity3D(
  SuperParticleSystem3D<T, PARTICLETYPE>& sParticleSystem, bool print)
  : SuperParticleSystemF3D<T, PARTICLETYPE>(sParticleSystem, 3), _print(print)
{
  this->getName() = "physVelocity";

  const int maxC = this->_sParticleSystem.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(
      new BlockParticleSystemPhysVelocity3D<T, PARTICLETYPE>(
        this->_sParticleSystem.getParticleSystem(iC),
        this->_sParticleSystem.getOverlap(),
        _print)
    );
  }
}



template<typename T, typename PARTICLETYPE>
BlockParticleSystemPhysVelocity3D<T, PARTICLETYPE>::BlockParticleSystemPhysVelocity3D(
  ParticleSystem3D<T, PARTICLETYPE>& particleSystem,
  int overlap,
  bool print)
  : BlockParticleSystemF3D<T, PARTICLETYPE>(particleSystem, 3),
    _overlap(overlap),
    _print(print)
{
  this->getName() = "physVelocity";
}

template<typename T, typename PARTICLETYPE>
bool BlockParticleSystemPhysVelocity3D<T, PARTICLETYPE>::operator()(T output[], const int input[])
{
  if (_print) {
    std::cout << input[0] << " | " << singleton::mpi().getRank() << std::endl;
  }

  //Check wether iterator exceded particle number for an individual core
  if (input[0] < this->_particleSystem.size()){ 
    Vector<T,3> velocity = this->_particleSystem[input[0]].getMobility().getVelocity(); 
    output[0] = velocity[0];
    output[1] = velocity[1];
    output[2] = velocity[2];
    return true;
  } else {
    return false;   //Return false if input[0] is too large for particular core
  }
}

}
#endif
