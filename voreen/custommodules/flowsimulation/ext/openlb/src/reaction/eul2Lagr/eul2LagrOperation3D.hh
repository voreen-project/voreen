/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2016 Thomas Henn, Davide Dapelo
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

#ifndef EUL2LAGR_OPERATION_HH
#define EUL2LAGR_OPERATION_HH

#include <string>
#include <iostream>
#include <set>
#include <vector>
#include <list>
#include <deque>

namespace olb {

template<typename T, template<typename U> class PARTICLETYPE, typename DESCRIPTOR>
Eul2LagrOperation3D<T,PARTICLETYPE,DESCRIPTOR>::Eul2LagrOperation3D (SuperLattice<T,DESCRIPTOR>& sLattice, int nParticles)
  : _nParticles(nParticles),
    _sLattice(sLattice)
{ }

template<typename T, template<typename U> class PARTICLETYPE, typename DESCRIPTOR>
void Eul2LagrOperation3D<T,PARTICLETYPE,DESCRIPTOR>::applyParticleOperation ( typename std::deque<PARTICLETYPE<T> >::iterator p,
    ParticleSystem3D<T, PARTICLETYPE>& pSys )
{
  int globic = pSys.getIGeometry();
  int locIC = this->_sLattice.getLoadBalancer().loc(globic);

  // particle's physical position
  T physPosP[] { p->getPos()[0], p->getPos()[1], p->getPos()[2] };

  // particle's dimensionless position, rounded at neighbouring voxel
  int latticeRoundedPosP[] {0, 0, 0};
  this->_sLattice.getCuboidGeometry().get(globic).getLatticeR( latticeRoundedPosP, physPosP );

  auto eul2LagrRho = _sLattice.getBlock(locIC).get (
                       latticeRoundedPosP[0], latticeRoundedPosP[1], latticeRoundedPosP[2] ).template getFieldPointer<descriptors::EUL2LAGR>();
  eul2LagrRho[0] += 1. / _nParticles;
}


}


#endif
