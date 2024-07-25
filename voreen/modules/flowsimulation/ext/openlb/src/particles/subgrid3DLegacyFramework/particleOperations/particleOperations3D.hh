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

#ifndef PARTICLE_OPERATIONS_3D_HH
#define PARTICLE_OPERATIONS_3D_HH

#include <string>
#include <iostream>
#include <set>
#include <vector>
#include <list>
#include <deque>

namespace olb {

template<typename T, template<typename U> class PARTICLETYPE>
RandomWalkOperation3D<T,PARTICLETYPE>::RandomWalkOperation3D(T dT, T diffusivity)
  : _dT(dT),
    _randomNormal( 0, util::sqrt(2.*diffusivity*dT) )
{ }

template<typename T, template<typename U> class PARTICLETYPE>
void RandomWalkOperation3D<T,PARTICLETYPE>::applyParticleOperation ( typename std::deque<PARTICLETYPE<T> >::iterator p,
    ParticleSystem3D<T, PARTICLETYPE>& pSys )
{
  T posOldArray[] { p->getPos()[0], p->getPos()[1], p->getPos()[2] };
  T posIncrementArray[] {0., 0., 0.};
  _randomNormal(posIncrementArray, posOldArray);
  std::vector<T> posNew { posOldArray[0]+posIncrementArray[0], posOldArray[1]+posIncrementArray[1], posOldArray[2]+posIncrementArray[2] };
  p->setPos(posNew);
  p->setVel({ (posNew[0]-posOldArray[0])/_dT, (posNew[1]-posOldArray[1])/_dT, (posNew[2]-posOldArray[2])/_dT });
}


template<typename T, template<typename U> class PARTICLETYPE>
RandomTruncatedWalkOperation3D<T,PARTICLETYPE>::RandomTruncatedWalkOperation3D(T dT, T diffusivity, T n)
  : _dT(dT),
    _randomNormal( 0, util::sqrt(2.*diffusivity*dT), n )
{ }

template<typename T, template<typename U> class PARTICLETYPE>
void RandomTruncatedWalkOperation3D<T,PARTICLETYPE>::applyParticleOperation ( typename std::deque<PARTICLETYPE<T> >::iterator p,
    ParticleSystem3D<T, PARTICLETYPE>& pSys )
{
  T posOldArray[] { p->getPos()[0], p->getPos()[1], p->getPos()[2] };
  T posIncrementArray0[] {0., 0., 0.};
  T posIncrementArray1[] {0., 0., 0.};
  T posIncrementArray2[] {0., 0., 0.};
  _randomNormal(posIncrementArray0, posOldArray);
  _randomNormal(posIncrementArray1, posOldArray);
  _randomNormal(posIncrementArray2, posOldArray);
  std::vector<T> posNew { posOldArray[0]+posIncrementArray0[0], posOldArray[1]+posIncrementArray1[0], posOldArray[2]+posIncrementArray2[0] };
  p->setPos(posNew);
  p->setVel({ (posNew[0]-posOldArray[0])/_dT, (posNew[1]-posOldArray[1])/_dT, (posNew[2]-posOldArray[2])/_dT });
}


template<typename T, template<typename U> class PARTICLETYPE, typename DESCRIPTOR>
PassiveAdvectionOperation3D<T,PARTICLETYPE,DESCRIPTOR>::PassiveAdvectionOperation3D (
  SuperLattice<T,DESCRIPTOR>& sLattice, UnitConverter<T,DESCRIPTOR>& converter )
  : _dT(converter.getPhysDeltaT()),
    _getVel(sLattice, converter)
{ }

template<typename T, template<typename U> class PARTICLETYPE, typename DESCRIPTOR>
void PassiveAdvectionOperation3D<T,PARTICLETYPE,DESCRIPTOR>::applyParticleOperation ( typename std::deque<PARTICLETYPE<T> >::iterator p,
    ParticleSystem3D<T, PARTICLETYPE>& pSys )
{
  T fluidVel[] {0., 0., 0.};
  _getVel(fluidVel, &p->getPos()[0], p->getCuboid());

  T oldPos[] {p->getPos()[0], p->getPos()[1], p->getPos()[2]};
  p->setPos({ oldPos[0]+fluidVel[0]*_dT, oldPos[1]+fluidVel[1]*_dT, oldPos[2]+fluidVel[2]*_dT });
  std::vector<T> oldVel { p->getVel() };
  p->setVel({ oldVel[0]+fluidVel[0], oldVel[1]+fluidVel[1], oldVel[2]+fluidVel[2] });
}

}



#endif
