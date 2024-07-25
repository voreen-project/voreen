/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2016 Thomas Henn, Mathias J. Krause, Davide Dapelo
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

#ifndef PARTICLE_OPERATIONS_3D_H
#define PARTICLE_OPERATIONS_3D_H

#include <set>
#include <vector>
#include <list>
#include <deque>
#include <string>
#include <iostream>

namespace olb {

template<typename T> class Particle3D;
template<typename T, template<typename U> class PARTICLETYPE> class ParticleSystem3D;

template<typename T, template<typename U> class PARTICLETYPE>
class ParticleOperation3D {
protected:
  ParticleOperation3D() {}
public:
  virtual void applyParticleOperation(typename std::deque<PARTICLETYPE<T> >::iterator p, ParticleSystem3D<T, PARTICLETYPE>& pSys) =0;
};

template<typename T, template<typename U> class PARTICLETYPE>
class RandomWalkOperation3D final : public ParticleOperation3D<T,PARTICLETYPE> {
public:
  RandomWalkOperation3D(T dT, T diffusivity);
  virtual void applyParticleOperation(typename std::deque<PARTICLETYPE<T> >::iterator p, ParticleSystem3D<T, PARTICLETYPE>& pSys) override;
private:
  T _dT;
  AnalyticalRandomNormal<3,T,T> _randomNormal;
};

template<typename T, template<typename U> class PARTICLETYPE>
class RandomTruncatedWalkOperation3D final : public ParticleOperation3D<T,PARTICLETYPE> {
public:
  RandomTruncatedWalkOperation3D(T dT, T diffusivity, T n);
  virtual void applyParticleOperation(typename std::deque<PARTICLETYPE<T> >::iterator p, ParticleSystem3D<T, PARTICLETYPE>& pSys) override;
private:
  T _dT;
  AnalyticalRandomTruncatedNormal<3,T,T> _randomNormal;
};

template<typename T, template<typename U> class PARTICLETYPE, typename DESCRIPTOR>
class PassiveAdvectionOperation3D final : public ParticleOperation3D<T,PARTICLETYPE> {
public:
  PassiveAdvectionOperation3D(SuperLattice<T,DESCRIPTOR>& sLattice, UnitConverter<T,DESCRIPTOR>& converter);
  virtual void applyParticleOperation(typename std::deque<PARTICLETYPE<T> >::iterator p, ParticleSystem3D<T, PARTICLETYPE>& pSys) override;
private:
  T _dT;
  SuperLatticeInterpPhysVelocity3D<T,DESCRIPTOR> _getVel;
};

}


#endif

