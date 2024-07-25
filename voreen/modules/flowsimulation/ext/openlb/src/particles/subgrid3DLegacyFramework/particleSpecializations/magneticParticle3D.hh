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

#ifndef MAGNETIC_PARTICLE_3D_HH
#define MAGNETIC_PARTICLE_3D_HH

#include <string>
#include <iostream>
#include <set>
#include <vector>
#include <list>
#include <deque>

#include "magneticParticle3D.h"

namespace olb {


/////////////////////////////////////// MagneticParticle3D ///////////////////////////////////////

template<typename T>
MagneticParticle3D<T>::MagneticParticle3D()
  : Particle3D<T>::Particle3D(), _dMoment(3, T()), _aVel(3, T()), _torque(3, T()), _magnetisation(T())
{ }

template<typename T>
MagneticParticle3D<T>::MagneticParticle3D(const MagneticParticle3D<T>& p)
  : Particle3D<T>::Particle3D(p), _dMoment(p._dMoment), _aVel(p._aVel), _torque(p._torque), _magnetisation(p._magnetisation), _sActivity(p._sActivity)
{ }

template<typename T>
MagneticParticle3D<T>::MagneticParticle3D(std::vector<T> pos, std::vector<T> vel, T mas,
    T rad, int id, T masAdd)
  : Particle3D<T>::Particle3D(pos, vel, mas, rad), _dMoment(3, T()), _aVel(3, T()), _torque(3, T()), _magnetisation(T())
{ }

template<typename T>
MagneticParticle3D<T>::MagneticParticle3D(std::vector<T> pos, std::vector<T> vel, T mas, T rad, int id,
    std::vector<T> dMoment, std::vector<T> aVel, std::vector<T> torque, T magnetisation)
  : Particle3D<T>::Particle3D(pos, vel, mas, rad, id),
    _dMoment(dMoment), _aVel(aVel), _torque(torque), _magnetisation(magnetisation)
{ }

template<typename T>
MagneticParticle3D<T>::MagneticParticle3D(std::vector<T> pos, std::vector<T> vel, T mas, T rad, int id,
    std::vector<T> dMoment, std::vector<T> aVel, std::vector<T> torque, T magnetisation, int sActivity)
  : Particle3D<T>::Particle3D(pos, vel, mas, rad, id),
    _dMoment(dMoment), _aVel(aVel), _torque(torque), _magnetisation(magnetisation), _sActivity(sActivity)
{ }

template<typename T>
inline void MagneticParticle3D<T>::resetTorque()
{
  for (int i = 0; i < 3; i++) {
    _torque[i] = 0.;
  }
}

template<typename T>
inline std::vector<T>& MagneticParticle3D<T>::getMoment()
{
  return _dMoment;
}

template<typename T>
inline const std::vector<T>& MagneticParticle3D<T>::getMoment() const
{
  return _dMoment;
}

template<typename T>
inline void MagneticParticle3D<T>::setMoment(std::vector<T> moment)
{
  _dMoment = moment;
}

template<typename T>
inline void MagneticParticle3D<T>::setAVel(std::vector<T> aVel)
{
  _aVel = aVel;
}

template<typename T>
inline std::vector<T>& MagneticParticle3D<T>::getAVel()
{
  return _aVel;
}

template<typename T>
inline const std::vector<T>& MagneticParticle3D<T>::getAVel() const
{
  return _aVel;
}

template<typename T>
inline void MagneticParticle3D<T>::setTorque(std::vector<T> torque)
{
  _torque = torque;
}

template<typename T>
inline void MagneticParticle3D<T>::setMagnetisation(T magnetisation)
{
  _magnetisation = magnetisation;
  //std::cout<< "Setting magnetisation: "<< _magnetisation << std::endl;
}

template<typename T>
inline void MagneticParticle3D<T>::setSActivity(int sActivity)
{
  _sActivity = sActivity;
}

template<typename T>
inline typename std::deque<std::list<MagneticParticle3D<T>*>>::iterator& MagneticParticle3D<T>::getAggloItr()
{
  return _aggloItr;
}

template<typename T>
inline void MagneticParticle3D<T>::setAggloItr(typename std::deque<std::list<MagneticParticle3D<T>*>>::iterator aggloItr)
{
  _aggloItr = aggloItr;
}

template<typename T>
inline std::vector<T>& MagneticParticle3D<T>::getTorque()
{
  return _torque;
}

template<typename T>
inline const std::vector<T>& MagneticParticle3D<T>::getTorque() const
{
  return _torque;
}

template<typename T>
inline T& MagneticParticle3D<T>::getMagnetisation()
{
  return _magnetisation;
}

template<typename T>
inline const T& MagneticParticle3D<T>::getMagnetisation() const
{
  return _magnetisation;
}

template<typename T>
inline int& MagneticParticle3D<T>::getSActivity()
{
  return MagneticParticle3D<T>::_sActivity;
}

template<typename T>
void MagneticParticle3D<T>::serialize(T serial[])
{
  for (int i = 0; i < 3; i++) {
    serial[i] = this->_pos[i];
    serial[i + 3] = this->_vel[i];
    serial[i + 6] = this->_force[i];
  }
  serial[9] = this->_mas;
  serial[10] = this->_rad;
  serial[11] = (double) this->_cuboid;
  serial[12] = (double) this->_active;
  serial[13] = (double) this->_id;
  for (int i = 0; i < 3; i++) {
    serial[i + 14] = this->_storeForce[i];
    serial[i + 17] = _dMoment[i];
    serial[i + 20] = _aVel[i];
    serial[i + 23] = _torque[i];
  }
  serial[26] = _magnetisation;
  serial[27] = (double) _sActivity;
}

template<typename T>
void MagneticParticle3D<T>::unserialize(T* data)
{
  for (int i = 0; i < 3; i++) {
    this->_pos[i] = data[i];
    this->_vel[i] = data[i + 3];
    this->_force[i] = data[i + 6];
  }
  this->_mas = data[9];
  this->_rad = data[10];
  this->_cuboid = (int) data[11];
  this->_active = (bool) data[12];
  this->_id = (int) data[13];
  for (int i = 0; i < 3; i++) {
    this->_storeForce[i] = data[i + 14];
    _dMoment[i] = data[i + 17];
    _aVel[i] = data[i + 20];
    _torque[i] = data[i + 23];
  }
  _magnetisation = data[26];
  _sActivity = (int) data[27];
}


/////////////////////////////////////// SimulateParticles<T, MagneticParticle3D> ///////////////////////////////////////

template<typename T>
SimulateParticles<T,MagneticParticle3D>::SimulateParticles(ParticleSystem3D<T, MagneticParticle3D>* ps)
  : _pSys(ps)
{ }

template<typename T>
inline void SimulateParticles<T,MagneticParticle3D>::simulate(T dT, bool scale)
{
  _pSys->resetMag();
  _pSys->computeForce();
  _pSys->explicitEuler(dT, scale);
  _pSys->integrateTorqueMag(dT);

#ifdef CollisionModels
  _pSys->partialElasticImpact(0.67);
#endif
}

template<typename T>
inline void SimulateParticles<T,MagneticParticle3D>::simulateWithTwoWayCoupling_Mathias ( T dT,
    ForwardCouplingModel<T,MagneticParticle3D>& forwardCoupling,
    BackCouplingModel<T,MagneticParticle3D>& backCoupling,
    int material, int subSteps, bool scale )
{
  for (int iSubStep=1; iSubStep<=subSteps; iSubStep++) {
    if (! _pSys->executeForwardCoupling(forwardCoupling) ) {
      std::cout << " on substep " << iSubStep << std::endl;
      singleton::exit(1);
    }
    _pSys->computeForce();
    _pSys->explicitEuler(dT/(T)(subSteps), scale);
    //_pSys->rungeKutta4(dT/(T)(subSteps));
  }
  _pSys->executeBackwardCoupling(backCoupling, material);
}

template<typename T>
inline void SimulateParticles<T,MagneticParticle3D>::simulateWithTwoWayCoupling_Davide ( T dT,
    ForwardCouplingModel<T,MagneticParticle3D>& forwardCoupling,
    BackCouplingModel<T,MagneticParticle3D>& backCoupling,
    int material, int subSteps, bool scale )
{
  for (int iSubStep=1; iSubStep<=subSteps; iSubStep++) {
    if (! _pSys->executeForwardCoupling(forwardCoupling) ) {
      std::cout << " on substep " << iSubStep << std::endl;
      singleton::exit(1);
    }
    _pSys->computeForce();
    _pSys->explicitEuler(dT/(T)(subSteps), scale);
    //_pSys->rungeKutta4(dT/(T)(subSteps));
    _pSys->executeBackwardCoupling(backCoupling, material, subSteps);
  }
}

template<typename T>
inline void SimulateParticles<T,MagneticParticle3D>::simulate(T dT, std::set<int> sActivityOfParticle, bool scale)
{
  _pSys->resetMag(sActivityOfParticle);
  _pSys->computeForce(sActivityOfParticle);
  _pSys->explicitEuler(dT, sActivityOfParticle, scale);
  _pSys->integrateTorqueMag(dT, sActivityOfParticle);

#ifdef CollisionModels
  _pSys->partialElasticImpact(0.67);
#endif

#ifdef CollisionModelsCombindedWithMechContactForce
  _pSys->partialElasticImpactForCombinationWithMechContactForce(0.67);
#endif
}


}

#endif
