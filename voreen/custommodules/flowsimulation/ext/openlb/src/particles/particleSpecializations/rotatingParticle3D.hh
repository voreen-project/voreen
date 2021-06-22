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

#ifndef ROTATING_PARTICLE_3D_HH
#define ROTATING_PARTICLE_3D_HH

#include <string>
#include <iostream>
#include <set>
#include <vector>
#include <list>
#include <deque>

#include "rotatingParticle3D.h"

namespace olb {


template<typename T>
RotatingParticle3D<T>::RotatingParticle3D()
  : Particle3D<T>::Particle3D(), _aVel(3, T()), _torque(3, T())
{ }

template<typename T>
RotatingParticle3D<T>::RotatingParticle3D(std::vector<T> pos, T mas, T rad)
  : Particle3D<T>::Particle3D(pos, mas, rad), _aVel(3, T()), _torque(3, T())
{ }

template<typename T>
RotatingParticle3D<T>::RotatingParticle3D(const RotatingParticle3D<T>& p)
  : Particle3D<T>::Particle3D(p), _aVel(p.getAVel()), _torque(p._torque)
{ }

template<typename T>
RotatingParticle3D<T>::RotatingParticle3D(std::vector<T> pos, std::vector<T> vel, T mas,
    T rad)
  : Particle3D<T>::Particle3D(pos, vel, mas, rad), _aVel(3, T()), _torque(3, T())
{ }

template<typename T>
inline std::vector<T>& RotatingParticle3D<T>::getAVel()
{
  return _aVel;
}

template<typename T>
inline const std::vector<T>& RotatingParticle3D<T>::getAVel() const
{
  return _aVel;
}

template<typename T>
inline std::vector<T>& RotatingParticle3D<T>::getTorque()
{
  return _torque;
}

template<typename T>
inline const std::vector<T>& RotatingParticle3D<T>::getTorque() const
{
  return _torque;
}

template<typename T>
void RotatingParticle3D<T>::serialize(T serial[])
{
  for (int i = 0; i < 3; i++) {
    serial[i] = this->_pos[i];
    serial[i + 3] = this->_vel[i];
    serial[i + 6] = this->_force[i];
  }
  serial[9] = this->_mas;
  serial[10] = this->_rad;
  serial[11] = this->_cuboid;
  serial[12] = (double) this->_active;
  serial[13] = (double) _aVel[0];
  serial[14] = (double) _aVel[1];
  serial[15] = (double) _aVel[2];
  serial[16] = (double) _torque[0];
  serial[17] = (double) _torque[1];
  serial[18] = (double) _torque[2];
}

template<typename T>
void RotatingParticle3D<T>::unserialize(T* data)
{
  for (int i = 0; i < 3; i++) {
    this->_pos[i] = data[i];
    this->_vel[i] = data[i + 3];
    this->_force[i] = data[i + 6];
  }
  this->_mas = data[9];
  this->_rad = data[10];
  this->_cuboid = int(data[11]);
  this->_active = (bool) data[12];
  _aVel[0] = (bool) data[13];
  _aVel[1] = (bool) data[14];
  _aVel[2] = (bool) data[15];
  _torque[0] = (bool) data[16];
  _torque[1] = (bool) data[17];
  _torque[2] = (bool) data[18];
}


/////////////////////////////////////// SimulateParticles<T, RotatingParticle3D> ///////////////////////////////////////

template<typename T>
SimulateParticles<T,RotatingParticle3D>::SimulateParticles(ParticleSystem3D<T, RotatingParticle3D>* ps)
  : _pSys(ps)
{ }

template<typename T>
inline void SimulateParticles<T,RotatingParticle3D>::simulate(T dT, bool scale)
{
  _pSys->computeForce();
  _pSys->explicitEuler(dT, scale);
  _pSys->integrateTorque(dT);
}


}

#endif
