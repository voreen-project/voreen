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

#ifndef AGG_PARTICLE_3D_HH
#define AGG_PARTICLE_3D_HH

#include <string>
#include <iostream>
#include <set>
#include <vector>
#include <list>
#include <deque>

#include "aggParticle3D.h"

namespace olb {


template<typename T>
AggParticle3D<T>::AggParticle3D()
  : Particle3D<T>::Particle3D()
{
  _aggl = false;
}

template<typename T>
AggParticle3D<T>::AggParticle3D(std::vector<T> pos, T mas, T rad)
  : Particle3D<T>::Particle3D(pos, mas, rad)
{
  _aggl = false;
}

template<typename T>
AggParticle3D<T>::AggParticle3D(const Particle3D<T>& p)
  : Particle3D<T>::Particle3D(p)
{
  _aggl = false;
}

template<typename T>
AggParticle3D<T>::AggParticle3D(std::vector<T> pos, std::vector<T> vel, T mas,
                                T rad)
  : Particle3D<T>::Particle3D(pos, vel, mas, rad)
{
  _aggl = false;
}

template<typename T>
inline void AggParticle3D<T>::setMass(T mas)
{
  this->_mas = mas;
}

template<typename T>
inline void AggParticle3D<T>::setRad(T rad)
{
  this->_rad = rad;
}

template<typename T>
inline const bool& AggParticle3D<T>::getAggl()
{
  return _aggl;
}

template<typename T>
inline const bool& AggParticle3D<T>::getAggl() const
{
  return _aggl;
}

template<typename T>
inline void AggParticle3D<T>::setAggl(bool aggl)
{
  _aggl = aggl;
  if (aggl) {
    this->setActive(false);
  }
}

template<typename T>
void AggParticle3D<T>::serialize(T serial[])
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
  serial[13] = (double) _aggl;
}

template<typename T>
void AggParticle3D<T>::unserialize(T* data)
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
  _aggl = (bool) data[13];
}


}

#endif
