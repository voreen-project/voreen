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

#ifndef EL_PARTICLE_3D_HH
#define EL_PARTICLE_3D_HH

#include <string>
#include <iostream>
#include <set>
#include <vector>
#include <list>
#include <deque>

#include "elParticle3D.h"

namespace olb {


template<typename T>
ElParticle3D<T>::ElParticle3D()
  : Particle3D<T>(),
    _charge(1.)
{ }

template<typename T>
ElParticle3D<T>::ElParticle3D(std::vector<T> pos, T mas, T rad, T charge)
  : Particle3D<T>(pos, mas, rad),
    _charge(charge)
{ }

template<typename T>
ElParticle3D<T>::ElParticle3D(std::vector<T> pos, std::vector<T> vel, T mas,
                              T rad, T charge)
  : Particle3D<T>(pos, vel, mas, rad),
    _charge(charge)
{
}

template<typename T>
ElParticle3D<T>::ElParticle3D(const ElParticle3D<T>& p)
  : Particle3D<T>(p),
    _charge(p._charge)
{ }

template<typename T>
void ElParticle3D<T>::serialize(T serial[])
{
  serial[0] = this->_pos[0];
  serial[1] = this->_pos[1];
  serial[2] = this->_pos[2];
  serial[3] = this->_vel[0];
  serial[4] = this->_vel[1];
  serial[5] = this->_vel[2];
  serial[6] = this->_rad;
  serial[7] = this->_mas;
  serial[8] = (double) this->_active;
  serial[9] = _charge;
}

template<typename T>
void ElParticle3D<T>::unserialize(T* data)
{
  this->_pos[0] = data[0];
  this->_pos[1] = data[1];
  this->_pos[2] = data[2];
  this->_vel[0] = data[3];
  this->_vel[1] = data[4];
  this->_vel[2] = data[5];
  this->_rad = data[6];
  this->_mas = data[7];
  this->_active = (bool) data[8];
  _charge = data[9];
}


}

#endif
