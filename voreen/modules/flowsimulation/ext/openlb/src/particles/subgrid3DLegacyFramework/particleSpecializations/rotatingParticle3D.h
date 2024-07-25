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

#ifndef ROTATING_PARTICLE_3D_H
#define ROTATING_PARTICLE_3D_H

#include <set>
#include <vector>
#include <list>
#include <deque>
#include <string>
#include <iostream>
#include "particles/subgrid3DLegacyFramework/particle3D.h"

namespace olb {


/////////////////////////////////////////// RotatingParticle3D ///////////////////////////////////////////

/*
 * Rotating Particles
 */
template<typename T>
class RotatingParticle3D: public Particle3D<T> {
public:
  RotatingParticle3D();
  RotatingParticle3D(std::vector<T> pos, T mas = 1., T rad = 1.);
  RotatingParticle3D(std::vector<T> pos, std::vector<T> vel, T mas = 1., T rad = 1.);
  RotatingParticle3D(const RotatingParticle3D<T>& p);
  inline std::vector<T>& getAVel();
  inline const std::vector<T>& getAVel() const;
  inline std::vector<T>& getTorque();
  inline const std::vector<T>& getTorque() const;
  void serialize(T serial[]);
  void unserialize(T*);

  static const int serialPartSize = 19;

private:
  std::vector<T> _aVel;
  std::vector<T> _torque;
};

/////////////////////////////////////// SimulateParticles<T, RotatingParticle3D> ///////////////////////////////////////
template<typename T>
class SimulateParticles<T, RotatingParticle3D> {
public:
  SimulateParticles(ParticleSystem3D<T, RotatingParticle3D>* ps);
  inline void simulate(T dT, bool scale = false);

private:
  ParticleSystem3D<T, RotatingParticle3D>* _pSys;
};

}
#endif

