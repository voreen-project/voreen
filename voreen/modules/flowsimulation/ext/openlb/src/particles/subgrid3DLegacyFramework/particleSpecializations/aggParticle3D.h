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

#ifndef AGG_PARTICLE_3D_H
#define AGG_PARTICLE_3D_H

#include <set>
#include <vector>
#include <list>
#include <deque>
#include <string>
#include <iostream>
#include "particles/subgrid3DLegacyFramework/particle3D.h"

namespace olb {

/*
 * Particles for Agglomeration
 */
template<typename T>
class AggParticle3D: public Particle3D<T> {
public:
  AggParticle3D();
  AggParticle3D(std::vector<T> pos, T mas = 1., T rad = 1.);
  AggParticle3D(std::vector<T> pos, std::vector<T> vel, T mas = 1., T rad = 1.);
  AggParticle3D(const Particle3D<T>& p);
  inline void setMass(T mas);
  inline void setRad(T rad);
  inline const bool& getAggl();
  inline const bool& getAggl() const;
  inline void setAggl(bool aggl);
  void serialize(T serial[]);
  void unserialize(T*);

  static const int serialPartSize = 14;

private:
  bool _aggl;
};


}
#endif

