/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014 Robin Trunk, Mathias J. Krause
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

#ifndef BOUNDARYSIMPLEREFLECTION3D_H
#define BOUNDARYSIMPLEREFLECTION3D_H

#include <set>
#include <math.h>
#include "geometry/superGeometry.h"
#include "particles/subgrid3DLegacyFramework/particleSystem3D.h"
#include "boundary3D.h"

namespace olb {

template<typename T, template<typename U> class PARTICLETYPE>
class ParticleSystem3D;

template<typename T, template<typename U> class PARTICLETYPE>
class SimpleReflectBoundary3D : public Boundary3D<T, PARTICLETYPE> {
public:
  SimpleReflectBoundary3D(T dT, SuperGeometry<T,3>& sg, std::set<int> materials);
  ~SimpleReflectBoundary3D() override {};
  void applyBoundary(typename std::deque<PARTICLETYPE<T> >::iterator& p, ParticleSystem3D<T, PARTICLETYPE>& psSys) override;
private:
  T _dT;
  SuperGeometry<T,3>& _sg;
  std::set<int> _materials;
  std::set<int>::iterator _matIter;
};

template<typename T, template<typename U> class PARTICLETYPE>
class ReflectBoundary3D : public Boundary3D<T, PARTICLETYPE> {
public:
  ReflectBoundary3D(T dT, SuperGeometry<T,3>& sg, std::set<int> materials);
  ~ReflectBoundary3D() override {};
  void applyBoundary(typename std::deque<PARTICLETYPE<T> >::iterator& p, ParticleSystem3D<T, PARTICLETYPE>& psSys) override;
private:
  Vector<T,3> computeDiscreteNormal(int latticeR[]);
  T _dT;
  SuperGeometry<T,3>& _sg;
  std::set<int> _materials;
  std::set<int>::iterator _matIter;
};

}
#endif
