/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2017 Marie-Luise Maier, Mathias J. Krause, Sascha Janz
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

/// Schiller and Naumann drag force
#ifndef SCHILLERNAUMANNDRAGFORCE_3D_HH
#define SCHILLERNAUMANNDRAGFORCE_3D_HH

#include <cmath>
#include "schillerNaumannDragForce3D.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace olb {

template<typename T, template<typename U> class PARTICLETYPE, typename DESCRIPTOR>
SchillerNaumannDragForce3D<T, PARTICLETYPE, DESCRIPTOR>::SchillerNaumannDragForce3D(SuperLatticeInterpPhysVelocity3D<T, DESCRIPTOR>& getVel, T dynVisc, T physDensity):
  Force3D<T, PARTICLETYPE>(), _getVel(getVel), _dynVisc(dynVisc),  _physDensity(physDensity)
{
//  this->_name = "SchillerNaumannDragForce3D";
}

template<typename T, template<typename U> class PARTICLETYPE, typename DESCRIPTOR>
void SchillerNaumannDragForce3D<T, PARTICLETYPE, DESCRIPTOR>::applyForce(
  typename std::deque<PARTICLETYPE<T> >::iterator p, int pInt,
  ParticleSystem3D<T, PARTICLETYPE>& pSys)
{

  Vector<T, 3> force = {T(0), T(0), T(0)} ;
  
  T fluidVel[3] = {0., 0., 0.};
  _getVel(fluidVel, &p->getPos()[0], p->getCuboid());

  T fluidVelAbs = util::sqrt(std::pow(fluidVel[0], 2.) + util::pow(fluidVel[1], 2.) + util::pow(fluidVel[2], 2.));
  T partVelAbs = util::sqrt(std::pow(p->getVel()[0], 2.) + util::pow(p->getVel()[1], 2.) + util::pow(p->getVel()[2], 2.));
  T particleRe = 2. * p->getRad() * abs(fluidVelAbs - partVelAbs) * _physDensity / _dynVisc;
  T coeffSN = 1. + 0.15 * util::pow(particleRe, 0.687);

  for (int i = 0; i < 3; i++) {

    force[i] = -1. * 6 * M_PI * p->getRad() * (p->getVel()[i] - fluidVel[i]) * _dynVisc * coeffSN;
    p->getForce()[i] += force[i] ;
  }

}
}

#endif // SchillerNaumannDragForce3D
