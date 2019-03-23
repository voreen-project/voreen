/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2008 Jonas Latt
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

/** \file
 * BGK Dynamics with adjustable speed of sound -- generic implementation.
 */
#ifndef CHOPARD_DYNAMICS_HH
#define CHOPARD_DYNAMICS_HH

#include "chopardDynamics.h"
#include "core/util.h"

namespace olb {

////////////////////// Class ChopardDynamics //////////////////////////

/** \param vs2_ speed of sound
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 */
template<typename T, template<typename U> class Lattice>
ChopardDynamics<T,Lattice>::ChopardDynamics (
  T vs2_, T omega_, Momenta<T,Lattice>& momenta_ )
  : BasicDynamics<T,Lattice>(momenta_),
    vs2(vs2_),
    omega(omega_)
{ }

/** With this constructor, the speed of sound is vs2 = cs2
 *  \param omega_ relaxation parameter, related to the dynamic viscosity
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 */
template<typename T, template<typename U> class Lattice>
ChopardDynamics<T,Lattice>::ChopardDynamics (
  T omega_, Momenta<T,Lattice>& momenta_ )
  : BasicDynamics<T,Lattice>(momenta_),
    vs2((T)1/Lattice<T>::invCs2),
    omega(omega_)
{ }

template<typename T, template<typename U> class Lattice>
T ChopardDynamics<T,Lattice>::computeEquilibrium
(int iPop, T rho, const T u[Lattice<T>::d], T uSqr) const
{
  return chopardEquilibrium(iPop, rho, u, uSqr, vs2);
}

template<typename T, template<typename U> class Lattice>
void ChopardDynamics<T,Lattice>::collide (
  Cell<T,Lattice>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho, u[Lattice<T>::d];
  this->_momenta.computeRhoU(cell, rho, u);
  T uSqr = chopardBgkCollision(cell, rho, u, vs2, omega);
  statistics.incrementStats(rho, uSqr);
}

template<typename T, template<typename U> class Lattice>
void ChopardDynamics<T,Lattice>::staticCollide (
  Cell<T,Lattice>& cell,
  const T u[Lattice<T>::d],
  LatticeStatistics<T>& statistics )
{
  T rho;
  rho = this->_momenta.computeRho(cell);
  T uSqr = chopardBgkCollision(cell, rho, u, vs2, omega);
  statistics.incrementStats(rho, uSqr);
}

template<typename T, template<typename U> class Lattice>
T ChopardDynamics<T,Lattice>::getOmega() const
{
  return omega;
}

template<typename T, template<typename U> class Lattice>
void ChopardDynamics<T,Lattice>::setOmega(T omega_)
{
  omega = omega_;
}

template<typename T, template<typename U> class Lattice>
T ChopardDynamics<T,Lattice>::getVs2() const
{
  return vs2;
}

template<typename T, template<typename U> class Lattice>
void ChopardDynamics<T,Lattice>::setVs2(T vs2_)
{
  vs2 = vs2_;
}

template<typename T, template<typename U> class Lattice>
T ChopardDynamics<T,Lattice>::chopardBgkCollision (
  Cell<T,Lattice>& cell, T rho, const T u[Lattice<T>::d], T vs2, T omega)
{
  const T uSqr = util::normSqr<T,Lattice<T>::d>(u);
  for (int iPop=0; iPop < Lattice<T>::q; ++iPop) {
    cell[iPop] *= (T)1-omega;
    cell[iPop] += omega * chopardEquilibrium(iPop, rho, u, uSqr, vs2);
  }
  return uSqr;
}

template<typename T, template<typename U> class Lattice>
T ChopardDynamics<T,Lattice>::chopardEquilibrium (
  int iPop, T rho, const T u[Lattice<T>::d], T uSqr, T vs2)
{
  if (iPop==0) {
    return rho*( (T)1 -
                 vs2*Lattice<T>::invCs2*((T)1-Lattice<T>::t[0]) -
                 Lattice<T>::t[0]/(T)2*Lattice<T>::invCs2*uSqr ) - Lattice<T>::t[0];
  } else {
    T c_u = T();
    for (int iD=0; iD < Lattice<T>::d; ++iD) {
      c_u += Lattice<T>::c[iPop][iD]*u[iD];
    }
    return rho * Lattice<T>::t[iPop] * Lattice<T>::invCs2* (
             vs2 + c_u +
             Lattice<T>::invCs2 / (T)2 * c_u*c_u -
             uSqr / (T)2
           ) - Lattice<T>::t[iPop];
  }
}


}

#endif
