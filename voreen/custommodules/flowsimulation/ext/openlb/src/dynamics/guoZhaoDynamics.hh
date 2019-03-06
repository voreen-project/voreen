/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2016-2017 Davide Dapelo, Mathias J. Krause
 *  OpenLB e-mail contact: info@openlb.net
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
 * Specific dynamics classes for Guo and Zhao (2002) porous model, with
 * which a Cell object can be instantiated -- generic implementation.
 */
#ifndef LB_GUOZHAO_DYNAMICS_HH
#define LB_GUOZHAO_DYNAMICS_HH

#include <algorithm>
#include <limits>
#include "dynamics/dynamics.h"
#include "core/cell.h"
#include "dynamics/guoZhaoLbHelpers.h"
#include "dynamics/firstOrderLbHelpers.h"

namespace olb {

////////////////////// Class GuoZhaoBGKdynamics /////////////////////////

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 */
template<typename T, template<typename U> class Lattice>
GuoZhaoBGKdynamics<T,Lattice>::GuoZhaoBGKdynamics (
  T omega, Momenta<T,Lattice>& momenta )
  : BasicDynamics<T,Lattice>(momenta),
    _omega(omega)
{
  // This ensures both that the constant sizeOfForce is defined in
  // ExternalField and that it has the proper size
  OLB_PRECONDITION( Lattice<T>::d == Lattice<T>::ExternalField::sizeOfForce );

  _epsilon = (T)1.0; // This to avoid a NaN error at the first timestep.
}

template<typename T, template<typename U> class Lattice>
T GuoZhaoBGKdynamics<T,Lattice>::computeEquilibrium(int iPop, T rho, const T u[Lattice<T>::d], T uSqr) const
{
//  int *foo = (int*)-1; // Making a bad pointer
//  cout << *foo; // Crashing the program
//  cout << "computeEquilibrium function reached. Stopping." << endl;
//  exit(1); // exits nicely
  return GuoZhaoLbHelpers<T,Lattice>::equilibrium(iPop, _epsilon, rho, u, uSqr);
}

template<typename T, template<typename U> class Lattice>
void GuoZhaoBGKdynamics<T,Lattice>::computeU (Cell<T,Lattice> const& cell, T u[Lattice<T>::d] ) const
{
  T rho;
  this->_momenta.computeRhoU(cell, rho, u);
  for (int iVel=0; iVel<Lattice<T>::d; ++iVel) {
    u[iVel] += cell.getExternal(forceBeginsAt)[iVel] / (T)2.;
  }
}

template<typename T, template<typename U> class Lattice>
void GuoZhaoBGKdynamics<T,Lattice>::computeRhoU (Cell<T,Lattice> const& cell, T& rho, T u[Lattice<T>::d] ) const
{
  this->_momenta.computeRhoU(cell, rho, u);
  for (int iVel=0; iVel<Lattice<T>::d; ++iVel) {
    u[iVel] += cell.getExternal(forceBeginsAt)[iVel] / (T)2.;
  }
}

template<typename T, template<typename U> class Lattice>
void GuoZhaoBGKdynamics<T,Lattice>::updateEpsilon (Cell<T,Lattice>& cell)
{
  _epsilon = *cell.getExternal(Lattice<T>::ExternalField::epsilonAt); //Copying epsilon from
  // external to member variable to provide access for computeEquilibrium.
}


template<typename T, template<typename U> class Lattice>
void GuoZhaoBGKdynamics<T,Lattice>::collide (
  Cell<T,Lattice>& cell,
  LatticeStatistics<T>& statistics )
{
  // Copying epsilon from
  // external to member variable to provide access for computeEquilibrium.
  updateEpsilon(cell);
  T rho, u[Lattice<T>::d];
  this->_momenta.computeRhoU(cell, rho, u);
  T* force = cell.getExternal(forceBeginsAt);
  for (int iVel=0; iVel<Lattice<T>::d; ++iVel) {
    u[iVel] += force[iVel] / (T)2.;
  }
  T uSqr = GuoZhaoLbHelpers<T,Lattice>::bgkCollision(cell, _epsilon, rho, u, _omega);
  GuoZhaoLbHelpers<T,Lattice>::updateGuoZhaoForce(cell, u);
  lbHelpers<T,Lattice>::addExternalForce(cell, u, _omega, rho);
  statistics.incrementStats(rho, uSqr);
}

template<typename T, template<typename U> class Lattice>
void GuoZhaoBGKdynamics<T,Lattice>::staticCollide (
  Cell<T,Lattice>& cell,
  const T u[Lattice<T>::d],
  LatticeStatistics<T>& statistics )
{
  // Copying epsilon from
  // external to member variable to provide access for computeEquilibrium.
  updateEpsilon(cell);
  T rho, uDummy[Lattice<T>::d];
  this->_momenta.computeRhoU(cell, rho, uDummy);
  T uSqr =GuoZhaoLbHelpers<T,Lattice>::bgkCollision(cell, _epsilon, rho, u, _omega);
  GuoZhaoLbHelpers<T,Lattice>::updateGuoZhaoForce(cell, u);
  lbHelpers<T,Lattice>::addExternalForce(cell, u, _omega, rho);
  statistics.incrementStats(rho, uSqr);
}

template<typename T, template<typename U> class Lattice>
T GuoZhaoBGKdynamics<T,Lattice>::getOmega() const
{
  return _omega;
}

template<typename T, template<typename U> class Lattice>
T GuoZhaoBGKdynamics<T,Lattice>::getEpsilon()
{
  return _epsilon;
}

template<typename T, template<typename U> class Lattice>
void GuoZhaoBGKdynamics<T,Lattice>::setOmega(T omega)
{
  _omega = omega;
}

}

#endif
