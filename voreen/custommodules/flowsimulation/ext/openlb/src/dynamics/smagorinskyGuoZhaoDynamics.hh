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
 * Specific dynamics classes for Guo and Zhao (2002) porous model
 * with a Smagorinsky LES turbulence model, with
 * which a Cell object can be instantiated -- generic implementation.
 */
#ifndef LB_SMAGO_GUOZHAO_DYNAMICS_HH
#define LB_SMAGO_GUOZHAO_DYNAMICS_HH

#include <algorithm>
#include <limits>
#include "dynamics/dynamics.h"
#include "core/cell.h"
#include "dynamics/guoZhaoLbHelpers.h"
#include "dynamics/firstOrderLbHelpers.h"

namespace olb {

////////////////////// Class SmagorinskyGuoZhaoBGKdynamics /////////////////////////

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 */
template<typename T, template<typename U> class Lattice>
SmagorinskyGuoZhaoBGKdynamics<T,Lattice>::SmagorinskyGuoZhaoBGKdynamics (T omega_,
    Momenta<T,Lattice>& momenta_, T smagoConst_, T dx_, T dt_ )
  : GuoZhaoBGKdynamics<T,Lattice>(omega_,momenta_), smagoConst(smagoConst_),
    preFactor(computePreFactor(omega_,smagoConst_) )
{ }

template<typename T, template<typename U> class Lattice>
void SmagorinskyGuoZhaoBGKdynamics<T,Lattice>::collide (
  Cell<T,Lattice>& cell,
  LatticeStatistics<T>& statistics )
{
  // Copying epsilon from
  // external to member variable to provide access for computeEquilibrium.
  this->updateEpsilon(cell);
  T rho, u[Lattice<T>::d], pi[util::TensorVal<Lattice<T> >::n];
  this->_momenta.computeAllMomenta(cell, rho, u, pi);
  T newOmega = computeOmega(this->getOmega(), preFactor, rho, pi);
  T* force = cell.getExternal(forceBeginsAt);
  for (int iVel=0; iVel<Lattice<T>::d; ++iVel) {
    u[iVel] += force[iVel] / (T)2.;
  }
  T uSqr = GuoZhaoLbHelpers<T,Lattice>::bgkCollision(cell, this->getEpsilon(), rho, u, newOmega);
  GuoZhaoLbHelpers<T,Lattice>::updateGuoZhaoForce(cell, u);
  lbHelpers<T,Lattice>::addExternalForce(cell, u, newOmega, rho);
  statistics.incrementStats(rho, uSqr);
}

template<typename T, template<typename U> class Lattice>
void SmagorinskyGuoZhaoBGKdynamics<T,Lattice>::staticCollide (
  Cell<T,Lattice>& cell,
  const T u[Lattice<T>::d],
  LatticeStatistics<T>& statistics )
{
  // Copying epsilon from
  // external to member variable to provide access for computeEquilibrium.
  this->updateEpsilon(cell);
  T rho, uTemp[Lattice<T>::d], pi[util::TensorVal<Lattice<T> >::n];
  this->_momenta.computeAllMomenta(cell, rho, uTemp, pi);
  T newOmega = computeOmega(this->getOmega(), preFactor, rho, pi);
  T uSqr =GuoZhaoLbHelpers<T,Lattice>::bgkCollision(cell, this->getEpsilon(), rho, u, newOmega);
  GuoZhaoLbHelpers<T,Lattice>::updateGuoZhaoForce(cell, u);
  lbHelpers<T,Lattice>::addExternalForce(cell, u, newOmega, rho);
  statistics.incrementStats(rho, uSqr);
}

template<typename T, template<typename U> class Lattice>
T SmagorinskyGuoZhaoBGKdynamics<T,Lattice>::getSmagorinskyOmega(Cell<T,Lattice>& cell )
{
  T rho, uTemp[Lattice<T>::d], pi[util::TensorVal<Lattice<T> >::n];
  this->_momenta.computeAllMomenta(cell, rho, uTemp, pi);
  T newOmega = computeOmega(this->getOmega(), preFactor, rho, pi);
  return newOmega;
}

template<typename T, template<typename U> class Lattice>
T SmagorinskyGuoZhaoBGKdynamics<T,Lattice>::computePreFactor(T omega_, T smagoConst_)
{
  return (T)smagoConst_*smagoConst_*Lattice<T>::invCs2*Lattice<T>::invCs2*2*sqrt(2);
}

template<typename T, template<typename U> class Lattice>
void SmagorinskyGuoZhaoBGKdynamics<T,Lattice>::setOmega(T omega)
{
//  _omega = omega;
  GuoZhaoBGKdynamics<T,Lattice>::setOmega(omega);
  preFactor = computePreFactor(omega, smagoConst);
}

template<typename T, template<typename U> class Lattice>
T SmagorinskyGuoZhaoBGKdynamics<T,Lattice>::computeOmega(T omega0, T preFactor_, T rho,
    T pi[util::TensorVal<Lattice<T> >::n] )
{
  T PiNeqNormSqr = pi[0]*pi[0] + 2.0*pi[1]*pi[1] + pi[2]*pi[2];
  if (util::TensorVal<Lattice<T> >::n == 6) {
    PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2*pi[4]*pi[4] +pi[5]*pi[5];
  }
  T PiNeqNorm    = sqrt(PiNeqNormSqr);
  /// Molecular realaxation time
  T tau_mol = 1. /omega0;
  /// Turbulent realaxation time
  T tau_turb = 0.5*(sqrt(tau_mol*tau_mol + preFactor_/rho*PiNeqNorm) - tau_mol);
  /// Effective realaxation time
  tau_eff = tau_mol+tau_turb;
  T omega_new= 1./tau_eff;
  return omega_new;
}

}

#endif
