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
#include "dynamics/lbm.h"
#include "dynamics/guoZhaoDynamics.h"

namespace olb {

////////////////////// Class GuoZhaoBGKdynamics /////////////////////////

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 */
template<typename T, typename DESCRIPTOR, typename MOMENTA>
GuoZhaoBGKdynamics<T,DESCRIPTOR,MOMENTA>::GuoZhaoBGKdynamics (T omega)
  : legacy::BasicDynamics<T,DESCRIPTOR,MOMENTA>(),
    _omega(omega)
{
  this->getName() = "GuoZhaoBGKdynamics";
  OLB_PRECONDITION( DESCRIPTOR::template provides<descriptors::FORCE>() );

  _epsilon = (T)1.0; // This to avoid a NaN error at the first timestep.
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
T GuoZhaoBGKdynamics<T,DESCRIPTOR,MOMENTA>::computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d]) const
{
  T uSqr = util::normSqr<T,DESCRIPTOR::d>(u);
  return GuoZhaoLbHelpers<T,DESCRIPTOR>::equilibrium(iPop, _epsilon, rho, u, uSqr);
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
void GuoZhaoBGKdynamics<T,DESCRIPTOR,MOMENTA>::computeU (ConstCell<T,DESCRIPTOR>& cell, T u[DESCRIPTOR::d] ) const
{
  T rho;
  this->computeRhoU(cell, rho, u);
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
void GuoZhaoBGKdynamics<T,DESCRIPTOR,MOMENTA>::computeRhoU (ConstCell<T,DESCRIPTOR>& cell, T& rho, T u[DESCRIPTOR::d] ) const
{
  MomentaF().computeRhoU(cell, rho, u);

  const T epsilon = cell.template getField<descriptors::EPSILON>();
  const T nu      = cell.template getField<descriptors::NU>();
  const T k       = cell.template getField<descriptors::K>();

  auto bodyF = cell.template getFieldPointer<descriptors::BODY_FORCE>();

  for (int iDim=0; iDim<DESCRIPTOR::d; ++iDim) {
    u[iDim] += 0.5*epsilon*bodyF[iDim];
  }

  const T uMag = util::sqrt( util::normSqr<T,DESCRIPTOR::d>(u) );
  const T Fe = 0.;//1.75/util::sqrt(150.*util::pow(epsilon,3));

  const T c_0 = 0.5*(1 + 0.5*epsilon*nu/k);
  const T c_1 = 0.5*epsilon*Fe/util::sqrt(k);

  for (int iDim=0; iDim<DESCRIPTOR::d; ++iDim) {
    u[iDim] /= (c_0 + util::sqrt(c_0*c_0 + c_1*uMag));
  }
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
void GuoZhaoBGKdynamics<T,DESCRIPTOR,MOMENTA>::updateEpsilon (Cell<T,DESCRIPTOR>& cell)
{
  // Copying epsilon from external to member variable to provide access for computeEquilibrium.
  _epsilon = cell.template getField<descriptors::EPSILON>();
}


template<typename T, typename DESCRIPTOR, typename MOMENTA>
CellStatistic<T> GuoZhaoBGKdynamics<T,DESCRIPTOR,MOMENTA>::collide(Cell<T,DESCRIPTOR>& cell)
{
  // Copying epsilon from
  // external to member variable to provide access for computeEquilibrium.
  updateEpsilon(cell);
  T rho, u[DESCRIPTOR::d];
  this->computeRhoU(cell, rho, u);
  GuoZhaoLbHelpers<T,DESCRIPTOR>::updateGuoZhaoForce(cell, u);
  T uSqr = GuoZhaoLbHelpers<T,DESCRIPTOR>::bgkCollision(cell, _epsilon, rho, u, _omega);
  GuoZhaoLbHelpers<T,DESCRIPTOR>::addExternalForce(cell, u, _omega, rho, _epsilon);
  return {rho, uSqr};
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
T GuoZhaoBGKdynamics<T,DESCRIPTOR,MOMENTA>::getOmega() const
{
  return _omega;
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
T GuoZhaoBGKdynamics<T,DESCRIPTOR,MOMENTA>::getEpsilon()
{
  return _epsilon;
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
void GuoZhaoBGKdynamics<T,DESCRIPTOR,MOMENTA>::setOmega(T omega)
{
  _omega = omega;
}


////////////////////// Class PorousGuoSimpleBGKdynamics /////////////////////////

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 */
template<typename T, typename DESCRIPTOR, typename MOMENTA>
PorousGuoSimpleBGKdynamics<T,DESCRIPTOR,MOMENTA>::PorousGuoSimpleBGKdynamics (T omega)
  : legacy::BasicDynamics<T,DESCRIPTOR,MOMENTA>(),
    _omega(omega)
{
  OLB_PRECONDITION( DESCRIPTOR::template provides<descriptors::FORCE>() );
}


template<typename T, typename DESCRIPTOR, typename MOMENTA>
void PorousGuoSimpleBGKdynamics<T,DESCRIPTOR,MOMENTA>::computeU (ConstCell<T,DESCRIPTOR>& cell, T u[DESCRIPTOR::d] ) const
{
  T rho;
  this->computeRhoU(cell, rho, u);
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
void PorousGuoSimpleBGKdynamics<T,DESCRIPTOR,MOMENTA>::computeRhoU (ConstCell<T,DESCRIPTOR>& cell, T& rho, T u[DESCRIPTOR::d] ) const
{
  MOMENTA().computeRhoU(cell, rho, u);
  const T porosity = cell.template getField<descriptors::POROSITY>();
  for (int iDim=0; iDim<DESCRIPTOR::d; ++iDim) {
    u[iDim] *= porosity;
  }
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
CellStatistic<T> PorousGuoSimpleBGKdynamics<T,DESCRIPTOR,MOMENTA>::collide (
  Cell<T,DESCRIPTOR>& cell)
{
  T rho, u[DESCRIPTOR::d];
  this->computeRhoU(cell, rho, u);
  GuoZhaoLbHelpers<T,DESCRIPTOR>::updateGuoForceSimple(cell, u);
  T uSqr = GuoZhaoLbHelpers<T,DESCRIPTOR>::bgkCollisionAndForcing(cell, rho, u, _omega);
  return {rho, uSqr};
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
T PorousGuoSimpleBGKdynamics<T,DESCRIPTOR,MOMENTA>::getOmega() const
{
  return _omega;
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
void PorousGuoSimpleBGKdynamics<T,DESCRIPTOR,MOMENTA>::setOmega(T omega)
{
  _omega = omega;
}




}

#endif
