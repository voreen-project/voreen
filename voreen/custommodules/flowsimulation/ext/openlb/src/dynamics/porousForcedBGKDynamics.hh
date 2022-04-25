/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2016 Asher Zarth, Thomas Henn, Mathias J. Krause, Jonas Latt, Jan E. Marquardt
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
 * BGK Dynamics for porous -- generic implementation.
 */
#ifndef POROUS_FORCED_BGK_DYNAMICS_HH
#define POROUS_FORCED_BGK_DYNAMICS_HH

#include "porousForcedBGKDynamics.h"
#include "core/cell.h"
#include "dynamics.h"
#include "core/util.h"
#include "lbm.h"
#include "math.h"

namespace olb {

////////////////////// Class PorousForcedBGKdynamics //////////////////////////

template<typename T, typename DESCRIPTOR, typename MOMENTA>
PorousForcedBGKdynamics<T,DESCRIPTOR,MOMENTA>::PorousForcedBGKdynamics (T omega)
  : legacy::BasicDynamics<T,DESCRIPTOR,MOMENTA>(), _omega(omega)
{
  this->getName() = "PorousForcedBGKdynamics";
  OLB_PRECONDITION( DESCRIPTOR::template provides<descriptors::FORCE>() );
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
void PorousForcedBGKdynamics<T,DESCRIPTOR,MOMENTA>::computeU (ConstCell<T,DESCRIPTOR>& cell, T u[DESCRIPTOR::d] ) const
{
  T rho;
  MomentaF().computeRhoU(cell, rho, u);
  for (int iVel=0; iVel<DESCRIPTOR::d; ++iVel) {
    u[iVel] += cell.template getFieldPointer<descriptors::FORCE>()[iVel] / (T)2.;
  }
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
void PorousForcedBGKdynamics<T,DESCRIPTOR,MOMENTA>::computeRhoU (ConstCell<T,DESCRIPTOR>& cell, T& rho, T u[DESCRIPTOR::d] ) const
{
  MomentaF().computeRhoU(cell, rho, u);
  for (int iVel=0; iVel<DESCRIPTOR::d; ++iVel) {
    u[iVel] += cell.template getFieldPointer<descriptors::FORCE>()[iVel] / (T)2.;
  }
}


template<typename T, typename DESCRIPTOR, typename MOMENTA>
CellStatistic<T> PorousForcedBGKdynamics<T,DESCRIPTOR,MOMENTA>::collide (
  Cell<T,DESCRIPTOR>& cell)
{
  T rho, u[DESCRIPTOR::d];
  MomentaF().computeRhoU(cell, rho, u);
  auto force = cell.template getField<descriptors::FORCE>();
  for (int iVel=0; iVel<DESCRIPTOR::d; ++iVel) {
    u[iVel] += force[iVel] / (T)2.;
  }
  T porosity = cell.template getField<descriptors::POROSITY>();
  for (int i=0; i<DESCRIPTOR::d; i++)  {
    u[i] *= porosity;
  }
  T uSqr = lbm<DESCRIPTOR>::bgkCollision(cell, rho, u, _omega);
  lbm<DESCRIPTOR>::addExternalForce(cell, rho, u, _omega, force);
  return {rho, uSqr};
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
T PorousForcedBGKdynamics<T,DESCRIPTOR,MOMENTA>::getOmega() const
{
  return _omega;
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
void PorousForcedBGKdynamics<T,DESCRIPTOR,MOMENTA>::setOmega(T omega)
{
  _omega = omega;
}


//////////////////// Class PorousParticleForcedBGKdynamics ////////////////////

template<typename T, typename DESCRIPTOR, typename MOMENTA, bool isStatic>
PorousParticleForcedBGKdynamics<T,DESCRIPTOR,MOMENTA,isStatic>::PorousParticleForcedBGKdynamics (T omega_)
  : BGKdynamics<T,DESCRIPTOR,MOMENTA>(omega_)
{}

template<typename T, typename DESCRIPTOR, typename MOMENTA, bool isStatic>
void PorousParticleForcedBGKdynamics<T,DESCRIPTOR,MOMENTA,isStatic>::computeU (
  ConstCell<T,DESCRIPTOR>& cell, T u[DESCRIPTOR::d]) const
{
  T rho;
  MOMENTA().computeRhoU(cell, rho, u);
  T u_tmp[3] = {0., 0., 0.};
  if ( cell.template getFieldPointer<descriptors::VELOCITY_DENOMINATOR>()[0] > std::numeric_limits<T>::epsilon()) {
    for (int i=0; i<DESCRIPTOR::d; i++)  {
      u_tmp[i] = (1.-cell.template getField<descriptors::POROSITY>())
                 * (cell.template getFieldPointer<descriptors::VELOCITY_NUMERATOR>()[i]
                    / cell.template getField<descriptors::VELOCITY_DENOMINATOR>()
                    - u[i]);
      u[i] += rho * u_tmp[i] / (T)2.;
    }
  }
}

template<typename T, typename DESCRIPTOR, typename MOMENTA, bool isStatic>
void PorousParticleForcedBGKdynamics<T,DESCRIPTOR,MOMENTA,isStatic>::computeRhoU (
  ConstCell<T,DESCRIPTOR>& cell, T& rho, T u[DESCRIPTOR::d]) const
{
  MOMENTA().computeRhoU(cell, rho, u);
  T u_tmp[3] = {0., 0., 0.};
  if ( cell.template getFieldPointer<descriptors::VELOCITY_DENOMINATOR>()[0] > std::numeric_limits<T>::epsilon()) {
    for (int i=0; i<DESCRIPTOR::d; i++)  {
      u_tmp[i] = (1.-cell.template getField<descriptors::POROSITY>())
                 * (cell.template getFieldPointer<descriptors::VELOCITY_NUMERATOR>()[i]
                    / cell.template getField<descriptors::VELOCITY_DENOMINATOR>()
                    - u[i]);
      u[i] += rho * u_tmp[i] / (T)2.;
    }
  }
}

template<typename T, typename DESCRIPTOR, typename MOMENTA, bool isStatic>
CellStatistic<T> PorousParticleForcedBGKdynamics<T,DESCRIPTOR,MOMENTA,isStatic>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho, u[DESCRIPTOR::d];
  MOMENTA().computeRhoU(cell, rho, u);
  const T uSqr = this->porousParticleBgkCollision(cell, rho, u, this->getOmega());
  statistics.incrementStats(rho, uSqr);
}

template<typename T, typename DESCRIPTOR, typename MOMENTA, bool isStatic>
T PorousParticleForcedBGKdynamics<T,DESCRIPTOR,MOMENTA,isStatic>::porousParticleBgkCollision(
  Cell<T,DESCRIPTOR>& cell, T rho, T u[DESCRIPTOR::d], T omega)
{
#if defined(FEATURE_HLBM_GUO_FORCING)
  return GuoForcing(cell, rho, u, omega);
#elif defined(FEATURE_HLBM_SHANCHEN_FORCING)
  return ShanChenForcing(cell, rho, u, omega);
#else
  return KupershtokhForcing(cell, rho, u, omega);
#endif
}

template<typename T, typename DESCRIPTOR, typename MOMENTA, bool isStatic>
T PorousParticleForcedBGKdynamics<T,DESCRIPTOR,MOMENTA,isStatic>::KupershtokhForcing(
  Cell<T,DESCRIPTOR>& cell, const T rho, const T u[DESCRIPTOR::d], const T omega)
{
  // external force, e.g. gravity or pressure
  const auto force = cell.template getField<descriptors::FORCE>();
  const auto velDenominator = cell.template getFieldPointer<descriptors::VELOCITY_DENOMINATOR>();
  const T uSqr = lbm<DESCRIPTOR>::bgkCollision(cell, rho, u, omega);
  T uPlus[DESCRIPTOR::d] = { };

  for (int iDim=0; iDim<DESCRIPTOR::d; ++iDim) {
    uPlus[iDim] = u[iDim] + force[iDim];
  }
  if (velDenominator[0] > std::numeric_limits<T>::epsilon()) {
    this->calculate(cell, uPlus);
  }
  const T uPlusSqr = util::normSqr<T,DESCRIPTOR::d>(uPlus);
  for (int tmp_iPop=0; tmp_iPop<DESCRIPTOR::q; ++tmp_iPop) {
    cell[tmp_iPop] += lbm<DESCRIPTOR>::equilibrium(tmp_iPop, rho, uPlus, uPlusSqr)
                      - lbm<DESCRIPTOR>::equilibrium(tmp_iPop, rho, u, uSqr);
  }

  return computeUSqr(cell);
  //return uSqr;
}

template<typename T, typename DESCRIPTOR, typename MOMENTA, bool isStatic>
T PorousParticleForcedBGKdynamics<T,DESCRIPTOR,MOMENTA,isStatic>::GuoForcing(
  Cell<T,DESCRIPTOR>& cell, const T rho, T u[DESCRIPTOR::d], const T omega)
{
  const auto forceExt = cell.template getField<descriptors::FORCE>();
  const auto velDenominator = cell.template getFieldPointer<descriptors::VELOCITY_DENOMINATOR>();
  T  force[DESCRIPTOR::d] = { };
  T  uPlus[DESCRIPTOR::d] = { };

  for (int iDim=0; iDim<DESCRIPTOR::d; ++iDim) {
    uPlus[iDim] = u[iDim]+forceExt[iDim];
  }
  if (velDenominator[0] > std::numeric_limits<T>::epsilon()) {
    this->calculate(cell, uPlus);
  }
  for (int iDim=0; iDim<DESCRIPTOR::d; ++iDim) {
    force[iDim] = uPlus[iDim]-u[iDim];
    u[iDim] += force[iDim] * 0.5;
  }

  const T uSqr = lbm<DESCRIPTOR>::bgkCollision(cell, rho, u, omega);

  for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    T c_u = T();
    for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
      c_u += descriptors::c<DESCRIPTOR>(iPop,iD)*u[iD];
    }
    c_u *= descriptors::invCs2<T,DESCRIPTOR>()*descriptors::invCs2<T,DESCRIPTOR>();
    T forceTerm = T();
    for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
      forceTerm +=
        (   ((T)descriptors::c<DESCRIPTOR>(iPop,iD)-u[iD]) * descriptors::invCs2<T,DESCRIPTOR>()
            + c_u * descriptors::c<DESCRIPTOR>(iPop,iD)                                      )
        * force[iD];
    }
    forceTerm *= descriptors::t<T,DESCRIPTOR>(iPop);
    forceTerm *= T(1) - omega/T(2);
    forceTerm *= rho;
    cell[iPop] += forceTerm;
  }

  // Reset contact helper if utilized
  if constexpr (DESCRIPTOR::template provides<descriptors::COLLISION_DETECTION>()) {
    cell.template setField<descriptors::COLLISION_DETECTION>(0);
  }

  return computeUSqr(cell);
  //return uSqr;
}

template<typename T, typename DESCRIPTOR, typename MOMENTA, bool isStatic>
T PorousParticleForcedBGKdynamics<T,DESCRIPTOR,MOMENTA,isStatic>::ShanChenForcing(
  Cell<T,DESCRIPTOR>& cell, const T rho, T u[DESCRIPTOR::d], const T omega)
{
  const auto velDenominator = cell.template getFieldPointer<descriptors::VELOCITY_DENOMINATOR>();
  const auto forceExt = cell.template getField<descriptors::FORCE>();
  const T tau = 1.0/omega;
  for (int iDim=0; iDim<DESCRIPTOR::d; ++iDim) {
    u[iDim] += forceExt[iDim]*tau;
  }
  if (velDenominator[0] > std::numeric_limits<T>::epsilon()) {
    this->calculate(cell, u);
  }

  // Reset contact helper if utilized
  if constexpr (DESCRIPTOR::template provides<descriptors::COLLISION_DETECTION>()) {
    cell.template setField<descriptors::COLLISION_DETECTION>(0);
  }

  return lbm<DESCRIPTOR>::bgkCollision(cell, rho, u, omega);
}

template<typename T, typename DESCRIPTOR, typename MOMENTA, bool isStatic>
T PorousParticleForcedBGKdynamics<T,DESCRIPTOR,MOMENTA,isStatic>::computeUSqr(Cell<T,DESCRIPTOR>& cell)
{
  T u[DESCRIPTOR::d];
  MOMENTA().computeU(cell, u);
  return util::normSqr<T,DESCRIPTOR::d>(u);
}

} // olb

#endif
