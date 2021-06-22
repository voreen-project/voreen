/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006-2015 Jonas Latt, Mathias J. Krause
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
 * A collection of dynamics classes (e.g. BGK) with which a Cell object
 * can be instantiated -- generic implementation.
 */
#ifndef LB_DYNAMICS_HH
#define LB_DYNAMICS_HH

#include <type_traits>
#include <numeric>
#include "dynamics.h"
#include "core/cell.h"
#include "core/cellD.h"
#include "lbHelpers.h"
#include "firstOrderLbHelpers.h"
#include "d3q13Helpers.h"
#include "rtlbmDescriptors.h"

namespace olb {

////////////////////// Class Dynamics ////////////////////////

template<typename T, typename DESCRIPTOR>
T Dynamics<T,DESCRIPTOR>::computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr) const
{
  return lbHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho, u, uSqr);
}

template<typename T, typename DESCRIPTOR>
void Dynamics<T,DESCRIPTOR>::iniEquilibrium(Cell<T,DESCRIPTOR>& cell, T rho, const T u[DESCRIPTOR::d])
{
  T uSqr = util::normSqr<T,DESCRIPTOR::d>(u);
  for (int iPop=0; iPop<DESCRIPTOR::q; ++iPop) {
    cell[iPop] = computeEquilibrium(iPop, rho, u, uSqr);
  }
}

template<typename T, typename DESCRIPTOR>
void Dynamics<T,DESCRIPTOR>::iniRegularized(
  Cell<T,DESCRIPTOR>&cell,
  T rho,
  const T u[DESCRIPTOR::d],
  const T pi[util::TensorVal<DESCRIPTOR >::n])
{
  T uSqr = util::normSqr<T,DESCRIPTOR::d>(u);
  for (int iPop=0; iPop<DESCRIPTOR::q; ++iPop) {
    cell[iPop] = computeEquilibrium(iPop, rho, u, uSqr) + firstOrderLbHelpers<T,DESCRIPTOR>::fromPiToFneq(iPop, pi);
  }
}


template<typename T, typename DESCRIPTOR>
void Dynamics<T,DESCRIPTOR>::setBoundaryIntersection(int iPop, T distance)
{ }

template<typename T, typename DESCRIPTOR>
bool Dynamics<T,DESCRIPTOR>::getBoundaryIntersection(int iPop, T point[DESCRIPTOR::d])
{
  return 0;
}

template<typename T, typename DESCRIPTOR>
void Dynamics<T,DESCRIPTOR>::defineRho(int iPop, T rho)
{ }

template<typename T, typename DESCRIPTOR>
void Dynamics<T,DESCRIPTOR>::defineU(const T u[DESCRIPTOR::d])
{ }

template<typename T, typename DESCRIPTOR>
void Dynamics<T,DESCRIPTOR>::defineU(int iPop, const T u[DESCRIPTOR::d])
{ }

template<typename T, typename DESCRIPTOR>
T Dynamics<T,DESCRIPTOR>::getVelocityCoefficient(int iPop)
{
  return 0;
}

template <typename T, typename DESCRIPTOR>
std::string& Dynamics<T,DESCRIPTOR>::getName()
{
  return _name;
}

template <typename T, typename DESCRIPTOR>
std::string const& Dynamics<T,DESCRIPTOR>::getName() const
{
  return _name;
}

////////////////////// Class BasicDynamics ////////////////////////

template<typename T, typename DESCRIPTOR>
BasicDynamics<T,DESCRIPTOR>::BasicDynamics(Momenta<T,DESCRIPTOR>& momenta)
  : _momenta(momenta)
{
  this->getName() = "BasicDynamics";  
}

template<typename T, typename DESCRIPTOR>
T BasicDynamics<T,DESCRIPTOR>::computeRho(ConstCell<T,DESCRIPTOR>& cell) const
{
  return _momenta.computeRho(cell);
}

template<typename T, typename DESCRIPTOR>
void BasicDynamics<T,DESCRIPTOR>::computeU (
  ConstCell<T,DESCRIPTOR>& cell,
  T u[DESCRIPTOR::d]) const
{
  _momenta.computeU(cell, u);
}

template<typename T, typename DESCRIPTOR>
void BasicDynamics<T,DESCRIPTOR>::computeJ (
  ConstCell<T,DESCRIPTOR>& cell,
  T j[DESCRIPTOR::d]) const
{
  _momenta.computeJ(cell, j);
}

template<typename T, typename DESCRIPTOR>
void BasicDynamics<T,DESCRIPTOR>::computeStress (
  ConstCell<T,DESCRIPTOR>& cell,
  T rho, const T u[DESCRIPTOR::d],
  T pi[util::TensorVal<DESCRIPTOR >::n] ) const
{
  _momenta.computeStress(cell, rho, u, pi);
}

template<typename T, typename DESCRIPTOR>
void BasicDynamics<T,DESCRIPTOR>::computeRhoU (
  ConstCell<T,DESCRIPTOR>& cell,
  T& rho, T u[DESCRIPTOR::d]) const
{
  _momenta.computeRhoU(cell, rho, u);
}

template<typename T, typename DESCRIPTOR>
void BasicDynamics<T,DESCRIPTOR>::computeAllMomenta (
  ConstCell<T,DESCRIPTOR>& cell,
  T& rho, T u[DESCRIPTOR::d],
  T pi[util::TensorVal<DESCRIPTOR >::n] ) const
{
  this->computeRhoU(cell, rho, u);
  this->computeStress(cell, rho, u, pi);
}

template<typename T, typename DESCRIPTOR>
void BasicDynamics<T,DESCRIPTOR>::defineRho(Cell<T,DESCRIPTOR>& cell, T rho)
{
  _momenta.defineRho(cell, rho);
}

template<typename T, typename DESCRIPTOR>
void BasicDynamics<T,DESCRIPTOR>::defineU (
  Cell<T,DESCRIPTOR>& cell,
  const T u[DESCRIPTOR::d])
{
  _momenta.defineU(cell, u);
}

template<typename T, typename DESCRIPTOR>
void BasicDynamics<T,DESCRIPTOR>::defineRhoU (
  Cell<T,DESCRIPTOR>& cell,
  T rho, const T u[DESCRIPTOR::d])
{
  _momenta.defineRhoU(cell, rho, u);
}

template<typename T, typename DESCRIPTOR>
void BasicDynamics<T,DESCRIPTOR>::defineAllMomenta (
  Cell<T,DESCRIPTOR>& cell,
  T rho, const T u[DESCRIPTOR::d],
  const T pi[util::TensorVal<DESCRIPTOR >::n] )
{
  _momenta.defineAllMomenta(cell, rho, u, pi);
}


////////////////////// Class BGKdynamics //////////////////////////

/** \param omega relaxation parameter, related to the dynamic viscosity
 *  \param momenta a Momenta object to know how to compute velocity momenta
 */
template<typename T, typename DESCRIPTOR>
BGKdynamics<T,DESCRIPTOR>::BGKdynamics (
  T omega, Momenta<T,DESCRIPTOR>& momenta )
  : BasicDynamics<T,DESCRIPTOR>(momenta),
    _omega(omega)
{
  this->getName() = "BGKdynamics";  
}

template<typename T, typename DESCRIPTOR>
void BGKdynamics<T,DESCRIPTOR>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho, u[DESCRIPTOR::d];
  this->_momenta.computeRhoU(cell, rho, u);
  T uSqr = lbHelpers<T,DESCRIPTOR>::bgkCollision(cell, rho, u, _omega);
  statistics.incrementStats(rho, uSqr);
}

template<typename T, typename DESCRIPTOR>
T BGKdynamics<T,DESCRIPTOR>::getOmega() const
{
  return _omega;
}

template<typename T, typename DESCRIPTOR>
void BGKdynamics<T,DESCRIPTOR>::setOmega(T omega)
{
  _omega = omega;
}


////////////////////// Class TRTdynamics //////////////////////////

/** \param omega relaxation parameter, related to the dynamic viscosity
 *  \param momenta a Momenta object to know how to compute velocity momenta
 */
template<typename T, typename DESCRIPTOR>
TRTdynamics<T,DESCRIPTOR>::TRTdynamics (
  T omega, Momenta<T,DESCRIPTOR>& momenta, T magicParameter )
  : BasicDynamics<T,DESCRIPTOR>(momenta),
    _omega(omega), _omega2(1/(magicParameter/(1/omega-0.5)+0.5)), _magicParameter(magicParameter)
{
  this->getName() = "TRTdynamics";  
}

template<typename T, typename DESCRIPTOR>
void TRTdynamics<T,DESCRIPTOR>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho, u[DESCRIPTOR::d];
  this->_momenta.computeRhoU(cell, rho, u);
  // T uSqr = lbHelpers<T,DESCRIPTOR>::bgkCollision(cell, rho, u, _omega);
  const T uSqr = util::normSqr<T,DESCRIPTOR::d>(u);
  T fPlus[DESCRIPTOR::q], fMinus[DESCRIPTOR::q];
  T fEq[DESCRIPTOR::q], fEqPlus[DESCRIPTOR::q], fEqMinus[DESCRIPTOR::q];

  for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    fPlus[iPop] = 0.5 * ( cell[iPop] + cell[descriptors::opposite<DESCRIPTOR>(iPop)] );
    fMinus[iPop] = 0.5 * ( cell[iPop] - cell[descriptors::opposite<DESCRIPTOR>(iPop)] );
    fEq[iPop] = lbHelpers<T, DESCRIPTOR>::equilibrium(iPop, rho, u, uSqr);
  }
  for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    fEqPlus[iPop] = 0.5 * ( fEq[iPop] + fEq[descriptors::opposite<DESCRIPTOR>(iPop)] );
    fEqMinus[iPop] = 0.5 * ( fEq[iPop] - fEq[descriptors::opposite<DESCRIPTOR>(iPop)] );
  }
  for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    cell[iPop] -= _omega * (fPlus[iPop] - fEqPlus[iPop]) + _omega2 * (fMinus[iPop] - fEqMinus[iPop]);
  }
  statistics.incrementStats(rho, uSqr);
}

template<typename T, typename DESCRIPTOR>
T TRTdynamics<T,DESCRIPTOR>::getOmega() const
{
  return _omega;
}

template<typename T, typename DESCRIPTOR>
void TRTdynamics<T,DESCRIPTOR>::setOmega(T omega)
{
  _omega = omega;
}

////////////////////// Class ConstRhoBGKdynamics //////////////////////////

/** \param omega relaxation parameter, related to the dynamic viscosity
 *  \param momenta a Momenta object to know how to compute velocity momenta
 */
template<typename T, typename DESCRIPTOR>
ConstRhoBGKdynamics<T,DESCRIPTOR>::ConstRhoBGKdynamics (
  T omega, Momenta<T,DESCRIPTOR>& momenta )
  : BasicDynamics<T,DESCRIPTOR>(momenta),
    _omega(omega)
{
  this->getName() = "ConstRhoBGKdynamics";  
}

template<typename T, typename DESCRIPTOR>
void ConstRhoBGKdynamics<T,DESCRIPTOR>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho, u[DESCRIPTOR::d];
  this->_momenta.computeRhoU(cell, rho, u);

  T deltaRho = (T)1 - (statistics).getAverageRho();
  T ratioRho = (T)1 + deltaRho/rho;

  T uSqr = lbHelpers<T,DESCRIPTOR>::constRhoBgkCollision (
             cell, rho, u, ratioRho, _omega );
  statistics.incrementStats(rho+deltaRho, uSqr);
}

template<typename T, typename DESCRIPTOR>
T ConstRhoBGKdynamics<T,DESCRIPTOR>::getOmega() const
{
  return _omega;
}

template<typename T, typename DESCRIPTOR>
void ConstRhoBGKdynamics<T,DESCRIPTOR>::setOmega(T omega)
{
  _omega = omega;
}

////////////////////// Class IncBGKdynamics //////////////////////////

/** \param omega relaxation parameter, related to the dynamic viscosity
 *  \param momenta a Momenta object to know how to compute velocity momenta
 */
template<typename T, typename DESCRIPTOR>
IncBGKdynamics<T,DESCRIPTOR>::IncBGKdynamics (
  T omega, Momenta<T,DESCRIPTOR>& momenta )
  : BasicDynamics<T,DESCRIPTOR>(momenta), _omega(omega)
{
  this->getName() = "IncBGKdynamics";  
}

template<typename T, typename DESCRIPTOR>
void IncBGKdynamics<T,DESCRIPTOR>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho = this->_momenta.computeRho(cell);
  T p = rho / descriptors::invCs2<T,DESCRIPTOR>();
  T j[DESCRIPTOR::d];
  this->_momenta.computeJ(cell, j);
  T uSqr = lbHelpers<T,DESCRIPTOR>::incBgkCollision(cell, p, j, _omega);
  statistics.incrementStats(rho, uSqr);
}

template<typename T, typename DESCRIPTOR>
T IncBGKdynamics<T,DESCRIPTOR>::getOmega() const
{
  return _omega;
}

template<typename T, typename DESCRIPTOR>
void IncBGKdynamics<T,DESCRIPTOR>::setOmega(T omega)
{
  _omega = omega;
}



////////////////////// Class RLBdynamics /////////////////////////

/** \param omega relaxation parameter, related to the dynamic viscosity
 *  \param momenta a Momenta object to know how to compute velocity momenta
 */
template<typename T, typename DESCRIPTOR>
RLBdynamics<T,DESCRIPTOR>::RLBdynamics (
  T omega, Momenta<T,DESCRIPTOR>& momenta )
  : BasicDynamics<T,DESCRIPTOR>(momenta),
    _omega(omega)
{
  this->getName() = "RLBdynamics";  
}

template<typename T, typename DESCRIPTOR>
void RLBdynamics<T,DESCRIPTOR>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho, u[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR >::n];
  this->_momenta.computeAllMomenta(cell, rho, u, pi);
  T uSqr = rlbHelpers<T,DESCRIPTOR>::rlbCollision(cell, rho, u, pi, _omega);
  statistics.incrementStats(rho, uSqr);
}

template<typename T, typename DESCRIPTOR>
T RLBdynamics<T,DESCRIPTOR>::getOmega() const
{
  return _omega;
}

template<typename T, typename DESCRIPTOR>
void RLBdynamics<T,DESCRIPTOR>::setOmega(T omega)
{
  _omega = omega;
}

////////////////////// Class CombinedRLBdynamics /////////////////////////

template<typename T, typename DESCRIPTOR, typename Dynamics>
CombinedRLBdynamics<T,DESCRIPTOR,Dynamics>::CombinedRLBdynamics (
  T omega, Momenta<T,DESCRIPTOR>& momenta )
  : BasicDynamics<T,DESCRIPTOR>(momenta),
    _boundaryDynamics(omega, momenta)
{
  this->getName() = "CombinedRLBdynamics";  
}

template<typename T, typename DESCRIPTOR, typename Dynamics>
T CombinedRLBdynamics<T,DESCRIPTOR,Dynamics>::
computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr) const
{
  return _boundaryDynamics.computeEquilibrium(iPop, rho, u, uSqr);
}

template<typename T, typename DESCRIPTOR, typename Dynamics>
void CombinedRLBdynamics<T,DESCRIPTOR,Dynamics>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{
  typedef DESCRIPTOR L;

  T rho, u[L::d], pi[util::TensorVal<DESCRIPTOR >::n];
  this->_momenta.computeAllMomenta(cell,rho,u,pi);

  T uSqr = util::normSqr<T,L::d>(u);

  for (int iPop = 0; iPop < L::q; ++iPop) {
    cell[iPop] = computeEquilibrium(iPop,rho,u,uSqr) +
                 firstOrderLbHelpers<T,DESCRIPTOR>::fromPiToFneq(iPop, pi);
  }

  _boundaryDynamics.collide(cell, statistics);
}

template<typename T, typename DESCRIPTOR, typename Dynamics>
T CombinedRLBdynamics<T,DESCRIPTOR,Dynamics>::getOmega() const
{
  return _boundaryDynamics.getOmega();
}

template<typename T, typename DESCRIPTOR, typename Dynamics>
void CombinedRLBdynamics<T,DESCRIPTOR,Dynamics>::setOmega(T omega)
{
  _boundaryDynamics.setOmega(omega);
}


////////////////////// Class ForcedBGKdynamics /////////////////////////

/** \param omega relaxation parameter, related to the dynamic viscosity
 *  \param momenta Momenta object to know how to compute velocity momenta
 */
template<typename T, typename DESCRIPTOR>
ForcedBGKdynamics<T,DESCRIPTOR>::ForcedBGKdynamics (
  T omega, Momenta<T,DESCRIPTOR>& momenta )
  : BasicDynamics<T,DESCRIPTOR>(momenta), _omega(omega)
{
  this->getName() = "ForcedBGKdynamics";  
  OLB_PRECONDITION( DESCRIPTOR::template provides<descriptors::FORCE>() );
}

template<typename T, typename DESCRIPTOR>
void ForcedBGKdynamics<T,DESCRIPTOR>::computeU (ConstCell<T,DESCRIPTOR>& cell, T u[DESCRIPTOR::d] ) const
{
  T rho;
  this->_momenta.computeRhoU(cell, rho, u);
  auto force = cell.template getFieldPointer<descriptors::FORCE>();
  for (int iVel=0; iVel < DESCRIPTOR::d; ++iVel) {
    u[iVel] += force[iVel] / (T)2.;
  }
}

template<typename T, typename DESCRIPTOR>
void ForcedBGKdynamics<T,DESCRIPTOR>::computeRhoU (ConstCell<T,DESCRIPTOR>& cell, T& rho, T u[DESCRIPTOR::d] ) const
{
  this->_momenta.computeRhoU(cell, rho, u);
  auto force = cell.template getFieldPointer<descriptors::FORCE>();
  for (int iVel=0; iVel<DESCRIPTOR::d; ++iVel) {
    u[iVel] += force[iVel] / (T)2.;
  }
}


template<typename T, typename DESCRIPTOR>
void ForcedBGKdynamics<T,DESCRIPTOR>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho, u[DESCRIPTOR::d];
  this->_momenta.computeRhoU(cell, rho, u);
  auto force = cell.template getFieldPointer<descriptors::FORCE>();
  for (int iVel=0; iVel<DESCRIPTOR::d; ++iVel) {
    u[iVel] += force[iVel] / (T)2.;
  }
  T uSqr = lbHelpers<T,DESCRIPTOR>::bgkCollision(cell, rho, u, _omega);
  lbHelpers<T,DESCRIPTOR>::addExternalForce(cell, u, _omega, rho);
  statistics.incrementStats(rho, uSqr);
}

template<typename T, typename DESCRIPTOR>
T ForcedBGKdynamics<T,DESCRIPTOR>::getOmega() const
{
  return _omega;
}

template<typename T, typename DESCRIPTOR>
void ForcedBGKdynamics<T,DESCRIPTOR>::setOmega(T omega)
{
  _omega = omega;
}

////////////////////// Class ForcedKupershtokhBGKdynamics /////////////////////////

/** \param omega relaxation parameter, related to the dynamic viscosity
 *  \param momenta Momenta object to know how to compute velocity momenta
 */
template<typename T, typename DESCRIPTOR>
ForcedKupershtokhBGKdynamics<T,DESCRIPTOR>::ForcedKupershtokhBGKdynamics (
  T omega, Momenta<T,DESCRIPTOR>& momenta )
  : BasicDynamics<T,DESCRIPTOR>(momenta), _omega(omega)
{
  this->getName() = "ForcedKupershtokhBGKdynamics";  
  OLB_PRECONDITION( DESCRIPTOR::template provides<descriptors::FORCE>() );
}

template<typename T, typename DESCRIPTOR>
void ForcedKupershtokhBGKdynamics<T,DESCRIPTOR>::computeU (ConstCell<T,DESCRIPTOR>& cell, T u[DESCRIPTOR::d] ) const
{
  T rho;
  this->_momenta.computeRhoU(cell, rho, u);
  for (int iVel=0; iVel<DESCRIPTOR::d; ++iVel) {
    u[iVel] += cell.template getFieldPointer<descriptors::FORCE>()[iVel] / (T)2.;
  }
}

template<typename T, typename DESCRIPTOR>
void ForcedKupershtokhBGKdynamics<T,DESCRIPTOR>::computeRhoU (ConstCell<T,DESCRIPTOR>& cell, T& rho, T u[DESCRIPTOR::d] ) const
{
  this->_momenta.computeRhoU(cell, rho, u);
  for (int iVel=0; iVel<DESCRIPTOR::d; ++iVel) {
    u[iVel] += cell.template getFieldPointer<descriptors::FORCE>()[iVel] / (T)2.;
  }
}


template<typename T, typename DESCRIPTOR>
void ForcedKupershtokhBGKdynamics<T,DESCRIPTOR>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho, u[DESCRIPTOR::d], uPlusDeltaU[DESCRIPTOR::d];
  this->_momenta.computeRhoU(cell, rho, u);
  auto force = cell.template getFieldPointer<descriptors::FORCE>();
  T uSqr = lbHelpers<T,DESCRIPTOR>::bgkCollision(cell, rho, u, _omega);

  for (int iVel=0; iVel<DESCRIPTOR::d; ++iVel) {
    uPlusDeltaU[iVel] = u[iVel] + force[iVel];
  }
  const T uPlusDeltaUSqr = util::normSqr<T,DESCRIPTOR::d>(uPlusDeltaU);

  for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    cell[iPop] += lbHelpers<T, DESCRIPTOR>::equilibrium(iPop, rho, uPlusDeltaU, uPlusDeltaUSqr)
                  - lbHelpers<T, DESCRIPTOR>::equilibrium(iPop, rho, u, uSqr);
  }

  statistics.incrementStats(rho, uSqr);
}

template<typename T, typename DESCRIPTOR>
T ForcedKupershtokhBGKdynamics<T,DESCRIPTOR>::getOmega() const
{
  return _omega;
}

template<typename T, typename DESCRIPTOR>
void ForcedKupershtokhBGKdynamics<T,DESCRIPTOR>::setOmega(T omega)
{
  _omega = omega;
}


////////////////////// Class ForcedIncBGKdynamics //////////////////////////

/** \param omega relaxation parameter, related to the dynamic viscosity
 *  \param momenta a Momenta object to know how to compute velocity momenta
 */
template<typename T, typename DESCRIPTOR>
ForcedIncBGKdynamics<T,DESCRIPTOR>::ForcedIncBGKdynamics (
  T omega, Momenta<T,DESCRIPTOR>& momenta )
  : ForcedBGKdynamics<T,DESCRIPTOR>(omega, momenta)
{ }


template<typename T, typename DESCRIPTOR>
T ForcedIncBGKdynamics<T,DESCRIPTOR>::computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr) const
{
  T p = rho / descriptors::invCs2<T,DESCRIPTOR>();
  return lbDynamicsHelpers<T,DESCRIPTOR>::incEquilibrium(iPop, u, uSqr, p);
}

template<typename T, typename DESCRIPTOR>
void ForcedIncBGKdynamics<T,DESCRIPTOR>::computeU (ConstCell<T,DESCRIPTOR>& cell, T u[DESCRIPTOR::d] ) const
{
  this->_momenta.computeJ(cell, u);
  auto force = cell.template getFieldPointer<descriptors::FORCE>();
  for (int iVel=0; iVel < DESCRIPTOR::d; ++iVel) {
    u[iVel] += force[iVel] / (T)2.;
  }
}

template<typename T, typename DESCRIPTOR>
void ForcedIncBGKdynamics<T,DESCRIPTOR>::computeRhoU (ConstCell<T,DESCRIPTOR>& cell, T& rho, T u[DESCRIPTOR::d] ) const
{
  rho = this->_momenta.computeRho(cell);
  this->_momenta.computeJ(cell, u);
  auto force = cell.template getFieldPointer<descriptors::FORCE>();
  for (int iVel=0; iVel<DESCRIPTOR::d; ++iVel) {
    u[iVel] += force[iVel] / (T)2.;
  }
}

template<typename T, typename DESCRIPTOR>
void ForcedIncBGKdynamics<T,DESCRIPTOR>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho = this->_momenta.computeRho(cell);
  T p = rho / descriptors::invCs2<T,DESCRIPTOR>();
  T u[DESCRIPTOR::d];
  this->_momenta.computeJ(cell, u);
  auto force = cell.template getFieldPointer<descriptors::FORCE>();
  for (int iVel=0; iVel<DESCRIPTOR::d; ++iVel) {
    u[iVel] += force[iVel] / (T)2.;
  }
  T uSqr = lbHelpers<T,DESCRIPTOR>::incBgkCollision(cell, p, u, this->_omega);
  lbHelpers<T,DESCRIPTOR>::addExternalForce(cell, u, this->_omega, 1.0);
  statistics.incrementStats(p, uSqr);
}

template<typename T, typename DESCRIPTOR>
ExternalTauForcedIncBGKdynamics<T,DESCRIPTOR>::ExternalTauForcedIncBGKdynamics (
  T omega, Momenta<T,DESCRIPTOR>& momenta )
  : ForcedIncBGKdynamics<T,DESCRIPTOR>(omega, momenta)
{ }

template<typename T, typename DESCRIPTOR>
void ExternalTauForcedIncBGKdynamics<T,DESCRIPTOR>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{
  T omega = 1.0 / cell.template getField<descriptors::TAU_EFF>();

  T rho = this->_momenta.computeRho(cell);
  T p = rho / descriptors::invCs2<T,DESCRIPTOR>();
  T u[DESCRIPTOR::d];
  this->_momenta.computeJ(cell, u);
  auto force = cell.template getFieldPointer<descriptors::FORCE>();
  for (int iVel=0; iVel<DESCRIPTOR::d; ++iVel) {
    u[iVel] += force[iVel] / (T)2.;
  }
  T uSqr = lbHelpers<T,DESCRIPTOR>::incBgkCollision(cell, p, u, omega);
  lbHelpers<T,DESCRIPTOR>::addExternalForce(cell, u, omega, 1.0);
  statistics.incrementStats(p, uSqr);
}


/** \param omega relaxation parameter, related to the dynamic viscosity
 *  \param momenta Momenta object to know how to compute velocity momenta
 *  \param sink counterpart of a source term
 */
template<typename T, typename DESCRIPTOR>
PoissonDynamics<T,DESCRIPTOR>::PoissonDynamics (
  T omega, Momenta<T,DESCRIPTOR>& momenta, T sink)
  : BasicDynamics<T,DESCRIPTOR>(momenta), _omega(omega), _sink(sink)
{
  this->getName() = "PoissonDynamics";  
}

template<typename T, typename DESCRIPTOR>
T PoissonDynamics<T,DESCRIPTOR>::computeEquilibrium( int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr ) const
{
  return descriptors::t<T,DESCRIPTOR>(iPop) * rho - descriptors::t<T,DESCRIPTOR>(iPop);
}


template<typename T, typename DESCRIPTOR>
T PoissonDynamics<T,DESCRIPTOR>::computeRho(ConstCell<T,DESCRIPTOR>& cell) const
{
  return lbHelpers<T,DESCRIPTOR>::computeRho(cell);
}

template<typename T, typename DESCRIPTOR>
void PoissonDynamics<T,DESCRIPTOR>::computeU(ConstCell<T,DESCRIPTOR>& cell, T u[DESCRIPTOR::d]) const
{
  for ( int iDim = 0; iDim < DESCRIPTOR::d; iDim++ ) {
    u[iDim] = T();
  }
}

template<typename T, typename DESCRIPTOR>
void PoissonDynamics<T,DESCRIPTOR>::computeJ(ConstCell<T,DESCRIPTOR>& cell, T j[DESCRIPTOR::d]) const
{
  lbHelpers<T,DESCRIPTOR>::computeJ(cell, j);
}

template<typename T, typename DESCRIPTOR>
void PoissonDynamics<T,DESCRIPTOR>::computeStress(ConstCell<T,DESCRIPTOR>& cell, T rho,
    const T u[DESCRIPTOR::d], T pi[util::TensorVal<DESCRIPTOR>::n] ) const
{
  for ( int iDim = 0; iDim < DESCRIPTOR::d; iDim++ ) {
    pi[iDim] = T();
  }
}

template<typename T, typename DESCRIPTOR>
void PoissonDynamics<T,DESCRIPTOR>::computeRhoU( ConstCell<T,DESCRIPTOR>& cell, T& rho, T u[DESCRIPTOR::d]) const
{
  rho = computeRho(cell);
  computeU(cell, u);
}

template<typename T, typename DESCRIPTOR>
void PoissonDynamics<T,DESCRIPTOR>::computeAllMomenta( ConstCell<T,DESCRIPTOR>& cell, T &rho,
    T u[DESCRIPTOR::q], T pi[util::TensorVal<DESCRIPTOR>::n] ) const
{
  rho = computeRho(cell);
  computeU(cell, u);
  computeStress(cell, rho, u, pi);
}


template<typename T, typename DESCRIPTOR>
void PoissonDynamics<T,DESCRIPTOR>::collide( Cell<T,DESCRIPTOR>& cell, LatticeStatistics<T>& statistics )
{
  T rho = computeRho(cell);

  for ( int iPop = 0; iPop < descriptors::q<DESCRIPTOR>(); ++iPop ) {
    cell[iPop] =
      (cell[iPop] + descriptors::t<T,DESCRIPTOR>(iPop))
      - _omega*( (cell[iPop] + descriptors::t<T,DESCRIPTOR>(iPop)) - descriptors::t<T,DESCRIPTOR>(iPop) * rho )
      - _sink*(cell[iPop] + descriptors::t<T,DESCRIPTOR>(iPop))
      - descriptors::t<T,DESCRIPTOR>(iPop);
  }

  statistics.incrementStats(rho, 0);
}

template<typename T, typename DESCRIPTOR>
T PoissonDynamics<T,DESCRIPTOR>::getOmega() const
{
  return _omega;
}

template<typename T, typename DESCRIPTOR>
void PoissonDynamics<T,DESCRIPTOR>::setOmega(T omega)
{
  _omega = omega;
}


template<typename T, typename DESCRIPTOR>
P1Dynamics<T,DESCRIPTOR>::P1Dynamics (
  T omega, Momenta<T,DESCRIPTOR>& momenta, T absorption, T scattering)
  : BasicDynamics<T,DESCRIPTOR>(momenta), _omega(omega), _absorption(absorption), _scattering(scattering)
{
  this->getName() = "P1Dynamics";  
}

template<typename T, typename DESCRIPTOR>
T P1Dynamics<T,DESCRIPTOR>::computeEquilibrium( int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr ) const
{
  std::array<T,DESCRIPTOR::d> u_array;
  for (int iDim = 0; iDim < DESCRIPTOR::d; iDim++) {
    u_array[iDim] = u[iDim];
  }
  return lbHelpers<T,DESCRIPTOR>::equilibriumP1( iPop, rho, u_array );
}


template<typename T, typename DESCRIPTOR>
T P1Dynamics<T,DESCRIPTOR>::computeRho(ConstCell<T,DESCRIPTOR>& cell) const
{
  return lbHelpers<T,DESCRIPTOR>::computeRho(cell);
}

template<typename T, typename DESCRIPTOR>
void P1Dynamics<T,DESCRIPTOR>::computeU(ConstCell<T,DESCRIPTOR>& cell, T u[DESCRIPTOR::d]) const
{
  for ( int iDim = 0; iDim < DESCRIPTOR::d; iDim++ ) {
    u[iDim] = T();
  }
}

template<typename T, typename DESCRIPTOR>
void P1Dynamics<T,DESCRIPTOR>::computeJ(ConstCell<T,DESCRIPTOR>& cell, T j[DESCRIPTOR::d]) const
{
  std::array<T,DESCRIPTOR::q> cellShifted;
  for (int iPop = 0; iPop <DESCRIPTOR::q; ++iPop) {
    cellShifted[iPop] = cell[iPop] + descriptors::t<T,DESCRIPTOR>(iPop);
  }
  std::array<T,DESCRIPTOR::d> moment1;
  moment1.fill( T(0) );
  // sum_j v_j f_j
  for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
    for (int iDim = 0; iDim < DESCRIPTOR::d; ++iDim) {
      moment1[iDim] += descriptors::c<DESCRIPTOR>(iPop,iDim)*cellShifted[iPop];
    }
  }
  for (int iDim = 0; iDim < DESCRIPTOR::d; ++iDim) {
    j[iDim] = moment1[iDim];
  }
}

template<typename T, typename DESCRIPTOR>
void P1Dynamics<T,DESCRIPTOR>::computeStress(ConstCell<T,DESCRIPTOR>& cell, T rho,
    const T u[DESCRIPTOR::d], T pi[util::TensorVal<DESCRIPTOR>::n] ) const
{
  for ( int iDim = 0; iDim < DESCRIPTOR::d; iDim++ ) {
    pi[iDim] = T();
  }
}

template<typename T, typename DESCRIPTOR>
void P1Dynamics<T,DESCRIPTOR>::computeRhoU( ConstCell<T,DESCRIPTOR>& cell, T& rho, T u[DESCRIPTOR::d]) const
{
  rho = computeRho(cell);
  computeU(cell, u);
}

template<typename T, typename DESCRIPTOR>
void P1Dynamics<T,DESCRIPTOR>::computeAllMomenta( ConstCell<T,DESCRIPTOR>& cell, T &rho,
    T u[DESCRIPTOR::q], T pi[util::TensorVal<DESCRIPTOR>::n] ) const
{
  rho = computeRho(cell);
  computeU(cell, u);
  computeStress(cell, rho, u, pi);
}


template<typename T, typename DESCRIPTOR>
void P1Dynamics<T,DESCRIPTOR>::collide( Cell<T,DESCRIPTOR>& cell, LatticeStatistics<T>& statistics )
{
  T rho = computeRho(cell);

  // compute equilibrium f^eq_i eq. (5.32) from DOI:10.1002/9780470177013

  // give me a std::array, please
  std::array<T,DESCRIPTOR::q> cellShifted;
  for (int iPop = 0; iPop <DESCRIPTOR::q; ++iPop) {
    cellShifted[iPop] = cell[iPop] + descriptors::t<T,DESCRIPTOR>(iPop);
  }
  double moment0 = std::accumulate(cellShifted.begin(),cellShifted.end(), T(0));
  T mom1[DESCRIPTOR::d];
  computeJ(cell, mom1);
  std::array<T,DESCRIPTOR::d> moment1;
  for (int iDim = 0; iDim < DESCRIPTOR::d; ++iDim) {
    moment1[iDim] = mom1[iDim];
  }
  std::array<T,DESCRIPTOR::q> fEq;
  for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
    fEq[iPop] = lbHelpers<T,DESCRIPTOR>::equilibriumP1( iPop, moment0, moment1 ) +descriptors::t<T,DESCRIPTOR>(iPop);
  }


  std::array<T,DESCRIPTOR::q> relaxCell;
  for (int iPop = 0; iPop <DESCRIPTOR::q; ++iPop) {
    relaxCell[iPop] =
      cellShifted[iPop]
      - _scattering*descriptors::norm_c<T,DESCRIPTOR>(iPop)*(cellShifted[iPop] - fEq[iPop])
      - _absorption*descriptors::norm_c<T,DESCRIPTOR>(iPop)*cellShifted[iPop];
  }

  for (int iPop = 0; iPop <DESCRIPTOR::q; ++iPop) {
    cell[iPop] = relaxCell[iPop] - descriptors::t<T,DESCRIPTOR>(iPop);
  }

  statistics.incrementStats(rho, 0);
}

template<typename T, typename DESCRIPTOR>
T P1Dynamics<T,DESCRIPTOR>::getOmega() const
{
  return _omega;
}

template<typename T, typename DESCRIPTOR>
void P1Dynamics<T,DESCRIPTOR>::setOmega(T omega)
{
  _omega = omega;
}

////////////////////// Class ResettingForcedBGKdynamics /////////////////////////

/** \param omega relaxation parameter, related to the dynamic viscosity
 *  \param momenta a Momenta object to know how to compute velocity momenta
 */
template<typename T, typename DESCRIPTOR>
ResettingForcedBGKdynamics<T,DESCRIPTOR>::ResettingForcedBGKdynamics (
  T omega, Momenta<T,DESCRIPTOR>& momenta )
  : ForcedBGKdynamics<T,DESCRIPTOR>(omega, momenta)
{
  OLB_PRECONDITION( DESCRIPTOR::template provides<descriptors::FORCE>() );
}

template<typename T, typename DESCRIPTOR>
void ResettingForcedBGKdynamics<T,DESCRIPTOR>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho, u[DESCRIPTOR::d];
  this->_momenta.computeRhoU(cell, rho, u);
  auto force = cell.template getFieldPointer<descriptors::FORCE>();
  if ( !util::nearZero(force[0]) || !util::nearZero(force[1]) || !util::nearZero(force[2]) ) // TODO: unnecessary??
    for (int iVel=0; iVel<DESCRIPTOR::d; ++iVel) {
      u[iVel] += force[iVel] / (T)2.;
    }
//  if (force[2] != 0)
//  std::cout << force[0] << " " << force[1] << " " << force[2] << std::endl;
  T uSqr = lbHelpers<T,DESCRIPTOR>::bgkCollision(cell, rho, u, this->_omega);
  lbHelpers<T,DESCRIPTOR>::addExternalForce(cell, u, this->_omega, rho);
  statistics.incrementStats(rho, uSqr);

  force[0] = _frc[0];
  force[1] = _frc[1];
  force[2] = _frc[2];
//  force[0] = 0.;
//  force[1] = 0.;
//  force[2] = 0.;
}

////////////////////// Class ForcedShanChenBGKdynamics /////////////////////////

/** \param omega relaxation parameter, related to the dynamic viscosity
 *  \param momenta a Momenta object to know how to compute velocity momenta
 */
template<typename T, typename DESCRIPTOR>
ForcedShanChenBGKdynamics<T,DESCRIPTOR>::ForcedShanChenBGKdynamics (
  T omega, Momenta<T,DESCRIPTOR>& momenta )
  : ForcedBGKdynamics<T,DESCRIPTOR>(omega, momenta )
{
  OLB_PRECONDITION( DESCRIPTOR::template provides<descriptors::FORCE>() );
}

template<typename T, typename DESCRIPTOR>
void ForcedShanChenBGKdynamics<T,DESCRIPTOR>::computeU (ConstCell<T,DESCRIPTOR>& cell, T u[DESCRIPTOR::d] ) const
{
  T rho;
  this->_momenta.computeRhoU(cell, rho, u);
  for (int iVel=0; iVel<DESCRIPTOR::d; ++iVel) {
    u[iVel] += cell.template getFieldPointer<descriptors::FORCE>()[iVel] / (T)2.;
  }
}

template<typename T, typename DESCRIPTOR>
void ForcedShanChenBGKdynamics<T,DESCRIPTOR>::computeRhoU (ConstCell<T,DESCRIPTOR>& cell, T& rho, T u[DESCRIPTOR::d] ) const
{
  this->_momenta.computeRhoU(cell, rho, u);
  for (int iVel=0; iVel<DESCRIPTOR::d; ++iVel) {
    u[iVel] += cell.template getFieldPointer<descriptors::FORCE>()[iVel] / (T)2.;
  }
}

template<typename T, typename DESCRIPTOR>
void ForcedShanChenBGKdynamics<T,DESCRIPTOR>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho, u[DESCRIPTOR::d];
  this->_momenta.computeRhoU(cell, rho, u);
  auto force = cell.template getFieldPointer<descriptors::FORCE>();
  for (int iVel=0; iVel<DESCRIPTOR::d; ++iVel) {
    u[iVel] += force[iVel] /  this->getOmega();
  }
  T uSqr = lbHelpers<T,DESCRIPTOR>::bgkCollision(cell, rho, u, this->getOmega() );
  uSqr=0.;
  for (int iVel=0; iVel<DESCRIPTOR::d; ++iVel) {
    u[iVel] += force[iVel] / (T)2.;
    u[iVel] -= force[iVel] /  this->getOmega();
    uSqr    += u[iVel]*u[iVel];
  }
  statistics.incrementStats(rho, uSqr);
}

////////////////////// Class ForcedTRTdynamics /////////////////////////

/** \param omega relaxation parameter, related to the dynamic viscosity
 *  \param momenta Momenta object to know how to compute velocity momenta
 */
template<typename T, typename DESCRIPTOR>
ForcedTRTdynamics<T,DESCRIPTOR>::ForcedTRTdynamics (
  T omega, Momenta<T,DESCRIPTOR>& momenta, T magicParameter )
  : BasicDynamics<T,DESCRIPTOR>(momenta),
    _omega(omega), _omega2(1/(magicParameter/(1/omega-0.5)+0.5)), _magicParameter(magicParameter)
{
  this->getName() = "ForcedTRTdynamics";  
  OLB_PRECONDITION( DESCRIPTOR::template provides<descriptors::FORCE>() );
}

template<typename T, typename DESCRIPTOR>
void ForcedTRTdynamics<T,DESCRIPTOR>::computeU (ConstCell<T,DESCRIPTOR>& cell, T u[DESCRIPTOR::d] ) const
{
  T rho;
  this->_momenta.computeRhoU(cell, rho, u);
  for (int iVel=0; iVel<DESCRIPTOR::d; ++iVel) {
    u[iVel] += cell.template getFieldPointer<descriptors::FORCE>()[iVel] / (T)2.;
  }
}

template<typename T, typename DESCRIPTOR>
void ForcedTRTdynamics<T,DESCRIPTOR>::computeRhoU (ConstCell<T,DESCRIPTOR>& cell, T& rho, T u[DESCRIPTOR::d] ) const
{
  this->_momenta.computeRhoU(cell, rho, u);
  for (int iVel=0; iVel<DESCRIPTOR::d; ++iVel) {
    u[iVel] += cell.template getFieldPointer<descriptors::FORCE>()[iVel] / (T)2.;
  }
}


template<typename T, typename DESCRIPTOR>
void ForcedTRTdynamics<T,DESCRIPTOR>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho, u[DESCRIPTOR::d];
  this->_momenta.computeRhoU(cell, rho, u);
  auto force = cell.template getFieldPointer<descriptors::FORCE>();
  for (int iVel=0; iVel<DESCRIPTOR::d; ++iVel) {
    u[iVel] += force[iVel] / (T)2.;
  }
  const T uSqr = util::normSqr<T,DESCRIPTOR::d>(u);
  T fPlus[DESCRIPTOR::q], fMinus[DESCRIPTOR::q];
  T fEq[DESCRIPTOR::q], fEqPlus[DESCRIPTOR::q], fEqMinus[DESCRIPTOR::q];

  for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    fPlus[iPop] = 0.5 * ( cell[iPop] + cell[descriptors::opposite<DESCRIPTOR>(iPop)] );
    fMinus[iPop] = 0.5 * ( cell[iPop] - cell[descriptors::opposite<DESCRIPTOR>(iPop)] );
    fEq[iPop] = lbHelpers<T, DESCRIPTOR>::equilibrium(iPop, rho, u, uSqr);
  }
  for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    fEqPlus[iPop] = 0.5 * ( fEq[iPop] + fEq[descriptors::opposite<DESCRIPTOR>(iPop)] );
    fEqMinus[iPop] = 0.5 * ( fEq[iPop] - fEq[descriptors::opposite<DESCRIPTOR>(iPop)] );
  }
  for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    cell[iPop] -= _omega * (fPlus[iPop] - fEqPlus[iPop]) + _omega2 * (fMinus[iPop] - fEqMinus[iPop]);
  }
  lbHelpers<T,DESCRIPTOR>::addExternalForce(cell, u, _omega, rho);
  statistics.incrementStats(rho, uSqr);
}

template<typename T, typename DESCRIPTOR>
T ForcedTRTdynamics<T,DESCRIPTOR>::getOmega() const
{
  return _omega;
}

template<typename T, typename DESCRIPTOR>
void ForcedTRTdynamics<T,DESCRIPTOR>::setOmega(T omega)
{
  _omega = omega;
}
////////////////////// Class D3Q13dynamics /////////////////////////

/** \param omega relaxation parameter, related to the dynamic viscosity
 *  \param momenta a Momenta object to know how to compute velocity momenta
 */
template<typename T, typename DESCRIPTOR>
D3Q13dynamics<T,DESCRIPTOR>::D3Q13dynamics (
  T omega, Momenta<T,DESCRIPTOR>& momenta )
  : BasicDynamics<T,DESCRIPTOR>(momenta)
{
  this->getName() = "D3Q13dynamics";  
  setOmega(omega);
}

template<typename T, typename DESCRIPTOR>
T D3Q13dynamics<T,DESCRIPTOR>::computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr) const
{
  // To get at the equilibrium, execute collision with relaxation parameters 1
  CellD<T,DESCRIPTOR> tmp;
  for (int pop=0; pop<DESCRIPTOR::q; ++pop) {
    tmp[pop] = descriptors::t<T,DESCRIPTOR>(pop);
  }
  d3q13Helpers<T>::collision(tmp, rho, u, (T)1, (T)1);
  return tmp[iPop];
}

template<typename T, typename DESCRIPTOR>
void D3Q13dynamics<T,DESCRIPTOR>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho, u[DESCRIPTOR::d];
  this->_momenta.computeRhoU(cell, rho, u);
  T uSqr = d3q13Helpers<T>::collision (
             cell, rho, u, lambda_nu, lambda_nu_prime );
  statistics.incrementStats(rho, uSqr);
}

template<typename T, typename DESCRIPTOR>
T D3Q13dynamics<T,DESCRIPTOR>::getOmega() const
{
  return (T)4 / ( (T)3/lambda_nu + (T)1/(T)2 );
}

template<typename T, typename DESCRIPTOR>
void D3Q13dynamics<T,DESCRIPTOR>::setOmega(T omega)
{
  lambda_nu = (T)3 / ( (T)4/omega - (T)1/(T)2 );
  lambda_nu_prime = (T)3 / ( (T)2/omega + (T)1/(T)2 );
}

////////////////////// Class Momenta //////////////////////////////

template<typename T, typename DESCRIPTOR>
void Momenta<T,DESCRIPTOR>::computeRhoU (
  ConstCell<T,DESCRIPTOR>& cell,
  T& rho, T u[DESCRIPTOR::d]) const
{
  rho = this->computeRho(cell);
  this->computeU(cell, u);

}

template<typename T, typename DESCRIPTOR>
void Momenta<T,DESCRIPTOR>::computeAllMomenta (
  ConstCell<T,DESCRIPTOR>& cell,
  T& rho, T u[DESCRIPTOR::d],
  T pi[util::TensorVal<DESCRIPTOR >::n] ) const
{
  this->computeRhoU(cell, rho, u);
  this->computeStress(cell, rho, u, pi);
}

template<typename T, typename DESCRIPTOR>
void Momenta<T,DESCRIPTOR>::defineRhoU (
  Cell<T,DESCRIPTOR>& cell,
  T rho, const T u[DESCRIPTOR::d])
{
  this->defineRho(cell, rho);
  this->defineU(cell, u);

}

////////////////////// Class BulkMomenta //////////////////////////

template<typename T, typename DESCRIPTOR>
T BulkMomenta<T,DESCRIPTOR>::computeRho(ConstCell<T,DESCRIPTOR>& cell) const
{
  return lbHelpers<T,DESCRIPTOR>::computeRho(cell);
}

template<typename T, typename DESCRIPTOR>
void BulkMomenta<T,DESCRIPTOR>::computeU(ConstCell<T,DESCRIPTOR>& cell, T u[DESCRIPTOR::d]) const
{
  T dummyRho;
  lbHelpers<T,DESCRIPTOR>::computeRhoU(cell, dummyRho, u);
}

template<typename T, typename DESCRIPTOR>
void BulkMomenta<T,DESCRIPTOR>::computeJ(ConstCell<T,DESCRIPTOR>& cell, T j[DESCRIPTOR::d]) const
{
  lbHelpers<T,DESCRIPTOR>::computeJ(cell, j);
}

template<typename T, typename DESCRIPTOR>
void BulkMomenta<T,DESCRIPTOR>::computeStress (
  ConstCell<T,DESCRIPTOR>& cell,
  T rho, const T u[DESCRIPTOR::d],
  T pi[util::TensorVal<DESCRIPTOR >::n] ) const
{
  lbHelpers<T,DESCRIPTOR>::computeStress(cell, rho, u, pi);
}

template<typename T, typename DESCRIPTOR>
void BulkMomenta<T,DESCRIPTOR>::computeRhoU (
  ConstCell<T,DESCRIPTOR>& cell,
  T& rho, T u[DESCRIPTOR::d] ) const
{
  lbHelpers<T,DESCRIPTOR>::computeRhoU(cell, rho,u);
}

template<typename T, typename DESCRIPTOR>
void BulkMomenta<T,DESCRIPTOR>::computeAllMomenta (
  ConstCell<T,DESCRIPTOR>& cell,
  T& rho, T u[DESCRIPTOR::d],
  T pi[util::TensorVal<DESCRIPTOR >::n] ) const
{
  lbHelpers<T,DESCRIPTOR>::computeAllMomenta(cell, rho, u, pi);
}

template<typename T, typename DESCRIPTOR>
void BulkMomenta<T,DESCRIPTOR>::defineRho(Cell<T,DESCRIPTOR>& cell, T rho)
{
  T oldRho, u[DESCRIPTOR::d];
  computeRhoU(cell, oldRho, u);
  T uSqr = util::normSqr<T,DESCRIPTOR::d>(u);
  T fNeq[DESCRIPTOR::q];
  lbHelpers<T,DESCRIPTOR>::computeFneq(cell, fNeq, oldRho, u);
  for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    cell[iPop] = lbHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho, u, uSqr) +
                 fNeq[iPop];
  }
}

template<typename T, typename DESCRIPTOR>
void BulkMomenta<T,DESCRIPTOR>::defineU (
  Cell<T,DESCRIPTOR>& cell,
  const T u[DESCRIPTOR::d])
{
  T rho, oldU[DESCRIPTOR::d];
  computeRhoU(cell, rho, oldU);
  T uSqr = util::normSqr<T,DESCRIPTOR::d>(u);
  T fNeq[DESCRIPTOR::q];
  lbHelpers<T,DESCRIPTOR>::computeFneq(cell, fNeq, rho, oldU);
  for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    cell[iPop] = lbHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho, u, uSqr) +
                 fNeq[iPop];
  }

}

template<typename T, typename DESCRIPTOR>
void BulkMomenta<T,DESCRIPTOR>::defineRhoU (
  Cell<T,DESCRIPTOR>& cell,
  T rho, const T u[DESCRIPTOR::d])
{
  T oldRho, oldU[DESCRIPTOR::d];
  computeRhoU(cell, oldRho, oldU);
  T uSqr = util::normSqr<T,DESCRIPTOR::d>(u);
  T fNeq[DESCRIPTOR::q];
  lbHelpers<T,DESCRIPTOR>::computeFneq(cell, fNeq, oldRho, oldU);
  for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    cell[iPop] = lbHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho, u, uSqr) +
                 fNeq[iPop];
  }
}

template<typename T, typename DESCRIPTOR>
void BulkMomenta<T,DESCRIPTOR>::defineAllMomenta (
  Cell<T,DESCRIPTOR>& cell,
  T rho, const T u[DESCRIPTOR::d],
  const T pi[util::TensorVal<DESCRIPTOR >::n] )
{
  T uSqr = util::normSqr<T,DESCRIPTOR::d>(u);
  for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    cell[iPop] = lbHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho, u, uSqr) +
                 firstOrderLbHelpers<T,DESCRIPTOR>::fromPiToFneq(iPop, pi);
  }
}

////////////////////// Class ExternalVelocityMomenta //////////////////////////

template<typename T, typename DESCRIPTOR>
T ExternalVelocityMomenta<T,DESCRIPTOR>::computeRho(ConstCell<T,DESCRIPTOR>& cell) const
{
  return lbHelpers<T,DESCRIPTOR>::computeRho(cell);
}

template<typename T, typename DESCRIPTOR>
void ExternalVelocityMomenta<T,DESCRIPTOR>::computeU(ConstCell<T,DESCRIPTOR>& cell, T u[DESCRIPTOR::d]) const
{
  auto uExt = cell.template getFieldPointer<descriptors::VELOCITY>();
  for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
    u[iD] = uExt[iD];
  }
}

template<typename T, typename DESCRIPTOR>
void ExternalVelocityMomenta<T,DESCRIPTOR>::computeJ(ConstCell<T,DESCRIPTOR>& cell, T j[DESCRIPTOR::d]) const
{
  const auto rho = computeRho(cell);
  auto uExt = cell.template getFieldPointer<descriptors::VELOCITY>();
  for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
    j[iD] = uExt[iD] * rho;
  }
}

template<typename T, typename DESCRIPTOR>
void ExternalVelocityMomenta<T,DESCRIPTOR>::computeStress (
  ConstCell<T,DESCRIPTOR>& cell,
  T rho, const T u[DESCRIPTOR::d],
  T pi[util::TensorVal<DESCRIPTOR >::n] ) const
{
  lbHelpers<T,DESCRIPTOR>::computeStress(cell, rho, u, pi);
}

template<typename T, typename DESCRIPTOR>
void ExternalVelocityMomenta<T,DESCRIPTOR>::computeRhoU (
  ConstCell<T,DESCRIPTOR>& cell,
  T& rho, T u[DESCRIPTOR::d] ) const
{
  rho = computeRho(cell);
  computeU(cell,u);
}

template<typename T, typename DESCRIPTOR>
void ExternalVelocityMomenta<T,DESCRIPTOR>::computeAllMomenta (
  ConstCell<T,DESCRIPTOR>& cell,
  T& rho, T u[DESCRIPTOR::d],
  T pi[util::TensorVal<DESCRIPTOR >::n] ) const
{
  computeRhoU(cell, rho,u);
  computeStress(cell, rho, u, pi);
}

template<typename T, typename DESCRIPTOR>
void ExternalVelocityMomenta<T,DESCRIPTOR>::defineRho(Cell<T,DESCRIPTOR>& cell, T rho)
{
  T oldRho, u[DESCRIPTOR::d];
  computeRhoU(cell, oldRho, u);
  T uSqr = util::normSqr<T,DESCRIPTOR::d>(u);
  T fNeq[DESCRIPTOR::q];
  lbHelpers<T,DESCRIPTOR>::computeFneq(cell, fNeq, oldRho, u);
  for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    cell[iPop] = lbHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho, u, uSqr) +
                 fNeq[iPop];
  }
}

template<typename T, typename DESCRIPTOR>
void ExternalVelocityMomenta<T,DESCRIPTOR>::defineU (
  Cell<T,DESCRIPTOR>& cell,
  const T u[DESCRIPTOR::d])
{
  auto uExt = cell.template getFieldPointer<descriptors::VELOCITY>();
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    uExt[iD] = u[iD];
  }
}

template<typename T, typename DESCRIPTOR>
void ExternalVelocityMomenta<T,DESCRIPTOR>::defineRhoU (
  Cell<T,DESCRIPTOR>& cell,
  T rho, const T u[DESCRIPTOR::d])
{
  defineRho(cell, rho);
  defineU(cell, u);
}

template<typename T, typename DESCRIPTOR>
void ExternalVelocityMomenta<T,DESCRIPTOR>::defineAllMomenta (
  Cell<T,DESCRIPTOR>& cell,
  T rho, const T u[DESCRIPTOR::d],
  const T pi[util::TensorVal<DESCRIPTOR >::n] )
{
  defineU(cell, u);
  T uSqr = util::normSqr<T,DESCRIPTOR::d>(u);
  for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    cell[iPop] = lbHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho, u, uSqr) +
                 firstOrderLbHelpers<T,DESCRIPTOR>::fromPiToFneq(iPop, pi);
  }
}

////////////////////// Class BounceBack ///////////////////////////

template<typename T, typename DESCRIPTOR>
BounceBack<T,DESCRIPTOR>::BounceBack()
{
  this->getName() = "BounceBack";  
  _rhoFixed=false;
}

template<typename T, typename DESCRIPTOR>
BounceBack<T,DESCRIPTOR>::BounceBack(T rho)
  :_rho(rho)
{
  this->getName() = "BounceBack";    
  _rhoFixed=true;
}

template<typename T, typename DESCRIPTOR>
void BounceBack<T,DESCRIPTOR>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{
  // !do not touch element 0!
  for (int iPop=1; iPop <= DESCRIPTOR::q/2; ++iPop) {
    std::swap(cell[iPop], cell[descriptors::opposite<DESCRIPTOR>(iPop)]);
  }
}

template<typename T, typename DESCRIPTOR>
T BounceBack<T,DESCRIPTOR>::computeRho(ConstCell<T,DESCRIPTOR>& cell) const
{

  if (_rhoFixed) {
    return _rho;
  }
  return lbHelpers<T,DESCRIPTOR>::computeRho(cell);
}

template<typename T, typename DESCRIPTOR>
void BounceBack<T,DESCRIPTOR>::computeU (
  ConstCell<T,DESCRIPTOR>& cell,
  T u[DESCRIPTOR::d]) const
{
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    u[iD] = T();
  }
}

template<typename T, typename DESCRIPTOR>
void BounceBack<T,DESCRIPTOR>::computeJ (
  ConstCell<T,DESCRIPTOR>& cell,
  T j[DESCRIPTOR::d]) const
{
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    j[iD] = T();
  }
}

template<typename T, typename DESCRIPTOR>
void BounceBack<T,DESCRIPTOR>::computeStress (
  ConstCell<T,DESCRIPTOR>& cell,
  T rho, const T u[DESCRIPTOR::d],
  T pi[util::TensorVal<DESCRIPTOR >::n] ) const
{
  for (int iPi=0; iPi<util::TensorVal<DESCRIPTOR >::n; ++iPi) {
    pi[iPi] = T();//std::numeric_limits<T>::signaling_NaN();
  }
}

template<typename T, typename DESCRIPTOR>
void BounceBack<T,DESCRIPTOR>::computeRhoU (
  ConstCell<T,DESCRIPTOR>& cell,
  T& rho, T u[DESCRIPTOR::d]) const
{
  rho = computeRho(cell);
  computeU(cell, u);
}

template<typename T, typename DESCRIPTOR>
void BounceBack<T,DESCRIPTOR>::computeAllMomenta (
  ConstCell<T,DESCRIPTOR>& cell,
  T& rho, T u[DESCRIPTOR::d],
  T pi[util::TensorVal<DESCRIPTOR >::n] ) const
{
  computeRhoU(cell, rho, u);
  computeStress(cell, rho, u, pi);
}

template<typename T, typename DESCRIPTOR>
void BounceBack<T,DESCRIPTOR>::defineRho(Cell<T,DESCRIPTOR>& cell, T rho)
{ }

template<typename T, typename DESCRIPTOR>
void BounceBack<T,DESCRIPTOR>::defineU (
  Cell<T,DESCRIPTOR>& cell,
  const T u[DESCRIPTOR::d])
{ }

template<typename T, typename DESCRIPTOR>
void BounceBack<T,DESCRIPTOR>::defineRhoU (
  Cell<T,DESCRIPTOR>& cell,
  T rho, const T u[DESCRIPTOR::d])
{ }

template<typename T, typename DESCRIPTOR>
void BounceBack<T,DESCRIPTOR>::defineAllMomenta (
  Cell<T,DESCRIPTOR>& cell,
  T rho, const T u[DESCRIPTOR::d],
  const T pi[util::TensorVal<DESCRIPTOR >::n] )
{ }

template<typename T, typename DESCRIPTOR>
T BounceBack<T,DESCRIPTOR>::getOmega() const
{
  return T();//std::numeric_limits<T>::signaling_NaN();
}

template<typename T, typename DESCRIPTOR>
void BounceBack<T,DESCRIPTOR>::setOmega(T omega)
{ }


////////////////////// Class BounceBackVelocity ///////////////////////////

template<typename T, typename DESCRIPTOR>
BounceBackVelocity<T,DESCRIPTOR>::BounceBackVelocity(const T u[DESCRIPTOR::d])
{
  this->getName() = "BounceBackVelocity";  
  _rhoFixed=false;
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    _u[iD] = u[iD];
  }
}

template<typename T, typename DESCRIPTOR>
BounceBackVelocity<T,DESCRIPTOR>::BounceBackVelocity(const T rho, const T u[DESCRIPTOR::d])
  :_rho(rho)
{
  this->getName() = "BounceBackVelocity";  
  _rhoFixed=true;
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    _u[iD] = u[iD];
  }
}

template<typename T, typename DESCRIPTOR>
void BounceBackVelocity<T,DESCRIPTOR>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{
  for (int iPop=1; iPop <= DESCRIPTOR::q/2; ++iPop) {
    std::swap(cell[iPop], cell[iPop+DESCRIPTOR::q/2]);
  }
  for (int iPop=1; iPop < DESCRIPTOR::q; ++iPop) {
    for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
      cell[iPop] += computeRho(cell)*_u[iD]*descriptors::c<DESCRIPTOR>(iPop,iD)*descriptors::t<T,DESCRIPTOR>(iPop)*2*descriptors::invCs2<T,DESCRIPTOR>();
    }
  }
}

template<typename T, typename DESCRIPTOR>
T BounceBackVelocity<T,DESCRIPTOR>::computeRho(ConstCell<T,DESCRIPTOR>& cell) const
{
  if (_rhoFixed) {
    return _rho;
  }
  return lbHelpers<T,DESCRIPTOR>::computeRho(cell);
}

template<typename T, typename DESCRIPTOR>
void BounceBackVelocity<T,DESCRIPTOR>::computeU (
  ConstCell<T,DESCRIPTOR>& cell,
  T u[DESCRIPTOR::d]) const
{
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    u[iD] = _u[iD];
  }
}

template<typename T, typename DESCRIPTOR>
void BounceBackVelocity<T,DESCRIPTOR>::computeJ (
  ConstCell<T,DESCRIPTOR>& cell,
  T j[DESCRIPTOR::d]) const
{
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    j[iD] = computeRho(cell)*_u[iD];
  }
}

template<typename T, typename DESCRIPTOR>
void BounceBackVelocity<T,DESCRIPTOR>::computeStress (
  ConstCell<T,DESCRIPTOR>& cell,
  T rho, const T u[DESCRIPTOR::d],
  T pi[util::TensorVal<DESCRIPTOR >::n] ) const
{
  for (int iPi=0; iPi<util::TensorVal<DESCRIPTOR >::n; ++iPi) {
    pi[iPi] = T();//std::numeric_limits<T>::signaling_NaN();
  }
}

template<typename T, typename DESCRIPTOR>
void BounceBackVelocity<T,DESCRIPTOR>::computeRhoU (
  ConstCell<T,DESCRIPTOR>& cell,
  T& rho, T u[DESCRIPTOR::d]) const
{
  rho = computeRho(cell);
  computeU(cell, u);
}

template<typename T, typename DESCRIPTOR>
void BounceBackVelocity<T,DESCRIPTOR>::computeAllMomenta (
  ConstCell<T,DESCRIPTOR>& cell,
  T& rho, T u[DESCRIPTOR::d],
  T pi[util::TensorVal<DESCRIPTOR >::n] ) const
{
  computeRhoU(cell, rho, u);
  computeStress(cell, rho, u, pi);
}

template<typename T, typename DESCRIPTOR>
void BounceBackVelocity<T,DESCRIPTOR>::defineRho(Cell<T,DESCRIPTOR>& cell, T rho)
{ }

template<typename T, typename DESCRIPTOR>
void BounceBackVelocity<T,DESCRIPTOR>::defineU (
  Cell<T,DESCRIPTOR>& cell,
  const T u[DESCRIPTOR::d])
{
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    _u[iD] = u[iD];
  }
}

template<typename T, typename DESCRIPTOR>
void BounceBackVelocity<T,DESCRIPTOR>::defineRhoU (
  Cell<T,DESCRIPTOR>& cell,
  T rho, const T u[DESCRIPTOR::d])
{
  defineRho(cell,rho);
  defineU(cell,u);
}

template<typename T, typename DESCRIPTOR>
void BounceBackVelocity<T,DESCRIPTOR>::defineAllMomenta (
  Cell<T,DESCRIPTOR>& cell,
  T rho, const T u[DESCRIPTOR::d],
  const T pi[util::TensorVal<DESCRIPTOR >::n] )
{ }

template<typename T, typename DESCRIPTOR>
T BounceBackVelocity<T,DESCRIPTOR>::getOmega() const
{
  return T();//std::numeric_limits<T>::signaling_NaN();
}

template<typename T, typename DESCRIPTOR>
void BounceBackVelocity<T,DESCRIPTOR>::setOmega(T omega)
{ }

////////////////////// Class BounceBackAnti ///////////////////////////

template<typename T, typename DESCRIPTOR>
BounceBackAnti<T,DESCRIPTOR>::BounceBackAnti()
{
  this->getName() = "BounceBackAnti";  
  _rhoFixed = false;
  _rho = T(1);
}

template<typename T, typename DESCRIPTOR>
BounceBackAnti<T,DESCRIPTOR>::BounceBackAnti(const T rho)
  :_rho(rho)
{
  this->getName() = "BounceBackAnti";   
  _rhoFixed = true;
}

template<typename T, typename DESCRIPTOR>
void BounceBackAnti<T,DESCRIPTOR>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{
  /*
    for (int iPop=1; iPop <= DESCRIPTOR::q/2; ++iPop) {
      std::swap(cell[iPop], cell[iPop+DESCRIPTOR::q/2]);
    }
    for (int iPop=1; iPop < DESCRIPTOR::q; ++iPop) {
      if (descriptors::c<DESCRIPTOR>(iPop,0) == -1)
        cell[iPop] = -cell[descriptors::opposite<DESCRIPTOR>(iPop)] + (computeRho(cell) - T(1))*descriptors::t<T,DESCRIPTOR>(iPop)*2;
    }
  */
  //T rho, u[DESCRIPTOR::d];
  //computeRhoU(cell, rho, u);
  //T uSqr = lbHelpers<T,DESCRIPTOR>::bgkCollision(cell, rho, u, 1.78571);
  //statistics.incrementStats(rho, uSqr);

}

template<typename T, typename DESCRIPTOR>
T BounceBackAnti<T,DESCRIPTOR>::computeRho(ConstCell<T,DESCRIPTOR>& cell) const
{

  if (_rhoFixed) {
    return _rho;
  }
  return lbHelpers<T,DESCRIPTOR>::computeRho(cell);
}

template<typename T, typename DESCRIPTOR>
void BounceBackAnti<T,DESCRIPTOR>::computeU (
  ConstCell<T,DESCRIPTOR>& cell,
  T u[DESCRIPTOR::d]) const
{
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    u[iD] = T();
  }
}

template<typename T, typename DESCRIPTOR>
void BounceBackAnti<T,DESCRIPTOR>::computeJ (
  ConstCell<T,DESCRIPTOR>& cell,
  T j[DESCRIPTOR::d]) const
{
  computeU(cell, j);
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    j[iD]*=computeRho(cell);
  }
}

template<typename T, typename DESCRIPTOR>
void BounceBackAnti<T,DESCRIPTOR>::computeStress (
  ConstCell<T,DESCRIPTOR>& cell,
  T rho, const T u[DESCRIPTOR::d],
  T pi[util::TensorVal<DESCRIPTOR >::n] ) const
{
  for (int iPi=0; iPi<util::TensorVal<DESCRIPTOR >::n; ++iPi) {
    pi[iPi] = T();//std::numeric_limits<T>::signaling_NaN();
  }
}

template<typename T, typename DESCRIPTOR>
void BounceBackAnti<T,DESCRIPTOR>::computeRhoU (
  ConstCell<T,DESCRIPTOR>& cell,
  T& rho, T u[DESCRIPTOR::d]) const
{
  rho = computeRho(cell);
  computeU(cell, u);
}

template<typename T, typename DESCRIPTOR>
void BounceBackAnti<T,DESCRIPTOR>::computeAllMomenta (
  ConstCell<T,DESCRIPTOR>& cell,
  T& rho, T u[DESCRIPTOR::d],
  T pi[util::TensorVal<DESCRIPTOR >::n] ) const
{
  computeRhoU(cell, rho, u);
  computeStress(cell, rho, u, pi);
}

template<typename T, typename DESCRIPTOR>
void BounceBackAnti<T,DESCRIPTOR>::defineRho(Cell<T,DESCRIPTOR>& cell, T rho)
{
  _rho = rho;
}

template<typename T, typename DESCRIPTOR>
void BounceBackAnti<T,DESCRIPTOR>::defineU (
  Cell<T,DESCRIPTOR>& cell,
  const T u[DESCRIPTOR::d])
{
}

template<typename T, typename DESCRIPTOR>
void BounceBackAnti<T,DESCRIPTOR>::defineRhoU (
  Cell<T,DESCRIPTOR>& cell,
  T rho, const T u[DESCRIPTOR::d])
{
  defineRho(cell,rho);
  defineU(cell,u);
}

template<typename T, typename DESCRIPTOR>
void BounceBackAnti<T,DESCRIPTOR>::defineAllMomenta (
  Cell<T,DESCRIPTOR>& cell,
  T rho, const T u[DESCRIPTOR::d],
  const T pi[util::TensorVal<DESCRIPTOR >::n] )
{ }

template<typename T, typename DESCRIPTOR>
T BounceBackAnti<T,DESCRIPTOR>::getOmega() const
{
  return T();//std::numeric_limits<T>::signaling_NaN();
}

template<typename T, typename DESCRIPTOR>
void BounceBackAnti<T,DESCRIPTOR>::setOmega(T omega)
{ }


////////////////////// Class partialBounceBack ///////////////////////////

template<typename T, typename DESCRIPTOR>
PartialBounceBack<T,DESCRIPTOR>::PartialBounceBack(T rf) : _rf(rf)
{
  this->getName() = "PartialBounceBack";  
}

template<typename T, typename DESCRIPTOR>
T PartialBounceBack<T, DESCRIPTOR>::computeEquilibrium ( int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr ) const
{
  return lbHelpers<T, DESCRIPTOR>::equilibriumFirstOrder( iPop, rho, u );
}

template<typename T, typename DESCRIPTOR>
void PartialBounceBack<T,DESCRIPTOR>::collide( Cell<T,DESCRIPTOR>& cell, LatticeStatistics<T>& statistics )
{
  for (int iPop=1; iPop <= DESCRIPTOR::q/2; ++iPop) {
    std::swap(cell[iPop], cell[iPop+DESCRIPTOR::q/2]);
  }
  for (int iPop=1; iPop < DESCRIPTOR::q; ++iPop) {
    cell[iPop] = (_rf -1) * (cell[iPop] + descriptors::t<T,DESCRIPTOR>(iPop)) - descriptors::t<T,DESCRIPTOR>(iPop);
  }
}


////////////////////// Class EquilibirumBoundaryFirstOrder ///////////////////////////

template<typename T, typename DESCRIPTOR>
EquilibirumBoundaryFirstOrder<T,DESCRIPTOR>::EquilibirumBoundaryFirstOrder(const T rho, const T u[DESCRIPTOR::d])
  :_rho(rho)
{
  this->getName() = "EquilibirumBoundaryFirstOrder";  
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    _u[iD] = u[iD];
  }
}

template<typename T, typename DESCRIPTOR>
void EquilibirumBoundaryFirstOrder<T,DESCRIPTOR>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{
  const T uSqr = util::normSqr<T, DESCRIPTOR::d>(_u);

  for ( int iPop = 0; iPop < DESCRIPTOR::q; ++iPop ) {
    cell[iPop] = lbHelpers<T, DESCRIPTOR>::equilibriumFirstOrder(iPop, _rho, _u);
  }

  statistics.incrementStats(_rho, uSqr);
}

template<typename T, typename DESCRIPTOR>
T EquilibirumBoundaryFirstOrder<T,DESCRIPTOR>::computeRho(ConstCell<T,DESCRIPTOR>& cell) const
{
  return _rho;
}

template<typename T, typename DESCRIPTOR>
void EquilibirumBoundaryFirstOrder<T,DESCRIPTOR>::computeU (
  ConstCell<T,DESCRIPTOR>& cell,
  T u[DESCRIPTOR::d]) const
{
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    u[iD] = _u[iD];
  }
}

template<typename T, typename DESCRIPTOR>
void EquilibirumBoundaryFirstOrder<T,DESCRIPTOR>::computeJ (
  ConstCell<T,DESCRIPTOR>& cell,
  T j[DESCRIPTOR::d]) const
{
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    j[iD] = _rho * _u[iD];
  }
}

template<typename T, typename DESCRIPTOR>
void EquilibirumBoundaryFirstOrder<T,DESCRIPTOR>::computeStress (
  ConstCell<T,DESCRIPTOR>& cell,
  T rho, const T u[DESCRIPTOR::d],
  T pi[util::TensorVal<DESCRIPTOR >::n] ) const
{
  for (int iPi=0; iPi<util::TensorVal<DESCRIPTOR >::n; ++iPi) {
    pi[iPi] = T();//std::numeric_limits<T>::signaling_NaN();
  }
}

template<typename T, typename DESCRIPTOR>
void EquilibirumBoundaryFirstOrder<T,DESCRIPTOR>::computeRhoU (
  ConstCell<T,DESCRIPTOR>& cell,
  T& rho, T u[DESCRIPTOR::d]) const
{
  rho = computeRho(cell);
  computeU(cell, u);
}

template<typename T, typename DESCRIPTOR>
void EquilibirumBoundaryFirstOrder<T,DESCRIPTOR>::computeAllMomenta (
  ConstCell<T,DESCRIPTOR>& cell,
  T& rho, T u[DESCRIPTOR::d],
  T pi[util::TensorVal<DESCRIPTOR >::n] ) const
{
  computeRhoU(cell, rho, u);
  computeStress(cell, rho, u, pi);
}

template<typename T, typename DESCRIPTOR>
void EquilibirumBoundaryFirstOrder<T,DESCRIPTOR>::defineRho(Cell<T,DESCRIPTOR>& cell, T rho)
{
  _rho = rho;
}

template<typename T, typename DESCRIPTOR>
void EquilibirumBoundaryFirstOrder<T,DESCRIPTOR>::defineU (
  Cell<T,DESCRIPTOR>& cell,
  const T u[DESCRIPTOR::d])
{
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    _u[iD] = u[iD];
  }
}

template<typename T, typename DESCRIPTOR>
void EquilibirumBoundaryFirstOrder<T,DESCRIPTOR>::defineRhoU (
  Cell<T,DESCRIPTOR>& cell,
  T rho, const T u[DESCRIPTOR::d])
{
  defineRho(cell, rho);
  defineU(cell, u);
}

template<typename T, typename DESCRIPTOR>
void EquilibirumBoundaryFirstOrder<T,DESCRIPTOR>::defineAllMomenta (
  Cell<T,DESCRIPTOR>& cell,
  T rho, const T u[DESCRIPTOR::d],
  const T pi[util::TensorVal<DESCRIPTOR >::n] )
{ }

template<typename T, typename DESCRIPTOR>
T EquilibirumBoundaryFirstOrder<T,DESCRIPTOR>::getOmega() const
{
  return T();//std::numeric_limits<T>::signaling_NaN();
}

template<typename T, typename DESCRIPTOR>
void EquilibirumBoundaryFirstOrder<T,DESCRIPTOR>::setOmega(T omega)
{ }


////////////////////// Class EquilibirumBoundarySecondOrder ///////////////////////////

template<typename T, typename DESCRIPTOR>
EquilibirumBoundarySecondOrder<T,DESCRIPTOR>::EquilibirumBoundarySecondOrder(const T rho, const T u[DESCRIPTOR::d])
  :_rho(rho)
{
  this->getName() = "EquilibirumBoundarySecondOrder";  
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    _u[iD] = u[iD];
  }
}

template<typename T, typename DESCRIPTOR>
void EquilibirumBoundarySecondOrder<T,DESCRIPTOR>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{
  const T uSqr = util::normSqr<T, DESCRIPTOR::d>(_u);

  for ( int iPop = 0; iPop < DESCRIPTOR::q; ++iPop ) {
    cell[iPop] = lbHelpers<T, DESCRIPTOR>::equilibrium(iPop, _rho, _u, uSqr);
  }

  statistics.incrementStats(_rho, uSqr);
}

template<typename T, typename DESCRIPTOR>
T EquilibirumBoundarySecondOrder<T,DESCRIPTOR>::computeRho(ConstCell<T,DESCRIPTOR>& cell) const
{
  return _rho;
}

template<typename T, typename DESCRIPTOR>
void EquilibirumBoundarySecondOrder<T,DESCRIPTOR>::computeU (
  ConstCell<T,DESCRIPTOR>& cell,
  T u[DESCRIPTOR::d]) const
{
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    u[iD] = _u[iD];
  }
}

template<typename T, typename DESCRIPTOR>
void EquilibirumBoundarySecondOrder<T,DESCRIPTOR>::computeJ (
  ConstCell<T,DESCRIPTOR>& cell,
  T j[DESCRIPTOR::d]) const
{
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    j[iD] = _rho * _u[iD];
  }
}

template<typename T, typename DESCRIPTOR>
void EquilibirumBoundarySecondOrder<T,DESCRIPTOR>::computeStress (
  ConstCell<T,DESCRIPTOR>& cell,
  T rho, const T u[DESCRIPTOR::d],
  T pi[util::TensorVal<DESCRIPTOR >::n] ) const
{
  for (int iPi=0; iPi<util::TensorVal<DESCRIPTOR >::n; ++iPi) {
    pi[iPi] = T();//std::numeric_limits<T>::signaling_NaN();
  }
}

template<typename T, typename DESCRIPTOR>
void EquilibirumBoundarySecondOrder<T,DESCRIPTOR>::computeRhoU (
  ConstCell<T,DESCRIPTOR>& cell,
  T& rho, T u[DESCRIPTOR::d]) const
{
  rho = computeRho(cell);
  computeU(cell, u);
}

template<typename T, typename DESCRIPTOR>
void EquilibirumBoundarySecondOrder<T,DESCRIPTOR>::computeAllMomenta (
  ConstCell<T,DESCRIPTOR>& cell,
  T& rho, T u[DESCRIPTOR::d],
  T pi[util::TensorVal<DESCRIPTOR >::n] ) const
{
  computeRhoU(cell, rho, u);
  computeStress(cell, rho, u, pi);
}

template<typename T, typename DESCRIPTOR>
void EquilibirumBoundarySecondOrder<T,DESCRIPTOR>::defineRho(Cell<T,DESCRIPTOR>& cell, T rho)
{
  _rho = rho;
}

template<typename T, typename DESCRIPTOR>
void EquilibirumBoundarySecondOrder<T,DESCRIPTOR>::defineU (
  Cell<T,DESCRIPTOR>& cell,
  const T u[DESCRIPTOR::d])
{
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    _u[iD] = u[iD];
  }
}

template<typename T, typename DESCRIPTOR>
void EquilibirumBoundarySecondOrder<T,DESCRIPTOR>::defineRhoU (
  Cell<T,DESCRIPTOR>& cell,
  T rho, const T u[DESCRIPTOR::d])
{
  defineRho(cell, rho);
  defineU(cell, u);
}

template<typename T, typename DESCRIPTOR>
void EquilibirumBoundarySecondOrder<T,DESCRIPTOR>::defineAllMomenta (
  Cell<T,DESCRIPTOR>& cell,
  T rho, const T u[DESCRIPTOR::d],
  const T pi[util::TensorVal<DESCRIPTOR >::n] )
{ }

template<typename T, typename DESCRIPTOR>
T EquilibirumBoundarySecondOrder<T,DESCRIPTOR>::getOmega() const
{
  return T();//std::numeric_limits<T>::signaling_NaN();
}

template<typename T, typename DESCRIPTOR>
void EquilibirumBoundarySecondOrder<T,DESCRIPTOR>::setOmega(T omega)
{ }



////////////////////// Class NoDynamics ///////////////////////////

template<typename T, typename DESCRIPTOR>
NoDynamics<T,DESCRIPTOR>::NoDynamics(T rho) :_rho(rho)
{
  this->getName() = "NoDynamics";  
}

template<typename T, typename DESCRIPTOR>
T NoDynamics<T,DESCRIPTOR>::computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr) const
{
  return T();
}

template<typename T, typename DESCRIPTOR>
void NoDynamics<T,DESCRIPTOR>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{ }

template<typename T, typename DESCRIPTOR>
T NoDynamics<T,DESCRIPTOR>::computeRho(ConstCell<T,DESCRIPTOR>& cell) const
{
  return _rho;
}

template<typename T, typename DESCRIPTOR>
void NoDynamics<T,DESCRIPTOR>::computeU (
  ConstCell<T,DESCRIPTOR>& cell,
  T u[DESCRIPTOR::d]) const
{
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    u[iD] = T();
  }
}

template<typename T, typename DESCRIPTOR>
void NoDynamics<T,DESCRIPTOR>::computeJ (
  ConstCell<T,DESCRIPTOR>& cell,
  T j[DESCRIPTOR::d]) const
{
  for (int iD=0; iD<DESCRIPTOR::d; ++iD) {
    j[iD] = T();
  }
}

template<typename T, typename DESCRIPTOR>
void NoDynamics<T,DESCRIPTOR>::computeStress (
  ConstCell<T,DESCRIPTOR>& cell,
  T rho, const T u[DESCRIPTOR::d],
  T pi[util::TensorVal<DESCRIPTOR >::n] ) const
{
  for (int iPi=0; iPi<util::TensorVal<DESCRIPTOR >::n; ++iPi) {
    pi[iPi] = T();//std::numeric_limits<T>::signaling_NaN();
  }
}

template<typename T, typename DESCRIPTOR>
void NoDynamics<T,DESCRIPTOR>::computeRhoU (
  ConstCell<T,DESCRIPTOR>& cell,
  T& rho, T u[DESCRIPTOR::d]) const
{
  rho = computeRho(cell);
  computeU(cell, u);
}

template<typename T, typename DESCRIPTOR>
void NoDynamics<T,DESCRIPTOR>::computeAllMomenta (
  ConstCell<T,DESCRIPTOR>& cell,
  T& rho, T u[DESCRIPTOR::d],
  T pi[util::TensorVal<DESCRIPTOR >::n] ) const
{
  computeRhoU(cell, rho, u);
  computeStress(cell, rho, u, pi);
}

template<typename T, typename DESCRIPTOR>
void NoDynamics<T,DESCRIPTOR>::defineRho(Cell<T,DESCRIPTOR>& cell, T rho)
{ }

template<typename T, typename DESCRIPTOR>
void NoDynamics<T,DESCRIPTOR>::defineU (
  Cell<T,DESCRIPTOR>& cell,
  const T u[DESCRIPTOR::d])
{ }

template<typename T, typename DESCRIPTOR>
void NoDynamics<T,DESCRIPTOR>::defineRhoU (
  Cell<T,DESCRIPTOR>& cell,
  T rho, const T u[DESCRIPTOR::d])
{ }

template<typename T, typename DESCRIPTOR>
void NoDynamics<T,DESCRIPTOR>::defineAllMomenta (
  Cell<T,DESCRIPTOR>& cell,
  T rho, const T u[DESCRIPTOR::d],
  const T pi[util::TensorVal<DESCRIPTOR >::n] )
{ }

template<typename T, typename DESCRIPTOR>
T NoDynamics<T,DESCRIPTOR>::getOmega() const
{
  return T();//std::numeric_limits<T>::signaling_NaN();
}

template<typename T, typename DESCRIPTOR>
void NoDynamics<T,DESCRIPTOR>::setOmega(T omega)
{ }

////////////////////// Class offDynamics ///////////////////////////

template<typename T, typename DESCRIPTOR>
OffDynamics<T,DESCRIPTOR>::OffDynamics(const T _location[DESCRIPTOR::d])
{
  this->getName() = "OffDynamics";  
  typedef DESCRIPTOR L;
  for (int iD = 0; iD < L::d; iD++) {
    location[iD] = _location[iD];
  }
  for (int iPop = 0; iPop < L::q; iPop++) {
    distances[iPop] = -1;
    velocityCoefficient[iPop] = 0;
    for (int iD = 0; iD < L::d; iD++) {
      boundaryIntersection[iPop][iD] = _location[iD];
      _u[iPop][iD] = T();
    }
  }
  _rho=T(1);
}

template<typename T, typename DESCRIPTOR>
OffDynamics<T,DESCRIPTOR>::OffDynamics(const T _location[DESCRIPTOR::d], T _distances[DESCRIPTOR::q])
{
  this->getName() = "OffDynamics";  
  typedef DESCRIPTOR L;
  for (int iD = 0; iD < L::d; iD++) {
    location[iD] = _location[iD];
  }
  for (int iPop = 0; iPop < L::q; iPop++) {
    distances[iPop] = _distances[iPop];
    velocityCoefficient[iPop] = 0;
    for (int iD = 0; iD < L::d; iD++) {
      boundaryIntersection[iPop][iD] = _location[iD] - _distances[iPop]*descriptors::c<L>(iPop,iD);
      _u[iPop][iD] = T();
    }
  }
  _rho=T(1);
}

template<typename T, typename DESCRIPTOR>
T OffDynamics<T,DESCRIPTOR>::computeRho(ConstCell<T,DESCRIPTOR>& cell) const
{
  /*typedef DESCRIPTOR L;
  T rhoTmp = T();
  T counter = T();
  int counter2 = int();
  for (int iPop = 0; iPop < L::q; iPop++) {
    if (distances[iPop] != -1) {
      rhoTmp += (cell[iPop] + descriptors::t<T,L>(iPop))*descriptors::t<T,L>(iPop);
      counter += descriptors::t<T,L>(iPop);
      counter2++;
    }
  }
  //if (rhoTmp/counter + 1<0.1999) std::cout << rhoTmp/counter2 + 1 <<std::endl;
  //if (rhoTmp/counter + 1>1.001) std::cout << rhoTmp/counter2 + 1 <<std::endl;
  return rhoTmp/counter/counter;*/
  return _rho;
}

template<typename T, typename DESCRIPTOR>
void OffDynamics<T,DESCRIPTOR>::computeU(ConstCell<T,DESCRIPTOR>& cell, T u[DESCRIPTOR::d] ) const
{
  typedef DESCRIPTOR L;
  for (int iD = 0; iD < L::d; iD++) {
    u[iD] = T();
  }
  int counter = 0;
  for (int iPop = 0; iPop < L::q; iPop++) {
    if ( !util::nearZero(distances[iPop]+1) ) {
      for (int iD = 0; iD < L::d; iD++) {
        u[iD] += _u[iPop][iD];
      }
      counter++;
    }
  }
  if (counter!=0) {
    for (int iD = 0; iD < L::d; iD++) {
      u[iD] /= counter;
    }
  }
  return;
}

template<typename T, typename DESCRIPTOR>
void OffDynamics<T,DESCRIPTOR>::setBoundaryIntersection(int iPop, T distance)
{
  /// direction points from the fluid node into the solid domain
  /// distance is the distance from the fluid node to the solid wall
  typedef DESCRIPTOR L;
  distances[iPop] = distance;
  for (int iD = 0; iD < L::d; iD++) {
    boundaryIntersection[iPop][iD] = location[iD] - distance*descriptors::c<L>(iPop,iD);
  }
}

template<typename T, typename DESCRIPTOR>
bool OffDynamics<T,DESCRIPTOR>::getBoundaryIntersection(int iPop, T intersection[DESCRIPTOR::d])
{
  typedef DESCRIPTOR L;
  if ( !util::nearZero(distances[iPop]+1) ) {
    for (int iD = 0; iD < L::d; iD++) {
      intersection[iD] = boundaryIntersection[iPop][iD];
    }
    return true;
  }
  return false;
}

template<typename T, typename DESCRIPTOR>
void OffDynamics<T,DESCRIPTOR>::defineRho(Cell<T,DESCRIPTOR>& cell, T rho)
{
  _rho=rho;
}

template<typename T, typename DESCRIPTOR>
void OffDynamics<T,DESCRIPTOR>::defineRho(int iPop, T rho)
{
  _rho=rho;
}

template<typename T, typename DESCRIPTOR>
void OffDynamics<T,DESCRIPTOR>::defineU (
  Cell<T,DESCRIPTOR>& cell,
  const T u[DESCRIPTOR::d])
{
  defineU(u);
}

template<typename T, typename DESCRIPTOR>
void OffDynamics<T,DESCRIPTOR>::defineU(const T u[DESCRIPTOR::d])
{
  typedef DESCRIPTOR L;
  for (int iPop = 0; iPop < L::q; iPop++) {
    if ( !util::nearZero(distances[iPop]+1) ) {
      defineU(iPop, u);
    }
  }
}

/// Bouzidi velocity boundary condition formulas for the Coefficients:
/** 2*     invCs2*weight*(c,u)  for dist < 1/2
 *  1/dist*invCs2*weight*(c,u)  for dist >= 1/2
 */

template<typename T, typename DESCRIPTOR>
void OffDynamics<T,DESCRIPTOR>::defineU(
  int iPop, const T u[DESCRIPTOR::d])
{
  OLB_PRECONDITION(distances[iPop] != -1)
  typedef DESCRIPTOR L;
  velocityCoefficient[iPop] = 0;
  // scalar product of c(iPop) and u
  for (int sum = 0; sum < L::d; sum++) { // +/- problem because of first stream than postprocess
    velocityCoefficient[iPop] -= descriptors::c<L>(iPop,sum)*u[sum];
  }
  // compute summand for boundary condition
  velocityCoefficient[iPop] *= 2*descriptors::invCs2<T,L>() * descriptors::t<T,L>(iPop);

  for (int iD = 0; iD < L::d; iD++) {
    _u[iPop][iD] = u[iD];
  }
}

template<typename T, typename DESCRIPTOR>
T OffDynamics<T,DESCRIPTOR>::getVelocityCoefficient(int iPop)
{
  return velocityCoefficient[iPop];
}

////////////////////// Class ZeroDistributionDynamics ///////////////////////////

template<typename T, typename DESCRIPTOR>
ZeroDistributionDynamics<T,DESCRIPTOR>::ZeroDistributionDynamics()
{
  this->getName() = "ZeroDistributionDynamics";  
}

template<typename T, typename DESCRIPTOR>
void ZeroDistributionDynamics<T,DESCRIPTOR>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{
  for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    cell[iPop] = -descriptors::t<T,DESCRIPTOR>(iPop);
  }
}

template<typename T, typename DESCRIPTOR>
T ZeroDistributionDynamics<T,DESCRIPTOR>::computeRho(ConstCell<T,DESCRIPTOR>& cell) const
{
  return lbHelpers<T,DESCRIPTOR>::computeRho(cell);
}

/////////////// Singletons //////////////////////////////////

namespace instances {

template<typename T, typename DESCRIPTOR>
BulkMomenta<T,DESCRIPTOR>& getBulkMomenta()
{
  static BulkMomenta<T,DESCRIPTOR> bulkMomentaSingleton;
  return bulkMomentaSingleton;
}

template<typename T, typename DESCRIPTOR>
ExternalVelocityMomenta<T,DESCRIPTOR>& getExternalVelocityMomenta()
{
  static ExternalVelocityMomenta<T,DESCRIPTOR> externalVelocityMomentaSingleton;
  return externalVelocityMomentaSingleton;
}

template<typename T, typename DESCRIPTOR>
BounceBack<T,DESCRIPTOR>& getBounceBack()
{
  static BounceBack<T,DESCRIPTOR> bounceBackSingleton;
  return bounceBackSingleton;
}

template<typename T, typename DESCRIPTOR>
PartialBounceBack<T,DESCRIPTOR>& getPartialBounceBack(const double rf)
{
  static PartialBounceBack<T,DESCRIPTOR> partialBounceBackSingleton(rf);
  return partialBounceBackSingleton;
}

template<typename T, typename DESCRIPTOR>
BounceBackVelocity<T,DESCRIPTOR>& getBounceBackVelocity(const double rho, const double u[DESCRIPTOR::d])
{
  static BounceBackVelocity<T,DESCRIPTOR> bounceBackSingleton(rho,u);
  return bounceBackSingleton;
}

template<typename T, typename DESCRIPTOR>
BounceBackAnti<T,DESCRIPTOR>& getBounceBackAnti(const double rho)
{
  static BounceBackAnti<T,DESCRIPTOR> bounceBackSingleton(rho);
  return bounceBackSingleton;
}

template<typename T, typename DESCRIPTOR>
EquilibirumBoundaryFirstOrder<T,DESCRIPTOR>& getEquilibirumBoundaryFirstOrder(const double rho, const double u[DESCRIPTOR::d])
{
  static EquilibirumBoundaryFirstOrder<T,DESCRIPTOR> equilibirumBoundaryFirstOrderSingleton(rho, u);
  return equilibirumBoundaryFirstOrderSingleton;
}

template<typename T, typename DESCRIPTOR>
EquilibirumBoundarySecondOrder<T,DESCRIPTOR>& getEquilibirumBoundarySecondOrder(const double rho, const double u[DESCRIPTOR::d])
{
  static EquilibirumBoundarySecondOrder<T,DESCRIPTOR> equilibirumBoundarySecondOrderSingleton(rho, u);
  return equilibirumBoundarySecondOrderSingleton;
}

template<typename T, typename DESCRIPTOR>
NoDynamics<T,DESCRIPTOR>& getNoDynamics(T rho)
{
  static NoDynamics<T,DESCRIPTOR> noDynamicsSingleton(rho);
  return noDynamicsSingleton;
}

template<typename T, typename DESCRIPTOR>
ZeroDistributionDynamics<T,DESCRIPTOR>& getZeroDistributionDynamics()
{
  static ZeroDistributionDynamics<T,DESCRIPTOR> zeroSingleton;
  return zeroSingleton;
}

}

}

#endif
