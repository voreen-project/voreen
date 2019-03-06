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

#include <algorithm>
#include <limits>
#include "dynamics.h"
#include "core/cell.h"
#include "lbHelpers.h"
#include "firstOrderLbHelpers.h"
#include "d3q13Helpers.h"

namespace olb {

////////////////////// Class Dynamics ////////////////////////

template<typename T, template<typename U> class Lattice>
T Dynamics<T,Lattice>::computeEquilibrium(int iPop, T rho, const T u[Lattice<T>::d], T uSqr) const
{
  return lbHelpers<T,Lattice>::equilibrium(iPop, rho, u, uSqr);
}

template<typename T, template<typename U> class Lattice>
void Dynamics<T,Lattice>::iniEquilibrium(Cell<T,Lattice>& cell, T rho, const T u[Lattice<T>::d])
{
  T uSqr = util::normSqr<T,Lattice<T>::d>(u);
  for (int iPop=0; iPop<Lattice<T>::q; ++iPop) {
    cell[iPop] = computeEquilibrium(iPop, rho, u, uSqr);
  }
}

template<typename T, template<typename U> class Lattice>
void Dynamics<T,Lattice>::setBoundaryIntersection(int iPop, T distance)
{ }

template<typename T, template<typename U> class Lattice>
bool Dynamics<T,Lattice>::getBoundaryIntersection(int iPop, T point[Lattice<T>::d])
{
  return 0;
}

template<typename T, template<typename U> class Lattice>
void Dynamics<T,Lattice>::defineRho(int iPop, T rho)
{ }

template<typename T, template<typename U> class Lattice>
void Dynamics<T,Lattice>::defineU(const T u[Lattice<T>::d])
{ }

template<typename T, template<typename U> class Lattice>
void Dynamics<T,Lattice>::defineU(int iPop, const T u[Lattice<T>::d])
{ }

template<typename T, template<typename U> class Lattice>
T Dynamics<T,Lattice>::getVelocityCoefficient(int iPop)
{
  return 0;
}

////////////////////// Class BasicDynamics ////////////////////////

template<typename T, template<typename U> class Lattice>
BasicDynamics<T,Lattice>::BasicDynamics(Momenta<T,Lattice>& momenta)
  : _momenta(momenta)
{ }

template<typename T, template<typename U> class Lattice>
T BasicDynamics<T,Lattice>::computeRho(Cell<T,Lattice> const& cell) const
{
  return _momenta.computeRho(cell);
}

template<typename T, template<typename U> class Lattice>
void BasicDynamics<T,Lattice>::computeU (
  Cell<T,Lattice> const& cell,
  T u[Lattice<T>::d]) const
{
  _momenta.computeU(cell, u);
}

template<typename T, template<typename U> class Lattice>
void BasicDynamics<T,Lattice>::computeJ (
  Cell<T,Lattice> const& cell,
  T j[Lattice<T>::d]) const
{
  _momenta.computeJ(cell, j);
}

template<typename T, template<typename U> class Lattice>
void BasicDynamics<T,Lattice>::computeStress (
  Cell<T,Lattice> const& cell,
  T rho, const T u[Lattice<T>::d],
  T pi[util::TensorVal<Lattice<T> >::n] ) const
{
  _momenta.computeStress(cell, rho, u, pi);
}

template<typename T, template<typename U> class Lattice>
void BasicDynamics<T,Lattice>::computeRhoU (
  Cell<T,Lattice> const& cell,
  T& rho, T u[Lattice<T>::d]) const
{
  _momenta.computeRhoU(cell, rho, u);
}

template<typename T, template<typename U> class Lattice>
void BasicDynamics<T,Lattice>::computeAllMomenta (
  Cell<T,Lattice> const& cell,
  T& rho, T u[Lattice<T>::d],
  T pi[util::TensorVal<Lattice<T> >::n] ) const
{
  this->computeRhoU(cell, rho, u);
  this->computeStress(cell, rho, u, pi);
}

template<typename T, template<typename U> class Lattice>
void BasicDynamics<T,Lattice>::defineRho(Cell<T,Lattice>& cell, T rho)
{
  _momenta.defineRho(cell, rho);
}

template<typename T, template<typename U> class Lattice>
void BasicDynamics<T,Lattice>::defineU (
  Cell<T,Lattice>& cell,
  const T u[Lattice<T>::d])
{
  _momenta.defineU(cell, u);
}

template<typename T, template<typename U> class Lattice>
void BasicDynamics<T,Lattice>::defineRhoU (
  Cell<T,Lattice>& cell,
  T rho, const T u[Lattice<T>::d])
{
  _momenta.defineRhoU(cell, rho, u);
}

template<typename T, template<typename U> class Lattice>
void BasicDynamics<T,Lattice>::defineAllMomenta (
  Cell<T,Lattice>& cell,
  T rho, const T u[Lattice<T>::d],
  const T pi[util::TensorVal<Lattice<T> >::n] )
{
  _momenta.defineAllMomenta(cell, rho, u, pi);
}


////////////////////// Class BGKdynamics //////////////////////////

/** \param omega relaxation parameter, related to the dynamic viscosity
 *  \param momenta a Momenta object to know how to compute velocity momenta
 */
template<typename T, template<typename U> class Lattice>
BGKdynamics<T,Lattice>::BGKdynamics (
  T omega, Momenta<T,Lattice>& momenta )
  : BasicDynamics<T,Lattice>(momenta),
    _omega(omega)
{ }

template<typename T, template<typename U> class Lattice>
void BGKdynamics<T,Lattice>::collide (
  Cell<T,Lattice>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho, u[Lattice<T>::d];
  this->_momenta.computeRhoU(cell, rho, u);
  T uSqr = lbHelpers<T,Lattice>::bgkCollision(cell, rho, u, _omega);
  statistics.incrementStats(rho, uSqr);
}


template<typename T, template<typename U> class Lattice>
void BGKdynamics<T,Lattice>::staticCollide (
  Cell<T,Lattice>& cell,
  const T u[Lattice<T>::d],
  LatticeStatistics<T>& statistics )
{
  T rho = this->_momenta.computeRho(cell);
  T uSqr = lbHelpers<T,Lattice>::bgkCollision(cell, rho, u, _omega);
  statistics.incrementStats(rho, uSqr);
}

template<typename T, template<typename U> class Lattice>
T BGKdynamics<T,Lattice>::getOmega() const
{
  return _omega;
}

template<typename T, template<typename U> class Lattice>
void BGKdynamics<T,Lattice>::setOmega(T omega)
{
  _omega = omega;
}


////////////////////// Class TRTdynamics //////////////////////////

/** \param omega relaxation parameter, related to the dynamic viscosity
 *  \param momenta a Momenta object to know how to compute velocity momenta
 */
template<typename T, template<typename U> class Lattice>
TRTdynamics<T,Lattice>::TRTdynamics (
  T omega, Momenta<T,Lattice>& momenta, T magicParameter )
  : BasicDynamics<T,Lattice>(momenta),
    _omega(omega), _omega2(1/(magicParameter/(1/omega-0.5)+0.5)), _magicParameter(magicParameter)
{ }

template<typename T, template<typename U> class Lattice>
void TRTdynamics<T,Lattice>::collide (
  Cell<T,Lattice>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho, u[Lattice<T>::d];
  this->_momenta.computeRhoU(cell, rho, u);
  // T uSqr = lbHelpers<T,Lattice>::bgkCollision(cell, rho, u, _omega);
  const T uSqr = util::normSqr<T,Lattice<T>::d>(u);
  T fPlus[Lattice<T>::q], fMinus[Lattice<T>::q];
  T fEq[Lattice<T>::q], fEqPlus[Lattice<T>::q], fEqMinus[Lattice<T>::q];

  for (int iPop=0; iPop < Lattice<T>::q; ++iPop) {
    fPlus[iPop] = 0.5 * ( cell[iPop] + cell[Lattice<T>::opposite[iPop]] );
    fMinus[iPop] = 0.5 * ( cell[iPop] - cell[Lattice<T>::opposite[iPop]] );
    fEq[iPop] = lbHelpers<T, Lattice>::equilibrium(iPop, rho, u, uSqr);
  }
  for (int iPop=0; iPop < Lattice<T>::q; ++iPop) {
    fEqPlus[iPop] = 0.5 * ( fEq[iPop] + fEq[Lattice<T>::opposite[iPop]] );
    fEqMinus[iPop] = 0.5 * ( fEq[iPop] - fEq[Lattice<T>::opposite[iPop]] );
  }
  for (int iPop=0; iPop < Lattice<T>::q; ++iPop) {
    cell[iPop] -= _omega * (fPlus[iPop] - fEqPlus[iPop]) + _omega2 * (fMinus[iPop] - fEqMinus[iPop]);
  }
  statistics.incrementStats(rho, uSqr);
}


template<typename T, template<typename U> class Lattice>
void TRTdynamics<T,Lattice>::staticCollide (
  Cell<T,Lattice>& cell,
  const T u[Lattice<T>::d],
  LatticeStatistics<T>& statistics )
{
  T rho = this->_momenta.computeRho(cell);
  T uSqr = lbHelpers<T,Lattice>::bgkCollision(cell, rho, u, _omega);
  statistics.incrementStats(rho, uSqr);
}

template<typename T, template<typename U> class Lattice>
T TRTdynamics<T,Lattice>::getOmega() const
{
  return _omega;
}

template<typename T, template<typename U> class Lattice>
void TRTdynamics<T,Lattice>::setOmega(T omega)
{
  _omega = omega;
}

////////////////////// Class ConstRhoBGKdynamics //////////////////////////

/** \param omega relaxation parameter, related to the dynamic viscosity
 *  \param momenta a Momenta object to know how to compute velocity momenta
 */
template<typename T, template<typename U> class Lattice>
ConstRhoBGKdynamics<T,Lattice>::ConstRhoBGKdynamics (
  T omega, Momenta<T,Lattice>& momenta )
  : BasicDynamics<T,Lattice>(momenta),
    _omega(omega)
{ }

template<typename T, template<typename U> class Lattice>
void ConstRhoBGKdynamics<T,Lattice>::collide (
  Cell<T,Lattice>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho, u[Lattice<T>::d];
  this->_momenta.computeRhoU(cell, rho, u);

  T deltaRho = (T)1 - (statistics).getAverageRho();
  T ratioRho = (T)1 + deltaRho/rho;

  T uSqr = lbHelpers<T,Lattice>::constRhoBgkCollision (
             cell, rho, u, ratioRho, _omega );
  statistics.incrementStats(rho+deltaRho, uSqr);
}

template<typename T, template<typename U> class Lattice>
void ConstRhoBGKdynamics<T,Lattice>::staticCollide (
  Cell<T,Lattice>& cell,
  const T u[Lattice<T>::d],
  LatticeStatistics<T>& statistics )
{
  T rho = this->_momenta.computeRho(cell);
  T uSqr = lbHelpers<T,Lattice>::bgkCollision(cell, rho, u, _omega);
  statistics.incrementStats(rho, uSqr);
}

template<typename T, template<typename U> class Lattice>
T ConstRhoBGKdynamics<T,Lattice>::getOmega() const
{
  return _omega;
}

template<typename T, template<typename U> class Lattice>
void ConstRhoBGKdynamics<T,Lattice>::setOmega(T omega)
{
  _omega = omega;
}

////////////////////// Class IncBGKdynamics //////////////////////////

/** \param omega relaxation parameter, related to the dynamic viscosity
 *  \param momenta a Momenta object to know how to compute velocity momenta
 */
template<typename T, template<typename U> class Lattice>
IncBGKdynamics<T,Lattice>::IncBGKdynamics (
  T omega, Momenta<T,Lattice>& momenta )
  : BasicDynamics<T,Lattice>(momenta), _omega(omega)
{ }

template<typename T, template<typename U> class Lattice>
void IncBGKdynamics<T,Lattice>::collide (
  Cell<T,Lattice>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho = this->_momenta.computeRho(cell);
  T p = rho / Lattice<T>::invCs2;
  T j[Lattice<T>::d];
  this->_momenta.computeJ(cell, j);
  T uSqr = lbHelpers<T,Lattice>::incBgkCollision(cell, p, j, _omega);
  statistics.incrementStats(rho, uSqr);
}

template<typename T, template<typename U> class Lattice>
void IncBGKdynamics<T,Lattice>::staticCollide (
  Cell<T,Lattice>& cell,
  const T j[Lattice<T>::d],
  LatticeStatistics<T>& statistics )
{
  T rho = this->_momenta.computeRho(cell);
  T p = rho / Lattice<T>::invCs2;
  T uSqr = lbHelpers<T,Lattice>::incBgkCollision(cell, p, j, _omega);
  statistics.incrementStats(rho, uSqr);
}

template<typename T, template<typename U> class Lattice>
T IncBGKdynamics<T,Lattice>::getOmega() const
{
  return _omega;
}

template<typename T, template<typename U> class Lattice>
void IncBGKdynamics<T,Lattice>::setOmega(T omega)
{
  _omega = omega;
}



////////////////////// Class RLBdynamics /////////////////////////

/** \param omega relaxation parameter, related to the dynamic viscosity
 *  \param momenta a Momenta object to know how to compute velocity momenta
 */
template<typename T, template<typename U> class Lattice>
RLBdynamics<T,Lattice>::RLBdynamics (
  T omega, Momenta<T,Lattice>& momenta )
  : BasicDynamics<T,Lattice>(momenta),
    _omega(omega)
{ }

template<typename T, template<typename U> class Lattice>
void RLBdynamics<T,Lattice>::collide (
  Cell<T,Lattice>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho, u[Lattice<T>::d], pi[util::TensorVal<Lattice<T> >::n];
  this->_momenta.computeAllMomenta(cell, rho, u, pi);
  T uSqr = rlbHelpers<T,Lattice>::rlbCollision(cell, rho, u, pi, _omega);
  statistics.incrementStats(rho, uSqr);
}

template<typename T, template<typename U> class Lattice>
void RLBdynamics<T,Lattice>::staticCollide (
  Cell<T,Lattice>& cell,
  const T u[Lattice<T>::d],
  LatticeStatistics<T>& statistics )
{
  T rho, uDummy[Lattice<T>::d], pi[util::TensorVal<Lattice<T> >::n];
  this->_momenta.computeAllMomenta(cell, rho, uDummy, pi);
  T uSqr = rlbHelpers<T,Lattice>::rlbCollision(cell, rho, u, pi, _omega);
  statistics.incrementStats(rho, uSqr);
}

template<typename T, template<typename U> class Lattice>
T RLBdynamics<T,Lattice>::getOmega() const
{
  return _omega;
}

template<typename T, template<typename U> class Lattice>
void RLBdynamics<T,Lattice>::setOmega(T omega)
{
  _omega = omega;
}

////////////////////// Class CombinedRLBdynamics /////////////////////////

template<typename T, template<typename U> class Lattice, typename Dynamics>
CombinedRLBdynamics<T,Lattice,Dynamics>::CombinedRLBdynamics (
  T omega, Momenta<T,Lattice>& momenta )
  : BasicDynamics<T,Lattice>(momenta),
    _boundaryDynamics(omega, momenta)
{ }

template<typename T, template<typename U> class Lattice, typename Dynamics>
T CombinedRLBdynamics<T,Lattice,Dynamics>::
computeEquilibrium(int iPop, T rho, const T u[Lattice<T>::d], T uSqr) const
{
  return _boundaryDynamics.computeEquilibrium(iPop, rho, u, uSqr);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void CombinedRLBdynamics<T,Lattice,Dynamics>::collide (
  Cell<T,Lattice>& cell,
  LatticeStatistics<T>& statistics )
{
  typedef Lattice<T> L;

  T rho, u[L::d], pi[util::TensorVal<Lattice<T> >::n];
  this->_momenta.computeAllMomenta(cell,rho,u,pi);

  T uSqr = util::normSqr<T,L::d>(u);

  for (int iPop = 0; iPop < L::q; ++iPop) {
    cell[iPop] = computeEquilibrium(iPop,rho,u,uSqr) +
                 firstOrderLbHelpers<T,Lattice>::fromPiToFneq(iPop, pi);
  }

  _boundaryDynamics.collide(cell, statistics);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void CombinedRLBdynamics<T,Lattice,Dynamics>::staticCollide (
  Cell<T,Lattice>& cell,
  const T u[Lattice<T>::d],
  LatticeStatistics<T>& statistics )
{
  typedef Lattice<T> L;

  T rho, falseU[L::d], pi[util::TensorVal<Lattice<T> >::n];
  this->_momenta.computeAllMomenta(cell, rho, falseU, pi);

  T uSqr = util::normSqr<T,L::d>(u);

  for (int iPop = 0; iPop < L::q; ++iPop) {
    cell[iPop] = computeEquilibrium(iPop,rho,u,uSqr) +
                 firstOrderLbHelpers<T,Lattice>::fromPiToFneq(iPop, pi);
  }

  _boundaryDynamics.staticCollide(cell, u, statistics);
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
T CombinedRLBdynamics<T,Lattice,Dynamics>::getOmega() const
{
  return _boundaryDynamics.getOmega();
}

template<typename T, template<typename U> class Lattice, typename Dynamics>
void CombinedRLBdynamics<T,Lattice,Dynamics>::setOmega(T omega)
{
  _boundaryDynamics.setOmega(omega);
}


////////////////////// Class ForcedBGKdynamics /////////////////////////

/** \param omega relaxation parameter, related to the dynamic viscosity
 *  \param momenta Momenta object to know how to compute velocity momenta
 */
template<typename T, template<typename U> class Lattice>
ForcedBGKdynamics<T,Lattice>::ForcedBGKdynamics (
  T omega, Momenta<T,Lattice>& momenta )
  : BasicDynamics<T,Lattice>(momenta), _omega(omega)
{
  // This ensures both that the constant sizeOfForce is defined in
  // ExternalField and that it has the proper size
  OLB_PRECONDITION( Lattice<T>::d == Lattice<T>::ExternalField::sizeOfForce );
}

template<typename T, template<typename U> class Lattice>
void ForcedBGKdynamics<T,Lattice>::computeU (Cell<T,Lattice> const& cell, T u[Lattice<T>::d] ) const
{
  T rho;
  this->_momenta.computeRhoU(cell, rho, u);
  for (int iVel=0; iVel<Lattice<T>::d; ++iVel) {
    u[iVel] += cell.getExternal(forceBeginsAt)[iVel] / (T)2.;
  }
}

template<typename T, template<typename U> class Lattice>
void ForcedBGKdynamics<T,Lattice>::computeRhoU (Cell<T,Lattice> const& cell, T& rho, T u[Lattice<T>::d] ) const
{
  this->_momenta.computeRhoU(cell, rho, u);
  for (int iVel=0; iVel<Lattice<T>::d; ++iVel) {
    u[iVel] += cell.getExternal(forceBeginsAt)[iVel] / (T)2.;
  }
}


template<typename T, template<typename U> class Lattice>
void ForcedBGKdynamics<T,Lattice>::collide (
  Cell<T,Lattice>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho, u[Lattice<T>::d];
  this->_momenta.computeRhoU(cell, rho, u);
  T* force = cell.getExternal(forceBeginsAt);
  if ( !util::nearZero(force[0]) || !util::nearZero(force[1]) || !util::nearZero(force[2]) ) // TODO: unnecessary??
    for (int iVel=0; iVel<Lattice<T>::d; ++iVel) {
      u[iVel] += force[iVel] / (T)2.;
    }
//  if (force[2] != 0)
//  std::cout << force[0] << " " << force[1] << " " << force[2] << std::endl;
  T uSqr = lbHelpers<T,Lattice>::bgkCollision(cell, rho, u, _omega);
  lbHelpers<T,Lattice>::addExternalForce(cell, u, _omega, rho);
  statistics.incrementStats(rho, uSqr);
}

template<typename T, template<typename U> class Lattice>
void ForcedBGKdynamics<T,Lattice>::staticCollide (
  Cell<T,Lattice>& cell,
  const T u[Lattice<T>::d],
  LatticeStatistics<T>& statistics )
{
  T rho, uDummy[Lattice<T>::d];
  this->_momenta.computeRhoU(cell, rho, uDummy);

  T uSqr = lbHelpers<T,Lattice>::bgkCollision(cell, rho, u, _omega);
  lbHelpers<T,Lattice>::addExternalForce(cell, u, _omega, rho);
  statistics.incrementStats(rho, uSqr);
}

template<typename T, template<typename U> class Lattice>
T ForcedBGKdynamics<T,Lattice>::getOmega() const
{
  return _omega;
}

template<typename T, template<typename U> class Lattice>
void ForcedBGKdynamics<T,Lattice>::setOmega(T omega)
{
  _omega = omega;
}

////////////////////// Class ForcedKupershtockhBGKdynamics /////////////////////////

/** \param omega relaxation parameter, related to the dynamic viscosity
 *  \param momenta Momenta object to know how to compute velocity momenta
 */
template<typename T, template<typename U> class Lattice>
ForcedKupershtokhBGKdynamics<T,Lattice>::ForcedKupershtokhBGKdynamics (
  T omega, Momenta<T,Lattice>& momenta )
  : BasicDynamics<T,Lattice>(momenta), _omega(omega)
{
  // This ensures both that the constant sizeOfForce is defined in
  // ExternalField and that it has the proper size
  OLB_PRECONDITION( Lattice<T>::d == Lattice<T>::ExternalField::sizeOfForce );
}

template<typename T, template<typename U> class Lattice>
void ForcedKupershtokhBGKdynamics<T,Lattice>::computeU (Cell<T,Lattice> const& cell, T u[Lattice<T>::d] ) const
{
  T rho;
  this->_momenta.computeRhoU(cell, rho, u);
  for (int iVel=0; iVel<Lattice<T>::d; ++iVel) {
    u[iVel] += cell.getExternal(forceBeginsAt)[iVel] / (T)2.;
  }
}

template<typename T, template<typename U> class Lattice>
void ForcedKupershtokhBGKdynamics<T,Lattice>::computeRhoU (Cell<T,Lattice> const& cell, T& rho, T u[Lattice<T>::d] ) const
{
  this->_momenta.computeRhoU(cell, rho, u);
  for (int iVel=0; iVel<Lattice<T>::d; ++iVel) {
    u[iVel] += cell.getExternal(forceBeginsAt)[iVel] / (T)2.;
  }
}


template<typename T, template<typename U> class Lattice>
void ForcedKupershtokhBGKdynamics<T,Lattice>::collide (
  Cell<T,Lattice>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho, u[Lattice<T>::d], uPlusDeltaU[Lattice<T>::d];
  this->_momenta.computeRhoU(cell, rho, u);
  T* force = cell.getExternal(forceBeginsAt);
  T uSqr = lbHelpers<T,Lattice>::bgkCollision(cell, rho, u, _omega);

  for (int iVel=0; iVel<Lattice<T>::d; ++iVel) {
    uPlusDeltaU[iVel] = u[iVel] + force[iVel];
  }
  const T uPlusDeltaUSqr = util::normSqr<T,Lattice<T>::d>(uPlusDeltaU);

  for (int iPop=0; iPop < Lattice<T>::q; ++iPop) {
    cell[iPop] += lbHelpers<T, Lattice>::equilibrium(iPop, rho, uPlusDeltaU, uPlusDeltaUSqr)
                - lbHelpers<T, Lattice>::equilibrium(iPop, rho, u, uSqr);
  }

  statistics.incrementStats(rho, uSqr);
}

template<typename T, template<typename U> class Lattice>
void ForcedKupershtokhBGKdynamics<T,Lattice>::staticCollide (
  Cell<T,Lattice>& cell,
  const T u[Lattice<T>::d],
  LatticeStatistics<T>& statistics )
{
  T rho, uDummy[Lattice<T>::d];
  this->_momenta.computeRhoU(cell, rho, uDummy);

  T uSqr = lbHelpers<T,Lattice>::bgkCollision(cell, rho, u, _omega);
  lbHelpers<T,Lattice>::addExternalForce(cell, u, _omega, rho);
  statistics.incrementStats(rho, uSqr);
}

template<typename T, template<typename U> class Lattice>
T ForcedKupershtokhBGKdynamics<T,Lattice>::getOmega() const
{
  return _omega;
}

template<typename T, template<typename U> class Lattice>
void ForcedKupershtokhBGKdynamics<T,Lattice>::setOmega(T omega)
{
  _omega = omega;
}

////////////////////// Class ResettingForcedBGKdynamics /////////////////////////

/** \param omega relaxation parameter, related to the dynamic viscosity
 *  \param momenta a Momenta object to know how to compute velocity momenta
 */
template<typename T, template<typename U> class Lattice>
ResettingForcedBGKdynamics<T,Lattice>::ResettingForcedBGKdynamics (
  T omega, Momenta<T,Lattice>& momenta )
  : ForcedBGKdynamics<T,Lattice>(omega, momenta)
{
  // This ensures both that the constant sizeOfForce is defined in
  // ExternalField and that it has the proper size
  OLB_PRECONDITION( Lattice<T>::d == Lattice<T>::ExternalField::sizeOfForce );
}

template<typename T, template<typename U> class Lattice>
void ResettingForcedBGKdynamics<T,Lattice>::collide (
  Cell<T,Lattice>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho, u[Lattice<T>::d];
  this->_momenta.computeRhoU(cell, rho, u);
  T* force = cell.getExternal(this->forceBeginsAt);
  if ( !util::nearZero(force[0]) || !util::nearZero(force[1]) || !util::nearZero(force[2]) ) // TODO: unnecessary??
    for (int iVel=0; iVel<Lattice<T>::d; ++iVel) {
      u[iVel] += force[iVel] / (T)2.;
    }
//  if (force[2] != 0)
//  std::cout << force[0] << " " << force[1] << " " << force[2] << std::endl;
  T uSqr = lbHelpers<T,Lattice>::bgkCollision(cell, rho, u, this->_omega);
  lbHelpers<T,Lattice>::addExternalForce(cell, u, this->_omega, rho);
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
template<typename T, template<typename U> class Lattice>
ForcedShanChenBGKdynamics<T,Lattice>::ForcedShanChenBGKdynamics (
  T omega, Momenta<T,Lattice>& momenta )
  : ForcedBGKdynamics<T,Lattice>(omega, momenta )
{
  // This ensures both that the constant sizeOfForce is defined in
  // ExternalField and that it has the proper size
  OLB_PRECONDITION( Lattice<T>::d == Lattice<T>::ExternalField::sizeOfForce );
}

template<typename T, template<typename U> class Lattice>
void ForcedShanChenBGKdynamics<T,Lattice>::computeU (Cell<T,Lattice> const& cell, T u[Lattice<T>::d] ) const
{
  T rho;
  this->_momenta.computeRhoU(cell, rho, u);
  for (int iVel=0; iVel<Lattice<T>::d; ++iVel) {
    u[iVel] += cell.getExternal(this->forceBeginsAt)[iVel] / (T)2.;
  }
}

template<typename T, template<typename U> class Lattice>
void ForcedShanChenBGKdynamics<T,Lattice>::computeRhoU (Cell<T,Lattice> const& cell, T& rho, T u[Lattice<T>::d] ) const
{
  this->_momenta.computeRhoU(cell, rho, u);
  for (int iVel=0; iVel<Lattice<T>::d; ++iVel) {
    u[iVel] += cell.getExternal(this->forceBeginsAt)[iVel] / (T)2.;
  }
}

template<typename T, template<typename U> class Lattice>
void ForcedShanChenBGKdynamics<T,Lattice>::collide (
  Cell<T,Lattice>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho, u[Lattice<T>::d];
  this->_momenta.computeRhoU(cell, rho, u);
  T* force = cell.getExternal(this->forceBeginsAt);
  for (int iVel=0; iVel<Lattice<T>::d; ++iVel) {
    u[iVel] += force[iVel] /  this->getOmega();
  }
  T uSqr = lbHelpers<T,Lattice>::bgkCollision(cell, rho, u, this->getOmega() );
  uSqr=0.;
  for (int iVel=0; iVel<Lattice<T>::d; ++iVel) {
    u[iVel] += force[iVel] / (T)2.;
    u[iVel] -= force[iVel] /  this->getOmega();
    uSqr    += pow(u[iVel],Lattice<T>::d);
  }
  statistics.incrementStats(rho, uSqr);
}

////////////////////// Class D3Q13dynamics /////////////////////////

/** \param omega relaxation parameter, related to the dynamic viscosity
 *  \param momenta a Momenta object to know how to compute velocity momenta
 */
template<typename T, template<typename U> class Lattice>
D3Q13dynamics<T,Lattice>::D3Q13dynamics (
  T omega, Momenta<T,Lattice>& momenta )
  : BasicDynamics<T,Lattice>(momenta)
{
  setOmega(omega);
}

template<typename T, template<typename U> class Lattice>
T D3Q13dynamics<T,Lattice>::computeEquilibrium(int iPop, T rho, const T u[Lattice<T>::d], T uSqr) const
{
  // To get at the equilibrium, execute collision with relaxation parameters 1
  Cell<T,Lattice> tmp;
  for (int pop=0; pop<Lattice<T>::q; ++pop) {
    tmp[pop] = Lattice<T>::t[pop];
  }
  d3q13Helpers<T>::collision(tmp, rho, u, (T)1, (T)1);
  return tmp[iPop];
}

template<typename T, template<typename U> class Lattice>
void D3Q13dynamics<T,Lattice>::collide (
  Cell<T,Lattice>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho, u[Lattice<T>::d];
  this->_momenta.computeRhoU(cell, rho, u);
  T uSqr = d3q13Helpers<T>::collision (
             cell, rho, u, lambda_nu, lambda_nu_prime );
  statistics.incrementStats(rho, uSqr);
}

template<typename T, template<typename U> class Lattice>
void D3Q13dynamics<T,Lattice>::staticCollide (
  Cell<T,Lattice>& cell,
  const T u[Lattice<T>::d],
  LatticeStatistics<T>& statistics )
{
  T rho = this->_momenta.computeRho(cell);
  T uSqr = d3q13Helpers<T>::collision (
             cell, rho, u, lambda_nu, lambda_nu_prime );
  statistics.incrementStats(rho, uSqr);
}

template<typename T, template<typename U> class Lattice>
T D3Q13dynamics<T,Lattice>::getOmega() const
{
  return (T)4 / ( (T)3/lambda_nu + (T)1/(T)2 );
}

template<typename T, template<typename U> class Lattice>
void D3Q13dynamics<T,Lattice>::setOmega(T omega)
{
  lambda_nu = (T)3 / ( (T)4/omega - (T)1/(T)2 );
  lambda_nu_prime = (T)3 / ( (T)2/omega + (T)1/(T)2 );
}

////////////////////// Class Momenta //////////////////////////////

template<typename T, template<typename U> class Lattice>
void Momenta<T,Lattice>::computeRhoU (
  Cell<T,Lattice> const& cell,
  T& rho, T u[Lattice<T>::d]) const
{
  rho = this->computeRho(cell);
  this->computeU(cell, u);

}

template<typename T, template<typename U> class Lattice>
void Momenta<T,Lattice>::computeAllMomenta (
  Cell<T,Lattice> const& cell,
  T& rho, T u[Lattice<T>::d],
  T pi[util::TensorVal<Lattice<T> >::n] ) const
{
  this->computeRhoU(cell, rho, u);
  this->computeStress(cell, rho, u, pi);
}

template<typename T, template<typename U> class Lattice>
void Momenta<T,Lattice>::defineRhoU (
  Cell<T,Lattice>& cell,
  T rho, const T u[Lattice<T>::d])
{
  this->defineRho(cell, rho);
  this->defineU(cell, u);

}

////////////////////// Class BulkMomenta //////////////////////////

template<typename T, template<typename U> class Lattice>
T BulkMomenta<T,Lattice>::computeRho(Cell<T,Lattice> const& cell) const
{
  return lbHelpers<T,Lattice>::computeRho(cell);
}

template<typename T, template<typename U> class Lattice>
void BulkMomenta<T,Lattice>::computeU(Cell<T,Lattice> const& cell, T u[Lattice<T>::d]) const
{
  T dummyRho;
  lbHelpers<T,Lattice>::computeRhoU(cell, dummyRho, u);
}

template<typename T, template<typename U> class Lattice>
void BulkMomenta<T,Lattice>::computeJ(Cell<T,Lattice> const& cell, T j[Lattice<T>::d]) const
{
  lbHelpers<T,Lattice>::computeJ(cell, j);
}

template<typename T, template<typename U> class Lattice>
void BulkMomenta<T,Lattice>::computeStress (
  Cell<T,Lattice> const& cell,
  T rho, const T u[Lattice<T>::d],
  T pi[util::TensorVal<Lattice<T> >::n] ) const
{
  lbHelpers<T,Lattice>::computeStress(cell, rho, u, pi);
}

template<typename T, template<typename U> class Lattice>
void BulkMomenta<T,Lattice>::computeRhoU (
  Cell<T,Lattice> const& cell,
  T& rho, T u[Lattice<T>::d] ) const
{
  lbHelpers<T,Lattice>::computeRhoU(cell, rho,u);
}

template<typename T, template<typename U> class Lattice>
void BulkMomenta<T,Lattice>::computeAllMomenta (
  Cell<T,Lattice> const& cell,
  T& rho, T u[Lattice<T>::d],
  T pi[util::TensorVal<Lattice<T> >::n] ) const
{
  this->computeRhoU(cell, rho, u);
  this->computeStress(cell, rho, u, pi);
}

template<typename T, template<typename U> class Lattice>
void BulkMomenta<T,Lattice>::defineRho(Cell<T,Lattice>& cell, T rho)
{
  T oldRho, u[Lattice<T>::d];
  computeRhoU(cell, oldRho, u);
  T uSqr = util::normSqr<T,Lattice<T>::d>(u);
  T fNeq[Lattice<T>::q];
  lbHelpers<T,Lattice>::computeFneq(cell, fNeq, oldRho, u);
  for (int iPop=0; iPop < Lattice<T>::q; ++iPop) {
    cell[iPop] = lbHelpers<T,Lattice>::equilibrium(iPop, rho, u, uSqr) +
                 fNeq[iPop];
  }
}

template<typename T, template<typename U> class Lattice>
void BulkMomenta<T,Lattice>::defineU (
  Cell<T,Lattice>& cell,
  const T u[Lattice<T>::d])
{
  T rho, oldU[Lattice<T>::d];
  computeRhoU(cell, rho, oldU);
  T uSqr = util::normSqr<T,Lattice<T>::d>(u);
  T fNeq[Lattice<T>::q];
  lbHelpers<T,Lattice>::computeFneq(cell, fNeq, rho, oldU);
  for (int iPop=0; iPop < Lattice<T>::q; ++iPop) {
    cell[iPop] = lbHelpers<T,Lattice>::equilibrium(iPop, rho, u, uSqr) +
                 fNeq[iPop];
  }

}

template<typename T, template<typename U> class Lattice>
void BulkMomenta<T,Lattice>::defineRhoU (
  Cell<T,Lattice>& cell,
  T rho, const T u[Lattice<T>::d])
{
  T oldRho, oldU[Lattice<T>::d];
  computeRhoU(cell, oldRho, oldU);
  T uSqr = util::normSqr<T,Lattice<T>::d>(u);
  T fNeq[Lattice<T>::q];
  lbHelpers<T,Lattice>::computeFneq(cell, fNeq, oldRho, oldU);
  for (int iPop=0; iPop < Lattice<T>::q; ++iPop) {
    cell[iPop] = lbHelpers<T,Lattice>::equilibrium(iPop, rho, u, uSqr) +
                 fNeq[iPop];
  }
}

template<typename T, template<typename U> class Lattice>
void BulkMomenta<T,Lattice>::defineAllMomenta (
  Cell<T,Lattice>& cell,
  T rho, const T u[Lattice<T>::d],
  const T pi[util::TensorVal<Lattice<T> >::n] )
{
  T uSqr = util::normSqr<T,Lattice<T>::d>(u);
  for (int iPop=0; iPop < Lattice<T>::q; ++iPop) {
    cell[iPop] = lbHelpers<T,Lattice>::equilibrium(iPop, rho, u, uSqr) +
                 firstOrderLbHelpers<T,Lattice>::fromPiToFneq(iPop, pi);
  }
}

////////////////////// Class ExternalVelocityMomenta //////////////////////////

template<typename T, template<typename U> class Lattice>
T ExternalVelocityMomenta<T,Lattice>::computeRho(Cell<T,Lattice> const& cell) const
{
  return lbHelpers<T,Lattice>::computeRho(cell);
}

template<typename T, template<typename U> class Lattice>
void ExternalVelocityMomenta<T,Lattice>::computeU(Cell<T,Lattice> const& cell, T u[Lattice<T>::d]) const
{
  T const* uExt = cell.getExternal(Lattice<T>::ExternalField::velocityBeginsAt);
  for (int iD=0; iD<Lattice<T>::d; ++iD) {
    u[iD] = uExt[iD];
  }
}

template<typename T, template<typename U> class Lattice>
void ExternalVelocityMomenta<T,Lattice>::computeJ(Cell<T,Lattice> const& cell, T j[Lattice<T>::d]) const
{
  T rho = computeRho(cell);
  T const* uExt = cell.getExternal(Lattice<T>::ExternalField::velocityBeginsAt);
  for (int iD=0; iD<Lattice<T>::d; ++iD) {
    j[iD] = uExt[iD]*rho;
  }
}

template<typename T, template<typename U> class Lattice>
void ExternalVelocityMomenta<T,Lattice>::computeStress (
  Cell<T,Lattice> const& cell,
  T rho, const T u[Lattice<T>::d],
  T pi[util::TensorVal<Lattice<T> >::n] ) const
{
  lbHelpers<T,Lattice>::computeStress(cell, rho, u, pi);
}

template<typename T, template<typename U> class Lattice>
void ExternalVelocityMomenta<T,Lattice>::computeRhoU (
  Cell<T,Lattice> const& cell,
  T& rho, T u[Lattice<T>::d] ) const
{
  rho = computeRho(cell);
  computeU(cell,u);
}

template<typename T, template<typename U> class Lattice>
void ExternalVelocityMomenta<T,Lattice>::computeAllMomenta (
  Cell<T,Lattice> const& cell,
  T& rho, T u[Lattice<T>::d],
  T pi[util::TensorVal<Lattice<T> >::n] ) const
{
  computeRhoU(cell, rho,u);
  computeStress(cell, rho, u, pi);
}

template<typename T, template<typename U> class Lattice>
void ExternalVelocityMomenta<T,Lattice>::defineRho(Cell<T,Lattice>& cell, T rho)
{
  T oldRho, u[Lattice<T>::d];
  computeRhoU(cell, oldRho, u);
  T uSqr = util::normSqr<T,Lattice<T>::d>(u);
  T fNeq[Lattice<T>::q];
  lbHelpers<T,Lattice>::computeFneq(cell, fNeq, oldRho, u);
  for (int iPop=0; iPop < Lattice<T>::q; ++iPop) {
    cell[iPop] = lbHelpers<T,Lattice>::equilibrium(iPop, rho, u, uSqr) +
                 fNeq[iPop];
  }
}

template<typename T, template<typename U> class Lattice>
void ExternalVelocityMomenta<T,Lattice>::defineU (
  Cell<T,Lattice>& cell,
  const T u[Lattice<T>::d])
{
  T* uExt = cell.getExternal(Lattice<T>::ExternalField::velocityBeginsAt);
  for (int iD=0; iD<Lattice<T>::d; ++iD) {
    uExt[iD] = u[iD];
  }
}

template<typename T, template<typename U> class Lattice>
void ExternalVelocityMomenta<T,Lattice>::defineRhoU (
  Cell<T,Lattice>& cell,
  T rho, const T u[Lattice<T>::d])
{
  defineRho(cell, rho);
  defineU(cell, u);
}

template<typename T, template<typename U> class Lattice>
void ExternalVelocityMomenta<T,Lattice>::defineAllMomenta (
  Cell<T,Lattice>& cell,
  T rho, const T u[Lattice<T>::d],
  const T pi[util::TensorVal<Lattice<T> >::n] )
{
  defineU(cell, u);
  T uSqr = util::normSqr<T,Lattice<T>::d>(u);
  for (int iPop=0; iPop < Lattice<T>::q; ++iPop) {
    cell[iPop] = lbHelpers<T,Lattice>::equilibrium(iPop, rho, u, uSqr) +
                 firstOrderLbHelpers<T,Lattice>::fromPiToFneq(iPop, pi);
  }
}

////////////////////// Class BounceBack ///////////////////////////

template<typename T, template<typename U> class Lattice>
BounceBack<T,Lattice>::BounceBack()
{
  _rhoFixed=false;
}

template<typename T, template<typename U> class Lattice>
BounceBack<T,Lattice>::BounceBack(T rho)
  :_rho(rho)
{
  _rhoFixed=true;
}

template<typename T, template<typename U> class Lattice>
void BounceBack<T,Lattice>::collide (
  Cell<T,Lattice>& cell,
  LatticeStatistics<T>& statistics )
{
  // !do not touch element 0!
  for (int iPop=1; iPop <= Lattice<T>::q/2; ++iPop) {
    std::swap(cell[iPop], cell[Lattice<T>::opposite[iPop]]);
  }
}

template<typename T, template<typename U> class Lattice>
void BounceBack<T,Lattice>::staticCollide (
  Cell<T,Lattice>& cell,
  const T u[Lattice<T>::d],
  LatticeStatistics<T>& statistics )
{
  this->collide(cell, statistics);
}

template<typename T, template<typename U> class Lattice>
T BounceBack<T,Lattice>::computeRho(Cell<T,Lattice> const& cell) const
{

  if (_rhoFixed) {
    return _rho;
  }
  return lbHelpers<T,Lattice>::computeRho(cell);
}

template<typename T, template<typename U> class Lattice>
void BounceBack<T,Lattice>::computeU (
  Cell<T,Lattice> const& cell,
  T u[Lattice<T>::d]) const
{
  for (int iD=0; iD<Lattice<T>::d; ++iD) {
    u[iD] = T();
  }
}

template<typename T, template<typename U> class Lattice>
void BounceBack<T,Lattice>::computeJ (
  Cell<T,Lattice> const& cell,
  T j[Lattice<T>::d]) const
{
  for (int iD=0; iD<Lattice<T>::d; ++iD) {
    j[iD] = T();
  }
}

template<typename T, template<typename U> class Lattice>
void BounceBack<T,Lattice>::computeStress (
  Cell<T,Lattice> const& cell,
  T rho, const T u[Lattice<T>::d],
  T pi[util::TensorVal<Lattice<T> >::n] ) const
{
  for (int iPi=0; iPi<util::TensorVal<Lattice<T> >::n; ++iPi) {
    pi[iPi] = T();//std::numeric_limits<T>::signaling_NaN();
  }
}

template<typename T, template<typename U> class Lattice>
void BounceBack<T,Lattice>::computeRhoU (
  Cell<T,Lattice> const& cell,
  T& rho, T u[Lattice<T>::d]) const
{
  rho = computeRho(cell);
  computeU(cell, u);
}

template<typename T, template<typename U> class Lattice>
void BounceBack<T,Lattice>::computeAllMomenta (
  Cell<T,Lattice> const& cell,
  T& rho, T u[Lattice<T>::d],
  T pi[util::TensorVal<Lattice<T> >::n] ) const
{
  computeRhoU(cell, rho, u);
  computeStress(cell, rho, u, pi);
}

template<typename T, template<typename U> class Lattice>
void BounceBack<T,Lattice>::defineRho(Cell<T,Lattice>& cell, T rho)
{ }

template<typename T, template<typename U> class Lattice>
void BounceBack<T,Lattice>::defineU (
  Cell<T,Lattice>& cell,
  const T u[Lattice<T>::d])
{ }

template<typename T, template<typename U> class Lattice>
void BounceBack<T,Lattice>::defineRhoU (
  Cell<T,Lattice>& cell,
  T rho, const T u[Lattice<T>::d])
{ }

template<typename T, template<typename U> class Lattice>
void BounceBack<T,Lattice>::defineAllMomenta (
  Cell<T,Lattice>& cell,
  T rho, const T u[Lattice<T>::d],
  const T pi[util::TensorVal<Lattice<T> >::n] )
{ }

template<typename T, template<typename U> class Lattice>
T BounceBack<T,Lattice>::getOmega() const
{
  return T();//std::numeric_limits<T>::signaling_NaN();
}

template<typename T, template<typename U> class Lattice>
void BounceBack<T,Lattice>::setOmega(T omega)
{ }


////////////////////// Class BounceBackVelocity ///////////////////////////

template<typename T, template<typename U> class Lattice>
BounceBackVelocity<T,Lattice>::BounceBackVelocity(const T u[Lattice<T>::d])
{
  _rhoFixed=false;
  for (int iD=0; iD<Lattice<T>::d; ++iD) {
    _u[iD] = u[iD];
  }
}

template<typename T, template<typename U> class Lattice>
BounceBackVelocity<T,Lattice>::BounceBackVelocity(const T rho, const T u[Lattice<T>::d])
  :_rho(rho)
{
  _rhoFixed=true;
  for (int iD=0; iD<Lattice<T>::d; ++iD) {
    _u[iD] = u[iD];
  }
}

template<typename T, template<typename U> class Lattice>
void BounceBackVelocity<T,Lattice>::collide (
  Cell<T,Lattice>& cell,
  LatticeStatistics<T>& statistics )
{
  for (int iPop=1; iPop <= Lattice<T>::q/2; ++iPop) {
    std::swap(cell[iPop], cell[iPop+Lattice<T>::q/2]);
  }
  for (int iPop=1; iPop < Lattice<T>::q; ++iPop) {
    for (int iD=0; iD<Lattice<T>::d; ++iD) {
      cell[iPop] += computeRho(cell)*_u[iD]*Lattice<T>::c[iPop][iD]*Lattice<T>::t[iPop]*2*Lattice<T>::invCs2;
    }
  }
}

template<typename T, template<typename U> class Lattice>
void BounceBackVelocity<T,Lattice>::staticCollide (
  Cell<T,Lattice>& cell,
  const T u[Lattice<T>::d],
  LatticeStatistics<T>& statistics )
{
  this->collide(cell, statistics);
}

template<typename T, template<typename U> class Lattice>
T BounceBackVelocity<T,Lattice>::computeRho(Cell<T,Lattice> const& cell) const
{
  if (_rhoFixed) {
    return _rho;
  }
  return lbHelpers<T,Lattice>::computeRho(cell);
}

template<typename T, template<typename U> class Lattice>
void BounceBackVelocity<T,Lattice>::computeU (
  Cell<T,Lattice> const& cell,
  T u[Lattice<T>::d]) const
{
  for (int iD=0; iD<Lattice<T>::d; ++iD) {
    u[iD] = _u[iD];
  }
}

template<typename T, template<typename U> class Lattice>
void BounceBackVelocity<T,Lattice>::computeJ (
  Cell<T,Lattice> const& cell,
  T j[Lattice<T>::d]) const
{
  for (int iD=0; iD<Lattice<T>::d; ++iD) {
    j[iD] = computeRho(cell)*_u[iD];
  }
}

template<typename T, template<typename U> class Lattice>
void BounceBackVelocity<T,Lattice>::computeStress (
  Cell<T,Lattice> const& cell,
  T rho, const T u[Lattice<T>::d],
  T pi[util::TensorVal<Lattice<T> >::n] ) const
{
  for (int iPi=0; iPi<util::TensorVal<Lattice<T> >::n; ++iPi) {
    pi[iPi] = T();//std::numeric_limits<T>::signaling_NaN();
  }
}

template<typename T, template<typename U> class Lattice>
void BounceBackVelocity<T,Lattice>::computeRhoU (
  Cell<T,Lattice> const& cell,
  T& rho, T u[Lattice<T>::d]) const
{
  rho = computeRho(cell);
  computeU(cell, u);
}

template<typename T, template<typename U> class Lattice>
void BounceBackVelocity<T,Lattice>::computeAllMomenta (
  Cell<T,Lattice> const& cell,
  T& rho, T u[Lattice<T>::d],
  T pi[util::TensorVal<Lattice<T> >::n] ) const
{
  computeRhoU(cell, rho, u);
  computeStress(cell, rho, u, pi);
}

template<typename T, template<typename U> class Lattice>
void BounceBackVelocity<T,Lattice>::defineRho(Cell<T,Lattice>& cell, T rho)
{ }

template<typename T, template<typename U> class Lattice>
void BounceBackVelocity<T,Lattice>::defineU (
  Cell<T,Lattice>& cell,
  const T u[Lattice<T>::d])
{
  for (int iD=0; iD<Lattice<T>::d; ++iD) {
    _u[iD] = u[iD];
  }
}

template<typename T, template<typename U> class Lattice>
void BounceBackVelocity<T,Lattice>::defineRhoU (
  Cell<T,Lattice>& cell,
  T rho, const T u[Lattice<T>::d])
{
  defineRho(cell,rho);
  defineU(cell,u);
}

template<typename T, template<typename U> class Lattice>
void BounceBackVelocity<T,Lattice>::defineAllMomenta (
  Cell<T,Lattice>& cell,
  T rho, const T u[Lattice<T>::d],
  const T pi[util::TensorVal<Lattice<T> >::n] )
{ }

template<typename T, template<typename U> class Lattice>
T BounceBackVelocity<T,Lattice>::getOmega() const
{
  return T();//std::numeric_limits<T>::signaling_NaN();
}

template<typename T, template<typename U> class Lattice>
void BounceBackVelocity<T,Lattice>::setOmega(T omega)
{ }

////////////////////// Class BounceBackAnti ///////////////////////////

template<typename T, template<typename U> class Lattice>
BounceBackAnti<T,Lattice>::BounceBackAnti()
{
  _rhoFixed = false;
  _rho = T(1);
}

template<typename T, template<typename U> class Lattice>
BounceBackAnti<T,Lattice>::BounceBackAnti(const T rho)
  :_rho(rho)
{
  _rhoFixed = true;
}

template<typename T, template<typename U> class Lattice>
void BounceBackAnti<T,Lattice>::collide (
  Cell<T,Lattice>& cell,
  LatticeStatistics<T>& statistics )
{
  /*
    for (int iPop=1; iPop <= Lattice<T>::q/2; ++iPop) {
      std::swap(cell[iPop], cell[iPop+Lattice<T>::q/2]);
    }
    for (int iPop=1; iPop < Lattice<T>::q; ++iPop) {
      if (Lattice<T>::c[iPop][0] == -1)
        cell[iPop] = -cell[Lattice<T>::opposite[iPop]] + (computeRho(cell) - T(1))*Lattice<T>::t[iPop]*2;
    }
  */
  //T rho, u[Lattice<T>::d];
  //computeRhoU(cell, rho, u);
  //T uSqr = lbHelpers<T,Lattice>::bgkCollision(cell, rho, u, 1.78571);
  //statistics.incrementStats(rho, uSqr);

}

template<typename T, template<typename U> class Lattice>
void BounceBackAnti<T,Lattice>::staticCollide (
  Cell<T,Lattice>& cell,
  const T u[Lattice<T>::d],
  LatticeStatistics<T>& statistics )
{
  this->collide(cell, statistics);
}

template<typename T, template<typename U> class Lattice>
T BounceBackAnti<T,Lattice>::computeRho(Cell<T,Lattice> const& cell) const
{

  if (_rhoFixed) {
    return _rho;
  }
  return lbHelpers<T,Lattice>::computeRho(cell);
}

template<typename T, template<typename U> class Lattice>
void BounceBackAnti<T,Lattice>::computeU (
  Cell<T,Lattice> const& cell,
  T u[Lattice<T>::d]) const
{
  for (int iD=0; iD<Lattice<T>::d; ++iD) {
    u[iD] = T();
  }
}

template<typename T, template<typename U> class Lattice>
void BounceBackAnti<T,Lattice>::computeJ (
  Cell<T,Lattice> const& cell,
  T j[Lattice<T>::d]) const
{
  computeU(cell, j);
  for (int iD=0; iD<Lattice<T>::d; ++iD) {
    j[iD]*=computeRho(cell);
  }
}

template<typename T, template<typename U> class Lattice>
void BounceBackAnti<T,Lattice>::computeStress (
  Cell<T,Lattice> const& cell,
  T rho, const T u[Lattice<T>::d],
  T pi[util::TensorVal<Lattice<T> >::n] ) const
{
  for (int iPi=0; iPi<util::TensorVal<Lattice<T> >::n; ++iPi) {
    pi[iPi] = T();//std::numeric_limits<T>::signaling_NaN();
  }
}

template<typename T, template<typename U> class Lattice>
void BounceBackAnti<T,Lattice>::computeRhoU (
  Cell<T,Lattice> const& cell,
  T& rho, T u[Lattice<T>::d]) const
{
  rho = computeRho(cell);
  computeU(cell, u);
}

template<typename T, template<typename U> class Lattice>
void BounceBackAnti<T,Lattice>::computeAllMomenta (
  Cell<T,Lattice> const& cell,
  T& rho, T u[Lattice<T>::d],
  T pi[util::TensorVal<Lattice<T> >::n] ) const
{
  computeRhoU(cell, rho, u);
  computeStress(cell, rho, u, pi);
}

template<typename T, template<typename U> class Lattice>
void BounceBackAnti<T,Lattice>::defineRho(Cell<T,Lattice>& cell, T rho)
{
  _rho = rho;
}

template<typename T, template<typename U> class Lattice>
void BounceBackAnti<T,Lattice>::defineU (
  Cell<T,Lattice>& cell,
  const T u[Lattice<T>::d])
{
}

template<typename T, template<typename U> class Lattice>
void BounceBackAnti<T,Lattice>::defineRhoU (
  Cell<T,Lattice>& cell,
  T rho, const T u[Lattice<T>::d])
{
  defineRho(cell,rho);
  defineU(cell,u);
}

template<typename T, template<typename U> class Lattice>
void BounceBackAnti<T,Lattice>::defineAllMomenta (
  Cell<T,Lattice>& cell,
  T rho, const T u[Lattice<T>::d],
  const T pi[util::TensorVal<Lattice<T> >::n] )
{ }

template<typename T, template<typename U> class Lattice>
T BounceBackAnti<T,Lattice>::getOmega() const
{
  return T();//std::numeric_limits<T>::signaling_NaN();
}

template<typename T, template<typename U> class Lattice>
void BounceBackAnti<T,Lattice>::setOmega(T omega)
{ }


////////////////////// Class partialBounceBack ///////////////////////////

template<typename T, template<typename U> class Lattice>
PartialBounceBack<T,Lattice>::PartialBounceBack(T rf) : _rf(rf)
{
}

template<typename T, template<typename U> class Lattice>
T PartialBounceBack<T, Lattice>::computeEquilibrium ( int iPop, T rho, const T u[Lattice<T>::d], T uSqr ) const
{
  return lbHelpers<T, Lattice>::equilibriumFirstOrder( iPop, rho, u );
}

template<typename T, template<typename U> class Lattice>
void PartialBounceBack<T,Lattice>::collide( Cell<T,Lattice>& cell, LatticeStatistics<T>& statistics )
{
  for (int iPop=1; iPop <= Lattice<T>::q/2; ++iPop) {
    std::swap(cell[iPop], cell[iPop+Lattice<T>::q/2]);
  }
  for (int iPop=1; iPop < Lattice<T>::q; ++iPop) {
    cell[iPop] = (_rf -1) * (cell[iPop] + Lattice<T>::t[iPop]) - Lattice<T>::t[iPop];
  }
}



////////////////////// Class NoDynamics ///////////////////////////

template<typename T, template<typename U> class Lattice>
NoDynamics<T,Lattice>::NoDynamics(T rho) :_rho(rho)
{
}

template<typename T, template<typename U> class Lattice>
T NoDynamics<T,Lattice>::computeEquilibrium(int iPop, T rho, const T u[Lattice<T>::d], T uSqr) const
{
  return T();
}

template<typename T, template<typename U> class Lattice>
void NoDynamics<T,Lattice>::collide (
  Cell<T,Lattice>& cell,
  LatticeStatistics<T>& statistics )
{ }

template<typename T, template<typename U> class Lattice>
void NoDynamics<T,Lattice>::staticCollide (
  Cell<T,Lattice>& cell,
  const T u[Lattice<T>::d],
  LatticeStatistics<T>& statistics )
{ }

template<typename T, template<typename U> class Lattice>
T NoDynamics<T,Lattice>::computeRho(Cell<T,Lattice> const& cell) const
{
  return _rho;
}

template<typename T, template<typename U> class Lattice>
void NoDynamics<T,Lattice>::computeU (
  Cell<T,Lattice> const& cell,
  T u[Lattice<T>::d]) const
{
  for (int iD=0; iD<Lattice<T>::d; ++iD) {
    u[iD] = T();
  }
}

template<typename T, template<typename U> class Lattice>
void NoDynamics<T,Lattice>::computeJ (
  Cell<T,Lattice> const& cell,
  T j[Lattice<T>::d]) const
{
  for (int iD=0; iD<Lattice<T>::d; ++iD) {
    j[iD] = T();
  }
}

template<typename T, template<typename U> class Lattice>
void NoDynamics<T,Lattice>::computeStress (
  Cell<T,Lattice> const& cell,
  T rho, const T u[Lattice<T>::d],
  T pi[util::TensorVal<Lattice<T> >::n] ) const
{
  for (int iPi=0; iPi<util::TensorVal<Lattice<T> >::n; ++iPi) {
    pi[iPi] = T();//std::numeric_limits<T>::signaling_NaN();
  }
}

template<typename T, template<typename U> class Lattice>
void NoDynamics<T,Lattice>::computeRhoU (
  Cell<T,Lattice> const& cell,
  T& rho, T u[Lattice<T>::d]) const
{
  rho = computeRho(cell);
  computeU(cell, u);
}

template<typename T, template<typename U> class Lattice>
void NoDynamics<T,Lattice>::computeAllMomenta (
  Cell<T,Lattice> const& cell,
  T& rho, T u[Lattice<T>::d],
  T pi[util::TensorVal<Lattice<T> >::n] ) const
{
  computeRhoU(cell, rho, u);
  computeStress(cell, rho, u, pi);
}

template<typename T, template<typename U> class Lattice>
void NoDynamics<T,Lattice>::defineRho(Cell<T,Lattice>& cell, T rho)
{ }

template<typename T, template<typename U> class Lattice>
void NoDynamics<T,Lattice>::defineU (
  Cell<T,Lattice>& cell,
  const T u[Lattice<T>::d])
{ }

template<typename T, template<typename U> class Lattice>
void NoDynamics<T,Lattice>::defineRhoU (
  Cell<T,Lattice>& cell,
  T rho, const T u[Lattice<T>::d])
{ }

template<typename T, template<typename U> class Lattice>
void NoDynamics<T,Lattice>::defineAllMomenta (
  Cell<T,Lattice>& cell,
  T rho, const T u[Lattice<T>::d],
  const T pi[util::TensorVal<Lattice<T> >::n] )
{ }

template<typename T, template<typename U> class Lattice>
T NoDynamics<T,Lattice>::getOmega() const
{
  return T();//std::numeric_limits<T>::signaling_NaN();
}

template<typename T, template<typename U> class Lattice>
void NoDynamics<T,Lattice>::setOmega(T omega)
{ }

////////////////////// Class offDynamics ///////////////////////////

template<typename T, template<typename U> class Lattice>
OffDynamics<T,Lattice>::OffDynamics(const T _location[Lattice<T>::d])
{
  typedef Lattice<T> L;
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

template<typename T, template<typename U> class Lattice>
OffDynamics<T,Lattice>::OffDynamics(const T _location[Lattice<T>::d], T _distances[Lattice<T>::q])
{
  typedef Lattice<T> L;
  for (int iD = 0; iD < L::d; iD++) {
    location[iD] = _location[iD];
  }
  for (int iPop = 0; iPop < L::q; iPop++) {
    distances[iPop] = _distances[iPop];
    velocityCoefficient[iPop] = 0;
    const int* c = L::c[iPop];
    for (int iD = 0; iD < L::d; iD++) {
      boundaryIntersection[iPop][iD] = _location[iD] - _distances[iPop]*c[iD];
      _u[iPop][iD] = T();
    }
  }
  _rho=T(1);
}

template<typename T, template<typename U> class Lattice>
T OffDynamics<T,Lattice>::computeRho(Cell<T,Lattice> const& cell) const
{
  /*typedef Lattice<T> L;
  T rhoTmp = T();
  T counter = T();
  int counter2 = int();
  for (int iPop = 0; iPop < L::q; iPop++) {
    if (distances[iPop] != -1) {
      rhoTmp += (cell[iPop] + L::t[iPop])*L::t[iPop];
      counter += L::t[iPop];
      counter2++;
    }
  }
  //if (rhoTmp/counter + 1<0.1999) std::cout << rhoTmp/counter2 + 1 <<std::endl;
  //if (rhoTmp/counter + 1>1.001) std::cout << rhoTmp/counter2 + 1 <<std::endl;
  return rhoTmp/counter/counter;*/
  return _rho;
}

template<typename T, template<typename U> class Lattice>
void OffDynamics<T,Lattice>::computeU(Cell<T,Lattice> const& cell, T u[Lattice<T>::d] ) const
{
  typedef Lattice<T> L;
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

template<typename T, template<typename U> class Lattice>
void OffDynamics<T,Lattice>::setBoundaryIntersection(int iPop, T distance)
{
  /// direction points from the fluid node into the solid domain
  /// distance is the distance from the fluid node to the solid wall
  typedef Lattice<T> L;
  distances[iPop] = distance;
  const int* c = L::c[iPop];
  for (int iD = 0; iD < L::d; iD++) {
    boundaryIntersection[iPop][iD] = location[iD] - distance*c[iD];
  }
}

template<typename T, template<typename U> class Lattice>
bool OffDynamics<T,Lattice>::getBoundaryIntersection(int iPop, T intersection[Lattice<T>::d])
{
  typedef Lattice<T> L;
  if ( !util::nearZero(distances[iPop]+1) ) {
    for (int iD = 0; iD < L::d; iD++) {
      intersection[iD] = boundaryIntersection[iPop][iD];
    }
    return true;
  }
  return false;
}

template<typename T, template<typename U> class Lattice>
void OffDynamics<T,Lattice>::defineRho(Cell<T,Lattice>& cell, T rho)
{
  _rho=rho;
}

template<typename T, template<typename U> class Lattice>
void OffDynamics<T,Lattice>::defineRho(int iPop, T rho)
{
  _rho=rho;
}

template<typename T, template<typename U> class Lattice>
void OffDynamics<T,Lattice>::defineU (
  Cell<T,Lattice>& cell,
  const T u[Lattice<T>::d])
{
  defineU(u);
}

template<typename T, template<typename U> class Lattice>
void OffDynamics<T,Lattice>::defineU(const T u[Lattice<T>::d])
{
  typedef Lattice<T> L;
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

template<typename T, template<typename U> class Lattice>
void OffDynamics<T,Lattice>::defineU(
  int iPop, const T u[Lattice<T>::d])
{
  OLB_PRECONDITION(distances[iPop] != -1)
  typedef Lattice<T> L;
  const int* c = L::c[iPop];
  velocityCoefficient[iPop] = 0;
  // scalar product of c(iPop) and u
  for (int sum = 0; sum < L::d; sum++) { // +/- problem because of first stream than postprocess
    velocityCoefficient[iPop] -= c[sum]*u[sum];
  }
  // compute summand for boundary condition
  velocityCoefficient[iPop] *= 2*L::invCs2 * L::t[iPop];

  for (int iD = 0; iD < L::d; iD++) {
    _u[iPop][iD] = u[iD];
  }
}

template<typename T, template<typename U> class Lattice>
T OffDynamics<T,Lattice>::getVelocityCoefficient(int iPop)
{
  return velocityCoefficient[iPop];
}

////////////////////// Class ZeroDistributionDynamics ///////////////////////////

template<typename T, template<typename U> class Lattice>
ZeroDistributionDynamics<T,Lattice>::ZeroDistributionDynamics()
{
}

template<typename T, template<typename U> class Lattice>
void ZeroDistributionDynamics<T,Lattice>::collide (
  Cell<T,Lattice>& cell,
  LatticeStatistics<T>& statistics )
{
  for (int iPop=0; iPop < Lattice<T>::q; ++iPop) {
    cell[iPop] = -Lattice<T>::t[iPop];
  }
}

template<typename T, template<typename U> class Lattice>
T ZeroDistributionDynamics<T,Lattice>::computeRho(Cell<T,Lattice> const& cell) const
{
  return lbHelpers<T,Lattice>::computeRho(cell);
}

/////////////// Singletons //////////////////////////////////

namespace instances {

template<typename T, template<typename U> class Lattice>
BulkMomenta<T,Lattice>& getBulkMomenta()
{
  static BulkMomenta<T,Lattice> bulkMomentaSingleton;
  return bulkMomentaSingleton;
}

template<typename T, template<typename U> class Lattice>
ExternalVelocityMomenta<T,Lattice>& getExternalVelocityMomenta()
{
  static ExternalVelocityMomenta<T,Lattice> externalVelocityMomentaSingleton;
  return externalVelocityMomentaSingleton;
}

template<typename T, template<typename U> class Lattice>
BounceBack<T,Lattice>& getBounceBack()
{
  static BounceBack<T,Lattice> bounceBackSingleton;
  return bounceBackSingleton;
}

template<typename T, template<typename U> class Lattice>
PartialBounceBack<T,Lattice>& getPartialBounceBack(const double rf)
{
  static PartialBounceBack<T,Lattice> partialBounceBackSingleton(rf);
  return partialBounceBackSingleton;
}


template<typename T, template<typename U> class Lattice>
BounceBackVelocity<T,Lattice>& getBounceBackVelocity(const double rho, const double u[Lattice<T>::d])
{
  static BounceBackVelocity<T,Lattice> bounceBackSingleton(rho,u);
  return bounceBackSingleton;
}

template<typename T, template<typename U> class Lattice>
BounceBackAnti<T,Lattice>& getBounceBackAnti(const double rho)
{
  static BounceBackAnti<T,Lattice> bounceBackSingleton(rho);
  return bounceBackSingleton;
}

template<typename T, template<typename U> class Lattice>
NoDynamics<T,Lattice>& getNoDynamics(T rho)
{
  static NoDynamics<T,Lattice> noDynamicsSingleton(rho);
  return noDynamicsSingleton;
}

template<typename T, template<typename U> class Lattice>
ZeroDistributionDynamics<T,Lattice>& getZeroDistributionDynamics()
{
  static ZeroDistributionDynamics<T,Lattice> zeroSingleton;
  return zeroSingleton;
}

}

}

#endif
