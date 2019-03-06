/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2016 Thomas Henn, Mathias J. Krause, Jonas Latt
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
#ifndef POROUS_BGK_DYNAMICS_HH
#define POROUS_BGK_DYNAMICS_HH

#include "porousBGKdynamics.h"
#include "core/cell.h"
#include "dynamics.h"
#include "core/util.h"
#include "lbHelpers.h"
#include "math.h"

namespace olb {

////////////////////// Class PorousBGKdynamics //////////////////////////

template<typename T, template<typename U> class Lattice>
PorousBGKdynamics<T,Lattice>::PorousBGKdynamics (
  T omega_, Momenta<T,Lattice>& momenta_)
  : BGKdynamics<T,Lattice>(omega_,momenta_),
    omega(omega_)
{ }

template<typename T, template<typename U> class Lattice>
void PorousBGKdynamics<T,Lattice>::collide (
  Cell<T,Lattice>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho, u[Lattice<T>::d];
  this->_momenta.computeRhoU(cell, rho, u);
  T* porosity = cell.getExternal(porosityIsAt);
  for (int i=0; i<Lattice<T>::d; i++)  {
    u[i] *= porosity[0];
  }
  T uSqr = lbHelpers<T,Lattice>::bgkCollision(cell, rho, u, omega);
  statistics.incrementStats(rho, uSqr);
}

template<typename T, template<typename U> class Lattice>
T PorousBGKdynamics<T,Lattice>::getOmega() const
{
  return omega;
}

template<typename T, template<typename U> class Lattice>
void PorousBGKdynamics<T,Lattice>::setOmega(T omega_)
{
  omega = omega_;
}


//////////////////// Class ExtendedPorousBGKdynamics ////////////////////

template<typename T, template<typename U> class Lattice>
ExtendedPorousBGKdynamics<T,Lattice>::ExtendedPorousBGKdynamics (
  T omega_, Momenta<T,Lattice>& momenta_)
  : BGKdynamics<T,Lattice>(omega_,momenta_),
    omega(omega_)
{
}

template<typename T, template<typename U> class Lattice>
void ExtendedPorousBGKdynamics<T,Lattice>::collide (
  Cell<T,Lattice>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho, u[Lattice<T>::d];
  this->_momenta.computeRhoU(cell, rho, u);
  T* porosity = cell.getExternal(porosityIsAt);
  T* localVelocity = cell.getExternal(localDragBeginsAt);

  cell.defineExternalField(localDragBeginsAt,Lattice<T>::d, u);

  for (int i=0; i<Lattice<T>::d; i++)  {
    u[i] *= porosity[0];
    u[i] += (1.-porosity[0]) * localVelocity[i];
  }
  T uSqr = lbHelpers<T,Lattice>::bgkCollision(cell, rho, u, omega);
  statistics.incrementStats(rho, uSqr);
}

template<typename T, template<typename U> class Lattice>
T ExtendedPorousBGKdynamics<T,Lattice>::getOmega() const
{
  return omega;
}

template<typename T, template<typename U> class Lattice>
void ExtendedPorousBGKdynamics<T,Lattice>::setOmega(T omega_)
{
  omega = omega_;
}

//////////////////// Class SubgridParticleBGKdynamics ////////////////////

template<typename T, template<typename U> class Lattice>
SubgridParticleBGKdynamics<T,Lattice>::SubgridParticleBGKdynamics (
  T omega_, Momenta<T,Lattice>& momenta_)
  : BGKdynamics<T,Lattice>(omega_,momenta_),
    omega(omega_)
{
  _fieldTmp[0] = T();
  _fieldTmp[1] = T();
  _fieldTmp[2] = T();
  _fieldTmp[3] = T();
}

template<typename T, template<typename U> class Lattice>
void SubgridParticleBGKdynamics<T,Lattice>::collide (
  Cell<T,Lattice>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho, u[Lattice<T>::d];
  this->_momenta.computeRhoU(cell, rho, u);
  T* porosity = cell.getExternal(porosityIsAt);
  T* extVelocity = cell.getExternal(localDragBeginsAt);
//  if (porosity[0] != 0) {
//    cout << "extVelocity: " << extVelocity[0] << " " <<  extVelocity[1] << " " <<  extVelocity[2] << " " << std::endl;
//    cout << "porosity: " << porosity[0] << std::endl;
//  }
  for (int i=0; i<Lattice<T>::d; i++)  {
    u[i] *= (1.-porosity[0]);
    u[i] += extVelocity[i];
  }
  T uSqr = lbHelpers<T,Lattice>::bgkCollision(cell, rho, u, omega);

  statistics.incrementStats(rho, uSqr);
  for (int i=0; i < 4; ++i) {
    cell.getExternal(0)[i] = 0; //_fieldTmp[i];
  }
}

template<typename T, template<typename U> class Lattice>
T SubgridParticleBGKdynamics<T,Lattice>::getOmega() const
{
  return omega;
}

template<typename T, template<typename U> class Lattice>
void SubgridParticleBGKdynamics<T,Lattice>::setOmega(T omega_)
{
  omega = omega_;
}

//////////////////// Class PorousParticleBGKdynamics ////////////////////

template<typename T, template<typename U> class Lattice>
PorousParticleBGKdynamics<T,Lattice>::PorousParticleBGKdynamics (
  T omega_, Momenta<T,Lattice>& momenta_)
  : BGKdynamics<T,Lattice>(omega_,momenta_),
    omega(omega_)
{}

template<typename T, template<typename U> class Lattice>
void PorousParticleBGKdynamics<T,Lattice>::collide (
  Cell<T,Lattice>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho, u[Lattice<T>::d];
  this->_momenta.computeRhoU(cell, rho, u);
  T* external = cell.getExternal(0);
  if (external[velDenominator] > std::numeric_limits<T>::epsilon()) {
    external[porosityIsAt] = 1.-external[porosityIsAt]; // 1-prod(1-smoothInd)
    for (int i=0; i<Lattice<T>::d; i++)  {
      u[i] += external[porosityIsAt] * (external[velNumerator+i] / external[velDenominator]- u[i]);
    }
  }
  T uSqr = lbHelpers<T,Lattice>::bgkCollision(cell, rho, u, omega);
  statistics.incrementStats(rho, uSqr);

  external[0] = 1.;
  for (int i=1; i < Lattice<T>::d+2; ++i) {
    external[i] = 0.;
  }

}

template<typename T, template<typename U> class Lattice>
T PorousParticleBGKdynamics<T,Lattice>::getOmega() const
{
  return omega;
}

template<typename T, template<typename U> class Lattice>
void PorousParticleBGKdynamics<T,Lattice>::setOmega(T omega_)
{
  omega = omega_;
}

//////////////////// Class KrauseHBGKdynamics ////////////////////

template<typename T, template<typename U> class Lattice>
KrauseHBGKdynamics<T,Lattice>::KrauseHBGKdynamics (T omega_,
    Momenta<T,Lattice>& momenta_, T smagoConst_, T dx_, T dt_ )
  : BGKdynamics<T,Lattice>(omega_,momenta_), omega(omega_), smagoConst(smagoConst_),
    preFactor(computePreFactor(omega_,smagoConst_) )
{
  _fieldTmp[0] = T(1);
  _fieldTmp[1] = T();
  _fieldTmp[2] = T();
  _fieldTmp[3] = T();
}

template<typename T, template<typename U> class Lattice>
void KrauseHBGKdynamics<T,Lattice>::collide (
  Cell<T,Lattice>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho, u[Lattice<T>::d];
  T newOmega[Lattice<T>::q];
  this->_momenta.computeRhoU(cell, rho, u);
  computeOmega(this->getOmega(), cell, preFactor, rho, u, newOmega);

  T vel_denom = *cell.getExternal(velDenominator);
  if (vel_denom > std::numeric_limits<T>::epsilon()) {
    T porosity = *cell.getExternal(porosityIsAt); // prod(1-smoothInd)
    T* vel_num = cell.getExternal(velNumerator);
    porosity = 1.-porosity; // 1-prod(1-smoothInd)
    for (int i=0; i<Lattice<T>::d; i++)  {
      u[i] += porosity * (vel_num[i] / vel_denom - u[i]);
    }
  }
  T uSqr = lbHelpers<T,Lattice>::bgkCollision(cell, rho, u, newOmega);
  statistics.incrementStats(rho, uSqr);

  for (int i=0; i < 4; ++i) {
    cell.getExternal(0)[i] = _fieldTmp[i];
  }
}

template<typename T, template<typename U> class Lattice>
T KrauseHBGKdynamics<T,Lattice>::getOmega() const
{
  return omega;
}

template<typename T, template<typename U> class Lattice>
void KrauseHBGKdynamics<T,Lattice>::setOmega(T omega)
{
  this->setOmega(omega);
  preFactor = computePreFactor(omega, smagoConst);
}

template<typename T, template<typename U> class Lattice>
T KrauseHBGKdynamics<T,Lattice>::computePreFactor(T omega, T smagoConst)
{
  return (T)smagoConst*smagoConst*Lattice<T>::invCs2*Lattice<T>::invCs2*2*sqrt(2);
}


template<typename T, template<typename U> class Lattice>
void KrauseHBGKdynamics<T,Lattice>::computeOmega(T omega0, Cell<T,Lattice>& cell, T preFactor, T rho,
    T u[Lattice<T>::d], T newOmega[Lattice<T>::q])
{
  T uSqr = u[0]*u[0];
  for (int iDim=0; iDim<Lattice<T>::d; iDim++) {
    uSqr += u[iDim]*u[iDim];
  }
  /// Molecular realaxation time
  T tau_mol = 1./omega0;

  for (int iPop=0; iPop<Lattice<T>::q; iPop++) {
    T fNeq = std::fabs(cell[iPop] - lbHelpers<T,Lattice>::equilibrium(iPop, rho, u, uSqr));
    /// Turbulent realaxation time
    T tau_turb = 0.5*(sqrt(tau_mol*tau_mol+(preFactor*fNeq))-tau_mol);
    /// Effective realaxation time
    tau_eff = tau_mol + tau_turb;
    newOmega[iPop] = 1./tau_eff;
  }
}


//////////////////// Class ParticlePorousBGKdynamics ////////////////////
/*
template<typename T, template<typename U> class Lattice>
ParticlePorousBGKdynamics<T,Lattice>::ParticlePorousBGKdynamics (
  T omega_, Momenta<T,Lattice>& momenta_)
  : BGKdynamics<T,Lattice>(omega_,momenta_),
    omega(omega_)
{ }

template<typename T, template<typename U> class Lattice>
void ParticlePorousBGKdynamics<T,Lattice>::collide (
  Cell<T,Lattice>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho, u[Lattice<T>::d];
  this->_momenta.computeRhoU(cell, rho, u);
  T* porosity = cell.getExternal(porosityIsAt);
  T* localVelocity = cell.getExternal(localDragBeginsAt);
  for (int i=0; i<Lattice<T>::d; i++)  {
    u[i] *= porosity[0];
    u[i] += localVelocity[i];
  }
  T uSqr = lbHelpers<T,Lattice>::bgkCollision(cell, rho, u, omega);
  statistics.incrementStats(rho, uSqr);

//  *cell.getExternal(porosityIsAt) = 1;
//  *cell.getExternal(localDragBeginsAt) = 0.;
//  *(cell.getExternal(localDragBeginsAt)+1) = 0.;
}

template<typename T, template<typename U> class Lattice>
T ParticlePorousBGKdynamics<T,Lattice>::getOmega() const {
  return omega;
}

template<typename T, template<typename U> class Lattice>
void ParticlePorousBGKdynamics<T,Lattice>::setOmega(T omega_) {
  omega = omega_;
}
*/

//////////////////// Class SmallParticleBGKdynamics ////////////////////

template<typename T, template<typename U> class Lattice>
SmallParticleBGKdynamics<T,Lattice>::SmallParticleBGKdynamics (
  T omega_, Momenta<T,Lattice>& momenta_)
  : BGKdynamics<T,Lattice>(omega_,momenta_),
    omega(omega_)
{ }

template<typename T, template<typename U> class Lattice>
void SmallParticleBGKdynamics<T,Lattice>::collide (
  Cell<T,Lattice>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho, u[Lattice<T>::d];
  this->_momenta.computeRhoU(cell, rho, u);
  T* porosity = cell.getExternal(porosityIsAt);
  T* localVelocity = cell.getExternal(localDragBeginsAt);

  //cout << porosity[0]  << " " <<   localVelocity[0]<< " " <<   localVelocity[1]<< " " <<   localVelocity[2]<<std::endl;
  for (int i=0; i<Lattice<T>::d; i++)  {
    u[i] *= porosity[0];
    u[i] += localVelocity[i];
  }
  T uSqr = lbHelpers<T,Lattice>::bgkCollision(cell, rho, u, omega);
  statistics.incrementStats(rho, uSqr);
}

template<typename T, template<typename U> class Lattice>
T SmallParticleBGKdynamics<T,Lattice>::getOmega() const
{
  return omega;
}

template<typename T, template<typename U> class Lattice>
void SmallParticleBGKdynamics<T,Lattice>::setOmega(T omega_)
{
  omega = omega_;
}

} // olb

#endif
