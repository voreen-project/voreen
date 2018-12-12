/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2008 Orestis Malaspinas, Andrea Parmigiani
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
#ifndef ADVECTION_DIFFUSION_DYNAMICS_HH
#define ADVECTION_DIFFUSION_DYNAMICS_HH

#include <algorithm>
#include <limits>
#include "advectionDiffusionDynamics.h"

namespace olb {


////////////////////// Class AdvectionDiffusionRLBdynamics //////////////////////////

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 */
//==================================================================//
//============= Regularized Model for Advection diffusion===========//
//==================================================================//

template<typename T, template<typename U> class Lattice>
AdvectionDiffusionRLBdynamics<T, Lattice>::AdvectionDiffusionRLBdynamics (
  T omega_, Momenta<T, Lattice>& momenta_ )
  : BasicDynamics<T, Lattice>( momenta_ ),
    omega( omega_ )
{ }

template<typename T, template<typename U> class Lattice>
T AdvectionDiffusionRLBdynamics<T, Lattice>::computeEquilibrium( int iPop, T rho,
    const T u[Lattice<T>::d], T uSqr ) const
{
  return lbHelpers<T, Lattice>::equilibriumFirstOrder( iPop, rho, u );
}


template<typename T, template<typename U> class Lattice>
void AdvectionDiffusionRLBdynamics<T, Lattice>::collide( Cell<T, Lattice>& cell,
    LatticeStatistics<T>& statistics )
{
  T temperature = this->_momenta.computeRho( cell );

  const T* u = cell.getExternal( Lattice<T>::ExternalField::velocityBeginsAt );

  T uSqr = lbHelpers<T, Lattice>::
           rlbCollision( cell, temperature, u, omega );

  statistics.incrementStats( temperature, uSqr );
}


template<typename T, template<typename U> class Lattice>
void AdvectionDiffusionRLBdynamics<T, Lattice>::staticCollide( Cell<T, Lattice>& cell,
    const T u[Lattice<T>::d] ,
    LatticeStatistics<T>& statistics )
{
  assert( false );
}

template<typename T, template<typename U> class Lattice>
T AdvectionDiffusionRLBdynamics<T, Lattice>::getOmega() const
{
  return omega;
}

template<typename T, template<typename U> class Lattice>
void AdvectionDiffusionRLBdynamics<T, Lattice>::setOmega( T omega_ )
{
  omega = omega_;
}

//==================================================================//
//============= BGK Model for Advection diffusion===========//
//==================================================================//

template<typename T, template<typename U> class Lattice>
AdvectionDiffusionBGKdynamics<T, Lattice>::AdvectionDiffusionBGKdynamics (
  T omega, Momenta<T, Lattice>& momenta )
  : BasicDynamics<T, Lattice>( momenta ),
    _omega(omega)
{ }

template<typename T, template<typename U> class Lattice>
AdvectionDiffusionBGKdynamics<T, Lattice>::AdvectionDiffusionBGKdynamics (
  const UnitConverter<T,Lattice>& converter, Momenta<T, Lattice>& momenta )
  : BasicDynamics<T, Lattice>( momenta ),
    _omega(converter.getLatticeRelaxationFrequency())
{ }

template<typename T, template<typename U> class Lattice>
T AdvectionDiffusionBGKdynamics<T, Lattice>::computeEquilibrium( int iPop, T rho,
    const T u[Lattice<T>::d], T uSqr ) const
{
  return lbHelpers<T, Lattice>::equilibriumFirstOrder( iPop, rho, u );
}


template<typename T, template<typename U> class Lattice>
void AdvectionDiffusionBGKdynamics<T, Lattice>::collide( Cell<T, Lattice>& cell,
    LatticeStatistics<T>& statistics )
{
  T temperature = this->_momenta.computeRho( cell );
  const T* u = cell.getExternal(Lattice<T>::ExternalField::velocityBeginsAt);

  T uSqr = lbHelpers<T, Lattice>::
           bgkCollision( cell, temperature, u, _omega );

  statistics.incrementStats( temperature, uSqr );
}


template<typename T, template<typename U> class Lattice>
void AdvectionDiffusionBGKdynamics<T, Lattice>::staticCollide( Cell<T, Lattice>& cell,
    const T u[Lattice<T>::d] ,
    LatticeStatistics<T>& statistics )
{
  assert( false );
}

template<typename T, template<typename U> class Lattice>
T AdvectionDiffusionBGKdynamics<T, Lattice>::getOmega() const
{
  return _omega;
}

template<typename T, template<typename U> class Lattice>
void AdvectionDiffusionBGKdynamics<T, Lattice>::setOmega( T omega )
{
  _omega = omega;
}


//==================================================================================//
//=========== BGK Model for Advection diffusion with Stokes drag and Smagorinsky====//
//==================================================================================//

template<typename T, template<typename U> class Lattice>
SmagorinskyParticleAdvectionDiffusionBGKdynamics<T, Lattice>::SmagorinskyParticleAdvectionDiffusionBGKdynamics (
  T omega_, Momenta<T,Lattice>& momenta_, T smagoConst_, T dx_, T dt_)
  : AdvectionDiffusionBGKdynamics<T,Lattice>(omega_,momenta_), smagoConst(smagoConst_), preFactor(computePreFactor(omega_,smagoConst_, dx_, dt_) )
{ }

template<typename T, template<typename U> class Lattice>
void SmagorinskyParticleAdvectionDiffusionBGKdynamics<T, Lattice>::collide(Cell<T,Lattice>& cell, LatticeStatistics<T>& statistics )
{
  T temperature, uad[Lattice<T>::d], pi[util::TensorVal<Lattice<T> >::n];
  this->_momenta.computeAllMomenta(cell, temperature, uad, pi);
  int offset = (statistics.getTime() % 2 == 0) ? Lattice<T>::ExternalField::velocityBeginsAt : Lattice<T>::ExternalField::velocity2BeginsAt;
  const T* u = cell.getExternal(offset);
  T newOmega = computeOmega(this->getOmega(), preFactor, temperature, pi);
  T uSqr = lbHelpers<T,Lattice>::bgkCollision(cell, temperature, u, newOmega);
  statistics.incrementStats(temperature, uSqr);
}

template<typename T, template<typename U> class Lattice>
void SmagorinskyParticleAdvectionDiffusionBGKdynamics<T, Lattice>::staticCollide(Cell<T,Lattice>& cell, const T u[Lattice<T>::d],
    LatticeStatistics<T>& statistics )
{
  assert( false );
}

template<typename T, template<typename U> class Lattice>
T SmagorinskyParticleAdvectionDiffusionBGKdynamics<T, Lattice>::getSmagorinskyOmega(Cell<T,Lattice>& cell)
{
  T temperature, uTemp[Lattice<T>::d], pi[util::TensorVal<Lattice<T> >::n];
  this->_momenta.computeAllMomenta(cell, temperature, uTemp, pi);
  T newOmega = computeOmega(this->getOmega(), preFactor, temperature, pi);
  return newOmega;
}

template<typename T, template<typename U> class Lattice>
void SmagorinskyParticleAdvectionDiffusionBGKdynamics<T, Lattice>::setOmega(T omega_)
{
  preFactor = computePreFactor(omega_, smagoConst, dx, dt);
}

template<typename T, template<typename U> class Lattice>
T SmagorinskyParticleAdvectionDiffusionBGKdynamics<T, Lattice>::computePreFactor(T omega_, T smagoConst_, T dx_, T dt_)
{
  return (T)(smagoConst_*smagoConst_*dx_*dx_)*Lattice<T>::invCs2/dt_*4*sqrt(2);
}

template<typename T, template<typename U> class Lattice>
T SmagorinskyParticleAdvectionDiffusionBGKdynamics<T, Lattice>::computeOmega(T omega0, T preFactor_, T rho, T pi[util::TensorVal<Lattice<T> >::n] )
{
  T PiNeqNormSqr = pi[0]*pi[0] + 2.0*pi[1]*pi[1] + pi[2]*pi[2];
  if (util::TensorVal<Lattice<T> >::n == 6) {
    PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2*pi[4]*pi[4] +pi[5]*pi[5];
  }
  T PiNeqNorm    = sqrt(PiNeqNormSqr);
  /// Molecular realaxation time
  T tau_mol = 1. /omega0;
  /// Turbulent realaxation time
  T tau_turb = 0.5*(sqrt(tau_mol*tau_mol+(preFactor_*tau_eff*PiNeqNorm))-tau_mol);
  /// Effective realaxation time
  tau_eff = tau_mol+tau_turb;
  T omega_new= 1./tau_eff;
  return omega_new;
}

//==================================================================//
//=========== BGK Model for Advection diffusion with Stokes Drag ====//
//==================================================================//

template<typename T, template<typename U> class Lattice>
ParticleAdvectionDiffusionBGKdynamics<T, Lattice>::ParticleAdvectionDiffusionBGKdynamics (
  T omega_, Momenta<T, Lattice>& momenta_ )
  : AdvectionDiffusionBGKdynamics<T,Lattice>(omega_,momenta_), omega( omega_ )
{ }

template<typename T, template<typename U> class Lattice>
void ParticleAdvectionDiffusionBGKdynamics<T, Lattice>::collide( Cell<T, Lattice>& cell,
    LatticeStatistics<T>& statistics )
{
  T temperature = this->_momenta.computeRho( cell );
  int offset = (statistics.getTime() % 2 == 0) ? Lattice<T>::ExternalField::velocityBeginsAt : Lattice<T>::ExternalField::velocity2BeginsAt;
  const T* u = cell.getExternal(offset);
  T uSqr = lbHelpers<T, Lattice>::
           bgkCollision( cell, temperature, u, omega );
  statistics.incrementStats( temperature, uSqr );
}


//==================================================================//
//================= MRT Model for Advection diffusion ==============//
//==================================================================//

template<typename T, template<typename U> class Lattice>
AdvectionDiffusionMRTdynamics<T, Lattice>::AdvectionDiffusionMRTdynamics(
    T omega, Momenta<T, Lattice>& momenta) :
    BasicDynamics<T, Lattice>(momenta), _omega(omega) {
  T rt[Lattice<T>::q]; // relaxation times vector.
  for (int iPop = 0; iPop < Lattice<T>::q; ++iPop) {
    rt[iPop] = Lattice<T>::S[iPop];
  }
  for (int iPop = 0; iPop < Lattice<T>::shearIndexes; ++iPop) {
    rt[Lattice<T>::shearViscIndexes[iPop]] = omega;
  }
  for (int iPop = 0; iPop < Lattice<T>::q; ++iPop) {
    for (int jPop = 0; jPop < Lattice<T>::q; ++jPop) {
      invM_S[iPop][jPop] = T();
      for (int kPop = 0; kPop < Lattice<T>::q; ++kPop) {
        if (kPop == jPop) {
          invM_S[iPop][jPop] += Lattice<T>::invM[iPop][kPop] * rt[kPop];
        }
      }
    }
  }

}

template<typename T, template<typename U> class Lattice>
T AdvectionDiffusionMRTdynamics<T, Lattice>::computeEquilibrium(int iPop, T rho,
    const T u[Lattice<T>::d], T uSqr) const {
  return lbHelpers<T, Lattice>::equilibrium(iPop, rho, u, uSqr);
}

template<typename T, template<typename U> class Lattice>
void AdvectionDiffusionMRTdynamics<T, Lattice>::collide(Cell<T, Lattice>& cell,
    LatticeStatistics<T>& statistics) {
  T temperature = this->_momenta.computeRho(cell);
  const T* u = cell.getExternal(Lattice<T>::ExternalField::velocityBeginsAt);

  T uSqr = lbHelpers<T, Lattice>::mrtCollision(cell, temperature, u, invM_S);

  statistics.incrementStats(temperature, uSqr);
}

template<typename T, template<typename U> class Lattice>
void AdvectionDiffusionMRTdynamics<T, Lattice>::staticCollide(
    Cell<T, Lattice>& cell, const T u[Lattice<T>::d],
    LatticeStatistics<T>& statistics) {
  assert(false);
}

template<typename T, template<typename U> class Lattice>
T AdvectionDiffusionMRTdynamics<T, Lattice>::getOmega() const {
  return _omega;
}

template<typename T, template<typename U> class Lattice>
void AdvectionDiffusionMRTdynamics<T, Lattice>::setOmega(T omega) {
  _omega = omega;
}



} // namespace olb



#endif
