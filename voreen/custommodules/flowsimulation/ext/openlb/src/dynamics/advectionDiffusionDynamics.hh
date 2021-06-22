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
#include "dynamics/lbHelpers.h"

namespace olb {


////////////////////// Class AdvectionDiffusionRLBdynamics //////////////////////////

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 */
//==================================================================//
//============= Regularized Model for Advection diffusion===========//
//==================================================================//

template<typename T, typename DESCRIPTOR>
AdvectionDiffusionRLBdynamics<T, DESCRIPTOR>::AdvectionDiffusionRLBdynamics (
  T omega_, Momenta<T, DESCRIPTOR>& momenta_ )
  : BasicDynamics<T, DESCRIPTOR>( momenta_ ),
    omega( omega_ )
{
  this->getName() = "AdvectionDiffusionRLBdynamics";  
}

template<typename T, typename DESCRIPTOR>
T AdvectionDiffusionRLBdynamics<T, DESCRIPTOR>::computeEquilibrium( int iPop, T rho,
    const T u[DESCRIPTOR::d], T uSqr ) const
{
  return lbHelpers<T, DESCRIPTOR>::equilibriumFirstOrder( iPop, rho, u );
}


template<typename T, typename DESCRIPTOR>
void AdvectionDiffusionRLBdynamics<T, DESCRIPTOR>::collide( Cell<T,DESCRIPTOR>& cell,
    LatticeStatistics<T>& statistics )
{
  T temperature = this->_momenta.computeRho( cell );

  const auto u = cell.template getField<descriptors::VELOCITY>();

  T uSqr = lbHelpers<T, DESCRIPTOR>::
           rlbCollision( cell, temperature, u.data(), omega );

  statistics.incrementStats( temperature, uSqr );
}

template<typename T, typename DESCRIPTOR>
T AdvectionDiffusionRLBdynamics<T, DESCRIPTOR>::getOmega() const
{
  return omega;
}

template<typename T, typename DESCRIPTOR>
void AdvectionDiffusionRLBdynamics<T, DESCRIPTOR>::setOmega( T omega_ )
{
  omega = omega_;
}


////////////////////// Class CombinedAdvectionDiffusionRLBdynamics /////////////////////////

template<typename T, typename DESCRIPTOR, typename Dynamics>
CombinedAdvectionDiffusionRLBdynamics<T,DESCRIPTOR,Dynamics>::CombinedAdvectionDiffusionRLBdynamics (
  T omega, Momenta<T,DESCRIPTOR>& momenta )
  : BasicDynamics<T,DESCRIPTOR>(momenta),
    _boundaryDynamics(omega, momenta)
{
  this->getName() = "CombinedAdvectionDiffusionRLBdynamics";  
}

template<typename T, typename DESCRIPTOR, typename Dynamics>
T CombinedAdvectionDiffusionRLBdynamics<T,DESCRIPTOR,Dynamics>::
computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr) const
{
  return lbHelpers<T, DESCRIPTOR>::equilibriumFirstOrder( iPop, rho, u );
}

template<typename T, typename DESCRIPTOR, typename Dynamics>
void CombinedAdvectionDiffusionRLBdynamics<T,DESCRIPTOR,Dynamics>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{
  typedef DESCRIPTOR l;

  T temperature = this->_momenta.computeRho( cell );
  auto u = cell.template getField<descriptors::VELOCITY>();

  T jNeq[DESCRIPTOR::d];
  // this->_momenta.computeJ( cell, jNeq );
  dynamic_cast<AdvectionDiffusionBoundaryMomenta<T,DESCRIPTOR>&>(this->_momenta).computeJneq( cell, jNeq );
  // cout << jNeq[0] << " " << jNeq[1] << " " << u[0] << " " << u[1] << endl;
  // stripe of equilibrium part u * T
  // for ( int iD = 0; iD < DESCRIPTOR::d; ++iD ) {
  //   jNeq[iD] -= u[iD] * temperature;
  // }
  // cout << jNeq[0] << " " << jNeq[1] << " " << u[0] << " " << u[1] << endl;

  for (int iPop = 0; iPop < l::q; ++iPop) {
    cell[iPop] = lbHelpers<T, DESCRIPTOR>::equilibriumFirstOrder( iPop, temperature, u.data() )
                 + firstOrderLbHelpers<T,DESCRIPTOR>::fromJneqToFneq(iPop, jNeq);
    // cout << firstOrderLbHelpers<T,DESCRIPTOR>::fromJneqToFneq(iPop,jNeq) << " ";
  }
  // cout << endl;
  // lbHelpers<T,DESCRIPTOR>::computeJ(cell, jNeq);
  // cout << jNeq[0] << " " << jNeq[1] << " " << u[0] << " " << u[1] << endl;

  _boundaryDynamics.collide(cell, statistics);
  // lbHelpers<T,DESCRIPTOR>::computeJ(cell, jNeq);
  // cout << jNeq[0] << " " << jNeq[1] << " " << u[0] << " " << u[1] << endl;
}

template<typename T, typename DESCRIPTOR, typename Dynamics>
T CombinedAdvectionDiffusionRLBdynamics<T,DESCRIPTOR,Dynamics>::getOmega() const
{
  return _boundaryDynamics.getOmega();
}

template<typename T, typename DESCRIPTOR, typename Dynamics>
void CombinedAdvectionDiffusionRLBdynamics<T,DESCRIPTOR,Dynamics>::setOmega(T omega)
{
  _boundaryDynamics.setOmega(omega);
}

//==================================================================//
//============= BGK Model for Advection diffusion===========//
//==================================================================//

template<typename T, typename DESCRIPTOR>
AdvectionDiffusionBGKdynamics<T, DESCRIPTOR>::AdvectionDiffusionBGKdynamics (
  T omega, Momenta<T, DESCRIPTOR>& momenta )
  : BasicDynamics<T, DESCRIPTOR>( momenta ),
    _omega(omega)
{
  this->getName() = "AdvectionDiffusionBGKdynamics";  
}

template<typename T, typename DESCRIPTOR>
AdvectionDiffusionBGKdynamics<T, DESCRIPTOR>::AdvectionDiffusionBGKdynamics (
  const UnitConverter<T,DESCRIPTOR>& converter, Momenta<T, DESCRIPTOR>& momenta )
  : BasicDynamics<T, DESCRIPTOR>( momenta ),
    _omega(converter.getLatticeRelaxationFrequency())
{
  this->getName() = "AdvectionDiffusionBGKdynamics";  
}

template<typename T, typename DESCRIPTOR>
T AdvectionDiffusionBGKdynamics<T, DESCRIPTOR>::computeEquilibrium( int iPop, T rho,
    const T u[DESCRIPTOR::d], T uSqr ) const
{
  return lbHelpers<T, DESCRIPTOR>::equilibriumFirstOrder( iPop, rho, u );
}


template<typename T, typename DESCRIPTOR>
void AdvectionDiffusionBGKdynamics<T, DESCRIPTOR>::collide( Cell<T,DESCRIPTOR>& cell,
    LatticeStatistics<T>& statistics )
{
  T temperature = this->_momenta.computeRho( cell );
  const auto u = cell.template getField<descriptors::VELOCITY>();

  T uSqr = lbHelpers<T, DESCRIPTOR>::
           bgkCollision( cell, temperature, u.data(), _omega );

  statistics.incrementStats( temperature, uSqr );
}

template<typename T, typename DESCRIPTOR>
T AdvectionDiffusionBGKdynamics<T, DESCRIPTOR>::getOmega() const
{
  return _omega;
}

template<typename T, typename DESCRIPTOR>
void AdvectionDiffusionBGKdynamics<T, DESCRIPTOR>::setOmega( T omega )
{
  _omega = omega;
}


//==================================================================//
//============= TRT Model for Advection diffusion===========//
//==================================================================//

template<typename T, typename DESCRIPTOR>
AdvectionDiffusionTRTdynamics<T, DESCRIPTOR>::AdvectionDiffusionTRTdynamics (
  T omega, Momenta<T, DESCRIPTOR>& momenta, T magicParameter )
  : AdvectionDiffusionBGKdynamics<T, DESCRIPTOR>( omega, momenta ),
    _omega2(1/(magicParameter/(1/omega-0.5)+0.5)), _magicParameter(magicParameter)
{
  this->getName() = "AdvectionDiffusionTRTdynamics";  
}

template<typename T, typename DESCRIPTOR>
void AdvectionDiffusionTRTdynamics<T, DESCRIPTOR>::collide( Cell<T,DESCRIPTOR>& cell,
    LatticeStatistics<T>& statistics )
{
  T temperature = this->_momenta.computeRho( cell );
  const auto u = cell.template getField<descriptors::VELOCITY>();
  const T uSqr = util::normSqr<T,DESCRIPTOR::d>(u);

  T fPlus[DESCRIPTOR::q], fMinus[DESCRIPTOR::q];
  T fEq[DESCRIPTOR::q], fEqPlus[DESCRIPTOR::q], fEqMinus[DESCRIPTOR::q];

  for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    fPlus[iPop] = 0.5 * ( cell[iPop] + cell[descriptors::opposite<DESCRIPTOR>(iPop)] );
    fMinus[iPop] = 0.5 * ( cell[iPop] - cell[descriptors::opposite<DESCRIPTOR>(iPop)] );
    fEq[iPop] = lbHelpers<T, DESCRIPTOR>::equilibrium(iPop, temperature, u.data(), uSqr);
  }
  for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    fEqPlus[iPop] = 0.5 * ( fEq[iPop] + fEq[descriptors::opposite<DESCRIPTOR>(iPop)] );
    fEqMinus[iPop] = 0.5 * ( fEq[iPop] - fEq[descriptors::opposite<DESCRIPTOR>(iPop)] );
  }
  for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    cell[iPop] -= _omega2 * (fPlus[iPop] - fEqPlus[iPop]) + this->_omega * (fMinus[iPop] - fEqMinus[iPop]);
  }

  statistics.incrementStats( temperature, 0. );
}



//==================================================================//
//============= BGK Model for Advection diffusion with source ===========//
//==================================================================//

template<typename T, typename DESCRIPTOR>
SourcedAdvectionDiffusionBGKdynamics<T, DESCRIPTOR>::SourcedAdvectionDiffusionBGKdynamics (
  T omega, Momenta<T, DESCRIPTOR>& momenta )
  : AdvectionDiffusionBGKdynamics<T, DESCRIPTOR>( omega, momenta ), _omegaMod(1. - 0.5 * omega)
{
  this->getName() = "SourcedAdvectionDiffusionBGKdynamics";  
  OLB_PRECONDITION( DESCRIPTOR::template provides<descriptors::SOURCE>() );
}


template<typename T, typename DESCRIPTOR>
void SourcedAdvectionDiffusionBGKdynamics<T, DESCRIPTOR>::collide( Cell<T,DESCRIPTOR>& cell,
    LatticeStatistics<T>& statistics )
{
  const auto u = cell.template getField<descriptors::VELOCITY>();
  const T source = cell.template getField<descriptors::SOURCE>();
  const T temperature = this->_momenta.computeRho( cell ) + 0.5 * source;

  T uSqr = lbHelpers<T, DESCRIPTOR>::
           bgkCollision( cell, temperature, u.data(), this->_omega );

  // Q / t_i = ( 1 - 0.5*_omega) * q
  const T sourceMod = source * _omegaMod;
  for ( int iPop = 0; iPop < DESCRIPTOR::q; iPop++ ) {
    cell[iPop] += sourceMod * descriptors::t<T,DESCRIPTOR>(iPop);
  }

  statistics.incrementStats( temperature, uSqr );
}

template<typename T, typename DESCRIPTOR>
T SourcedAdvectionDiffusionBGKdynamics<T,DESCRIPTOR>::computeRho(ConstCell<T,DESCRIPTOR>& cell) const
{
  const T source = cell.template getField<descriptors::SOURCE>();
  return this->_momenta.computeRho( cell ) + 0.5 * source;
}

template<typename T, typename DESCRIPTOR>
void SourcedAdvectionDiffusionBGKdynamics<T,DESCRIPTOR>::computeRhoU (
  ConstCell<T,DESCRIPTOR>& cell,
  T& rho, T u[DESCRIPTOR::d]) const
{
  this->_momenta.computeRhoU( cell, rho, u );
  const T source = cell.template getField<descriptors::SOURCE>();
  rho += 0.5 * source;
}

//==================================================================//
//============= BGK Model for Advection diffusion with phase change  ===========//
//==================================================================//

template<typename T, typename DESCRIPTOR>
TotalEnthalpyAdvectionDiffusionBGKdynamics<T, DESCRIPTOR>::TotalEnthalpyAdvectionDiffusionBGKdynamics (
  T omega, Momenta<T, DESCRIPTOR>& momenta, T T_s_, T T_l_, T cp_s_, T cp_l_, T lambda_s_, T lambda_l_, T l_)
  : AdvectionDiffusionBGKdynamics<T, DESCRIPTOR>( omega, momenta ),
    _T_s(T_s_), _T_l(T_l_), _cp_s(cp_s_), _cp_l(cp_l_), _lambda_s(lambda_s_), _lambda_l(lambda_l_), _l(l_), _H_s(cp_s_ * T_s_), _H_l(cp_l_ * T_l_ + l_), _cp_ref(2. * cp_s_ * cp_l_ / (cp_s_ + cp_l_))
{
  this->getName() = "TotalEnthalpyAdvectionDiffusionBGKdynamics";  
  // OLB_PRECONDITION( DESCRIPTOR::template provides<descriptors::SOURCE>() );
}

template<typename T, typename DESCRIPTOR>
T TotalEnthalpyAdvectionDiffusionBGKdynamics<T, DESCRIPTOR>::computeEquilibrium( int iPop, T rho,
    const T u[DESCRIPTOR::d], T uSqr ) const
{
  const T temperature = computeTemperature( rho );
  const T liquid_fraction = computeLiquidFraction( rho );
  const T cp = (1. - liquid_fraction) * _cp_s + liquid_fraction * _cp_l;

  T c_u = T();
  for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
    c_u += descriptors::c<DESCRIPTOR>(iPop,iD)*u[iD];
  }

  T f_eq;
  if (iPop == 0) {
    f_eq = rho - _cp_ref * temperature
           + cp * temperature * descriptors::t<T,DESCRIPTOR>(0) * ( _cp_ref / cp
               - descriptors::invCs2<T,DESCRIPTOR>() * (T)0.5 * uSqr )
           - descriptors::t<T,DESCRIPTOR>(0);
  }
  else {
    f_eq = cp * temperature * descriptors::t<T,DESCRIPTOR>(iPop) * ( _cp_ref / cp + descriptors::invCs2<T,DESCRIPTOR>() * c_u
           + descriptors::invCs2<T,DESCRIPTOR>() * descriptors::invCs2<T,DESCRIPTOR>() * (T)0.5 * c_u *c_u
           - descriptors::invCs2<T,DESCRIPTOR>() * (T)0.5 * uSqr )
           - descriptors::t<T,DESCRIPTOR>(iPop);
  }

  return f_eq;
}


template<typename T, typename DESCRIPTOR>
void TotalEnthalpyAdvectionDiffusionBGKdynamics<T, DESCRIPTOR>::collide( Cell<T,DESCRIPTOR>& cell,
    LatticeStatistics<T>& statistics )
{
  const T enthalpy = this->_momenta.computeRho( cell );
  const T temperature = computeTemperature( enthalpy );
  const T liquid_fraction = computeLiquidFraction( enthalpy );
  const T lambda = (1. - liquid_fraction) * _lambda_s + liquid_fraction * _lambda_l;
  const T cp = (1. - liquid_fraction) * _cp_s + liquid_fraction * _cp_l;
  const T omega = 1. / ( lambda / _cp_ref * descriptors::invCs2<T,DESCRIPTOR>() + 0.5 );

  auto u = cell.template getFieldPointer<descriptors::VELOCITY>();

  const T uSqr = util::normSqr<T,DESCRIPTOR::d>(u);

  const T f_eq = enthalpy - _cp_ref * temperature
                 + cp * temperature * descriptors::t<T,DESCRIPTOR>(0) * ( _cp_ref / cp
                     - descriptors::invCs2<T,DESCRIPTOR>() * (T)0.5 * uSqr )
                 - descriptors::t<T,DESCRIPTOR>(0);
  cell[0] *= (T)1 - omega;
  cell[0] += omega * f_eq;
  for (int iPop=1; iPop < DESCRIPTOR::q; ++iPop) {
    T c_u = T();
    for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
      c_u += descriptors::c<DESCRIPTOR>(iPop,iD)*u[iD];
    }
    const T f_eq = cp * temperature * descriptors::t<T,DESCRIPTOR>(iPop) * ( _cp_ref / cp + descriptors::invCs2<T,DESCRIPTOR>() * c_u
                   + descriptors::invCs2<T,DESCRIPTOR>() * descriptors::invCs2<T,DESCRIPTOR>() * (T)0.5 * c_u *c_u
                   - descriptors::invCs2<T,DESCRIPTOR>() * (T)0.5 * uSqr )
                   - descriptors::t<T,DESCRIPTOR>(iPop);
    cell[iPop] *= (T)1 - omega;
    cell[iPop] += omega * f_eq;
  }

  statistics.incrementStats( enthalpy, uSqr );
}

template<typename T, typename DESCRIPTOR>
T TotalEnthalpyAdvectionDiffusionBGKdynamics<T, DESCRIPTOR>::computeTemperature( T enthalpy ) const
{
  T temperature = T();
  if (enthalpy <= _H_s) {
    temperature = _T_s - (_H_s - enthalpy) / _cp_s;
  }
  else if (enthalpy >= _H_l) {
    temperature = _T_l + (enthalpy - _H_l) / _cp_l;
  }
  else {
    temperature = (_H_l - enthalpy) / (_H_l - _H_s) * _T_s + (enthalpy - _H_s) / (_H_l - _H_s) * _T_l;
  }
  return temperature;
}

template<typename T, typename DESCRIPTOR>
T TotalEnthalpyAdvectionDiffusionBGKdynamics<T, DESCRIPTOR>::computeLiquidFraction( T enthalpy ) const
{
  T liquid_fraction = T();
  if (enthalpy <= _H_s) {
    liquid_fraction = 0.;
  }
  else if (enthalpy >= _H_l) {
    liquid_fraction = 1.;
  }
  else {
    liquid_fraction = (enthalpy - _H_s) / _l;
  }
  return liquid_fraction;
}

//==================================================================//
//============= TRT Model for Advection diffusion with phase change  ===========//
//==================================================================//

template<typename T, typename DESCRIPTOR>
TotalEnthalpyAdvectionDiffusionTRTdynamics<T, DESCRIPTOR>::TotalEnthalpyAdvectionDiffusionTRTdynamics (
  T omega, Momenta<T, DESCRIPTOR>& momenta, T magicParameter_, T T_s_, T T_l_, T cp_s_, T cp_l_, T lambda_s_, T lambda_l_, T l_ )
  : TotalEnthalpyAdvectionDiffusionBGKdynamics<T, DESCRIPTOR>( omega, momenta, T_s_, T_l_, cp_s_, cp_l_, lambda_s_, lambda_l_, l_ ), _magicParameter(magicParameter_)
{
}

template<typename T, typename DESCRIPTOR>
void TotalEnthalpyAdvectionDiffusionTRTdynamics<T, DESCRIPTOR>::collide( Cell<T,DESCRIPTOR>& cell,
    LatticeStatistics<T>& statistics )
{
  const T enthalpy = this->_momenta.computeRho( cell );
  const T liquid_fraction = this->computeLiquidFraction( enthalpy );
  const T lambda = (1. - liquid_fraction) * this->_lambda_s + liquid_fraction * this->_lambda_l;
  const T omega = 1. / ( lambda / this->_cp_ref * descriptors::invCs2<T,DESCRIPTOR>() + 0.5 );
  const T omega2 = 1./(_magicParameter/(1./omega-0.5)+0.5);

  const auto u = cell.template getField<descriptors::VELOCITY>();

  const T uSqr = util::normSqr<T,DESCRIPTOR::d>(u);

  T fPlus[DESCRIPTOR::q], fMinus[DESCRIPTOR::q];
  T fEq[DESCRIPTOR::q], fEqPlus[DESCRIPTOR::q], fEqMinus[DESCRIPTOR::q];

  for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    fPlus[iPop] = 0.5 * ( cell[iPop] + cell[descriptors::opposite<DESCRIPTOR>(iPop)] );
    fMinus[iPop] = 0.5 * ( cell[iPop] - cell[descriptors::opposite<DESCRIPTOR>(iPop)] );
    fEq[iPop] = this->computeEquilibrium( iPop, enthalpy, u.data(), uSqr );
  }

  for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    fEqPlus[iPop] = 0.5 * ( fEq[iPop] + fEq[descriptors::opposite<DESCRIPTOR>(iPop)] );
    fEqMinus[iPop] = 0.5 * ( fEq[iPop] - fEq[descriptors::opposite<DESCRIPTOR>(iPop)] );
  }

  for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    cell[iPop] -= omega2 * (fPlus[iPop] - fEqPlus[iPop]) + omega * (fMinus[iPop] - fEqMinus[iPop]);
  }

  statistics.incrementStats( enthalpy, uSqr );
}

//==================================================================//
//============= BGK Model for advection diffusion phase field equation  ===========//
//==================================================================//

template<typename T, typename DESCRIPTOR>
PhaseFieldAdvectionDiffusionBGKdynamics<T, DESCRIPTOR>::PhaseFieldAdvectionDiffusionBGKdynamics (
  T omega, Momenta<T, DESCRIPTOR>& momenta, T interface_thickness )
  : AdvectionDiffusionBGKdynamics<T, DESCRIPTOR>( omega, momenta ),
    _mobility( (1.0/omega-0.5)/descriptors::invCs2<T,DESCRIPTOR>() ),
    _interface_thickness( interface_thickness )
{
  this->getName() = "PhaseFieldAdvectionDiffusionBGKdynamics";  
  OLB_PRECONDITION( DESCRIPTOR::template provides<descriptors::SOURCE>() );
}


template<typename T, typename DESCRIPTOR>
void PhaseFieldAdvectionDiffusionBGKdynamics<T, DESCRIPTOR>::collide( Cell<T,DESCRIPTOR>& cell,
    LatticeStatistics<T>& statistics )
{
  const T phi = this->_momenta.computeRho( cell );

  const auto u = cell.template getField<descriptors::VELOCITY>();
  const T uSqr = util::normSqr<T,DESCRIPTOR::d>(u.data());

  auto n = cell.template getFieldPointer<descriptors::INTERPHASE_NORMAL>();
  for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    T c_u = T();
    T c_n = T();
    for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
      c_u += descriptors::c<DESCRIPTOR>(iPop,iD)*u[iD];
      c_n += descriptors::c<DESCRIPTOR>(iPop,iD)*n[iD];
    }
    T f_eq = this->computeEquilibrium( iPop, phi, u.data(), uSqr );
    f_eq += descriptors::t<T,DESCRIPTOR>(iPop) * _mobility * descriptors::invCs2<T,DESCRIPTOR>() *
            ( 4.0 * phi * (1.0 - phi) / _interface_thickness ) * c_n;

    cell[iPop] *= (T)1 - this->_omega;
    cell[iPop] += this->_omega * f_eq;
  }

  statistics.incrementStats( phi, uSqr );
}

//==================================================================================//
//=========== BGK Model for Advection diffusion with Stokes drag and Smagorinsky====//
//==================================================================================//

template<typename T, typename DESCRIPTOR>
SmagorinskyParticleAdvectionDiffusionBGKdynamics<T, DESCRIPTOR>::SmagorinskyParticleAdvectionDiffusionBGKdynamics (
  T omega_, Momenta<T,DESCRIPTOR>& momenta_, T smagoConst_, T dx_, T dt_)
  : AdvectionDiffusionBGKdynamics<T,DESCRIPTOR>(omega_,momenta_), smagoConst(smagoConst_), preFactor(computePreFactor(omega_,smagoConst_, dx_, dt_) )
{
  this->getName() = "SmagorinskyParticleAdvectionDiffusionBGKdynamics";  
}

template<typename T, typename DESCRIPTOR>
void SmagorinskyParticleAdvectionDiffusionBGKdynamics<T, DESCRIPTOR>::collide(Cell<T,DESCRIPTOR>& cell, LatticeStatistics<T>& statistics )
{
  T temperature, uad[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR >::n];
  this->_momenta.computeAllMomenta(cell, temperature, uad, pi);
  auto u = (statistics.getTime() % 2 == 0) ? cell.template getField<descriptors::VELOCITY>() : cell.template getField<descriptors::VELOCITY2>();
  T newOmega = computeOmega(this->getOmega(), preFactor, temperature, pi);
  T uSqr = lbHelpers<T,DESCRIPTOR>::bgkCollision(cell, temperature, u.data(), newOmega);
  statistics.incrementStats(temperature, uSqr);
}

template<typename T, typename DESCRIPTOR>
T SmagorinskyParticleAdvectionDiffusionBGKdynamics<T, DESCRIPTOR>::getSmagorinskyOmega(Cell<T,DESCRIPTOR>& cell)
{
  T temperature, uTemp[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR >::n];
  this->_momenta.computeAllMomenta(cell, temperature, uTemp, pi);
  T newOmega = computeOmega(this->getOmega(), preFactor, temperature, pi);
  return newOmega;
}

template<typename T, typename DESCRIPTOR>
void SmagorinskyParticleAdvectionDiffusionBGKdynamics<T, DESCRIPTOR>::setOmega(T omega_)
{
  preFactor = computePreFactor(omega_, smagoConst, dx, dt);
}

template<typename T, typename DESCRIPTOR>
T SmagorinskyParticleAdvectionDiffusionBGKdynamics<T, DESCRIPTOR>::computePreFactor(T omega_, T smagoConst_, T dx_, T dt_)
{
  return (T)(smagoConst_*smagoConst_*dx_*dx_)*descriptors::invCs2<T,DESCRIPTOR>()/dt_*4*sqrt(2);
}

template<typename T, typename DESCRIPTOR>
T SmagorinskyParticleAdvectionDiffusionBGKdynamics<T, DESCRIPTOR>::computeOmega(T omega0, T preFactor_, T rho, T pi[util::TensorVal<DESCRIPTOR >::n] )
{
  T PiNeqNormSqr = pi[0]*pi[0] + 2.0*pi[1]*pi[1] + pi[2]*pi[2];
  if (util::TensorVal<DESCRIPTOR >::n == 6) {
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

template<typename T, typename DESCRIPTOR>
ParticleAdvectionDiffusionBGKdynamics<T, DESCRIPTOR>::ParticleAdvectionDiffusionBGKdynamics (
  T omega_, Momenta<T, DESCRIPTOR>& momenta_ )
  : AdvectionDiffusionBGKdynamics<T,DESCRIPTOR>(omega_,momenta_), omega( omega_ )
{
  this->getName() = "ParticleAdvectionDiffusionBGKdynamics";  
}

template<typename T, typename DESCRIPTOR>
void ParticleAdvectionDiffusionBGKdynamics<T, DESCRIPTOR>::collide( Cell<T,DESCRIPTOR>& cell,
    LatticeStatistics<T>& statistics )
{
  T temperature = this->_momenta.computeRho( cell );
  const auto u = (statistics.getTime() % 2 == 0) ? cell.template getField<descriptors::VELOCITY>() : cell.template getField<descriptors::VELOCITY2>();
  T uSqr = lbHelpers<T, DESCRIPTOR>::
           bgkCollision( cell, temperature, u.data(), omega );
  statistics.incrementStats( temperature, uSqr );
}


//==================================================================//
//================= MRT Model for Advection diffusion ==============//
//==================================================================//

template<typename T, typename DESCRIPTOR>
AdvectionDiffusionMRTdynamics<T, DESCRIPTOR>::AdvectionDiffusionMRTdynamics(
  T omega, Momenta<T, DESCRIPTOR>& momenta) :
  BasicDynamics<T, DESCRIPTOR>(momenta), _omega(omega)
{
  this->getName() = "AdvectionDiffusionMRTdynamics";  
  T rt[DESCRIPTOR::q]; // relaxation times vector.
  for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
    rt[iPop] = descriptors::s<T,DESCRIPTOR>(iPop);
  }
  for (int iPop = 0; iPop < descriptors::shearIndexes<DESCRIPTOR>(); ++iPop) {
    rt[descriptors::shearViscIndexes<DESCRIPTOR>(iPop)] = omega;
  }
  for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
    for (int jPop = 0; jPop < DESCRIPTOR::q; ++jPop) {
      invM_S[iPop][jPop] = T();
      for (int kPop = 0; kPop < DESCRIPTOR::q; ++kPop) {
        if (kPop == jPop) {
          invM_S[iPop][jPop] += descriptors::invM<T,DESCRIPTOR>(iPop,kPop) * rt[kPop];
        }
      }
    }
  }

}

template<typename T, typename DESCRIPTOR>
T AdvectionDiffusionMRTdynamics<T, DESCRIPTOR>::computeEquilibrium(int iPop, T rho,
    const T u[DESCRIPTOR::d], T uSqr) const
{
  return lbHelpers<T, DESCRIPTOR>::equilibrium(iPop, rho, u, uSqr);
}

template<typename T, typename DESCRIPTOR>
void AdvectionDiffusionMRTdynamics<T, DESCRIPTOR>::collide(Cell<T,DESCRIPTOR>& cell,
    LatticeStatistics<T>& statistics)
{
  T temperature = this->_momenta.computeRho(cell);
  const auto u = cell.template getField<descriptors::VELOCITY>();

  T uSqr = lbHelpers<T, DESCRIPTOR>::mrtCollision(cell, temperature, u.data(), invM_S);

  statistics.incrementStats(temperature, uSqr);
}

template<typename T, typename DESCRIPTOR>
T AdvectionDiffusionMRTdynamics<T, DESCRIPTOR>::getOmega() const
{
  return _omega;
}

template<typename T, typename DESCRIPTOR>
void AdvectionDiffusionMRTdynamics<T, DESCRIPTOR>::setOmega(T omega)
{
  _omega = omega;
}



} // namespace olb



#endif
