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

template<typename T, typename DESCRIPTOR>
PorousBGKdynamics<T,DESCRIPTOR>::PorousBGKdynamics (
  T omega_, Momenta<T,DESCRIPTOR>& momenta_)
  : BGKdynamics<T,DESCRIPTOR>(omega_,momenta_),
    omega(omega_)
{ }

template<typename T, typename DESCRIPTOR>
void PorousBGKdynamics<T,DESCRIPTOR>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho, u[DESCRIPTOR::d];
  this->_momenta.computeRhoU(cell, rho, u);
  T porosity = cell.template getField<descriptors::POROSITY>();
  for (int i=0; i<DESCRIPTOR::d; i++)  {
    u[i] *= porosity;
  }
  T uSqr = lbHelpers<T,DESCRIPTOR>::bgkCollision(cell, rho, u, omega);
  statistics.incrementStats(rho, uSqr);
}

template<typename T, typename DESCRIPTOR>
T PorousBGKdynamics<T,DESCRIPTOR>::getOmega() const
{
  return omega;
}

template<typename T, typename DESCRIPTOR>
void PorousBGKdynamics<T,DESCRIPTOR>::setOmega(T omega_)
{
  omega = omega_;
}


//////////////////// Class ExtendedPorousBGKdynamics ////////////////////

template<typename T, typename DESCRIPTOR>
ExtendedPorousBGKdynamics<T,DESCRIPTOR>::ExtendedPorousBGKdynamics (
  T omega_, Momenta<T,DESCRIPTOR>& momenta_)
  : BGKdynamics<T,DESCRIPTOR>(omega_,momenta_),
    omega(omega_)
{
  this->getName() = "ExtendedPorousBGKdynamics";
}

template<typename T, typename DESCRIPTOR>
void ExtendedPorousBGKdynamics<T,DESCRIPTOR>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho, u[DESCRIPTOR::d];
  this->_momenta.computeRhoU(cell, rho, u);
  T porosity = cell.template getField<descriptors::POROSITY>();
  auto localVelocity = cell.template getFieldPointer<descriptors::LOCAL_DRAG>();

  cell.template setField<descriptors::LOCAL_DRAG>(u);

  for (int i=0; i<DESCRIPTOR::d; i++)  {
    u[i] *= porosity[0];
    u[i] += (1.-porosity[0]) * localVelocity[i];
  }
  T uSqr = lbHelpers<T,DESCRIPTOR>::bgkCollision(cell, rho, u, omega);
  statistics.incrementStats(rho, uSqr);
}

template<typename T, typename DESCRIPTOR>
T ExtendedPorousBGKdynamics<T,DESCRIPTOR>::getOmega() const
{
  return omega;
}

template<typename T, typename DESCRIPTOR>
void ExtendedPorousBGKdynamics<T,DESCRIPTOR>::setOmega(T omega_)
{
  omega = omega_;
}

//////////////////// Class SubgridParticleBGKdynamics ////////////////////

template<typename T, typename DESCRIPTOR>
SubgridParticleBGKdynamics<T,DESCRIPTOR>::SubgridParticleBGKdynamics (
  T omega_, Momenta<T,DESCRIPTOR>& momenta_)
  : BGKdynamics<T,DESCRIPTOR>(omega_,momenta_),
    omega(omega_)
{
  _fieldTmp[0] = T();
  _fieldTmp[1] = T();
  _fieldTmp[2] = T();
  _fieldTmp[3] = T();
}

template<typename T, typename DESCRIPTOR>
void SubgridParticleBGKdynamics<T,DESCRIPTOR>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho, u[DESCRIPTOR::d];
  this->_momenta.computeRhoU(cell, rho, u);
  T porosity = cell.template getField<descriptors::POROSITY>();
  auto extVelocity = cell.template getFieldPointer<descriptors::LOCAL_DRAG>();
//  if (porosity[0] != 0) {
//    cout << "extVelocity: " << extVelocity[0] << " " <<  extVelocity[1] << " " <<  extVelocity[2] << " " << std::endl;
//    cout << "porosity: " << porosity[0] << std::endl;
//  }
  for (int i=0; i<DESCRIPTOR::d; i++)  {
    u[i] *= (1.-porosity);
    u[i] += extVelocity[i];
  }
  T uSqr = lbHelpers<T,DESCRIPTOR>::bgkCollision(cell, rho, u, omega);

  statistics.incrementStats(rho, uSqr);
  cell.template setField<descriptors::POROSITY>(0);
  cell.template setField<descriptors::VELOCITY_NUMERATOR>(0);
  cell.template setField<descriptors::VELOCITY_DENOMINATOR>(0);
}

template<typename T, typename DESCRIPTOR>
T SubgridParticleBGKdynamics<T,DESCRIPTOR>::getOmega() const
{
  return omega;
}

template<typename T, typename DESCRIPTOR>
void SubgridParticleBGKdynamics<T,DESCRIPTOR>::setOmega(T omega_)
{
  omega = omega_;
}

//////////////////// Class PorousParticleDynamics ////////////////////

template<typename T, typename DESCRIPTOR, bool isStatic>
template <bool isStatic_>
std::enable_if_t<isStatic_>
PorousParticleDynamics<T,DESCRIPTOR,isStatic>::calculate(ConstCell<T,DESCRIPTOR>& cell, T* pVelocity) {
  for (int i=0; i<DESCRIPTOR::d; i++)  {
    pVelocity[i] -= (1.-(cell.template getField<descriptors::POROSITY>())) * pVelocity[i];
  }
}

template<typename T, typename DESCRIPTOR, bool isStatic>
template <bool isStatic_>
std::enable_if_t<!isStatic_>
PorousParticleDynamics<T,DESCRIPTOR,isStatic>::calculate(Cell<T,DESCRIPTOR>& cell, T* pVelocity) {
  for (int i=0; i<DESCRIPTOR::d; i++)  {
    pVelocity[i] += (1.-cell.template getField<descriptors::POROSITY>())
        * (cell.template getFieldPointer<descriptors::VELOCITY_NUMERATOR>()[i]
          / cell.template getField<descriptors::VELOCITY_DENOMINATOR>() - pVelocity[i]);
  }
  cell.template setField<descriptors::POROSITY>(1.);
  cell.template setField<descriptors::VELOCITY_DENOMINATOR>(0.);
}

//////////////////// Class PorousParticleBGKdynamics ////////////////////

template<typename T, typename DESCRIPTOR, bool isStatic>
PorousParticleBGKdynamics<T,DESCRIPTOR,isStatic>::PorousParticleBGKdynamics (
  T omega_, Momenta<T,DESCRIPTOR>& momenta_)
  : BGKdynamics<T,DESCRIPTOR>(omega_,momenta_)
{}

template<typename T, typename DESCRIPTOR, bool isStatic>
void PorousParticleBGKdynamics<T,DESCRIPTOR,isStatic>::computeU (
  ConstCell<T,DESCRIPTOR>& cell, T u[DESCRIPTOR::d]) const
{
  T rho;
  this->_momenta.computeRhoU(cell, rho, u);
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

template<typename T, typename DESCRIPTOR, bool isStatic>
void PorousParticleBGKdynamics<T,DESCRIPTOR,isStatic>::computeRhoU (
  ConstCell<T,DESCRIPTOR>& cell, T& rho, T u[DESCRIPTOR::d]) const
{
  this->_momenta.computeRhoU(cell, rho, u);
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

template<typename T, typename DESCRIPTOR, bool isStatic>
void PorousParticleBGKdynamics<T,DESCRIPTOR,isStatic>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho, u[DESCRIPTOR::d];
  this->_momenta.computeRhoU(cell, rho, u);
  T uSqr = this->porousParticleBgkCollision(cell, rho, u, this->getOmega());
  statistics.incrementStats(rho, uSqr);
}

template<typename T, typename DESCRIPTOR, bool isStatic>
T PorousParticleBGKdynamics<T,DESCRIPTOR,isStatic>::porousParticleBgkCollision(Cell<T,DESCRIPTOR>& cell, T rho, T u[DESCRIPTOR::d], T omega)
{
  auto velDenominator = cell.template getFieldPointer<descriptors::VELOCITY_DENOMINATOR>();

#if defined(FEATURE_HLBM_SHANCHEN_FORCING)
  if (velDenominator[0] > std::numeric_limits<T>::epsilon()) {
    T  u_tmp[DESCRIPTOR::d] = { };
    for(int iDim=0; iDim<DESCRIPTOR::d; iDim++)
      u_tmp[iDim] = u[iDim];
    this->calculate(cell, u);
#ifdef FEATURE_HLBM_MLA
    auto velNumerator   = cell.template getFieldPointer<descriptors::VELOCITY_NUMERATOR>();
    const T tmp_uSqr = util::normSqr<T,DESCRIPTOR::d>(u);
    for(int iDim=0; iDim<DESCRIPTOR::d; iDim++) {
      for(int iPop=0; iPop<DESCRIPTOR::q; iPop++) {
        velNumerator[iDim] -= descriptors::c<DESCRIPTOR>(iPop,iDim) * omega
                * (lbDynamicsHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho, u, tmp_uSqr)
                   - lbDynamicsHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho, u_tmp, tmp_uSqr_2) );
      }
    }
#endif
  }
  T uSqr = lbHelpers<T,DESCRIPTOR>::bgkCollision(cell, rho, u, omega);
#elif defined(FEATURE_HLBM_GUO_FORCING)
  T  force[DESCRIPTOR::d] = { };
#ifdef FEATURE_HLBM_MLA
  T  u_saved[DESCRIPTOR::d] = { };
  for(int iDim=0; iDim<DESCRIPTOR::d; iDim++)
    u_saved[iDim] = u[iDim];
#endif
  bool particle = velDenominator[0] > std::numeric_limits<T>::epsilon();
  if (particle) {
    T  uPlus[DESCRIPTOR::d] = { };
    for(int iDim=0; iDim<DESCRIPTOR::d; iDim++)
      uPlus[iDim] = u[iDim];
    this->calculate(cell, uPlus);
    for(int iDim=0; iDim<DESCRIPTOR::d; iDim++) {
      force[iDim] = uPlus[iDim]-u[iDim];
      u[iDim] += force[iDim] * 0.5;
    }
  }
  T uSqr = lbHelpers<T,DESCRIPTOR>::bgkCollision(cell, rho, u, omega);
#ifdef FEATURE_HLBM_MLA
  for(int iDim=0; iDim<DESCRIPTOR::d; iDim++) {
    velNumerator[iDim] = 0.;
  }
#endif
  // lbHelpers::addExternalForce is restricted to a force stored in the descriptor field <FORCE>
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
#ifdef FEATURE_HLBM_MLA
    if (particle) {
      const T tmp_uSqr_mla = util::normSqr<T,DESCRIPTOR::d>(u_saved);
      const T tmp_uSqr_mla_2 = util::normSqr<T,DESCRIPTOR::d>(u);
      for(int iDim=0; iDim<DESCRIPTOR::d; iDim++) {
        velNumerator[iDim] += descriptors::c<DESCRIPTOR>(iPop,iDim) * ( omega
                                * (lbDynamicsHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho, u_saved, tmp_uSqr_mla)
                                   - lbDynamicsHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho, u, tmp_uSqr_mla_2) )
                                - forceTerm );
      }
    }
#endif
  }
#ifdef FEATURE_HLBM_MLA
  // should yield same result as above, however this tests the forcing
  // and is for now more consistent with the idea of mla
  //for(int iDim=0; iDim<DESCRIPTOR::d; iDim++)
  //  velNumerator[iDim] = -force[iDim]*rho;
#endif
#else
  // use Kuperstokh forcing by default
  T uSqr = lbHelpers<T,DESCRIPTOR>::bgkCollision(cell, rho, u, omega);
  T uPlus[DESCRIPTOR::d] = { };
  T diff[DESCRIPTOR::q] = {};
  if (velDenominator[0] > std::numeric_limits<T>::epsilon()) {
    for(int iDim=0; iDim<DESCRIPTOR::d; iDim++)
      uPlus[iDim] = u[iDim];
    this->calculate(cell, uPlus);
    const T uPlusSqr = util::normSqr<T,DESCRIPTOR::d>(uPlus);
    for(int tmp_iPop=0; tmp_iPop<DESCRIPTOR::q; tmp_iPop++) {
      diff[tmp_iPop] += lbHelpers<T,DESCRIPTOR>::equilibrium(tmp_iPop, rho, uPlus, uPlusSqr)
                      - lbHelpers<T,DESCRIPTOR>::equilibrium(tmp_iPop, rho, u, uSqr);
      cell[tmp_iPop] += diff[tmp_iPop];
    }
#ifdef FEATURE_HLBM_MLA
    for(int iDim=0; iDim<DESCRIPTOR::d; iDim++) {
      velNumerator[iDim] = 0.;
      for(int tmp_iPop=0; tmp_iPop<DESCRIPTOR::q; tmp_iPop++)
        velNumerator[iDim] -= descriptors::c<DESCRIPTOR>(tmp_iPop,iDim) * diff[tmp_iPop];
    }
    // should yield same result as above, however this tests the forcing
    // and is for now more consistent with the idea of mla
    //for(int iDim=0; iDim<DESCRIPTOR::d; iDim++)
      //*(velNumerator+iDim) = -rho*(uPlus[iDim]-u[iDim]);
#endif
  }
#endif

  return uSqr;
}


//////////////////// Class DBBParticleBGKdynamics ////////////////////

template<typename T, typename DESCRIPTOR, bool isStatic>
DBBParticleBGKdynamics<T,DESCRIPTOR,isStatic>::DBBParticleBGKdynamics (
  T omega_, Momenta<T,DESCRIPTOR>& momenta_)
  : BGKdynamics<T,DESCRIPTOR>(omega_,momenta_)
{
  this->getName() = "DBBParticleBGKdynamics";
}

template<typename T, typename DESCRIPTOR, bool isStatic>
void DBBParticleBGKdynamics<T,DESCRIPTOR,isStatic>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho, u[DESCRIPTOR::d],eta[DESCRIPTOR::q],uPlus[DESCRIPTOR::d],tmp_cell[(DESCRIPTOR::q+1)/2];
  this->_momenta.computeRhoU(cell, rho, u);
  T uSqr = this->dbbParticleBgkCollision(cell, rho, u, eta, uPlus, tmp_cell, this->getOmega());
  statistics.incrementStats(rho, uSqr);
}

template<typename T, typename DESCRIPTOR, bool isStatic>
T DBBParticleBGKdynamics<T,DESCRIPTOR,isStatic>::dbbParticleBgkCollision(Cell<T,DESCRIPTOR>& cell, T rho, T u[DESCRIPTOR::d], T eta[DESCRIPTOR::d], T uPlus[DESCRIPTOR::d],T tmp_cell[(DESCRIPTOR::q+1)/2], T omega)
{

  T tmpMomentumLoss[DESCRIPTOR::d] = { };

  auto velNumerator   = cell.template getFieldPointer<descriptors::VELOCITY_NUMERATOR>();
  auto zeta = cell.template getFieldPointer<descriptors::ZETA>();
  auto velDenominator = cell.template getFieldPointer<descriptors::VELOCITY_DENOMINATOR>();

  if(*(velDenominator)>1)
    rho/=T(*(velDenominator));
  for(int tmp_iPop=1; 2*tmp_iPop<DESCRIPTOR::q; tmp_iPop++) {
    eta[tmp_iPop]=6.*descriptors::t<T,DESCRIPTOR>(tmp_iPop)*rho*(descriptors::c<DESCRIPTOR>(tmp_iPop,0)*(*velNumerator)+descriptors::c<DESCRIPTOR>(tmp_iPop,1)*(*(velNumerator+1)));
    tmp_cell[tmp_iPop]=(*(zeta+tmp_iPop))*(-cell[tmp_iPop]+cell[descriptors::opposite<DESCRIPTOR>(tmp_iPop)]+eta[tmp_iPop]);
    cell[tmp_iPop]+=tmp_cell[tmp_iPop]/(1.+2.*(*(zeta+tmp_iPop)));
    cell[descriptors::opposite<DESCRIPTOR>(tmp_iPop)]-=tmp_cell[tmp_iPop]/(1.+2.*(*(zeta+tmp_iPop)));
    *(zeta+tmp_iPop) = 0.;
    *(zeta+descriptors::opposite<DESCRIPTOR>(tmp_iPop)) = 0.;
  }

  cell.template setField<descriptors::POROSITY>(1.);
  *velDenominator=0.;

  this->_momenta.computeRhoU(cell, rho, uPlus);

  T uPlusSqr = lbHelpers<T,DESCRIPTOR>::bgkCollision(cell, rho, uPlus, omega);

  T uSqr = util::normSqr<T,DESCRIPTOR::d>(u);

  T diff[DESCRIPTOR::q] = {};
  for(int tmp_iPop=0; tmp_iPop<DESCRIPTOR::q; tmp_iPop++) {
    diff[tmp_iPop] += lbHelpers<T,DESCRIPTOR>::equilibrium(tmp_iPop, rho, uPlus, uPlusSqr)
      - lbHelpers<T,DESCRIPTOR>::equilibrium(tmp_iPop, rho, u, uSqr);

    for(int iDim=0; iDim<DESCRIPTOR::d; iDim++) {
      tmpMomentumLoss[iDim] -= descriptors::c<DESCRIPTOR>(tmp_iPop,iDim) * diff[tmp_iPop];
    }
  }

  for(int i_dim=0; i_dim<DESCRIPTOR::d; i_dim++) {
    *(velNumerator+i_dim) = tmpMomentumLoss[i_dim];
  }

  return uSqr;
}




//////////////////// Class KrauseHBGKdynamics ////////////////////

template<typename T, typename DESCRIPTOR>
KrauseHBGKdynamics<T,DESCRIPTOR>::KrauseHBGKdynamics (T omega_,
    Momenta<T,DESCRIPTOR>& momenta_, T smagoConst_, T dx_, T dt_ )
  : BGKdynamics<T,DESCRIPTOR>(omega_,momenta_), omega(omega_), smagoConst(smagoConst_),
    preFactor(computePreFactor(omega_,smagoConst_) )
{
  _fieldTmp[0] = T(1);
  _fieldTmp[1] = T();
  _fieldTmp[2] = T();
  _fieldTmp[3] = T();
}

template<typename T, typename DESCRIPTOR>
void KrauseHBGKdynamics<T,DESCRIPTOR>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho, u[DESCRIPTOR::d];
  T newOmega[DESCRIPTOR::q];
  this->_momenta.computeRhoU(cell, rho, u);
  computeOmega(this->getOmega(), cell, preFactor, rho, u, newOmega);

  T vel_denom = cell.template getField<descriptors::VELOCITY_DENOMINATOR>();
  if (vel_denom > std::numeric_limits<T>::epsilon()) {
    T porosity = cell.template getField<descriptors::POROSITY>(); // prod(1-smoothInd)
    auto vel_num = cell.template getFieldPointer<descriptors::VELOCITY_NUMERATOR>();
    porosity = 1.-porosity; // 1-prod(1-smoothInd)
    for (int i=0; i<DESCRIPTOR::d; i++)  {
      u[i] += porosity * (vel_num[i] / vel_denom - u[i]);
    }
  }
  T uSqr = lbHelpers<T,DESCRIPTOR>::bgkCollision(cell, rho, u, newOmega);
  statistics.incrementStats(rho, uSqr);

  cell.template setField<descriptors::POROSITY>(_fieldTmp[0]);
  cell.template setField<descriptors::VELOCITY_NUMERATOR>({_fieldTmp[1], _fieldTmp[2]});
  cell.template setField<descriptors::VELOCITY_DENOMINATOR>(_fieldTmp[3]);
}

template<typename T, typename DESCRIPTOR>
T KrauseHBGKdynamics<T,DESCRIPTOR>::getOmega() const
{
  return omega;
}

template<typename T, typename DESCRIPTOR>
void KrauseHBGKdynamics<T,DESCRIPTOR>::setOmega(T omega)
{
  this->setOmega(omega);
  preFactor = computePreFactor(omega, smagoConst);
}

template<typename T, typename DESCRIPTOR>
T KrauseHBGKdynamics<T,DESCRIPTOR>::computePreFactor(T omega, T smagoConst)
{
  return (T)smagoConst*smagoConst*descriptors::invCs2<T,DESCRIPTOR>()*descriptors::invCs2<T,DESCRIPTOR>()*2*sqrt(2);
}


template<typename T, typename DESCRIPTOR>
void KrauseHBGKdynamics<T,DESCRIPTOR>::computeOmega(T omega0, Cell<T,DESCRIPTOR>& cell, T preFactor, T rho,
    T u[DESCRIPTOR::d], T newOmega[DESCRIPTOR::q])
{
  T uSqr = u[0]*u[0];
  for (int iDim=0; iDim<DESCRIPTOR::d; iDim++) {
    uSqr += u[iDim]*u[iDim];
  }
  /// Molecular realaxation time
  T tau_mol = 1./omega0;

  for (int iPop=0; iPop<DESCRIPTOR::q; iPop++) {
    T fNeq = std::fabs(cell[iPop] - lbHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho, u, uSqr));
    /// Turbulent realaxation time
    T tau_turb = 0.5*(sqrt(tau_mol*tau_mol+(preFactor*fNeq))-tau_mol);
    /// Effective realaxation time
    tau_eff = tau_mol + tau_turb;
    newOmega[iPop] = 1./tau_eff;
  }
}


//////////////////// Class ParticlePorousBGKdynamics ////////////////////
/*
template<typename T, typename DESCRIPTOR>
ParticlePorousBGKdynamics<T,DESCRIPTOR>::ParticlePorousBGKdynamics (
  T omega_, Momenta<T,DESCRIPTOR>& momenta_)
  : BGKdynamics<T,DESCRIPTOR>(omega_,momenta_),
    omega(omega_)
{ }

template<typename T, typename DESCRIPTOR>
void ParticlePorousBGKdynamics<T,DESCRIPTOR>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho, u[DESCRIPTOR::d];
  this->_momenta.computeRhoU(cell, rho, u);
  auto porosity = cell.template getFieldPointer<descriptors::POROSITY>();
  auto localVelocity = cell.template getFieldPointer<descriptors::LOCAL_DRAG>();
  for (int i=0; i<DESCRIPTOR::d; i++)  {
    u[i] *= porosity[0];
    u[i] += localVelocity[i];
  }
  T uSqr = lbHelpers<T,DESCRIPTOR>::bgkCollision(cell, rho, u, omega);
  statistics.incrementStats(rho, uSqr);

//  *cell.template getFieldPointer<descriptors::POROSITY>() = 1;
//  *cell.template getFieldPointer<descriptors::LOCAL_DRAG>() = 0.;
//  *(cell.template getFieldPointer<descriptors::LOCAL_DRAG>()+1) = 0.;
}

template<typename T, typename DESCRIPTOR>
T ParticlePorousBGKdynamics<T,DESCRIPTOR>::getOmega() const {
  return omega;
}

template<typename T, typename DESCRIPTOR>
void ParticlePorousBGKdynamics<T,DESCRIPTOR>::setOmega(T omega_) {
  omega = omega_;
}
*/

//////////////////// Class SmallParticleBGKdynamics ////////////////////

template<typename T, typename DESCRIPTOR>
SmallParticleBGKdynamics<T,DESCRIPTOR>::SmallParticleBGKdynamics (
  T omega_, Momenta<T,DESCRIPTOR>& momenta_)
  : BGKdynamics<T,DESCRIPTOR>(omega_,momenta_),
    omega(omega_)
{ }

template<typename T, typename DESCRIPTOR>
void SmallParticleBGKdynamics<T,DESCRIPTOR>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho, u[DESCRIPTOR::d];
  this->_momenta.computeRhoU(cell, rho, u);
  T porosity = cell.template getField<descriptors::POROSITY>();
  auto localVelocity = cell.template getFieldPointer<descriptors::LOCAL_DRAG>();

  //cout << porosity[0]  << " " <<   localVelocity[0]<< " " <<   localVelocity[1]<< " " <<   localVelocity[2]<<std::endl;
  for (int i=0; i<DESCRIPTOR::d; i++)  {
    u[i] *= porosity;
    u[i] += localVelocity[i];
  }
  T uSqr = lbHelpers<T,DESCRIPTOR>::bgkCollision(cell, rho, u, omega);
  statistics.incrementStats(rho, uSqr);
}

template<typename T, typename DESCRIPTOR>
T SmallParticleBGKdynamics<T,DESCRIPTOR>::getOmega() const
{
  return omega;
}

template<typename T, typename DESCRIPTOR>
void SmallParticleBGKdynamics<T,DESCRIPTOR>::setOmega(T omega_)
{
  omega = omega_;
}

////////////////////// Class PSMBGKdynamics //////////////////////////

/** \param omega relaxation parameter, related to the dynamic viscosity
 *  \param momenta a Momenta object to know how to compute velocity momenta
 */
template<typename T, typename DESCRIPTOR>
PSMBGKdynamics<T,DESCRIPTOR>::PSMBGKdynamics (
  T omega_, Momenta<T,DESCRIPTOR>& momenta_, int mode_ )
  : BGKdynamics<T,DESCRIPTOR>(omega_, momenta_),
    omega(omega_), paramA(1. / omega_ - 0.5)
{
  mode = (Mode) mode_;
}

template<typename T, typename DESCRIPTOR>
void PSMBGKdynamics<T,DESCRIPTOR>::computeU (ConstCell<T,DESCRIPTOR>& cell, T u[DESCRIPTOR::d] ) const
{
  T rho;
  this->_momenta.computeRhoU(cell, rho, u);
//  T epsilon = 1. - *(cell.template getFieldPointer<descriptors::POROSITY>());
//  // speed up paramB
//  T paramB = (epsilon * paramA) / ((1. - epsilon) + paramA);
//  // speed up paramC
//  T paramC = (1. - paramB);
//  for (int iVel=0; iVel<DESCRIPTOR::d; ++iVel) {
//    u[iVel] = paramC * u[iVel] +
//              paramB * cell.template getFieldPointer<descriptors::VELOCITY_SOLID>()[iVel];
//  }
}

template<typename T, typename DESCRIPTOR>
void PSMBGKdynamics<T,DESCRIPTOR>::computeRhoU (ConstCell<T,DESCRIPTOR>& cell, T& rho, T u[DESCRIPTOR::d] ) const
{
    this->_momenta.computeRhoU(cell, rho, u);
//  T epsilon = 1. - *(cell.template getFieldPointer<descriptors::POROSITY>());
//  // speed up paramB
//  T paramB = (epsilon * paramA) / ((1. - epsilon) + paramA);
//  // speed up paramC
//  T paramC = (1. - paramB);
//  for (int iVel=0; iVel<DESCRIPTOR::d; ++iVel) {
//    u[iVel] = paramC * u[iVel] +
//              paramB * cell.template getFieldPointer<descriptors::VELOCITY_SOLID>()[iVel];
//  }
}

template<typename T, typename DESCRIPTOR>
void PSMBGKdynamics<T,DESCRIPTOR>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho, u[DESCRIPTOR::d], uSqr;
  T epsilon = 1. - cell.template getField<descriptors::POROSITY>();

  this->_momenta.computeRhoU(cell, rho, u);
  // velocity at the boundary
  auto u_s = cell.template getField<descriptors::VELOCITY_SOLID>();

  if (epsilon < 1e-5) {
    uSqr = lbHelpers<T,DESCRIPTOR>::bgkCollision(cell, rho, u, omega);
  } else {
    // speed up paramB
    T paramB = (epsilon * paramA) / ((1. - epsilon) + paramA);
    // speed up paramC
    T paramC = (1. - paramB);

    T omega_s[DESCRIPTOR::q];
    T cell_tmp[DESCRIPTOR::q];

    const T uSqr_s = util::normSqr<T,DESCRIPTOR::d>(u_s.data());

    uSqr = util::normSqr<T,DESCRIPTOR::d>(u);
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      cell_tmp[iPop] = cell[iPop];
      switch(mode){
        case M2: omega_s[iPop] = (lbDynamicsHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho, u_s.data(), uSqr_s ) - cell[iPop])
                         + (1 - omega) * (cell[iPop] - lbDynamicsHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho, u, uSqr )); break;
        case M3: omega_s[iPop] = (cell[descriptors::opposite<DESCRIPTOR>(iPop)] - lbDynamicsHelpers<T,DESCRIPTOR>::equilibrium(descriptors::opposite<DESCRIPTOR>(iPop), rho, u_s.data(), uSqr_s ))
          - (cell[iPop] - lbDynamicsHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho, u_s.data(), uSqr_s ));
      }

    }

    uSqr = lbHelpers<T,DESCRIPTOR>::bgkCollision(cell, rho, u, omega);

    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      cell[iPop] = cell_tmp[iPop] + paramC * (cell[iPop] - cell_tmp[iPop]);
      cell[iPop] += paramB * omega_s[iPop];
    }
    for (int iVel=0; iVel<DESCRIPTOR::d; ++iVel) {
      u[iVel] = paramC * u[iVel] + paramB * u_s[iVel];
    }
    uSqr = util::normSqr<T,DESCRIPTOR::d>(u);
  }
  statistics.incrementStats(rho, uSqr);
}

template<typename T, typename DESCRIPTOR>
T PSMBGKdynamics<T,DESCRIPTOR>::getOmega() const
{
  return omega;
}

template<typename T, typename DESCRIPTOR>
void PSMBGKdynamics<T,DESCRIPTOR>::setOmega(T omega_)
{
  paramA = (1. / omega_ - 0.5);
  omega = omega_;

}

////////////////////// Class ForcedPSMBGKdynamics //////////////////////////

/** \param omega relaxation parameter, related to the dynamic viscosity
 *  \param momenta a Momenta object to know how to compute velocity momenta
 */
template<typename T, typename DESCRIPTOR>
ForcedPSMBGKdynamics<T,DESCRIPTOR>::ForcedPSMBGKdynamics (
  T omega_, Momenta<T,DESCRIPTOR>& momenta_, int mode_ )
  : ForcedBGKdynamics<T,DESCRIPTOR>(omega_, momenta_),
    omega(omega_), paramA(1. / omega_ - 0.5)
{
  mode = (Mode) mode_;
}


template<typename T, typename DESCRIPTOR>
void ForcedPSMBGKdynamics<T,DESCRIPTOR>::computeU (ConstCell<T,DESCRIPTOR>& cell, T u[DESCRIPTOR::d] ) const
{
  T rho;
  this->_momenta.computeRhoU(cell, rho, u);
  T epsilon = 1. - cell.template getField<descriptors::POROSITY>();
  // speed up paramB
  T paramB = (epsilon * paramA) / ((1. - epsilon) + paramA);
  // speed up paramC
  T paramC = (1. - paramB);
  for (int iVel=0; iVel<DESCRIPTOR::d; ++iVel) {
    u[iVel] = paramC * (u[iVel] + cell.template getFieldPointer<descriptors::FORCE>()[iVel] / (T)2.) +
              paramB * cell.template getFieldPointer<descriptors::VELOCITY_SOLID>()[iVel];
  }
}

template<typename T, typename DESCRIPTOR>
void ForcedPSMBGKdynamics<T,DESCRIPTOR>::computeRhoU (ConstCell<T,DESCRIPTOR>& cell, T& rho, T u[DESCRIPTOR::d] ) const
{
  this->_momenta.computeRhoU(cell, rho, u);
  T epsilon = 1. - cell.template getField<descriptors::POROSITY>();
  // speed up paramB
  T paramB = (epsilon * paramA) / ((1. - epsilon) + paramA);
  // speed up paramC
  T paramC = (1. - paramB);
  for (int iVel=0; iVel<DESCRIPTOR::d; ++iVel) {
    u[iVel] = paramC * (u[iVel] + cell.template getFieldPointer<descriptors::FORCE>()[iVel] / (T)2.) +
              paramB * cell.template getFieldPointer<descriptors::VELOCITY_SOLID>()[iVel];
  }
}

template<typename T, typename DESCRIPTOR>
void ForcedPSMBGKdynamics<T,DESCRIPTOR>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho, u[DESCRIPTOR::d], uSqr=0;
  T epsilon = 1. - cell.template getField<descriptors::POROSITY>();

  this->_momenta.computeRhoU(cell, rho, u);

  auto force = cell.template getFieldPointer<descriptors::FORCE>();
  for (int iVel=0; iVel<DESCRIPTOR::d; ++iVel) {
    u[iVel] += force[iVel] / (T)2.;
  }
  // velocity at the boundary
  auto u_s = cell.template getField<descriptors::VELOCITY_SOLID>();

  if (epsilon < 1e-5) {
    uSqr = lbHelpers<T,DESCRIPTOR>::bgkCollision(cell, rho, u, omega);
    lbHelpers<T,DESCRIPTOR>::addExternalForce(cell, u, omega, rho);
  } else {
    // speed up paramB
    T paramB = (epsilon * paramA) / ((1. - epsilon) + paramA);
    // speed up paramC
    T paramC = (1. - paramB);

    T omega_s[DESCRIPTOR::q];
    T cell_tmp[DESCRIPTOR::q];

    const T uSqr_s = util::normSqr<T,DESCRIPTOR::d>(u_s.data());

    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      cell_tmp[iPop] = cell[iPop];
      switch(mode){
        case M2: omega_s[iPop] = (lbDynamicsHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho, u_s.data(), uSqr_s ) - cell[iPop])
                         + (1 - omega) * (cell[iPop] - lbDynamicsHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho, u, uSqr )); break;
        case M3: omega_s[iPop] = (cell[descriptors::opposite<DESCRIPTOR>(iPop)] - lbDynamicsHelpers<T,DESCRIPTOR>::equilibrium(descriptors::opposite<DESCRIPTOR>(iPop), rho, u_s.data(), uSqr_s ))
               - (cell[iPop] - lbDynamicsHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho, u_s.data(), uSqr_s ));
      }
    }

    uSqr = lbHelpers<T,DESCRIPTOR>::bgkCollision(cell, rho, u, omega);
    lbHelpers<T,DESCRIPTOR>::addExternalForce(cell, u, omega, rho);

    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      cell[iPop] = cell_tmp[iPop] + paramC * (cell[iPop] - cell_tmp[iPop]);
      cell[iPop] += paramB * omega_s[iPop];
    }
    for (int iVel=0; iVel<DESCRIPTOR::d; ++iVel) {
      u[iVel] = paramC * u[iVel] + paramB * u_s[iVel];
    }
    uSqr = util::normSqr<T,DESCRIPTOR::d>(u);
  }
  statistics.incrementStats(rho, uSqr);
}

template<typename T, typename DESCRIPTOR>
T ForcedPSMBGKdynamics<T,DESCRIPTOR>::getOmega() const
{
  return omega;
}

template<typename T, typename DESCRIPTOR>
void ForcedPSMBGKdynamics<T,DESCRIPTOR>::setOmega(T omega_)
{
  paramA = (1. / omega_ - 0.5);
  omega = omega_;

}


} // olb

#endif
