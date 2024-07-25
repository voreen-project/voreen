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

#ifndef LEGACY_POROUS_BGK_DYNAMICS_HH
#define LEGACY_POROUS_BGK_DYNAMICS_HH

#include "porousBGKdynamics.h"
#include "core/cell.h"
#include "dynamics.h"
#include "core/util.h"
#include "dynamics/lbm.h"
#include "math.h"

namespace olb {

//////////////////// Class SubgridParticleBGKdynamics ////////////////////

template<typename T, typename DESCRIPTOR, typename MOMENTA>
SubgridParticleBGKdynamics<T,DESCRIPTOR,MOMENTA>::SubgridParticleBGKdynamics (T omega_)
  : legacy::BGKdynamics<T,DESCRIPTOR,MOMENTA>(omega_), omega(omega_)
{
  _fieldTmp[0] = T();
  _fieldTmp[1] = T();
  _fieldTmp[2] = T();
  _fieldTmp[3] = T();
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
CellStatistic<T> SubgridParticleBGKdynamics<T,DESCRIPTOR,MOMENTA>::collide (
  Cell<T,DESCRIPTOR>& cell )
{
  T rho, u[DESCRIPTOR::d];
  MOMENTA().computeRhoU(cell, rho, u);
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
  T uSqr = lbm<DESCRIPTOR>::bgkCollision(cell, rho, u, omega);

  //statistics.incrementStats(rho, uSqr);
  cell.template setField<descriptors::POROSITY>(0);
  cell.template setField<descriptors::VELOCITY_NUMERATOR>(0);
  cell.template setField<descriptors::VELOCITY_DENOMINATOR>(0);
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
T SubgridParticleBGKdynamics<T,DESCRIPTOR,MOMENTA>::getOmega() const
{
  return omega;
}

//////////////////// Class PorousParticleDynamics ////////////////////

template<typename T, typename DESCRIPTOR, bool isStatic>
template <bool isStatic_>
std::enable_if_t<isStatic_>
PorousParticleDynamics<T,DESCRIPTOR,isStatic>::calculate(ConstCell<T,DESCRIPTOR>& cell, T* pVelocity)
{
  for (int i=0; i<DESCRIPTOR::d; i++)  {
    pVelocity[i] -= (1.-(cell.template getField<descriptors::POROSITY>())) * pVelocity[i];
  }
}

template<typename T, typename DESCRIPTOR, bool isStatic>
template <bool isStatic_>
std::enable_if_t<!isStatic_>
PorousParticleDynamics<T,DESCRIPTOR,isStatic>::calculate(Cell<T,DESCRIPTOR>& cell, T* pVelocity)
{
  for (int i=0; i<DESCRIPTOR::d; i++)  {
    pVelocity[i] += (1.-cell.template getField<descriptors::POROSITY>())
                    * (cell.template getFieldPointer<descriptors::VELOCITY_NUMERATOR>()[i]
                       / cell.template getField<descriptors::VELOCITY_DENOMINATOR>() - pVelocity[i]);
  }
  // reset external field for next timestep
  cell.template setField<descriptors::POROSITY>(1.);
  cell.template setField<descriptors::VELOCITY_DENOMINATOR>(0.);
  cell.template setField<descriptors::VELOCITY_NUMERATOR>({0.,0.,0.});
}

////////////// Temporary Helper Class PorousParticleBGK //////////////////

template<typename T, typename DESCRIPTOR, typename MOMENTA, bool isStatic>
PorousParticleBGK<T,DESCRIPTOR,MOMENTA,isStatic>::PorousParticleBGK(T omega_)
  : legacy::BGKdynamics<T,DESCRIPTOR,MOMENTA>(omega_)
{}

template<typename T, typename DESCRIPTOR, typename MOMENTA, bool isStatic>
T PorousParticleBGK<T,DESCRIPTOR,MOMENTA,isStatic>::porousParticleBgkCollision(Cell<T,DESCRIPTOR>& cell, T rho, T u[DESCRIPTOR::d], T omega)
{
  auto velNumerator   = cell.template getFieldPointer<descriptors::VELOCITY_NUMERATOR>();
  auto velDenominator = cell.template getFieldPointer<descriptors::VELOCITY_DENOMINATOR>();

#if defined(FEATURE_HLBM_SHANCHEN_FORCING)
  if (velDenominator[0] > std::numeric_limits<T>::epsilon()) {
    T  u_tmp[DESCRIPTOR::d] = { };
    for (int iDim=0; iDim<DESCRIPTOR::d; iDim++) {
      u_tmp[iDim] = u[iDim];
    }
    this->calculate(cell, u);
#ifdef FEATURE_HLBM_MLA
    const T tmp_uSqr = util::normSqr<T,DESCRIPTOR::d>(u);
    for (int iDim=0; iDim<DESCRIPTOR::d; iDim++) {
      for (int iPop=0; iPop<DESCRIPTOR::q; iPop++) {
        velNumerator[iDim] -= descriptors::c<DESCRIPTOR>(iPop,iDim) * omega
                              * (equilibrium<DESCRIPTOR>::secondOrder(iPop, rho, u, tmp_uSqr)
                                 - equilibrium<DESCRIPTOR>::secondOrder(iPop, rho, u_tmp, tmp_uSqr_2) );
      }
    }
#endif
  }
  T uSqr = lbm<DESCRIPTOR>::bgkCollision(cell, rho, u, omega);
#elif defined(FEATURE_HLBM_GUO_FORCING)
  T  force[DESCRIPTOR::d] = { };
#ifdef FEATURE_HLBM_MLA
  T  u_saved[DESCRIPTOR::d] = { };
  for (int iDim=0; iDim<DESCRIPTOR::d; iDim++) {
    u_saved[iDim] = u[iDim];
  }
#endif
  bool particle = velDenominator[0] > std::numeric_limits<T>::epsilon();
  if (particle) {
    T  uPlus[DESCRIPTOR::d] = { };
    for (int iDim=0; iDim<DESCRIPTOR::d; iDim++) {
      uPlus[iDim] = u[iDim];
    }
    this->calculate(cell, uPlus);
    for (int iDim=0; iDim<DESCRIPTOR::d; iDim++) {
      force[iDim] = uPlus[iDim]-u[iDim];
      u[iDim] += force[iDim] * 0.5;
    }
  }
  T uSqr = lbm<DESCRIPTOR>::bgkCollision(cell, rho, u, omega);
#ifdef FEATURE_HLBM_MLA
  for (int iDim=0; iDim<DESCRIPTOR::d; iDim++) {
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
      for (int iDim=0; iDim<DESCRIPTOR::d; iDim++) {
        velNumerator[iDim] += descriptors::c<DESCRIPTOR>(iPop,iDim) * ( omega
                              * (lbm<DESCRIPTOR>::equilibrium(iPop, rho, u_saved, tmp_uSqr_mla)
                                 - lbm<DESCRIPTOR>::equilibrium(iPop, rho, u, tmp_uSqr_mla_2) )
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
  T uSqr = lbm<DESCRIPTOR>::bgkCollision(cell, rho, u, omega);
  T uPlus[DESCRIPTOR::d] = { };
  T diff[DESCRIPTOR::q] = {};
  if (velDenominator[0] > std::numeric_limits<T>::epsilon()) {
    for (int iDim=0; iDim<DESCRIPTOR::d; iDim++) {
      uPlus[iDim] = u[iDim];
    }
    this->calculate(cell, uPlus);
    const T uPlusSqr = util::normSqr<T,DESCRIPTOR::d>(uPlus);
    for (int tmp_iPop=0; tmp_iPop<DESCRIPTOR::q; tmp_iPop++) {
      diff[tmp_iPop] += equilibrium<DESCRIPTOR>::secondOrder(tmp_iPop, rho, uPlus, uPlusSqr)
                        - equilibrium<DESCRIPTOR>::secondOrder(tmp_iPop, rho, u, uSqr);
      cell[tmp_iPop] += diff[tmp_iPop];
    }
#ifdef FEATURE_HLBM_MLA
    for (int iDim=0; iDim<DESCRIPTOR::d; iDim++) {
      velNumerator[iDim] = 0.;
      for (int tmp_iPop=0; tmp_iPop<DESCRIPTOR::q; tmp_iPop++) {
        velNumerator[iDim] -= descriptors::c<DESCRIPTOR>(tmp_iPop,iDim) * diff[tmp_iPop];
      }
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

template<typename T, typename DESCRIPTOR, typename MOMENTA, bool isStatic>
DBBParticleBGKdynamics<T,DESCRIPTOR,MOMENTA,isStatic>::DBBParticleBGKdynamics (T omega_)
  : BGKdynamics<T,DESCRIPTOR,MOMENTA>(omega_)
{
  this->getName() = "DBBParticleBGKdynamics";
}

template<typename T, typename DESCRIPTOR, typename MOMENTA, bool isStatic>
CellStatistic<T> DBBParticleBGKdynamics<T,DESCRIPTOR,MOMENTA,isStatic>::collide (
  Cell<T,DESCRIPTOR>& cell )
{
  T rho, u[DESCRIPTOR::d],eta[DESCRIPTOR::q],uPlus[DESCRIPTOR::d],tmp_cell[(DESCRIPTOR::q+1)/2];
  MOMENTA().computeRhoU(cell, rho, u);
  T uSqr = this->dbbParticleBgkCollision(cell, rho, u, eta, uPlus, tmp_cell, this->getOmega());
  //statistics.incrementStats(rho, uSqr);
}

template<typename T, typename DESCRIPTOR, typename MOMENTA, bool isStatic>
T DBBParticleBGKdynamics<T,DESCRIPTOR,MOMENTA,isStatic>::dbbParticleBgkCollision(Cell<T,DESCRIPTOR>& cell, T rho, T u[DESCRIPTOR::d], T eta[DESCRIPTOR::d], T uPlus[DESCRIPTOR::d],T tmp_cell[(DESCRIPTOR::q+1)/2], T omega)
{

  T tmpMomentumLoss[DESCRIPTOR::d] = { };

  auto velNumerator   = cell.template getFieldPointer<descriptors::VELOCITY_NUMERATOR>();
  auto zeta = cell.template getFieldPointer<descriptors::ZETA>();
  auto velDenominator = cell.template getFieldPointer<descriptors::VELOCITY_DENOMINATOR>();

  if (*(velDenominator)>1) {
    rho/=T(*(velDenominator));
  }
  for (int tmp_iPop=1; 2*tmp_iPop<DESCRIPTOR::q; tmp_iPop++) {
    eta[tmp_iPop]=6.*descriptors::t<T,DESCRIPTOR>(tmp_iPop)*rho*(descriptors::c<DESCRIPTOR>(tmp_iPop,0)*(*velNumerator)+descriptors::c<DESCRIPTOR>(tmp_iPop,1)*(*(velNumerator+1)));
    tmp_cell[tmp_iPop]=(*(zeta+tmp_iPop))*(-cell[tmp_iPop]+cell[descriptors::opposite<DESCRIPTOR>(tmp_iPop)]+eta[tmp_iPop]);
    cell[tmp_iPop]+=tmp_cell[tmp_iPop]/(1.+2.*(*(zeta+tmp_iPop)));
    cell[descriptors::opposite<DESCRIPTOR>(tmp_iPop)]-=tmp_cell[tmp_iPop]/(1.+2.*(*(zeta+tmp_iPop)));
    *(zeta+tmp_iPop) = 0.;
    *(zeta+descriptors::opposite<DESCRIPTOR>(tmp_iPop)) = 0.;
  }

  cell.template setField<descriptors::POROSITY>(1.);
  *velDenominator=0.;

  MOMENTA().computeRhoU(cell, rho, uPlus);

  T uPlusSqr = lbm<DESCRIPTOR>::bgkCollision(cell, rho, uPlus, omega);

  T uSqr = util::normSqr<T,DESCRIPTOR::d>(u);

  T diff[DESCRIPTOR::q] = {};
  for (int tmp_iPop=0; tmp_iPop<DESCRIPTOR::q; tmp_iPop++) {
    diff[tmp_iPop] += lbm<DESCRIPTOR>::equilibrium(tmp_iPop, rho, uPlus, uPlusSqr)
                      - lbm<DESCRIPTOR>::equilibrium(tmp_iPop, rho, u, uSqr);

    for (int iDim=0; iDim<DESCRIPTOR::d; iDim++) {
      tmpMomentumLoss[iDim] -= descriptors::c<DESCRIPTOR>(tmp_iPop,iDim) * diff[tmp_iPop];
    }
  }

  for (int i_dim=0; i_dim<DESCRIPTOR::d; i_dim++) {
    *(velNumerator+i_dim) = tmpMomentumLoss[i_dim];
  }

  return uSqr;
}




//////////////////// Class KrauseHBGKdynamics ////////////////////

template<typename T, typename DESCRIPTOR, typename MOMENTA>
KrauseHBGKdynamics<T,DESCRIPTOR,MOMENTA>::KrauseHBGKdynamics (T omega_, T smagoConst_, T dx_, T dt_ )
  : BGKdynamics<T,DESCRIPTOR,MOMENTA>(omega_), omega(omega_), smagoConst(smagoConst_),
    preFactor(computePreFactor(omega_,smagoConst_) )
{
  _fieldTmp[0] = T(1);
  _fieldTmp[1] = T();
  _fieldTmp[2] = T();
  _fieldTmp[3] = T();
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
CellStatistic<T> KrauseHBGKdynamics<T,DESCRIPTOR,MOMENTA>::collide (
  Cell<T,DESCRIPTOR>& cell )
{
  T rho, u[DESCRIPTOR::d];
  T newOmega[DESCRIPTOR::q];
  MOMENTA().computeRhoU(cell, rho, u);
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
  T uSqr = lbm<DESCRIPTOR>::bgkCollision(cell, rho, u, newOmega);
  //statistics.incrementStats(rho, uSqr);

  cell.template setField<descriptors::POROSITY>(_fieldTmp[0]);
  cell.template setField<descriptors::VELOCITY_NUMERATOR>({_fieldTmp[1], _fieldTmp[2]});
  cell.template setField<descriptors::VELOCITY_DENOMINATOR>(_fieldTmp[3]);
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
T KrauseHBGKdynamics<T,DESCRIPTOR,MOMENTA>::getOmega() const
{
  return omega;
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
T KrauseHBGKdynamics<T,DESCRIPTOR,MOMENTA>::computePreFactor(T omega, T smagoConst)
{
  return (T)smagoConst*smagoConst*descriptors::invCs2<T,DESCRIPTOR>()*descriptors::invCs2<T,DESCRIPTOR>()*2*util::sqrt(2);
}


template<typename T, typename DESCRIPTOR, typename MOMENTA>
void KrauseHBGKdynamics<T,DESCRIPTOR,MOMENTA>::computeOmega(T omega0, Cell<T,DESCRIPTOR>& cell, T preFactor, T rho,
    T u[DESCRIPTOR::d], T newOmega[DESCRIPTOR::q])
{
  T uSqr = u[0]*u[0];
  for (int iDim=0; iDim<DESCRIPTOR::d; iDim++) {
    uSqr += u[iDim]*u[iDim];
  }
  /// Molecular realaxation time
  T tau_mol = 1./omega0;

  for (int iPop=0; iPop<DESCRIPTOR::q; iPop++) {
    T fNeq = util::fabs(cell[iPop] - lbm<DESCRIPTOR>::equilibrium(iPop, rho, u, uSqr));
    /// Turbulent realaxation time
    T tau_turb = 0.5*(util::sqrt(tau_mol*tau_mol+(preFactor*fNeq))-tau_mol);
    /// Effective realaxation time
    tau_eff = tau_mol + tau_turb;
    newOmega[iPop] = 1./tau_eff;
  }
}


//////////////////// Class SmallParticleBGKdynamics ////////////////////

template<typename T, typename DESCRIPTOR, typename MOMENTA>
SmallParticleBGKdynamics<T,DESCRIPTOR,MOMENTA>::SmallParticleBGKdynamics (T omega_)
  : legacy::BGKdynamics<T,DESCRIPTOR,MOMENTA>(omega_),
    omega(omega_)
{ }

template<typename T, typename DESCRIPTOR, typename MOMENTA>
CellStatistic<T> SmallParticleBGKdynamics<T,DESCRIPTOR,MOMENTA>::collide (
  Cell<T,DESCRIPTOR>& cell )
{
  T rho, u[DESCRIPTOR::d];
  MOMENTA().computeRhoU(cell, rho, u);
  T porosity = cell.template getField<descriptors::POROSITY>();
  auto localVelocity = cell.template getFieldPointer<descriptors::LOCAL_DRAG>();

  //cout << porosity[0]  << " " <<   localVelocity[0]<< " " <<   localVelocity[1]<< " " <<   localVelocity[2]<<std::endl;
  for (int i=0; i<DESCRIPTOR::d; i++)  {
    u[i] *= porosity;
    u[i] += localVelocity[i];
  }
  T uSqr = lbm<DESCRIPTOR>::bgkCollision(cell, rho, u, omega);
  //statistics.incrementStats(rho, uSqr);
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
T SmallParticleBGKdynamics<T,DESCRIPTOR,MOMENTA>::getOmega() const
{
  return omega;
}


} // olb

#endif
