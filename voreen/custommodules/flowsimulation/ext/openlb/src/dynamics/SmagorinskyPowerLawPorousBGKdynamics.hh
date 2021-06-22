/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012, 2015 Mathias J. Krause, Vojtech Cvrcekt, Davide Dapelo
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
 * Porous-particle BGK Dynamics with adjusted omega
 * and Smagorinsky turbulence model -- generic implementation.
 * Strain rate similar to "J.Boyd, J. Buick and S.Green: A second-order accurate lattice Boltzmann non-Newtonian flow model"
 * Power Law similar to "Huidan Yu, Sharath S. Girimaji, Li-Shi Luo - DNS and LES of decaying isotropic turbulence with and without frame rotation using lattice Boltzmann method"
 */
#ifndef SMAGORINSKY_POWER_LAW_POROUS_BGK_DYNAMICS_HH
#define SMAGORINSKY_POWER_LAW_POROUS_BGK_DYNAMICS_HH

#include "SmagorinskyPowerLawPorousBGKdynamics.h"
#include "SmagorinskyPorousParticleBGKdynamics.hh"
#include "SmagorinskyPowerLawBGKdynamics.hh"
#include "math.h"

namespace olb {

////////////////////// Class SmagorinskyPowerLawPorousParticleBGKdynamics //////////////////////////

/** \param vs2_ speed of sound
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 */
template<typename T, typename DESCRIPTOR>
SmagorinskyPowerLawPorousParticleBGKdynamics<T,DESCRIPTOR>::SmagorinskyPowerLawPorousParticleBGKdynamics (
  T omega_, Momenta<T,DESCRIPTOR>& momenta_, T m_, T n_ , T nuMin, T nuMax, T smagoConst_)
  : SmagorinskyPorousParticleBGKdynamics<T,DESCRIPTOR>(omega_, momenta_, smagoConst_),
    PowerLawDynamics<T,DESCRIPTOR>(m_, n_, nuMin, nuMax)
{ }

template<typename T, typename DESCRIPTOR>
void SmagorinskyPowerLawPorousParticleBGKdynamics<T,DESCRIPTOR>::collide (
  Cell<T,DESCRIPTOR>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho, u[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR >::n];
  this->_momenta.computeAllMomenta(cell, rho, u, pi);
//  T oldOmega = this->getOmega(); //compute with constant omega
  T oldOmega = cell.template getField<descriptors::OMEGA>(); //compute with dynamic omega
  T omegaPL = this->computeOmegaPL(cell, oldOmega, rho, pi);
  T newOmega = computeEffectiveOmega(cell, omegaPL);
  T uSqr = this->porousParticleBgkCollision(cell, rho, u, newOmega);
  statistics.incrementStats(rho, uSqr);
  cell.template setField<descriptors::OMEGA>(omegaPL); //compute with dynamic omega
}

template<typename T, typename DESCRIPTOR>
T SmagorinskyPowerLawPorousParticleBGKdynamics<T,DESCRIPTOR>::computeEffectiveOmega(Cell<T,DESCRIPTOR>& cell, T omega)
{
  T rho, u[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR >::n];
  this->_momenta.computeAllMomenta(cell, rho, u, pi);
  T PiNeqNormSqr = pi[0]*pi[0] + 2.0*pi[1]*pi[1] + pi[2]*pi[2];
  if (util::TensorVal<DESCRIPTOR >::n == 6) {
    PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2*pi[4]*pi[4] +pi[5]*pi[5];
  }
  T PiNeqNorm    = sqrt(PiNeqNormSqr);
  /// Molecular realaxation time
  T tau_mol = 1. /omega;
  /// Turbulent realaxation time
  T tau_turb = 0.5*(sqrt(tau_mol*tau_mol + this->getPreFactor()/rho*PiNeqNorm) - tau_mol);
  /// Effective realaxation time
  T tau_eff = tau_mol+tau_turb;
  T omega_new= 1./tau_eff;
  return omega_new;
}

}

#endif
