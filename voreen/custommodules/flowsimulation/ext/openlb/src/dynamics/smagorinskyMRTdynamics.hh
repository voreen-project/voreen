/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Patrick Nathen, Mathias J. Krause
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
 * MRT Dynamics with adjusted omega -- generic implementation.
 */
#ifndef SMAGORINSKY_MRT_DYNAMICS_HH
#define SMAGORINSKY_MRT_DYNAMICS_HH

#include <algorithm>
#include <limits>
#include "smagorinskyMRTdynamics.h"
#include "mrtDynamics.h"
#include "mrtHelpers.h"
#include "core/cell.h"
#include "core/util.h"
#include "math.h"


using namespace std;
namespace olb {

////////////////////// Class SmagorinskyMRTdynamics //////////////////////////

/** \param vs2_ speed of sound
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 */



template<typename T, template<typename U> class Lattice>
SmagorinskyMRTdynamics<T,Lattice>::SmagorinskyMRTdynamics (
  T omega_, Momenta<T,Lattice>& momenta_, T smagoConst_, T dx_, T dt_ )
  : MRTdynamics<T,Lattice>(omega_, momenta_/*, smagoConst_*/),
    smagoConst(smagoConst_),
    preFactor(computePreFactor(omega_,smagoConst_) )
{

  T rtSGS[Lattice<T>::q]; // relaxation times vector for SGS approach.
  for (int iPop  = 0; iPop < Lattice<T>::q; ++iPop) {
    rtSGS[iPop] = Lattice<T>::S[iPop];
  }
  for (int iPop  = 0; iPop < Lattice<T>::shearIndexes; ++iPop) {
    rtSGS[Lattice<T>::shearViscIndexes[iPop]] = omega;
  }
  for (int iPop = 0; iPop < Lattice<T>::q; ++iPop) {
    for (int jPop = 0; jPop < Lattice<T>::q; ++jPop) {
      invM_S_SGS[iPop][jPop] = T();
      for (int kPop = 0; kPop < Lattice<T>::q; ++kPop) {
        if (kPop == jPop) {
          invM_S_SGS[iPop][jPop] += Lattice<T>::invM[iPop][kPop] * rtSGS[kPop];
        }
      }
    }
  }
}


template<typename T, template<typename U> class Lattice>
void SmagorinskyMRTdynamics<T,Lattice>::collide(
  Cell<T,Lattice>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho, u[Lattice<T>::d], pi[util::TensorVal<Lattice<T> >::n];
  this->_momenta.computeAllMomenta(cell, rho, u, pi);
  T newOmega = computeOmega(this->getOmega(), preFactor, rho, pi);
  for (int iPop = 0; iPop < Lattice<T>::q; ++iPop) {
    for (int jPop = 0; jPop < Lattice<T>::shearIndexes; ++jPop) {
      invM_S_SGS[iPop][Lattice<T>::shearViscIndexes[jPop]] = Lattice<T>::invM[iPop][Lattice<T>::shearViscIndexes[jPop]] * newOmega;
    }
  }

  T uSqr = mrtHelpers<T,Lattice>::mrtSGSCollision(cell, rho, u, newOmega, invM_S_SGS);
  statistics.incrementStats(rho, uSqr);
}


template<typename T, template<typename U> class Lattice>
void SmagorinskyMRTdynamics<T,Lattice>::staticCollide(
  Cell<T,Lattice>& cell,
  const T u[Lattice<T>::d],
  LatticeStatistics<T>& statistics )
{
  T rho, uTemp[Lattice<T>::d], pi[util::TensorVal<Lattice<T> >::n];
  this->_momenta.computeAllMomenta(cell, rho, uTemp, pi);
  T newOmega = computeOmega(this->getOmega(), preFactor, rho, pi);
  for (int iPop = 0; iPop < Lattice<T>::q; ++iPop) {
    for (int jPop = 0; jPop < Lattice<T>::shearIndexes; ++jPop) {
      invM_S_SGS[iPop][Lattice<T>::shearViscIndexes[jPop]] = Lattice<T>::invM[iPop][Lattice<T>::shearViscIndexes[jPop]] * newOmega;
    }
  }

  T uSqr = mrtHelpers<T,Lattice>::mrtSGSCollision(cell, rho, u, newOmega, invM_S_SGS);
  statistics.incrementStats(rho, uSqr);
}

template<typename T, template<typename U> class Lattice>
void SmagorinskyMRTdynamics<T,Lattice>::setOmega(T omega)
{
  this->setOmega(omega);
  preFactor = computePreFactor(omega, smagoConst);
}

template<typename T, template<typename U> class Lattice>
T SmagorinskyMRTdynamics<T,Lattice>::getSmagorinskyOmega(Cell<T,Lattice>& cell )
{
  T rho, uTemp[Lattice<T>::d], pi[util::TensorVal<Lattice<T> >::n];
  this->_momenta.computeAllMomenta(cell, rho, uTemp, pi);
  T newOmega = computeOmega(this->getOmega(), preFactor, rho, pi);
  return newOmega;
}

template<typename T, template<typename U> class Lattice>
T SmagorinskyMRTdynamics<T,Lattice>::computePreFactor(T omega, T smagoConst)
{
  return (T)smagoConst*smagoConst*Lattice<T>::invCs2*Lattice<T>::invCs2*2*sqrt(2);
}

template<typename T, template<typename U> class Lattice>
T SmagorinskyMRTdynamics<T,Lattice>::computeOmega(T omega0, T preFactor, T rho, T pi[util::TensorVal<Lattice<T> >::n] )
{

  T PiNeqNormSqr = pi[0]*pi[0] + 2.0*pi[1]*pi[1] + pi[2]*pi[2];
  if (util::TensorVal<Lattice<T> >::n == 6) {
    PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2*pi[4]*pi[4] +pi[5]*pi[5];
  }
  T PiNeqNorm    = sqrt(PiNeqNormSqr);
  /// Molecular realaxation time
  T tau_mol = 1. /omega0;
  /// Turbulent realaxation time
  T tau_turb = 0.5*(sqrt(tau_mol*tau_mol + preFactor/rho*PiNeqNorm) - tau_mol);
  /// Effective realaxation time
  tau_eff = tau_mol+tau_turb;
  T omega_new= 1./tau_eff;
  return omega_new;
}



template<typename T, template<typename U> class Lattice>
void SmagorinskyForcedMRTdynamics<T,Lattice>::collide(
  Cell<T,Lattice>& cell,
  LatticeStatistics<T>& statistics )
{
  T rho, u[Lattice<T>::d], pi[util::TensorVal<Lattice<T> >::n];
  this->_momenta.computeAllMomenta(cell, rho, u, pi);
  T newOmega = computeOmega(this->getOmega(), this->preFactor, rho, pi);
  for (int iPop = 0; iPop < Lattice<T>::q; ++iPop) {
    for (int jPop = 0; jPop < Lattice<T>::shearIndexes; ++jPop) {
      this->invM_S_SGS[iPop][Lattice<T>::shearViscIndexes[jPop]] = Lattice<T>::invM[iPop][Lattice<T>::shearViscIndexes[jPop]] * newOmega;
    }
  }

  T uSqr = mrtHelpers<T,Lattice>::mrtSGSCollision(cell, rho, u, newOmega, this->invM_S_SGS);

  mrtHelpers<T,Lattice>::addExternalForce(cell, rho, u, this->invM_S_SGS);

  statistics.incrementStats(rho, uSqr);
}

}

#endif
