/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007, 2017 Orestis Malaspinas, Jonas Latt, Mathias J. Krause
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
 * A collection of entropic dynamics classes (e.g. Entropic, ForcedEntropic, Entropic, ForcedEntropic) with which a Cell object
 * can be instantiated -- generic implementation.
 */
#ifndef ENTROPIC_LB_DYNAMICS_HH
#define ENTROPIC_LB_DYNAMICS_HH

#include <algorithm>
#include <limits>
#include "lbm.h"
#include "entropicDynamics.h"
#include "entropicLbHelpers.h"

namespace olb {

//==============================================================================//
/////////////////////////// Class EntropicEqDynamics ///////////////////////////////
//==============================================================================//
/** \param omega_ relaxation parameter, related to the dynamic viscosity
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 */
template<typename T, typename DESCRIPTOR, typename MOMENTA>
EntropicEqDynamics<T,DESCRIPTOR,MOMENTA>::EntropicEqDynamics (T omega_)
  : legacy::BasicDynamics<T,DESCRIPTOR,MOMENTA>(),
    omega(omega_)
{
  this->getName() = "EntropicEqDynamics";
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
T EntropicEqDynamics<T,DESCRIPTOR,MOMENTA>::computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr) const
{
  return entropicLbHelpers<T,DESCRIPTOR>::equilibrium(iPop,rho,u);
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
CellStatistic<T> EntropicEqDynamics<T,DESCRIPTOR,MOMENTA>::collide (
  Cell<T,DESCRIPTOR>& cell)
{
  typedef DESCRIPTOR L;
  typedef entropicLbHelpers<T,DESCRIPTOR> eLbH;

  T rho, u[DESCRIPTOR::d];
  MOMENTA().computeRhoU(cell, rho, u);
  T uSqr = util::normSqr<T,L::d>(u);

  for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    cell[iPop] += omega * (eLbH::equilibrium(iPop,rho,u) - cell[iPop]);
  }

  //statistics.incrementStats(rho, uSqr);
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
T EntropicEqDynamics<T,DESCRIPTOR,MOMENTA>::getOmega() const
{
  return omega;
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
void EntropicEqDynamics<T,DESCRIPTOR,MOMENTA>::setOmega(T omega_)
{
  omega = omega_;
}


//====================================================================//
//////////////////// Class ForcedEntropicEqDynamics //////////////////////
//====================================================================//

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template<typename T, typename DESCRIPTOR, typename MOMENTA>
ForcedEntropicEqDynamics<T,DESCRIPTOR,MOMENTA>::ForcedEntropicEqDynamics (T omega_)
  : legacy::BasicDynamics<T,DESCRIPTOR,MOMENTA>(),
    omega(omega_)
{
  this->getName() = "ForcedEntropicEqDynamics";
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
T ForcedEntropicEqDynamics<T,DESCRIPTOR,MOMENTA>::computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr) const
{
  return entropicLbHelpers<T,DESCRIPTOR>::equilibrium(iPop,rho,u);
}


template<typename T, typename DESCRIPTOR, typename MOMENTA>
CellStatistic<T> ForcedEntropicEqDynamics<T,DESCRIPTOR,MOMENTA>::collide (
  Cell<T,DESCRIPTOR>& cell)
{
  typedef DESCRIPTOR L;
  typedef entropicLbHelpers<T,DESCRIPTOR> eLbH;

  T rho, u[DESCRIPTOR::d];
  MOMENTA().computeRhoU(cell, rho, u);

  T* force = cell.template getFieldPointer<descriptors::FORCE>();
  for (int iDim=0; iDim<DESCRIPTOR::d; ++iDim) {
    u[iDim] += force[iDim] / (T)2.;
  }
  T uSqr = util::normSqr<T,L::d>(u);

  for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    cell[iPop] += omega * (eLbH::equilibrium(iPop,rho,u) - cell[iPop]);
  }

  lbm<DESCRIPTOR>::addExternalForce(cell, u, omega);

  //statistics.incrementStats(rho, uSqr);
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
T ForcedEntropicEqDynamics<T,DESCRIPTOR,MOMENTA>::getOmega() const
{
  return omega;
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
void ForcedEntropicEqDynamics<T,DESCRIPTOR,MOMENTA>::setOmega(T omega_)
{
  omega = omega_;
}


//==============================================================================//
/////////////////////////// Class EntropicDynamics ///////////////////////////////
//==============================================================================//
/** \param omega_ relaxation parameter, related to the dynamic viscosity
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 */
template<typename T, typename DESCRIPTOR, typename MOMENTA>
EntropicDynamics<T,DESCRIPTOR,MOMENTA>::EntropicDynamics (T omega_)
  : legacy::BasicDynamics<T,DESCRIPTOR,MOMENTA>(),
    omega(omega_)
{
  this->getName() = "EntropicDynamics";
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
T EntropicDynamics<T,DESCRIPTOR,MOMENTA>::computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr) const
{
  return entropicLbHelpers<T,DESCRIPTOR>::equilibrium(iPop,rho,u);
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
CellStatistic<T> EntropicDynamics<T,DESCRIPTOR,MOMENTA>::collide (
  Cell<T,DESCRIPTOR>& cell)
{
  typedef DESCRIPTOR L;
  typedef entropicLbHelpers<T,DESCRIPTOR> eLbH;

  T rho, u[DESCRIPTOR::d];
  MOMENTA().computeRhoU(cell, rho, u);
  T uSqr = util::normSqr<T,L::d>(u);

  T f[L::q], fEq[L::q], fNeq[L::q];
  for (int iPop = 0; iPop < L::q; ++iPop) {
    fEq[iPop]  = eLbH::equilibrium(iPop,rho,u);
    fNeq[iPop] = cell[iPop] - fEq[iPop];
    f[iPop]    = cell[iPop] + descriptors::t<T,L>(iPop);
    fEq[iPop] += descriptors::t<T,L>(iPop);
  }
  //==============================================================================//
  //============= Evaluation of alpha using a Newton Raphson algorithm ===========//
  //==============================================================================//

  T alpha = 2.0;
  bool converged = getAlpha(alpha,f,fNeq);
  if (!converged) {
    std::cout << "Newton-Raphson failed to converge.\n";
    exit(1);
  }

  OLB_ASSERT(converged,"Entropy growth failed to converge!");

  T omegaTot = omega / 2.0 * alpha;
  for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    cell[iPop] *= (T)1-omegaTot;
    cell[iPop] += omegaTot * (fEq[iPop]-descriptors::t<T,L>(iPop));
  }

  //statistics.incrementStats(rho, uSqr);
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
T EntropicDynamics<T,DESCRIPTOR,MOMENTA>::getOmega() const
{
  return omega;
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
void EntropicDynamics<T,DESCRIPTOR,MOMENTA>::setOmega(T omega_)
{
  omega = omega_;
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
T EntropicDynamics<T,DESCRIPTOR,MOMENTA>::computeEntropy(const T f[])
{
  typedef DESCRIPTOR L;
  T entropy = T();
  for (int iPop = 0; iPop < L::q; ++iPop) {
    OLB_ASSERT(f[iPop] > T(), "f[iPop] <= 0");
    entropy += f[iPop]*util::log(f[iPop]/descriptors::t<T,L>(iPop));
  }

  return entropy;
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
T EntropicDynamics<T,DESCRIPTOR,MOMENTA>::computeEntropyGrowth(const T f[], const T fNeq[], const T &alpha)
{
  typedef DESCRIPTOR L;

  T fAlphaFneq[L::q];
  for (int iPop = 0; iPop < L::q; ++iPop) {
    fAlphaFneq[iPop] = f[iPop] - alpha*fNeq[iPop];
  }

  return computeEntropy(f) - computeEntropy(fAlphaFneq);
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
T EntropicDynamics<T,DESCRIPTOR,MOMENTA>::computeEntropyGrowthDerivative(const T f[], const T fNeq[], const T &alpha)
{
  typedef DESCRIPTOR L;

  T entropyGrowthDerivative = T();
  for (int iPop = 0; iPop < L::q; ++iPop) {
    T tmp = f[iPop] - alpha*fNeq[iPop];
    OLB_ASSERT(tmp > T(), "f[iPop] - alpha*fNeq[iPop] <= 0");
    entropyGrowthDerivative += fNeq[iPop]*(util::log(tmp/descriptors::t<T,L>(iPop)));
  }

  return entropyGrowthDerivative;
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
bool EntropicDynamics<T,DESCRIPTOR,MOMENTA>::getAlpha(T &alpha, const T f[], const T fNeq[])
{
  const T epsilon = std::numeric_limits<T>::epsilon();

  T alphaGuess = T();
  const T var = 100.0;
  const T errorMax = epsilon*var;
  T error = 1.0;
  int count = 0;
  for (count = 0; count < 10000; ++count) {
    T entGrowth = computeEntropyGrowth(f,fNeq,alpha);
    T entGrowthDerivative = computeEntropyGrowthDerivative(f,fNeq,alpha);
    if ((error < errorMax) || (util::fabs(entGrowth) < var*epsilon)) {
      return true;
    }
    alphaGuess = alpha - entGrowth /
                 entGrowthDerivative;
    error = util::fabs(alpha-alphaGuess);
    alpha = alphaGuess;
  }
  return false;
}

//====================================================================//
//////////////////// Class ForcedEntropicDynamics //////////////////////
//====================================================================//

/** \param omega_ relaxation parameter, related to the dynamic viscosity
 */
template<typename T, typename DESCRIPTOR, typename MOMENTA>
ForcedEntropicDynamics<T,DESCRIPTOR,MOMENTA>::ForcedEntropicDynamics (T omega_)
  : legacy::BasicDynamics<T,DESCRIPTOR,MOMENTA>(),
    omega(omega_)
{
  this->getName() = "ForcedEntropicDynamics";
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
T ForcedEntropicDynamics<T,DESCRIPTOR,MOMENTA>::computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr) const
{
  return entropicLbHelpers<T,DESCRIPTOR>::equilibrium(iPop,rho,u);
}


template<typename T, typename DESCRIPTOR, typename MOMENTA>
CellStatistic<T> ForcedEntropicDynamics<T,DESCRIPTOR,MOMENTA>::collide (
  Cell<T,DESCRIPTOR>& cell)
{
  typedef DESCRIPTOR L;
  typedef entropicLbHelpers<T,DESCRIPTOR> eLbH;

  T rho, u[DESCRIPTOR::d];
  MOMENTA().computeRhoU(cell, rho, u);
  T uSqr = util::normSqr<T,L::d>(u);

  T f[L::q], fEq[L::q], fNeq[L::q];
  for (int iPop = 0; iPop < L::q; ++iPop) {
    fEq[iPop]  = eLbH::equilibrium(iPop,rho,u);
    fNeq[iPop] = cell[iPop] - fEq[iPop];
    f[iPop]    = cell[iPop] + descriptors::t<T,L>(iPop);
    fEq[iPop] += descriptors::t<T,L>(iPop);
  }
  //==============================================================================//
  //============= Evaluation of alpha using a Newton Raphson algorithm ===========//
  //==============================================================================//

  T alpha = 2.0;
  bool converged = getAlpha(alpha,f,fNeq);
  if (!converged) {
    std::cout << "Newton-Raphson failed to converge.\n";
    exit(1);
  }

  OLB_ASSERT(converged,"Entropy growth failed to converge!");

  T* force = cell.template getFieldPointer<descriptors::FORCE>();
  for (int iDim=0; iDim<DESCRIPTOR::d; ++iDim) {
    u[iDim] += force[iDim] / (T)2.;
  }
  uSqr = util::normSqr<T,L::d>(u);
  T omegaTot = omega / 2.0 * alpha;
  for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
    cell[iPop] *= (T)1-omegaTot;
    cell[iPop] += omegaTot * eLbH::equilibrium(iPop,rho,u);
  }
  lbm<DESCRIPTOR>::addExternalForce(cell, u, omegaTot);

  //statistics.incrementStats(rho, uSqr);
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
T ForcedEntropicDynamics<T,DESCRIPTOR,MOMENTA>::getOmega() const
{
  return omega;
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
void ForcedEntropicDynamics<T,DESCRIPTOR,MOMENTA>::setOmega(T omega_)
{
  omega = omega_;
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
T ForcedEntropicDynamics<T,DESCRIPTOR,MOMENTA>::computeEntropy(const T f[])
{
  typedef DESCRIPTOR L;
  T entropy = T();
  for (int iPop = 0; iPop < L::q; ++iPop) {
    OLB_ASSERT(f[iPop] > T(), "f[iPop] <= 0");
    entropy += f[iPop]*util::log(f[iPop]/descriptors::t<T,L>(iPop));
  }

  return entropy;
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
T ForcedEntropicDynamics<T,DESCRIPTOR,MOMENTA>::computeEntropyGrowth(const T f[], const T fNeq[], const T &alpha)
{
  typedef DESCRIPTOR L;

  T fAlphaFneq[L::q];
  for (int iPop = 0; iPop < L::q; ++iPop) {
    fAlphaFneq[iPop] = f[iPop] - alpha*fNeq[iPop];
  }

  return computeEntropy(f) - computeEntropy(fAlphaFneq);
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
T ForcedEntropicDynamics<T,DESCRIPTOR,MOMENTA>::computeEntropyGrowthDerivative(const T f[], const T fNeq[], const T &alpha)
{
  typedef DESCRIPTOR L;

  T entropyGrowthDerivative = T();
  for (int iPop = 0; iPop < L::q; ++iPop) {
    T tmp = f[iPop] - alpha*fNeq[iPop];
    OLB_ASSERT(tmp > T(), "f[iPop] - alpha*fNeq[iPop] <= 0");
    entropyGrowthDerivative += fNeq[iPop]*util::log(tmp/descriptors::t<T,L>(iPop));
  }

  return entropyGrowthDerivative;
}

template<typename T, typename DESCRIPTOR, typename MOMENTA>
bool ForcedEntropicDynamics<T,DESCRIPTOR,MOMENTA>::getAlpha(T &alpha, const T f[], const T fNeq[])
{
  const T epsilon = std::numeric_limits<T>::epsilon();

  T alphaGuess = T();
  const T var = 100.0;
  const T errorMax = epsilon*var;
  T error = 1.0;
  int count = 0;
  for (count = 0; count < 10000; ++count) {
    T entGrowth = computeEntropyGrowth(f,fNeq,alpha);
    T entGrowthDerivative = computeEntropyGrowthDerivative(f,fNeq,alpha);
    if ((error < errorMax) || (util::fabs(entGrowth) < var*epsilon)) {
      return true;
    }
    alphaGuess = alpha - entGrowth /
                 entGrowthDerivative;
    error = util::fabs(alpha-alphaGuess);
    alpha = alphaGuess;
  }
  return false;
}

}

#endif

