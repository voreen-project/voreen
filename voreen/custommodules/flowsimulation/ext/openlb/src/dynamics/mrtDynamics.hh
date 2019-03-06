/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Jonas Latt
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
#ifndef MRT_DYNAMICS_HH
#define MRT_DYNAMICS_HH

#include <algorithm>
#include <limits>
#include "mrtHelpers.h"

namespace olb {

//==============================================================================//
/////////////////////////// Class MRTdynamics ///////////////////////////////
//==============================================================================//
/** \param omega_ relaxation parameter, related to the dynamic viscosity
 *  \param momenta_ a Momenta object to know how to compute velocity momenta
 *  \param lambda_ will be used as an
 */

// Original implementation based on:
// D'Humieres et al., "Multiple-relaxation-time lattice Boltzmann models in three dimensions",
// Phil: Trans. R. soc. Lond. A (2002) 360, 437-451
// and
// Yu et al,, "LES of turbulent square jet flow using an MRT lattice Boltzmann model",
// Computers & Fluids 35 (2006), 957-965
template<typename T, template<typename U> class Lattice>
MRTdynamics<T,Lattice>::MRTdynamics (
  T omega_, Momenta<T,Lattice>& momenta_ )
  : BasicDynamics<T,Lattice>(momenta_), omega(omega_), lambda(omega_)
{
  T rt[Lattice<T>::q]; // relaxation times vector.
  for (int iPop  = 0; iPop < Lattice<T>::q; ++iPop) {
    rt[iPop] = Lattice<T>::S[iPop];
  }
  for (int iPop  = 0; iPop < Lattice<T>::shearIndexes; ++iPop) {
    rt[Lattice<T>::shearViscIndexes[iPop]] = omega;
  }
  for (int iPop = 0; iPop < Lattice<T>::q; ++iPop) {
    for (int jPop = 0; jPop < Lattice<T>::q; ++jPop) {
      invM_S[iPop][jPop] = T();
      for (int kPop = 0; kPop < Lattice<T>::q; ++kPop) {
        if (kPop == jPop) {
          invM_S[iPop][jPop] += Lattice<T>::invM[iPop][kPop] *
                                rt[kPop];
        }
      }
    }
  }

}

template<typename T, template<typename U> class Lattice>
T MRTdynamics<T,Lattice>::computeEquilibrium(int iPop, T rho, const T u[Lattice<T>::d], T uSqr) const
{
  return lbHelpers<T,Lattice>::equilibrium(iPop, rho, u, uSqr);
}

template<typename T, template<typename U> class Lattice>
void MRTdynamics<T,Lattice>::collide (
  Cell<T,Lattice>& cell,
  LatticeStatistics<T>& statistics )
{
  typedef Lattice<T> L;
  typedef mrtHelpers<T,Lattice> mrtH;

  T rho, u[L::d];
  this->_momenta.computeRhoU(cell, rho, u);

  T uSqr = mrtH::mrtCollision(cell,rho,u,invM_S);

  statistics.incrementStats(rho, uSqr);
}

template<typename T, template<typename U> class Lattice>
void MRTdynamics<T,Lattice>::staticCollide (
  Cell<T,Lattice>& cell,
  const T u[Lattice<T>::d],
  LatticeStatistics<T>& statistics )
{
  typedef mrtHelpers<T,Lattice> mrtH;

  T rho = T(1);
  T uSqr = mrtH::mrtCollision(cell, rho, u, invM_S);
  statistics.incrementStats(rho, uSqr);
}

template<typename T, template<typename U> class Lattice>
T MRTdynamics<T,Lattice>::getOmega() const
{
  return omega;
}

template<typename T, template<typename U> class Lattice>
void MRTdynamics<T,Lattice>::setOmega(T omega_)
{
  omega = omega_;
}

template<typename T, template<typename U> class Lattice>
T MRTdynamics<T,Lattice>::getLambda() const
{
  return lambda;
}

template<typename T, template<typename U> class Lattice>
void MRTdynamics<T,Lattice>::setLambda(T lambda_)
{
  lambda = lambda_;
}




template<typename T, template<typename U> class Lattice>
MRTdynamics2<T,Lattice>::MRTdynamics2 (
  T omega_, Momenta<T,Lattice>& momenta_ )
  : MRTdynamics<T,Lattice>(omega_, momenta_)
{
  T rt[Lattice<T>::q]; // relaxation times vector.
  for (int iPop  = 0; iPop < Lattice<T>::q; ++iPop) {
    rt[iPop] = Lattice<T>::S_2[iPop];
  }
  for (int iPop  = 0; iPop < Lattice<T>::shearIndexes; ++iPop) {
    rt[Lattice<T>::shearViscIndexes[iPop]] = omega;
  }
  for (int iPop = 0; iPop < Lattice<T>::q; ++iPop) {
    for (int jPop = 0; jPop < Lattice<T>::q; ++jPop) {
      invM_S_2[iPop][jPop] = T();
      for (int kPop = 0; kPop < Lattice<T>::q; ++kPop) {
        if (kPop == jPop) {
          invM_S_2[iPop][jPop] += Lattice<T>::invM[iPop][kPop] *
                                  rt[kPop];
        }
      }
    }
  }
}


// Stabalized MRT scheme with uniform relaxation times
template<typename T, template<typename U> class Lattice>
void MRTdynamics2<T,Lattice>::collide (
  Cell<T,Lattice>& cell,
  LatticeStatistics<T>& statistics )
{
  typedef Lattice<T> L;
  typedef mrtHelpers<T,Lattice> mrtH;

  T rho, u[L::d];
  this->_momenta.computeRhoU(cell, rho, u);

  T uSqr = mrtH::mrtCollision(cell,rho,u,invM_S_2);

  statistics.incrementStats(rho, uSqr);
}



template<typename T, template<typename U> class Lattice>
ForcedMRTdynamics<T,Lattice>::ForcedMRTdynamics (
  T omega_, Momenta<T,Lattice>& momenta_ )
  : MRTdynamics<T,Lattice>(omega_, momenta_)
{
}

template<typename T, template<typename U> class Lattice>
void ForcedMRTdynamics<T,Lattice>::collide (
  Cell<T,Lattice>& cell,
  LatticeStatistics<T>& statistics )
{
  typedef Lattice<T> L;
  typedef mrtHelpers<T,Lattice> mrtH;

  T rho, u[L::d];
  this->_momenta.computeRhoU(cell, rho, u);

  T uSqr = mrtH::mrtCollision(cell,rho,u,this->invM_S);
  mrtH::addExternalForce(cell, rho, u, this->invM_S);

  statistics.incrementStats(rho, uSqr);
}


} // end namespace

#endif

