/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2017 Albert Mink, Christopher McHardy
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
 * A collection of radiative transport dynamics classes -- header file.
 */

#ifndef RTLBM_DYNAMICS_H
#define RTLBM_DYNAMICS_H

#include "dynamics/dynamics.h"

namespace olb {

/**
 * Solves poisson equation, according to Mink et al 2016. First order scheme.
 *
 * \f[ D \Delta \Phi = \sigma * \Phi \f]
 *
 * \param _omega       is relaxation parameter and accounts the time scale of the problem
 * \param _sink        corresponds to \f$ sink = \sigma / 8 \f$
 */
template<typename T, template<typename U> class Lattice>
class RTLBMdynamicsMink : public BasicDynamics<T, Lattice> {
public:
  /// Constructor
  RTLBMdynamicsMink( T omega, Momenta<T, Lattice>& momenta, T latticeAbsorption, T latticeScattering );
  /// Compute equilibrium distribution function
  T computeEquilibrium( int iPop, T rho, const T u[Lattice<T>::d], T uSqr ) const override;
  /// Collision step
  void collide( Cell<T, Lattice>& cell, LatticeStatistics<T>& statistics ) override;
  /// Collide with fixed velocity
  void staticCollide( Cell<T, Lattice>& cell, const T u[Lattice<T>::d], LatticeStatistics<T>& statistics ) override;
  /// Get local relaxation parameter of the dynamics
  T getOmega() const override;
  /// Set local relaxation parameter of the dynamics
  void setOmega( T omega ) override;
  T getSink() const;
private:
  T _omega;
  T _sink;    // 3*latticeAbs*(latticeAbs+latticeScat) / 8
};

/**
 * Solves poisson equation, according to Mink et al 2016.
 * Second order scheme, due to no collision with fixed relaxation time.
 *
 * Target equation (Poisson problem)
 * \f[ \Delta \Phi = f(\Phi) \f]
 * for linear function \f$ f(\Phi) = \alpha \Phi \f$.
 *
 * Where parameter \f$ \alpha \f$ (nondimensional) is related to sink parameter of dynamics by
 * \f[ sink = \alpha / 8 \f]
 *

 * \param _sink     given by \f$ sink = \alpha / 8 \f$
 */
template<typename T, template<typename U> class Lattice>
class RTLBMconstDynamicsMink : public BasicDynamics<T, Lattice> {
public:
  /// Constructor
  RTLBMconstDynamicsMink( Momenta<T, Lattice>& momenta, T latticeAbsorption, T latticeScattering );
  /// Compute equilibrium distribution function
  T computeEquilibrium( int iPop, T rho, const T u[Lattice<T>::d], T uSqr ) const override;
  /// Collision step
  void collide( Cell<T, Lattice>& cell, LatticeStatistics<T>& statistics ) override;
  /// Collide with fixed velocity
  void staticCollide( Cell<T, Lattice>& cell, const T u[Lattice<T>::d], LatticeStatistics<T>& statistics ) override;
  /// Get local relaxation parameter of the dynamics
  T getOmega() const override;
  /// Set local relaxation parameter of the dynamics
  void setOmega( T omega ) override;
private:
  T _sink;
};


/**
 * Solves RTE according Christopher McHardy et al 2016.
 * absorption and scattering coefficient:
 * \f$ \sigma_a \f$ and \f$ \sigma_s \f$
 *
 * \param omega             change into beta the extinction coefficient
 * \param singleScatAlbedo  is the single scattering albedo, given by \f$ \frac{\sigma_s}{sigma_a + sigma_s} \f$
 */
template<typename T, template<typename U> class Lattice>
class RTLBMdynamicsMcHardy : public BasicDynamics<T, Lattice> {
public:
  /// Constructor
  RTLBMdynamicsMcHardy( Momenta<T,Lattice>& momenta, T latticeAbsorption, T latticeScattering );
  /// Compute equilibrium distribution function
  T computeEquilibrium( int iPop, T rho, const T u[Lattice<T>::d], T uSqr ) const override;
  /// Collision step
  void collide( Cell<T, Lattice>& cell, LatticeStatistics<T>& statistics ) override;
  /// Collide with fixed velocity
  void staticCollide( Cell<T, Lattice>& cell, const T u[Lattice<T>::d], LatticeStatistics<T>& statistics ) override;
  /// Get local relaxation parameter of the dynamics
  T getOmega() const override;
  /// Set local relaxation parameter of the dynamics
  void setOmega( T omega ) override;
  T getSink() const;

protected:
  T _absorption;
  T _scattering;
};

template<typename T, template<typename U> class Lattice>
class RTLBMdynamicsMcHardyWH : public RTLBMdynamicsMcHardy<T, Lattice> {
public:
  /// Constructor
  RTLBMdynamicsMcHardyWH( Momenta<T,Lattice>& momenta, T latticeAbsorption, T latticeScattering );
  /// Compute equilibrium distribution function
  T computeEquilibrium( int iPop, T rho, const T u[Lattice<T>::d], T uSqr ) const override;
  /// Collision step
  void collide( Cell<T, Lattice>& cell, LatticeStatistics<T>& statistics ) override;
  /// Collide with fixed velocity
  void staticCollide( Cell<T, Lattice>& cell, const T u[Lattice<T>::d], LatticeStatistics<T>& statistics ) override;
  /// Get local relaxation parameter of the dynamics
  T getOmega() const override;
  /// Set local relaxation parameter of the dynamics
  void setOmega( T omega ) override;
private:
//  T _absorption;
//  T _scattering;
};


}  // namespace olb

#endif
