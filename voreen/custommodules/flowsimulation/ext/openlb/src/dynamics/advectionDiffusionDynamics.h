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
 * can be instantiated -- header file.
 */
#ifndef ADVECTION_DIFFUSION_DYNAMICS_H
#define ADVECTION_DIFFUSION_DYNAMICS_H

#include "latticeDescriptors.h"
#include "dynamics/dynamics.h"

namespace olb {

// ========= the RLB advection diffusion dynamics ========//
/// it uses the regularized approximation that can be found in the thesis of J. Latt (2007).
template<typename T, template<typename U> class Lattice>
class AdvectionDiffusionRLBdynamics : public BasicDynamics<T, Lattice> {
public:
  /// Constructor
  AdvectionDiffusionRLBdynamics( T omega_, Momenta<T, Lattice>& momenta_ );
  /// Compute equilibrium distribution function
  T computeEquilibrium( int iPop, T rho, const T u[Lattice<T>::d], T uSqr ) const override;
  /// Collision step
  void collide( Cell<T, Lattice>& cell, LatticeStatistics<T>& statistics ) override;
  /// Collide with fixed velocity
  void staticCollide( Cell<T, Lattice>& cell, const T u[Lattice<T>::d],
                              LatticeStatistics<T>& statistics ) override;
  /// Get local relaxation parameter of the dynamics
  T getOmega() const override;
  /// Set local relaxation parameter of the dynamics
  void setOmega( T omega_ ) override;
private:
  T omega;  ///< relaxation parameter
};

// ========= the BGK advection diffusion dynamics ========//
/// This approach contains a slight error in the diffusion term.
template<typename T, template<typename U> class Lattice>
class AdvectionDiffusionBGKdynamics : public BasicDynamics<T, Lattice> {
public:
  /// Constructor
  AdvectionDiffusionBGKdynamics( T omega, Momenta<T, Lattice>& momenta );
  AdvectionDiffusionBGKdynamics( const UnitConverter<T,Lattice>& converter, Momenta<T, Lattice>& momenta );
  /// Compute equilibrium distribution function
  T computeEquilibrium( int iPop, T rho, const T u[Lattice<T>::d], T uSqr ) const override;
  /// Collision step
  void collide( Cell<T, Lattice>& cell, LatticeStatistics<T>& statistics ) override;
  /// Collide with fixed velocity
  void staticCollide( Cell<T, Lattice>& cell, const T u[Lattice<T>::d],
                              LatticeStatistics<T>& statistics ) override;
  /// Get local relaxation parameter of the dynamics
  T getOmega() const override;
  /// Set local relaxation parameter of the dynamics
  void setOmega( T omega ) override;
private:
  T _omega;  ///< relaxation parameter
};

// ========= the BGK advection diffusion Stokes drag dynamics with a Smagorinsky turbulence model ========//
/// This approach contains a slight error in the diffusion term.
template<typename T, template<typename U> class Lattice>
class SmagorinskyParticleAdvectionDiffusionBGKdynamics : public olb::AdvectionDiffusionBGKdynamics<T,Lattice> {
public:
  /// Constructor
  SmagorinskyParticleAdvectionDiffusionBGKdynamics(T omega_, Momenta<T,Lattice>& momenta_, T smagoConst_, T dx_, T dt_);
  /// Collision step
  virtual void collide(Cell<T,Lattice>& cell, LatticeStatistics<T>& statistics );
  /// Collide with fixed velocity
  virtual void staticCollide(Cell<T,Lattice>& cell, const T u[Lattice<T>::d],
                             LatticeStatistics<T>& statistics );
  /// Get local smagorinsky relaxation parameter of the dynamics
  virtual T getSmagorinskyOmega(Cell<T,Lattice>& cell);
  /// Set local relaxation parameter of the dynamics
  virtual void setOmega(T omega_);

private:
  /// Computes a constant prefactor in order to speed up the computation
  T computePreFactor(T omega_, T smagoConst_, T dx_, T dt_);
  /// Computes the local smagorinsky relaxation parameter
  T computeOmega(T omega0, T preFacto_r, T rho, T pi[util::TensorVal<Lattice<T> >::n] );

  /// effective collision time based upon Smagorisnky approach
  T tau_eff;
  /// Smagorinsky constant
  T smagoConst;
  /// Precomputed constant which speeeds up the computation
  T preFactor;
  T dx;
  T dt;
};

// ========= the BGK advection diffusion Stokes drag dynamics  ========//
/// This approach contains a slight error in the diffusion term.
template<typename T, template<typename U> class Lattice>
class ParticleAdvectionDiffusionBGKdynamics : public olb::AdvectionDiffusionBGKdynamics<T,Lattice> {
public:
  /// Constructor
  ParticleAdvectionDiffusionBGKdynamics(T omega_, Momenta<T,Lattice>& momenta_);
  /// Collision step
  void collide(Cell<T,Lattice>& cell, LatticeStatistics<T>& statistics ) override;
private:
  T omega;  ///< relaxation parameter
};


// ========= the MRT advection diffusion dynamics ========//
/// This approach is based on the multi-distribution LBM model.
/// The couplintg is done using the Boussinesq approximation
template<typename T, template<typename U> class Lattice>
class AdvectionDiffusionMRTdynamics : public BasicDynamics<T, Lattice> {
public:
  /// Constructor
  AdvectionDiffusionMRTdynamics( T omega, Momenta<T, Lattice>& momenta );
  /// Compute equilibrium distribution function
  T computeEquilibrium( int iPop, T rho, const T u[Lattice<T>::d], T uSqr ) const override;
  /// Collision step
  void collide( Cell<T, Lattice>& cell, LatticeStatistics<T>& statistics ) override;
  /// Collide with fixed velocity
  void staticCollide( Cell<T, Lattice>& cell, const T u[Lattice<T>::d],
                              LatticeStatistics<T>& statistics ) override;
  /// Get local relaxation parameter of the dynamics
  T getOmega() const override;
  /// Set local relaxation parameter of the dynamics
  void setOmega( T omega ) override;
private:
  T _omega;  ///< relaxation parameter
protected:
  T invM_S[Lattice<T>::q][Lattice<T>::q]; ///< inverse relaxation times matrix
};

} // namespace olb

#endif
