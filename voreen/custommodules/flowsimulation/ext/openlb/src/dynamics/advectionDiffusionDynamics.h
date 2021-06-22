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
#include "core/unitConverter.h"

namespace olb {

// ========= the RLB advection diffusion dynamics ========//
/// it uses the regularized approximation that can be found in the thesis of J. Latt (2007).
template<typename T, typename DESCRIPTOR>
class AdvectionDiffusionRLBdynamics : public BasicDynamics<T, DESCRIPTOR> {
public:
  /// Constructor
  AdvectionDiffusionRLBdynamics( T omega_, Momenta<T, DESCRIPTOR>& momenta_ );
  /// Compute equilibrium distribution function
  T computeEquilibrium( int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr ) const override;
  /// Collision step
  void collide( Cell<T,DESCRIPTOR>& cell, LatticeStatistics<T>& statistics ) override;
  /// Get local relaxation parameter of the dynamics
  T getOmega() const override;
  /// Set local relaxation parameter of the dynamics
  void setOmega( T omega_ ) override;
private:
  T omega;  ///< relaxation parameter
};

/// Implementation of Regularized BGK collision, followed by any Dynamics
template<typename T, typename DESCRIPTOR, typename Dynamics>
class CombinedAdvectionDiffusionRLBdynamics : public BasicDynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  CombinedAdvectionDiffusionRLBdynamics(T omega, Momenta<T,DESCRIPTOR>& momenta);
  /// Compute equilibrium distribution function
  T computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr) const override;
  /// Collision step
  void collide(Cell<T,DESCRIPTOR>& cell,
                       LatticeStatistics<T>& statistics_) override;
  /// Get local relaxation parameter of the dynamics
  T getOmega() const override;
  /// Set local relaxation parameter of the dynamics
  void setOmega(T omega) override;
private:
  Dynamics _boundaryDynamics;
};

// ========= the BGK advection diffusion dynamics ========//
/// This approach contains a slight error in the diffusion term.
template<typename T, typename DESCRIPTOR>
class AdvectionDiffusionBGKdynamics : public BasicDynamics<T, DESCRIPTOR> {
public:
  /// Constructor
  AdvectionDiffusionBGKdynamics( T omega, Momenta<T, DESCRIPTOR>& momenta );
  AdvectionDiffusionBGKdynamics( const UnitConverter<T,DESCRIPTOR>& converter, Momenta<T, DESCRIPTOR>& momenta );
  /// Compute equilibrium distribution function
  T computeEquilibrium( int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr ) const override;
  /// Collision step
  void collide( Cell<T,DESCRIPTOR>& cell, LatticeStatistics<T>& statistics ) override;
  /// Get local relaxation parameter of the dynamics
  T getOmega() const override;
  /// Set local relaxation parameter of the dynamics
  void setOmega( T omega ) override;
protected:
  T _omega;  ///< relaxation parameter
};



// ========= the TRT advection diffusion dynamics ========//
/// This approach contains a slight error in the diffusion term.
template<typename T, typename DESCRIPTOR>
class AdvectionDiffusionTRTdynamics : public AdvectionDiffusionBGKdynamics<T, DESCRIPTOR> {
public:
  /// Constructor
  AdvectionDiffusionTRTdynamics( T omega, Momenta<T, DESCRIPTOR>& momenta, T magicParameter );
  /// Collision step
  void collide( Cell<T,DESCRIPTOR>& cell, LatticeStatistics<T>& statistics ) override;
protected:
  T _omega2; /// relaxation parameter for odd moments
  T _magicParameter;
};


// ======= BGK advection diffusion dynamics with source term  ======//
// following Seta, T. (2013). Implicit temperature-correction-based
// immersed-boundary thermal lattice Boltzmann method for the simulation
// of natural convection. Physical Review E, 87(6), 063304.
template<typename T, typename DESCRIPTOR>
class SourcedAdvectionDiffusionBGKdynamics : public AdvectionDiffusionBGKdynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  SourcedAdvectionDiffusionBGKdynamics(T omega_, Momenta<T,DESCRIPTOR>& momenta_ );
  /// Collision step
  void collide(Cell<T,DESCRIPTOR>& cell, LatticeStatistics<T>& statistics ) override;
  /// Compute Density
  T computeRho(ConstCell<T,DESCRIPTOR>& cell) const override;
  /// Compute fluid velocity and particle density on the cell.
  void computeRhoU (
    ConstCell<T,DESCRIPTOR>& cell,
    T& rho, T u[DESCRIPTOR::d]) const override;
private:
  const T _omegaMod;
};

// ======= BGK advection diffusion dynamics for solid-liquid phase change  ======//
// following Huang, R. (2015). Phase interface effects in the total
// enthalpy-based lattice Boltzmann model for solid–liquid phase change.
// Journal of Computational Physics, 294, 345-362.
template<typename T, typename DESCRIPTOR>
class TotalEnthalpyAdvectionDiffusionBGKdynamics : public AdvectionDiffusionBGKdynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  TotalEnthalpyAdvectionDiffusionBGKdynamics(T omega_, Momenta<T,DESCRIPTOR>& momenta_,
      T T_s_, T T_l_, T cp_s_, T cp_l_, T lambda_s_, T lambda_l_, T l_);
  /// Compute equilibrium distribution function
  T computeEquilibrium( int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr ) const override;
  /// Collision step
  void collide(Cell<T,DESCRIPTOR>& cell, LatticeStatistics<T>& statistics ) override;
  T computeTemperature( T enthalpy ) const;
  T computeLiquidFraction( T enthalpy ) const;
protected:
  const T _T_s, _T_l, _cp_s, _cp_l, _lambda_s, _lambda_l, _l;
  const T _H_s, _H_l, _cp_ref;
};

// ======= TRT advection diffusion dynamics for solid-liquid phase change  ======//
// following Huang, R. (2015). Phase interface effects in the total
// enthalpy-based lattice Boltzmann model for solid–liquid phase change.
// Journal of Computational Physics, 294, 345-362.
template<typename T, typename DESCRIPTOR>
class TotalEnthalpyAdvectionDiffusionTRTdynamics : public TotalEnthalpyAdvectionDiffusionBGKdynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  TotalEnthalpyAdvectionDiffusionTRTdynamics(T omega_, Momenta<T,DESCRIPTOR>& momenta_,
      T magicParameter_,
      T T_s_, T T_l_, T cp_s_, T cp_l_, T lambda_s_, T lambda_l_, T l_);
  /// Collision step
  void collide(Cell<T,DESCRIPTOR>& cell, LatticeStatistics<T>& statistics ) override;
private:
  const T _magicParameter;
};

// ======= BGK advection diffusion dynamics for phase field equation  ======//
// following Fakhari, Abbas, et al. (2017). Improved locality of the phase-field
// lattice-Boltzmann model for immiscible fluids at high density ratios.
// Physical Review E 96.5, 053301.
template<typename T, typename DESCRIPTOR>
class PhaseFieldAdvectionDiffusionBGKdynamics : public AdvectionDiffusionBGKdynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  PhaseFieldAdvectionDiffusionBGKdynamics(T omega_, Momenta<T,DESCRIPTOR>& momenta_, T interface_thickness = 3.0 );
  /// Collision step
  void collide(Cell<T,DESCRIPTOR>& cell, LatticeStatistics<T>& statistics ) override;
protected:
  const T _mobility, _interface_thickness;
};

// ========= the BGK advection diffusion Stokes drag dynamics with a Smagorinsky turbulence model ========//
/// This approach contains a slight error in the diffusion term.
template<typename T, typename DESCRIPTOR>
class SmagorinskyParticleAdvectionDiffusionBGKdynamics : public olb::AdvectionDiffusionBGKdynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  SmagorinskyParticleAdvectionDiffusionBGKdynamics(T omega_, Momenta<T,DESCRIPTOR>& momenta_, T smagoConst_, T dx_, T dt_);
  /// Collision step
  virtual void collide(Cell<T,DESCRIPTOR>& cell, LatticeStatistics<T>& statistics );
  /// Get local smagorinsky relaxation parameter of the dynamics
  virtual T getSmagorinskyOmega(Cell<T,DESCRIPTOR>& cell);
  /// Set local relaxation parameter of the dynamics
  virtual void setOmega(T omega_);

private:
  /// Computes a constant prefactor in order to speed up the computation
  T computePreFactor(T omega_, T smagoConst_, T dx_, T dt_);
  /// Computes the local smagorinsky relaxation parameter
  T computeOmega(T omega0, T preFacto_r, T rho, T pi[util::TensorVal<DESCRIPTOR >::n] );

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
template<typename T, typename DESCRIPTOR>
class ParticleAdvectionDiffusionBGKdynamics : public olb::AdvectionDiffusionBGKdynamics<T,DESCRIPTOR> {
public:
  /// Constructor
  ParticleAdvectionDiffusionBGKdynamics(T omega_, Momenta<T,DESCRIPTOR>& momenta_);
  /// Collision step
  void collide(Cell<T,DESCRIPTOR>& cell, LatticeStatistics<T>& statistics ) override;
private:
  T omega;  ///< relaxation parameter
};


// ========= the MRT advection diffusion dynamics ========//
/// This approach is based on the multi-distribution LBM model.
/// The couplintg is done using the Boussinesq approximation
template<typename T, typename DESCRIPTOR>
class AdvectionDiffusionMRTdynamics : public BasicDynamics<T, DESCRIPTOR> {
public:
  /// Constructor
  AdvectionDiffusionMRTdynamics( T omega, Momenta<T, DESCRIPTOR>& momenta );
  /// Compute equilibrium distribution function
  T computeEquilibrium( int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr ) const override;
  /// Collision step
  void collide( Cell<T,DESCRIPTOR>& cell, LatticeStatistics<T>& statistics ) override;
  /// Get local relaxation parameter of the dynamics
  T getOmega() const override;
  /// Set local relaxation parameter of the dynamics
  void setOmega( T omega ) override;
private:
  T _omega;  ///< relaxation parameter
protected:
  T invM_S[DESCRIPTOR::q][DESCRIPTOR::q]; ///< inverse relaxation times matrix
};

} // namespace olb

#endif
