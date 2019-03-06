/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012, 2015 Mathias J. Krause, Vojtech Cvrcek
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
 * BGK Dynamics with adjusted omega -- header file.
 * Strain rate similar to "J.Boyd, J. Buick and S.Green: A second-order accurate lattice Boltzmann non-Newtonian flow model"
 * Power Law similar to "Huidan Yu, Sharath S. Girimaji, Li-Shi Luo - DNS and LES of decaying isotropic turbulence with and without frame rotation using lattice Boltzmann method"
 *
 */
#ifndef POWER_LAW_BGK_DYNAMICS_H
#define POWER_LAW_BGK_DYNAMICS_H

#include "dynamics/dynamics.h"
#include "core/cell.h"

namespace olb {

/// Implementation of the BGK collision step
template<typename T, template<typename U> class Lattice>
class PowerLawBGKdynamics : public BGKdynamics<T,Lattice> {
public:
  /// Constructor
  /// m,n...parameter in the power law model, dt...the explizit typed time step from LBM model
  PowerLawBGKdynamics(T omega_, Momenta<T,Lattice>& momenta_, T m_=0.1, T n_=.5, T dt_=T(0.0016));

  /// Collision step
  void collide(Cell<T,Lattice>& cell,
                       LatticeStatistics<T>& statistics_) override;
  /// Collide with fixed velocity
  void staticCollide(Cell<T,Lattice>& cell,
                             const T u[Lattice<T>::d],
                             LatticeStatistics<T>& statistics_) override;

  /// Set local relaxation parameter of the dynamics
  void setOmega(T omega_) override;

  /// Get local powerLaw relaxation parameter of the dynamics
  virtual T getPowerLawOmega(Cell<T,Lattice>& cell_);

private:
  /// Computes a constant prefactor in order to speed up the computation
  T computePreFactor(T omega_, T smagoConst_);

  /// Computes the local powerLaw relaxation parameter
  T computeOmega(T omega0_, T preFactor_, T rho_, T pi_[util::TensorVal<Lattice<T> >::n] );

private:
  T smagoConst;  ///< PowerLaw constant
  T preFactor; ///< Precomputed constant which speeeds up the computation
  T m;
  T n;
  T dt;
};

template<typename T, template<typename U> class Lattice>
class PowerLawForcedBGKdynamics : public BGKdynamics<T,Lattice> {
public:
  /// Constructor
  /// m,n...parameter in the power law model, dt...the explizit typed time step from LBM model
  PowerLawForcedBGKdynamics(T omega_, Momenta<T,Lattice>& momenta_, T m_=0.1, T n_=.5, T dt_=T(0.0016));

  /// Collision step
  virtual void collide(Cell<T,Lattice>& cell,
                       LatticeStatistics<T>& statistics_);
  /// Collide with fixed velocity
  virtual void staticCollide(Cell<T,Lattice>& cell,
                             const T u[Lattice<T>::d],
                             LatticeStatistics<T>& statistics_);

  /// Set local relaxation parameter of the dynamics
  virtual void setOmega(T omega_);

  /// Get local powerLaw relaxation parameter of the dynamics
  virtual T getPowerLawOmega(Cell<T,Lattice>& cell_);

private:
  /// Computes a constant prefactor in order to speed up the computation
  T computePreFactor(T omega_, T smagoConst_);

  /// Computes the local powerLaw relaxation parameter
  T computeOmega(T omega0_, T preFactor_, T rho_, T pi_[util::TensorVal<Lattice<T> >::n] );

private:
  T smagoConst;  ///< PowerLaw constant
  T preFactor; ///< Precomputed constant which speeeds up the computation
  T m;
  T n;
  T dt;
};
}

#endif
