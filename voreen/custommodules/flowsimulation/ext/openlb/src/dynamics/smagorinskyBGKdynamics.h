/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012-2015 Mathias J. Krause, Jonas Latt, Patrick Nathen
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
 */
#ifndef SMAGORINSKY_BGK_DYNAMICS_H
#define SMAGORINSKY_BGK_DYNAMICS_H

#include "dynamics/dynamics.h"
#include "core/cell.h"

#include<complex> // For shear kalman Smagorinsky - Populations

namespace olb {

/// Interface for the Large-Eddy-Simulation dynamics classes
template<typename T, template<typename U> class Lattice>
struct LESDynamics {
  /// Destructor: virtual to enable inheritance
  virtual ~LESDynamics() { }
  /// Get local effective relaxation parameter of the dynamics
  virtual T getEffectiveOmega(Cell<T,Lattice>& cell_) =0;

};

/// Implementation of Smagorinsky Dynamics
template<typename T, template<typename U> class Lattice>
class SmagorinskyDynamics : public LESDynamics<T,Lattice> {
public:
  /// Constructor
  SmagorinskyDynamics(T smagoConst_);

private:
  /// Smagorinsky constant
  T smagoConst;

protected:
  /// get the constant preFactor variable used to speed up calculations
  virtual T getPreFactor();
  /// get the Smagorinsky constant
  virtual T getSmagoConst();
  /// Compute constant prefactor variable in order to speed up the computation
  virtual T computePreFactor();
  /// Precomputed constant which speeeds up the computation
  T preFactor;
};

/// Implementation of the Smagorinsky BGK collision step
template<typename T, template<typename U> class Lattice>
class SmagorinskyBGKdynamics : public SmagorinskyDynamics<T,Lattice>, public BGKdynamics<T,Lattice> {
public:
  /// Constructor
  SmagorinskyBGKdynamics(T omega_, Momenta<T,Lattice>& momenta_,
                         T smagoConst_);
  /// Collision step
  void collide(Cell<T,Lattice>& cell, LatticeStatistics<T>& statistics_) override;
  /// Collide with fixed velocity
  void staticCollide(Cell<T,Lattice>& cell_, const T u[Lattice<T>::d],
                             LatticeStatistics<T>& statistics_) override;
  /// Set local relaxation parameter of the dynamics
  void setOmega(T omega_) override;
  /// Get local smagorinsky relaxation parameter of the dynamics
  T getEffectiveOmega(Cell<T,Lattice>& cell_) override;

protected:
  /// Computes the local smagorinsky relaxation parameter
  T computeEffectiveOmega(Cell<T,Lattice>& cell);
};

/// Implementation of the ForcedBGK collision step
template<typename T, template<typename U> class Lattice>
class SmagorinskyForcedBGKdynamics : public SmagorinskyDynamics<T,Lattice>, public ForcedBGKdynamics<T,Lattice> {
public:
  /// Constructor
  SmagorinskyForcedBGKdynamics(T omega_, Momenta<T,Lattice>& momenta_, T smagoConst_);
  /// Collision step
  virtual void collide(Cell<T,Lattice>& cell, LatticeStatistics<T>& statistics_);
  /// Collide with fixed velocity
  virtual void staticCollide(Cell<T,Lattice>& cell_, const T u[Lattice<T>::d],
                             LatticeStatistics<T>& statistics_);
  /// Set local relaxation parameter of the dynamics
  virtual void setOmega(T omega_);
  /// Get local smagorinsky relaxation parameter of the dynamics
  virtual T getEffectiveOmega(Cell<T,Lattice>& cell_);

protected:
  /// Computes the local smagorinsky relaxation parameter
  T computeEffectiveOmega(Cell<T,Lattice>& cell_);
};

/// Implementation of the consistent Strain Smagorinsky BGK collision step
///
/// Consistent subgrid scale modelling for lattice Boltzmann methods
/// Orestis Malaspinas and Pierre Sagaut
/// Journal of Fluid Mechanics / Volume / June 2012, pp 514-542
/// DOI: http://dx.doi.org/10.1017/jfm.2012.155

template<typename T, template<typename U> class Lattice>
class ConStrainSmagorinskyBGKdynamics : public SmagorinskyBGKdynamics<T,Lattice> {
public:
  /// Constructor
  ConStrainSmagorinskyBGKdynamics(T omega_, Momenta<T,Lattice>& momenta_,
                                  T smagoConst_=T(.1));
protected:
  /// Computes the local smagorinsky relaxation parameter
  T computeEffectiveOmega(Cell<T,Lattice>& cell_);
};

/// Implementation of the consistent Smagorinsky BGK collision step
///
/// Consistent subgrid scale modelling for lattice Boltzmann methods
/// Orestis Malaspinas and Pierre Sagaut
/// Journal of Fluid Mechanics / Volume / June 2012, pp 514-542
/// DOI: http://dx.doi.org/10.1017/jfm.2012.155

template<typename T, template<typename U> class Lattice>
class ConSmagorinskyBGKdynamics : public SmagorinskyBGKdynamics<T,Lattice> {
public:
  /// Constructor
  ConSmagorinskyBGKdynamics(T omega_, Momenta<T,Lattice>& momenta_, T smagoConst_);
protected:
  /// should be remove --> David
  T computeEffectiveOmega(Cell<T,Lattice>& cell_);

};

/// Implementation of a the dynamic Smarorinsky BGK collision step
template<typename T, template<typename U> class Lattice>
class DynSmagorinskyBGKdynamics : public SmagorinskyBGKdynamics<T,Lattice> {
public:
  /// Constructor
  DynSmagorinskyBGKdynamics(T omega_, Momenta<T,Lattice>& momenta_);

protected:
  /// Computes the local smagorinsky relaxation parameter
  T computeEffectiveOmega(Cell<T,Lattice>& cell);
  const static T smagoConstIsAt = Lattice<T>::ExternalField::smagoConstIsAt;
};

/// Implementation of the ADM BGK collision step

template<typename T, template<typename U> class Lattice>
class ADMBGKdynamics : public BGKdynamics<T,Lattice> {
public:
  /// Constructor
  ADMBGKdynamics(T omega_, Momenta<T,Lattice>& momenta_);
  /// Collision step
  virtual void collide(Cell<T,Lattice>& cell, LatticeStatistics<T>& statistics_);
  /// Collide with fixed velocity
  virtual void staticCollide(Cell<T,Lattice>& cell, const T u[Lattice<T>::d],
                             LatticeStatistics<T>& statistics_);
private:
  T omega;
  static const int rhoIsAt = Lattice<T>::ExternalField::rhoIsAt;
  static const int velocityBeginsAt = Lattice<T>::ExternalField::velocityBeginsAt;
  static const int sizeOfVelocity = Lattice<T>::ExternalField::sizeOfVelocity;
};

/// Implementation of the ForcedADMBGK collision step
template<typename T, template<typename U> class Lattice>
class ForcedADMBGKdynamics : public BGKdynamics<T,Lattice> {
public:
  /// Constructor
  ForcedADMBGKdynamics(T omega_, Momenta<T,Lattice>& momenta_);

  /// Collision step
  virtual void collide(Cell<T,Lattice>& cell,
                       LatticeStatistics<T>& statistics_);
  /// Collide with fixed velocity
  virtual void staticCollide(Cell<T,Lattice>& cell,
                             const T u[Lattice<T>::d],
                             LatticeStatistics<T>& statistics_);
private:
  T omega;
  static const int forceBeginsAt = Lattice<T>::ExternalField::forceBeginsAt;
  static const int sizeOfForce   = Lattice<T>::ExternalField::sizeOfForce;
  static const int filRhoIsAt = Lattice<T>::ExternalField::filRhoIsAt;
  static const int localFilVelXBeginsAt = Lattice<T>::ExternalField::localFilVelXBeginsAt;
  static const int localFilVelYBeginsAt = Lattice<T>::ExternalField::localFilVelYBeginsAt;
  static const int localFilVelZBeginsAt = Lattice<T>::ExternalField::localFilVelZBeginsAt;
};

/// Implementation of a Shear Smarorinsky BGK collision step
/// Shown good results for wall-bounded flows
/// Leveque et al.: Shear-Improved Smagorinsky Model for Large-Eddy Simulation
/// of Wall-Bounded Turbulent Flows
/// DOI: http://dx.doi.org/10.1017/S0022112006003429

template<typename T, template<typename U> class Lattice>
class ShearSmagorinskyBGKdynamics : public SmagorinskyBGKdynamics<T,Lattice> {
public:
  /// Constructor
  ShearSmagorinskyBGKdynamics(T omega_, Momenta<T,Lattice>& momenta_, T smagoConst_);
  /// Collision step
  virtual void collide(Cell<T,Lattice>& cell, LatticeStatistics<T>& statistics_);
  /// Collide with fixed velocity
  virtual void staticCollide(Cell<T,Lattice>& cell_, const T u[Lattice<T>::d],
                             LatticeStatistics<T>& statistics_);
  /// Get Effective Omega stored in a external field
  virtual T getEffectiveOmega(Cell<T,Lattice>& cell);
protected:
  /// Computes the local smagorinsky relaxation parameter
  T computeEffectiveOmega(Cell<T,Lattice>& cell, int iT);
  /// The external field variables' positions
  const static int avShearIsAt = Lattice<T>::ExternalField::avShearIsAt;

};

/// Implementation of the ForcedBGK collision step
template<typename T, template<typename U> class Lattice>
class ShearSmagorinskyForcedBGKdynamics : public SmagorinskyForcedBGKdynamics<T,Lattice> {
public:
  /// Constructor
  ShearSmagorinskyForcedBGKdynamics(T omega_, Momenta<T,Lattice>& momenta_, T smagoConst_);
  /// Collision step
  virtual void collide(Cell<T,Lattice>& cell, LatticeStatistics<T>& statistics_);
  /// Collide with fixed velocity
  virtual void staticCollide(Cell<T,Lattice>& cell_, const T u[Lattice<T>::d],
                             LatticeStatistics<T>& statistics_);
  /// Get Effective Omega stored in a external field
  virtual T getEffectiveOmega(Cell<T,Lattice>& cell);
protected:
  /// Computes the local smagorinsky relaxation parameter
  T computeEffectiveOmega(Cell<T,Lattice>& cell, int iT);
  // Define current time step
  /// Smagorinsky constant
  const static int avShearIsAt = Lattice<T>::ExternalField::avShearIsAt;
};

/// Implementation of the ForcedBGK collision step
template<typename T, template<typename U> class Lattice>
class SmagorinskyLinearVelocityForcedBGKdynamics : public SmagorinskyForcedBGKdynamics<T,Lattice> {
public:
  /// Constructor
  SmagorinskyLinearVelocityForcedBGKdynamics(T omega_, Momenta<T,Lattice>& momenta_,
      T smagoConst_);
  /// Collision step
  virtual void collide(Cell<T,Lattice>& cell, LatticeStatistics<T>& statistics_);
};

/// Implementation of the BGK collision step
template<typename T, template<typename U> class Lattice>
class KrauseBGKdynamics : public SmagorinskyBGKdynamics<T,Lattice> {
public:
  /// Constructor
  KrauseBGKdynamics(T omega_, Momenta<T,Lattice>& momenta_, T smagoConst_);
  /// Collision step
  void collide(Cell<T,Lattice>& cell, LatticeStatistics<T>& statistics_) override;
  /// Collide with fixed velocity
  void staticCollide(Cell<T,Lattice>& cell, const T u[Lattice<T>::d],
                             LatticeStatistics<T>& statistics_) override;
  /// Get local smagorinsky relaxation parameter of the dynamics
  T getEffectiveOmega(Cell<T,Lattice>& cell_) override;

private:
  /// Computes a constant prefactor in order to speed up the computation
  T computePreFactor() override;
  /// Computes the local smagorinsky relaxation parameter
  void computeEffectiveOmega(T omega0, Cell<T,Lattice>& cell, T preFactor_, T rho,
                             T u[Lattice<T>::d], T newOmega[Lattice<T>::q]);
  T preFactor;
};


/// Implementation of the BGK collision step
template<typename T, template<typename U> class Lattice>
class WALEBGKdynamics : public SmagorinskyBGKdynamics<T,Lattice> {
public:
  /// Constructor
  WALEBGKdynamics(T omega_, Momenta<T,Lattice>& momenta_, T smagoConst_);

protected:
  /// Computes a constant prefactor in order to speed up the computation
  T computePreFactor();
  /// Computes the local smagorinsky relaxation parameter
  T computeEffectiveOmega(Cell<T,Lattice>& cell_);
};

/// Implementation of the BGK collision step
template<typename T, template<typename U> class Lattice>
class FDKalmanShearSmagorinskyBGKdynamics : public SmagorinskyBGKdynamics<T,Lattice> {
public:
  /// Constructor
  FDKalmanShearSmagorinskyBGKdynamics(T omega_, Momenta<T,Lattice>& momenta_, T smagoConst_, T u_char_lat, T f_char_lat);
  /// Get local effective relaxation parameter of the dynamics
  virtual T getEffectiveOmega(Cell<T,Lattice>& cell_);

protected:
  /// Computes a constant prefactor in order to speed up the computation
  virtual T computePreFactor();
  /// Computes the local smagorinsky relaxation parameter
  virtual T computeOmega(Cell<T,Lattice>& cell_);

  // The variance of increment of kalman filtered velocity
  T VarInVelKal;
  T UCharLat;
private:
  void computeNormStrainRate(Cell<T,Lattice>& cell, int PosVelGrad, T& NormStrainRate);
  void KalmanStep(Cell<T,Lattice>& cell);
};



////////////////////////////////////////////////////////////////////////////////
/// Implementation of a Shear Smarorinsky BGK collision step with Kalman Filter
//
/// Leveque et al.: Shear-Improved Smagorinsky Model for Large-Eddy Simulation
/// of Wall-Bounded Turbulent Flows
///
/// Boudet et al. (2016) A Kalman filter adapted of the estimation of mean gradients
//   in the a large-eddy simulation of unsteady turbulent flows.

template<typename T, template<typename U> class Lattice>
class ShearKalmanSmagorinskyBGKdynamics : public SmagorinskyBGKdynamics<T,Lattice> {
public:
  /// Constructor
  ShearKalmanSmagorinskyBGKdynamics(T omega_, Momenta<T,Lattice>& momenta_,
                                     T smagoConst_, T u_char_lat, T f_char_lat);
  /// Get local effective relaxation parameter of the dynamics
  virtual T getEffectiveOmega(Cell<T,Lattice>& cell_);

protected:
  /// Computes the local smagorinsky relaxation parameter
  T computeEffectiveOmega(Cell<T,Lattice>& cell_);
private:
  /// Updates the filtered velocity with a Kalman procedure
  void KalmanStep(Cell<T,Lattice>& cell);
  /// Computes the kalman filtered velocity and strain rate using the filtered population stored in a externa field
  void computeKalmanUStress(Cell<T,Lattice>& cell, T (&KalmanU)[Lattice<T>::d], T (&KalmanPi)[util::TensorVal<Lattice<T> >::n]);
  /// Computes The Kalman filtered velocity using the filtered populations stored in a external field
  void computeKalmanU(Cell<T,Lattice>& cell, T (&KalmanU)[Lattice<T>::d]);
  /// Computes the Kalman filtered strain rate using the filtered populations stored in a external field
  void computeKalmanStress(Cell<T,Lattice>& cell, T (&KalmanU)[Lattice<T>::d], T (&KalmanPi)[util::TensorVal<Lattice<T> >::n]);
  /// Computes instantaneous tau_sgs and update kalman tau_sgs
  void computeAndupdateTauSgs(Cell<T,Lattice>& cell, T rho, T pi[util::TensorVal<Lattice<T> >::n],
                              T KalmanPiNeqN[util::TensorVal<Lattice<T> >::n], T KalmanPiNeqN1[util::TensorVal<Lattice<T> >::n],
                              T K, T &tau_sgs);
  /// Methods to compute the square Norm of second order moment non-quilibrium distribution function
  void computeNormSOM(T pi[util::TensorVal<Lattice<T> >::n], T &piNorm);
  void computeNormSOM(T pi1[util::TensorVal<Lattice<T> >::n], T pi2[util::TensorVal<Lattice<T> >::n], T rho, T &piNorm);
  void computeNormSOM(T pi[util::TensorVal<Lattice<T> >::n], T rho, T &piNorm);
  /// Compute the instantaneous tau_sgs
  void computeTauSgs(Cell<T,Lattice>& cell, T rho, T KalmanPiNeqNormSqr, T KalmanInstPiNeqNormSqr, T PiNeqNormSqr, T K, T &tau_sgs);
  void computeRoots4thPoly(T A, T B, T C, T D, T E, std::complex<T> (&Roots)[4]);
  // Update the local kalman tau_sgs stored in a external field
  void updateTauSgsKalman(Cell<T,Lattice>& cell, T NN, T Nn1, T n1n1, T N1N1, T K, T tau_sgs_n1);

  // The variance of increment of kalman filtered velocity
  T VarInVelKal;
  T UCharLat;
  /// The external field variables' positions
  static const int ErrorCovarianceIsAt = Lattice<T>::ExternalField::ErrorCovarianceIsAt;
  static const int VarianceIsAt = Lattice<T>::ExternalField::VarianceIsAt;
  static const int TauSgsIsAt = Lattice<T>::ExternalField::TauSgsIsAt;
  static const int FilteredPopulationIsAt = Lattice<T>::ExternalField::FilteredPopulationIsAt;
};



}

#endif
