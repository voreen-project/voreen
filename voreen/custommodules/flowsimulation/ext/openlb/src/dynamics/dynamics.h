/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Jonas Latt, Mathias J. Krause
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
#ifndef LB_DYNAMICS_H
#define LB_DYNAMICS_H

#include "latticeDescriptors.h"
#include "core/util.h"
#include "core/postProcessing.h"
#include "core/latticeStatistics.h"

namespace olb {


template<typename T, template<typename U> class Lattice> class Cell;

/// Interface for the dynamics classes
template<typename T, template<typename U> class Lattice>
struct Dynamics {
  /// Destructor: virtual to enable inheritance
  virtual ~Dynamics() { }
  /// Implementation of the collision step
  virtual void collide(Cell<T,Lattice>& cell,
                       LatticeStatistics<T>& statistics_) =0;
  /// Collide with fixed velocity
  virtual void staticCollide(Cell<T,Lattice>& cell,
                             const T u[Lattice<T>::d],
                             LatticeStatistics<T>& statistics_) =0;
  /// Compute equilibrium distribution function
  virtual T computeEquilibrium(int iPop, T rho, const T u[Lattice<T>::d], T uSqr) const;
  /// Initialize cell at equilibrium distribution
  void iniEquilibrium(Cell<T,Lattice>& cell, T rho, const T u[Lattice<T>::d]);
  /// Compute particle density on the cell.
  /** \return particle density
   */
  virtual T computeRho(Cell<T,Lattice> const& cell) const =0;
  /// Compute fluid velocity on the cell.
  /** \param u fluid velocity
   */
  virtual void computeU( Cell<T,Lattice> const& cell,
                         T u[Lattice<T>::d] ) const =0;
  /// Compute fluid momentum (j=rho*u) on the cell.
  /** \param j fluid momentum
   */
  virtual void computeJ( Cell<T,Lattice> const& cell,
                         T j[Lattice<T>::d] ) const =0;
  /// Compute components of the stress tensor on the cell.
  /** \param pi stress tensor */
  virtual void computeStress (
    Cell<T,Lattice> const& cell,
    T rho, const T u[Lattice<T>::d],
    T pi[util::TensorVal<Lattice<T> >::n] ) const =0;
  /// Compute fluid velocity and particle density on the cell.
  /** \param rho particle density
   *  \param u fluid velocity
   */
  virtual void computeRhoU (
    Cell<T,Lattice> const& cell,
    T& rho, T u[Lattice<T>::d]) const =0;
  /// Compute all momenta on the cell, up to second order.
  /** \param rho particle density
   *  \param u fluid velocity
   *  \param pi stress tensor
   */
  virtual void computeAllMomenta (
    Cell<T,Lattice> const& cell,
    T& rho, T u[Lattice<T>::d],
    T pi[util::TensorVal<Lattice<T> >::n] ) const =0;
  /// Set particle density on the cell.
  /** \param rho particle density
   */
  virtual void defineRho(Cell<T,Lattice>& cell, T rho) =0;
  virtual void defineRho(int iPop, T rho);
  /// Set fluid velocity on the cell.
  /** \param u fluid velocity
   */
  virtual void defineU(Cell<T,Lattice>& cell,
                       const T u[Lattice<T>::d]) =0;
  /// Functions for offLattice Velocity boundary conditions
  virtual void setBoundaryIntersection(int iPop, T distance);
  virtual bool getBoundaryIntersection(int iPop, T point[Lattice<T>::d]);
  virtual void defineU(const T u[Lattice<T>::d]);
  virtual void defineU(int iPop, const T u[Lattice<T>::d]);
  virtual T    getVelocityCoefficient(int iPop);

  /// Define fluid velocity and particle density on the cell.
  /** \param rho particle density
   *  \param u fluid velocity
   */
  virtual void defineRhoU (
    Cell<T,Lattice>& cell,
    T rho, const T u[Lattice<T>::d]) =0;
  /// Define all momenta on the cell, up to second order.
  /** \param rho particle density
   *  \param u fluid velocity
   *  \param pi stress tensor
   */
  virtual void defineAllMomenta (
    Cell<T,Lattice>& cell,
    T rho, const T u[Lattice<T>::d],
    const T pi[util::TensorVal<Lattice<T> >::n] ) =0;
  /// Get local relaxation parameter of the dynamics
  virtual T getOmega() const =0;
  /// Set local relaxation parameter of the dynamics
  virtual void setOmega(T omega) =0;
};

/// Interface for classes that compute velocity momenta
/** This class is useful for example to distinguish between bulk and
 * boundary nodes, given that on the boundaries, a particular strategy
 * must be applied to compute the velocity momenta.
 */
template<typename T, template<typename U> class Lattice>
struct Momenta {
  /// Destructor: virtual to enable inheritance
  virtual ~Momenta() { }
  /// Compute particle density on the cell.
  virtual T computeRho(Cell<T,Lattice> const& cell) const =0;
  /// Compute fluid velocity on the cell.
  virtual void computeU (
    Cell<T,Lattice> const& cell,
    T u[Lattice<T>::d] ) const =0;
  /// Compute fluid momentum on the cell.
  virtual void computeJ (
    Cell<T,Lattice> const& cell,
    T j[Lattice<T>::d] ) const =0;
  /// Compute components of the stress tensor on the cell.
  virtual void computeStress (
    Cell<T,Lattice> const& cell,
    T rho, const T u[Lattice<T>::d],
    T pi[util::TensorVal<Lattice<T> >::n] ) const =0;
  /// Compute fluid velocity and particle density on the cell.
  virtual void computeRhoU (
    Cell<T,Lattice> const& cell,
    T& rho, T u[Lattice<T>::d]) const;
  /// Compute all momenta on the cell, up to second order.
  virtual void computeAllMomenta (
    Cell<T,Lattice> const& cell,
    T& rho, T u[Lattice<T>::d],
    T pi[util::TensorVal<Lattice<T> >::n] ) const;
  /// Set particle density on the cell.
  virtual void defineRho(Cell<T,Lattice>& cell, T rho) =0;
  /// Set fluid velocity on the cell.
  virtual void defineU(Cell<T,Lattice>& cell,
                       const T u[Lattice<T>::d]) =0;
  /// Define fluid velocity and particle density on the cell.
  virtual void defineRhoU (
    Cell<T,Lattice>& cell,
    T rho, const T u[Lattice<T>::d]);
  /// Define all momenta on the cell, up to second order.
  virtual void defineAllMomenta (
    Cell<T,Lattice>& cell,
    T rho, const T u[Lattice<T>::d],
    const T pi[util::TensorVal<Lattice<T> >::n] ) =0;
};

/// Abstract base for dynamics classes
/** In this version of the Dynamics classes, the computation of the
 * velocity momenta is taken care of by an object of type Momenta.
 */
template<typename T, template<typename U> class Lattice>
class BasicDynamics : public Dynamics<T,Lattice> {
public:
  /// Must be contructed with an object of type Momenta
  BasicDynamics(Momenta<T,Lattice>& momenta);
  /// Implemented via the Momenta object
  T computeRho(Cell<T,Lattice> const& cell) const override;
  /// Implemented via the Momenta object
  void computeU (
    Cell<T,Lattice> const& cell,
    T u[Lattice<T>::d] ) const override;
  /// Implemented via the Momenta object
  void computeJ (
    Cell<T,Lattice> const& cell,
    T j[Lattice<T>::d] ) const override;
  /// Implemented via the Momenta object
  void computeStress (
    Cell<T,Lattice> const& cell,
    T rho, const T u[Lattice<T>::d],
    T pi[util::TensorVal<Lattice<T> >::n] ) const override;
  /// Implemented via the Momenta object
  void computeRhoU (
    Cell<T,Lattice> const& cell,
    T& rho, T u[Lattice<T>::d]) const override;
  /// Implemented via the Momenta object
  void computeAllMomenta (
    Cell<T,Lattice> const& cell,
    T& rho, T u[Lattice<T>::d],
    T pi[util::TensorVal<Lattice<T> >::n] ) const override;
  /// Implemented via the Momenta object
  void defineRho(Cell<T,Lattice>& cell, T rho) override;
  /// Implemented via the Momenta object
  void defineU(Cell<T,Lattice>& cell,
                       const T u[Lattice<T>::d]) override;
  /// Implemented via the Momenta object
  void defineRhoU (
    Cell<T,Lattice>& cell,
    T rho, const T u[Lattice<T>::d]) override;
  /// Implemented via the Momenta object
  void defineAllMomenta (
    Cell<T,Lattice>& cell,
    T rho, const T u[Lattice<T>::d],
    const T pi[util::TensorVal<Lattice<T> >::n] ) override;
protected:
  Momenta<T,Lattice>& _momenta;  ///< computation of velocity momenta
};

/// Implementation of a generic dynamics to realize a pressure drop at a periodic boundary
template<typename T, template<typename U> class Lattice, typename BaseDynamics>
class PeriodicPressureDynamics : public BaseDynamics {

public:
  /// Constructor
  PeriodicPressureDynamics(BaseDynamics& baseDynamics, T densityOffset, int nx, int ny, int nz=0) : BaseDynamics(baseDynamics), _densityOffset(densityOffset), _nx(nx), _ny(ny), _nz(nz) {};

  /// Implementation of the collision step
  void collide(Cell<T,Lattice>& cell,
                       LatticeStatistics<T>& statistics_) override{
    BaseDynamics::collide(cell,statistics_);
    for (int iPop=0; iPop < Lattice<T>::q; ++iPop) {
      if ( ((_nx==1 || _nx==-1) && Lattice<T>::c[iPop][0]==_nx) || ((_ny==1 || _ny==-1) && Lattice<T>::c[iPop][1]==_ny) || (Lattice<T>::d==3 && !_nz && Lattice<T>::c[iPop][2]==_nz) ) {
        cell[iPop] += (cell[iPop] + Lattice<T>::t[iPop])*_densityOffset;
      }
    }
};
  /// Collide with fixed velocity
  void staticCollide(Cell<T,Lattice>& cell,
                             const T u[Lattice<T>::d],
                             LatticeStatistics<T>& statistics_) override {
    BaseDynamics::staticCollide(cell,u,statistics_);
  for (int iPop=0; iPop < Lattice<T>::q; ++iPop) {
    if ( ((_nx==1 || _nx==-1) && Lattice<T>::c[iPop][0]==_nx) || ((_ny==1 || _ny==-1) && Lattice<T>::c[iPop][1]==_ny) || (Lattice<T>::d==3 && !_nz && Lattice<T>::c[iPop][2]==_nz) ) {
      cell[iPop] += (cell[iPop] + Lattice<T>::t[iPop])*_densityOffset;
    }
  }
};

private:
  T _densityOffset;
  int _nx, _ny, _nz;
};

/// Implementation of the BGK collision step
template<typename T, template<typename U> class Lattice>
class BGKdynamics : public BasicDynamics<T,Lattice> {
public:
  /// Constructor
  BGKdynamics(T omega, Momenta<T,Lattice>& momenta);
  /// Collision step
  void collide(Cell<T,Lattice>& cell,
                       LatticeStatistics<T>& statistics_) override;
  /// Collide with fixed velocity
  void staticCollide(Cell<T,Lattice>& cell,
                             const T u[Lattice<T>::d],
                             LatticeStatistics<T>& statistics_) override;
  /// Get local relaxation parameter of the dynamics
  T getOmega() const override;
  /// Set local relaxation parameter of the dynamics
  void setOmega(T omega) override;
private:
  T _omega;  ///< relaxation parameter
};

/// Implementation of the TRT collision step
template<typename T, template<typename U> class Lattice>
class TRTdynamics : public BasicDynamics<T,Lattice> {
public:
  /// Constructor
  TRTdynamics(T omega, Momenta<T,Lattice>& momenta, T magicParameter);
  /// Collision step
  void collide(Cell<T,Lattice>& cell,
                       LatticeStatistics<T>& statistics_) override;
  /// Collide with fixed velocity
  void staticCollide(Cell<T,Lattice>& cell,
                             const T u[Lattice<T>::d],
                             LatticeStatistics<T>& statistics_) override;
  /// Get local relaxation parameter of the dynamics
  T getOmega() const override;
  /// Set local relaxation parameter of the dynamics
  void setOmega(T omega) override;
private:
  T _omega;  ///< relaxation parameter
  T _omega2; /// relaxation parameter for odd moments
  T _magicParameter;
};

/// Implementation of the pressure-corrected BGK collision step
template<typename T, template<typename U> class Lattice>
class ConstRhoBGKdynamics : public BasicDynamics<T,Lattice> {
public:
  /// Constructor
  ConstRhoBGKdynamics(T omega, Momenta<T,Lattice>& momenta);
  /// Collision step
  void collide(Cell<T,Lattice>& cell,
                       LatticeStatistics<T>& statistics_) override;
  /// Collide with fixed velocity
  void staticCollide(Cell<T,Lattice>& cell,
                             const T u[Lattice<T>::d],
                             LatticeStatistics<T>& statistics_) override;
  /// Get local relaxation parameter of the dynamics
  T getOmega() const override;
  /// Set local relaxation parameter of the dynamics
  void setOmega(T omega) override;
private:
  T _omega;  ///< relaxation parameter
};

/// Implementation of the so-called incompressible collision step
template<typename T, template<typename U> class Lattice>
class IncBGKdynamics : public BasicDynamics<T,Lattice> {
public:
  /// Constructor
  IncBGKdynamics(T omega, Momenta<T,Lattice>& momenta);
  /// Collision step
  void collide(Cell<T,Lattice>& cell,
                       LatticeStatistics<T>& statistics_) override;
  /// Collide with fixed velocity
  void staticCollide(Cell<T,Lattice>& cell,
                             const T u[Lattice<T>::d],
                             LatticeStatistics<T>& statistics_) override;
  /// Get local relaxation parameter of the dynamics
  T getOmega() const override;
  /// Set local relaxation parameter of the dynamics
  void setOmega(T omega) override;
private:
  T _omega;  ///< relaxation parameter
};



/// Implementation of the Regularized BGK collision step
/** This model is substantially more stable than plain BGK, and has roughly
 * the same efficiency. However, it cuts out the modes at higher Knudsen
 * numbers and can not be used in the regime of rarefied gases.
 */
template<typename T, template<typename U> class Lattice>
class RLBdynamics : public BasicDynamics<T,Lattice> {
public:
  /// Constructor
  RLBdynamics(T omega, Momenta<T,Lattice>& momenta);
  /// Collision step
  void collide(Cell<T,Lattice>& cell,
                       LatticeStatistics<T>& statistics_) override;
  /// Collide with fixed velocity
  void staticCollide(Cell<T,Lattice>& cell,
                             const T u[Lattice<T>::d],
                             LatticeStatistics<T>& statistics_) override;
  /// Get local relaxation parameter of the dynamics
  T getOmega() const override;
  /// Set local relaxation parameter of the dynamics
  void setOmega(T omega) override;
private:
  T _omega;  ///< relaxation parameter
};

/// Implementation of Regularized BGK collision, followed by any Dynamics
template<typename T, template<typename U> class Lattice, typename Dynamics>
class CombinedRLBdynamics : public BasicDynamics<T,Lattice> {
public:
  /// Constructor
  CombinedRLBdynamics(T omega, Momenta<T,Lattice>& momenta);
  /// Compute equilibrium distribution function
  T computeEquilibrium(int iPop, T rho, const T u[Lattice<T>::d], T uSqr) const override;
  /// Collision step
  void collide(Cell<T,Lattice>& cell,
                       LatticeStatistics<T>& statistics_) override;
  /// Collide with fixed velocity
  void staticCollide(Cell<T,Lattice>& cell,
                             const T u[Lattice<T>::d],
                             LatticeStatistics<T>& statistics_) override;
  /// Get local relaxation parameter of the dynamics
  T getOmega() const override;
  /// Set local relaxation parameter of the dynamics
  void setOmega(T omega) override;
private:
  Dynamics _boundaryDynamics;
};

/// Implementation of the BGK collision step with external force
template<typename T, template<typename U> class Lattice>
class ForcedBGKdynamics : public BasicDynamics<T,Lattice> {
public:
  /// Constructor
  ForcedBGKdynamics(T omega, Momenta<T,Lattice>& momenta);
  ///  Compute fluid velocity on the cell.
  void computeU (
    Cell<T,Lattice> const& cell,
    T u[Lattice<T>::d] ) const override;
  /// Compute fluid velocity and particle density on the cell.
  void computeRhoU (
    Cell<T,Lattice> const& cell,
    T& rho, T u[Lattice<T>::d]) const override;
  /// Collision step
  void collide(Cell<T,Lattice>& cell,
                       LatticeStatistics<T>& statistics_) override;
  /// Collide with fixed velocity
  void staticCollide(Cell<T,Lattice>& cell,
                             const T u[Lattice<T>::d],
                             LatticeStatistics<T>& statistics_) override;
  /// Get local relaxation parameter of the dynamics
  T getOmega() const override;
  /// Set local relaxation parameter of the dynamics
  void setOmega(T omega) override;
protected:
  T _omega;  ///< relaxation parameter
  static const int forceBeginsAt = Lattice<T>::ExternalField::forceBeginsAt;
  static const int sizeOfForce   = Lattice<T>::ExternalField::sizeOfForce;
};

/// Implementation of the BGK collision step with external force
template<typename T, template<typename U> class Lattice>
class ForcedKupershtokhBGKdynamics : public BasicDynamics<T,Lattice> {
public:
  /// Constructor
  ForcedKupershtokhBGKdynamics(T omega, Momenta<T,Lattice>& momenta);
  ///  Compute fluid velocity on the cell.
  void computeU (
    Cell<T,Lattice> const& cell,
    T u[Lattice<T>::d] ) const override;
  /// Compute fluid velocity and particle density on the cell.
  void computeRhoU (
    Cell<T,Lattice> const& cell,
    T& rho, T u[Lattice<T>::d]) const override;
  /// Collision step
  void collide(Cell<T,Lattice>& cell,
                       LatticeStatistics<T>& statistics_) override;
  /// Collide with fixed velocity
  void staticCollide(Cell<T,Lattice>& cell,
                             const T u[Lattice<T>::d],
                             LatticeStatistics<T>& statistics_) override;
  /// Get local relaxation parameter of the dynamics
  T getOmega() const override;
  /// Set local relaxation parameter of the dynamics
  void setOmega(T omega) override;
protected:
  T _omega;  ///< relaxation parameter
  static const int forceBeginsAt = Lattice<T>::ExternalField::forceBeginsAt;
  static const int sizeOfForce   = Lattice<T>::ExternalField::sizeOfForce;
};

/// Implementation of the BGK collision step with external force
template<typename T, template<typename U> class Lattice>
class ResettingForcedBGKdynamics : public ForcedBGKdynamics<T,Lattice> {
public:
  ResettingForcedBGKdynamics(T omega, Momenta<T,Lattice>& momenta);
  /// Collision step
  void collide(Cell<T,Lattice>& cell,
                       LatticeStatistics<T>& statistics_) override;
  inline void setForce(T force[3])
  {
//    _frc[0] = force[0];
//    _frc[1] = force[1];
//    _frc[2] = force[2];
    _frc[0] = 0.0;
    _frc[1] = 0.0;
    _frc[2] = 0.0;
  }
private:
  T _frc[3];
};

/// Other Implementation of the BGK collision step with external force
template<typename T, template<typename U> class Lattice>
class ForcedShanChenBGKdynamics : public ForcedBGKdynamics<T,Lattice> {
public:
  /// Constructor
  ForcedShanChenBGKdynamics(T omega, Momenta<T,Lattice>& momenta);
  ///  Compute fluid velocity on the cell.
  void computeU (
    Cell<T,Lattice> const& cell,
    T u[Lattice<T>::d] ) const override;
  /// Compute fluid velocity and particle density on the cell.
  void computeRhoU (
    Cell<T,Lattice> const& cell,
    T& rho, T u[Lattice<T>::d]) const override;
  /// Collision step
  void collide(Cell<T,Lattice>& cell,
                       LatticeStatistics<T>& statistics_) override;
};

/// Implementation of the 3D D3Q13 dynamics
/** This is (so far) the minimal existing 3D model, with only 13
 * directions. Three different relaxation times are used to achieve
 * asymptotic hydrodynamics, isotropy and galilean invariance.
 */
template<typename T, template<typename U> class Lattice>
class D3Q13dynamics : public BasicDynamics<T,Lattice> {
public:
  /// Constructor
  D3Q13dynamics(T omega, Momenta<T,Lattice>& momenta);
  /// Compute equilibrium distribution function
  T computeEquilibrium(int iPop, T rho, const T u[Lattice<T>::d], T uSqr) const override;
  /// Collision step
  void collide(Cell<T,Lattice>& cell,
                       LatticeStatistics<T>& statistics_) override;
  /// Collide with fixed velocity
  void staticCollide(Cell<T,Lattice>& cell,
                             const T u[Lattice<T>::d],
                             LatticeStatistics<T>& statistics_) override;
  /// Get local relaxation parameter of the dynamics
  T getOmega() const override;
  /// Set local relaxation parameter of the dynamics
  void setOmega(T omega) override;
private:
  T lambda_nu;        ///< first relaxation parameter
  T lambda_nu_prime;  ///< second relaxation parameter
};

/// Standard computation of velocity momenta in the bulk
template<typename T, template<typename U> class Lattice>
struct BulkMomenta : public Momenta<T,Lattice> {
  /// Compute particle density on the cell.
  T computeRho(Cell<T,Lattice> const& cell) const override;
  /// Compute fluid velocity on the cell.
  void computeU (
    Cell<T,Lattice> const& cell,
    T u[Lattice<T>::d] ) const override;
  /// Compute fluid momentum on the cell.
  void computeJ (
    Cell<T,Lattice> const& cell,
    T j[Lattice<T>::d] ) const override;
  /// Compute components of the stress tensor on the cell.
  void computeStress (
    Cell<T,Lattice> const& cell,
    T rho, const T u[Lattice<T>::d],
    T pi[util::TensorVal<Lattice<T> >::n] ) const override;
  /// Compute fluid velocity and particle density on the cell.
  void computeRhoU (
    Cell<T,Lattice> const& cell,
    T& rho, T u[Lattice<T>::d]) const override;
  /// Compute all momenta on the cell, up to second order.
  void computeAllMomenta (
    Cell<T,Lattice> const& cell,
    T& rho, T u[Lattice<T>::d],
    T pi[util::TensorVal<Lattice<T> >::n] ) const override;
  /// Set particle density on the cell.
  void defineRho(Cell<T,Lattice>& cell, T rho) override;
  /// Set fluid velocity on the cell.
  void defineU(Cell<T,Lattice>& cell,
                       const T u[Lattice<T>::d]) override;
  /// Define fluid velocity and particle density on the cell.
  void defineRhoU (
    Cell<T,Lattice>& cell,
    T rho, const T u[Lattice<T>::d]) override;
  /// Define all momenta on the cell, up to second order.
  void defineAllMomenta (
    Cell<T,Lattice>& cell,
    T rho, const T u[Lattice<T>::d],
    const T pi[util::TensorVal<Lattice<T> >::n] ) override;
};

/// Velocity is stored in external scalar (and computed e.g. in a PostProcessor)
template<typename T, template<typename U> class Lattice>
struct ExternalVelocityMomenta : public Momenta<T,Lattice> {
  /// Compute particle density on the cell.
  T computeRho(Cell<T,Lattice> const& cell) const override;
  /// Compute fluid velocity on the cell.
  void computeU (
    Cell<T,Lattice> const& cell,
    T u[Lattice<T>::d] ) const override;
  /// Compute fluid momentum on the cell.
  void computeJ (
    Cell<T,Lattice> const& cell,
    T j[Lattice<T>::d] ) const override;
  /// Compute components of the stress tensor on the cell.
  void computeStress (
    Cell<T,Lattice> const& cell,
    T rho, const T u[Lattice<T>::d],
    T pi[util::TensorVal<Lattice<T> >::n] ) const override;
  /// Compute fluid velocity and particle density on the cell.
  void computeRhoU (
    Cell<T,Lattice> const& cell,
    T& rho, T u[Lattice<T>::d]) const override;
  /// Compute all momenta on the cell, up to second order.
  void computeAllMomenta (
    Cell<T,Lattice> const& cell,
    T& rho, T u[Lattice<T>::d],
    T pi[util::TensorVal<Lattice<T> >::n] ) const override;
  /// Set particle density on the cell.
  void defineRho(Cell<T,Lattice>& cell, T rho) override;
  /// Set fluid velocity on the cell.
  void defineU(Cell<T,Lattice>& cell,
                       const T u[Lattice<T>::d]) override;
  /// Define fluid velocity and particle density on the cell.
  void defineRhoU (
    Cell<T,Lattice>& cell,
    T rho, const T u[Lattice<T>::d]) override;
  /// Define all momenta on the cell, up to second order.
  void defineAllMomenta (
    Cell<T,Lattice>& cell,
    T rho, const T u[Lattice<T>::d],
    const T pi[util::TensorVal<Lattice<T> >::n] ) override;
};

/// Implementation of "bounce-back" dynamics
/** This is a very popular way to implement no-slip boundary conditions,
 * because the dynamics are independent of the orientation of the boundary.
 * It is a special case, because it implements no usual LB dynamics.
 * For that reason, it derives directly from the class Dynamics.
 *
 * The code works for both 2D and 3D lattices.
 */
template<typename T, template<typename U> class Lattice>
class BounceBack : public Dynamics<T,Lattice> {
public:
  /// A fictitious density value on bounce-back in not fixed on nodes via this constructor.
  BounceBack();
  /// You may fix a fictitious density value on bounce-back nodes via this constructor.
  BounceBack(T rho);
  /// Collision step
  void collide(Cell<T,Lattice>& cell,
                       LatticeStatistics<T>& statistics_) override;
  /// Collide with fixed velocity
  void staticCollide(Cell<T,Lattice>& cell,
                             const T u[Lattice<T>::d],
                             LatticeStatistics<T>& statistics_) override;
  /// Yields 1;
  T computeRho(Cell<T,Lattice> const& cell) const override;
  /// Yields 0;
  void computeU (
    Cell<T,Lattice> const& cell,
    T u[Lattice<T>::d] ) const override;
  /// Yields 0;
  void computeJ (
    Cell<T,Lattice> const& cell,
    T j[Lattice<T>::d] ) const override;
  /// Yields NaN
  void computeStress (
    Cell<T,Lattice> const& cell,
    T rho, const T u[Lattice<T>::d],
    T pi[util::TensorVal<Lattice<T> >::n] ) const override;
  void computeRhoU (
    Cell<T,Lattice> const& cell,
    T& rho, T u[Lattice<T>::d]) const override;
  void computeAllMomenta (
    Cell<T,Lattice> const& cell,
    T& rho, T u[Lattice<T>::d],
    T pi[util::TensorVal<Lattice<T> >::n] ) const override;
  /// Does nothing
  void defineRho(Cell<T,Lattice>& cell, T rho) override;
  /// Does nothing
  void defineU(Cell<T,Lattice>& cell,
                       const T u[Lattice<T>::d]) override;
  /// Does nothing
  void defineRhoU (
    Cell<T,Lattice>& cell,
    T rho, const T u[Lattice<T>::d]) override;
  /// Does nothing
  void defineAllMomenta (
    Cell<T,Lattice>& cell,
    T rho, const T u[Lattice<T>::d],
    const T pi[util::TensorVal<Lattice<T> >::n] ) override;
  /// Yields NaN
  T getOmega() const override;
  /// Does nothing
  void setOmega(T omega) override;
private:
  T _rho;
  bool _rhoFixed;
};

/// Implementation of "bounce-back velocity" dynamics
/** This is a very popular way to implement no-slip boundary conditions,
 * because the dynamics are independent of the orientation of the boundary.
 * It is a special case, because it implements no usual LB dynamics.
 * For that reason, it derives directly from the class Dynamics. It
 * fixes the velociy to a given velocity _u.
 *
 * The code works for both 2D and 3D lattices.
 */
template<typename T, template<typename U> class Lattice>
class BounceBackVelocity : public Dynamics<T,Lattice> {
public:
  /// A fictitious density value on bounce-back in not fixed on nodes via this constructor.
  BounceBackVelocity(const T u[Lattice<T>::d]);
  /// You may fix a fictitious density value on bounce-back nodes via this constructor.
  BounceBackVelocity(const T rho, const T u[Lattice<T>::d]);
  /// Collision step, bounce back with a fixed velocity _u
  void collide(Cell<T,Lattice>& cell,
                       LatticeStatistics<T>& statistics_) override;
  /// Collide with fixed velocity
  void staticCollide(Cell<T,Lattice>& cell,
                             const T u[Lattice<T>::d],
                             LatticeStatistics<T>& statistics_) override;
  /// Retuns rho (if defined else zero)
  T computeRho(Cell<T,Lattice> const& cell) const override;
  /// Retuns _u
  void computeU (
    Cell<T,Lattice> const& cell,
    T u[Lattice<T>::d] ) const override;
  /// Retuns rho (if defined else zero) times _u
  void computeJ (
    Cell<T,Lattice> const& cell,
    T j[Lattice<T>::d] ) const override;
  /// Yields NaN
  void computeStress (
    Cell<T,Lattice> const& cell,
    T rho, const T u[Lattice<T>::d],
    T pi[util::TensorVal<Lattice<T> >::n] ) const override;
  /// Retuns rho (if defined else zero) and _u
  void computeRhoU (
    Cell<T,Lattice> const& cell,
    T& rho, T u[Lattice<T>::d]) const override;
  void computeAllMomenta (
    Cell<T,Lattice> const& cell,
    T& rho, T u[Lattice<T>::d],
    T pi[util::TensorVal<Lattice<T> >::n] ) const override;
  /// Devines the velocity rho
  void defineRho(Cell<T,Lattice>& cell, T rho) override;
  /// Devines the velocity _u
  void defineU(Cell<T,Lattice>& cell,
                       const T u[Lattice<T>::d]) override;
  /// Devines rho and _u
  void defineRhoU (
    Cell<T,Lattice>& cell,
    T rho, const T u[Lattice<T>::d]) override;
  /// Does nothing
  void defineAllMomenta (
    Cell<T,Lattice>& cell,
    T rho, const T u[Lattice<T>::d],
    const T pi[util::TensorVal<Lattice<T> >::n] ) override;
  /// Yields NaN
  T getOmega() const override;
  /// Does nothing
  void setOmega(T omega) override;
private:
  T _rho;
  bool _rhoFixed;
  T _u[Lattice<T>::d];
};

/// Implementation of "bounce-back anti" dynamics
/** This is a way to implement a Dirichlet rho/pressure boundary conditions,
 * because the dynamics are independent of the orientation of the boundary.
 * It is a special case, because it implements no usual LB dynamics.
 * For that reason, it derives directly from the class Dynamics. It
 * fixes the rho to a given _rho.
 *
 * The code works for both 2D and 3D lattices.
 */
template<typename T, template<typename U> class Lattice>
class BounceBackAnti : public Dynamics<T,Lattice> {
public:
  /// A fictitious density value on bounce-back in not fixed on nodes via this constructor.
  BounceBackAnti();
  /// You may fix a fictitious density value on bounce-back nodes via this constructor.
  BounceBackAnti(T rho);
  /// Collision step, bounce back with a fixed velocity _u
  void collide(Cell<T,Lattice>& cell,
                       LatticeStatistics<T>& statistics_) override;
  /// Collide with fixed velocity
  void staticCollide(Cell<T,Lattice>& cell,
                             const T u[Lattice<T>::d],
                             LatticeStatistics<T>& statistics_) override;
  /// Retuns rho (if defined else zero)
  T computeRho(Cell<T,Lattice> const& cell) const override;
  /// Retuns _u
  void computeU (
    Cell<T,Lattice> const& cell,
    T u[Lattice<T>::d] ) const override;
  /// Retuns rho (if defined else zero) times _u
  void computeJ (
    Cell<T,Lattice> const& cell,
    T j[Lattice<T>::d] ) const override;
  /// Yields NaN
  void computeStress (
    Cell<T,Lattice> const& cell,
    T rho, const T u[Lattice<T>::d],
    T pi[util::TensorVal<Lattice<T> >::n] ) const override;
  /// Retuns rho (if defined else zero) and _u
  void computeRhoU (
    Cell<T,Lattice> const& cell,
    T& rho, T u[Lattice<T>::d]) const override;
  void computeAllMomenta (
    Cell<T,Lattice> const& cell,
    T& rho, T u[Lattice<T>::d],
    T pi[util::TensorVal<Lattice<T> >::n] ) const override;
  /// Devines the velocity rho
  void defineRho(Cell<T,Lattice>& cell, T rho) override;
  /// Devines the velocity _u
  void defineU(Cell<T,Lattice>& cell,
                       const T u[Lattice<T>::d]) override;
  /// Devines rho and _u
  void defineRhoU (
    Cell<T,Lattice>& cell,
    T rho, const T u[Lattice<T>::d]) override;
  /// Does nothing
  void defineAllMomenta (
    Cell<T,Lattice>& cell,
    T rho, const T u[Lattice<T>::d],
    const T pi[util::TensorVal<Lattice<T> >::n] ) override;
  /// Yields NaN
  T getOmega() const override;
  /// Does nothing
  void setOmega(T omega) override;
private:
  T _rho;
  bool _rhoFixed;
  T _u[Lattice<T>::d];
};


/** Corresponds to macro Robin boundary, micro Fresnel surface
 *  Motivated by Hiorth et al. 2008; doi 10.1002/fld.1822
 */
template<typename T, template<typename U> class Lattice>
class PartialBounceBack final: public BounceBack<T,Lattice> {
public:
  PartialBounceBack(T rf);
  T computeEquilibrium( int iPop, T rho, const T u[Lattice<T>::d], T uSqr ) const override;
  /// Collision step
  void collide(Cell<T,Lattice>& cell, LatticeStatistics<T>& statistics_) override;
private:
  T _rf;
};


/// Implementation of a "dead cell" that does nothing
template<typename T, template<typename U> class Lattice>
class NoDynamics : public Dynamics<T,Lattice> {
public:
  /// You may fix a fictitious density value on no dynamics node via this constructor.
  NoDynamics(T rho = T(1) );
  /// Yields 0;
  T computeEquilibrium(int iPop, T rho, const T u[Lattice<T>::d], T uSqr) const override;
  /// Collision step
  void collide(Cell<T,Lattice>& cell,
                       LatticeStatistics<T>& statistics_) override;
  /// Collide with fixed velocity
  void staticCollide(Cell<T,Lattice>& cell,
                             const T u[Lattice<T>::d],
                             LatticeStatistics<T>& statistics_) override;
  /// Yields 1;
  T computeRho(Cell<T,Lattice> const& cell) const override;
  /// Yields 0;
  void computeU (
    Cell<T,Lattice> const& cell,
    T u[Lattice<T>::d] ) const override;
  /// Yields 0;
  void computeJ (
    Cell<T,Lattice> const& cell,
    T j[Lattice<T>::d] ) const override;
  /// Yields NaN
  void computeStress (
    Cell<T,Lattice> const& cell,
    T rho, const T u[Lattice<T>::d],
    T pi[util::TensorVal<Lattice<T> >::n] ) const override;
  void computeRhoU (
    Cell<T,Lattice> const& cell,
    T& rho, T u[Lattice<T>::d]) const override;
  void computeAllMomenta (
    Cell<T,Lattice> const& cell,
    T& rho, T u[Lattice<T>::d],
    T pi[util::TensorVal<Lattice<T> >::n] ) const override;
  /// Does nothing
  void defineRho(Cell<T,Lattice>& cell, T rho) override;
  /// Does nothing
  void defineU(Cell<T,Lattice>& cell,
                       const T u[Lattice<T>::d]) override;
  /// Does nothing
  void defineRhoU (
    Cell<T,Lattice>& cell,
    T rho, const T u[Lattice<T>::d]) override;
  /// Does nothing
  void defineAllMomenta (
    Cell<T,Lattice>& cell,
    T rho, const T u[Lattice<T>::d],
    const T pi[util::TensorVal<Lattice<T> >::n] ) override;
  /// Yields NaN
  T getOmega() const override;
  /// Does nothing
  void setOmega(T omega) override;

private:
  /// Default rho=1
  T _rho;
};

/// Dynamics for offLattice boundary conditions
/// OffDynamics are basically NoDynamics with the additional functionality
/// to store given velocities exactly at boundary links.
template<typename T, template<typename U> class Lattice>
class OffDynamics : public NoDynamics<T,Lattice> {
public:
  /// Constructor
  OffDynamics(const T _location[Lattice<T>::d]);
  /// Constructor
  OffDynamics(const T _location[Lattice<T>::d], T _distances[Lattice<T>::q]);
  /// Returns local stored rho which is updated if the bc is used as velocity!=0 condition
  T computeRho(Cell<T,Lattice> const& cell) const override;
  /// Returns an average of the locally stored u
  void computeU(Cell<T,Lattice> const& cell, T u[Lattice<T>::d] ) const override;
  /// Set Intersection of the link and the boundary
  void setBoundaryIntersection(int iPop, T distance) override;
  /// Get Intersection of the link and the boundary
  bool getBoundaryIntersection(int iPop, T intersection[Lattice<T>::d]) override;
  /// Set particle density on the cell.
  void defineRho(Cell<T,Lattice>& cell, T rho) override;
  /// Set single velocity
  void defineRho(int iPop, T rho) override;
  /// Set fluid velocity on the cell.
  void defineU(Cell<T,Lattice>& cell, const T u[Lattice<T>::d]) override;
  /// Set constant velocity
  void defineU(const T u[Lattice<T>::d]) override;
  /// Set single velocity
  void defineU(int iPop, const T u[Lattice<T>::d]) override;
  /// Get VelocitySummand for Bouzidi-Boundary Condition
  T getVelocityCoefficient(int iPop) override;

private:
  T _rho;
  T _u[Lattice<T>::q][Lattice<T>::d];
  T location[Lattice<T>::d];
  T distances[Lattice<T>::q];
  T boundaryIntersection[Lattice<T>::q][Lattice<T>::d];
  T velocityCoefficient[Lattice<T>::q];
};

/// Implementation of density sink by setting a zero distribution on the cell
template<typename T, template<typename U> class Lattice>
class ZeroDistributionDynamics : public NoDynamics<T,Lattice> {
public:
  /// Constructor.
  ZeroDistributionDynamics();
  /// Collision step
  void collide(Cell<T,Lattice>& cell,
                       LatticeStatistics<T>& statistics_) override;
  /// Yields 1
  T computeRho(Cell<T,Lattice> const& cell) const override;
};


namespace instances {

template<typename T, template<typename U> class Lattice>
BulkMomenta<T,Lattice>& getBulkMomenta();

template<typename T, template<typename U> class Lattice>
ExternalVelocityMomenta<T,Lattice>& getExternalVelocityMomenta();

template<typename T, template<typename U> class Lattice>
BounceBack<T,Lattice>& getBounceBack();

template<typename T, template<typename U> class Lattice>
PartialBounceBack<T,Lattice>& getPartialBounceBack(const double rf);

template<typename T, template<typename U> class Lattice>
BounceBackVelocity<T,Lattice>& getBounceBackVelocity(const double rho, const double u[Lattice<T>::d]);

template<typename T, template<typename U> class Lattice>
BounceBackAnti<T,Lattice>& getBounceBackAnti(const double rho);

template<typename T, template<typename U> class Lattice>
NoDynamics<T,Lattice>& getNoDynamics(T rho = T(1) );

template<typename T, template<typename U> class Lattice>
ZeroDistributionDynamics<T,Lattice>& getZeroDistributionDynamics();

}

}

#endif
