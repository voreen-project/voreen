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
 * This object is a MRT LB dynamics as described in D.Yu et al. in
 * Progress in Aerospace Sciences 39 (2003) 329-367
 */
#ifndef MRT_DYNAMICS_H
#define MRT_DYNAMICS_H

#include "dynamics/dynamics.h"

namespace olb {

template<typename T, template<typename U> class Lattice> class Cell;


/// Implementation of the entropic collision step
template<typename T, template<typename U> class Lattice>
class MRTdynamics : public BasicDynamics<T,Lattice> {
public:
  /// Constructor
  MRTdynamics(T omega_, Momenta<T,Lattice>& momenta_);
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
  void setOmega(T omega_) override;
  /// Get local relaxation parameter of the dynamics
  T getLambda() const;
  /// Set local relaxation parameter of the dynamics
  void setLambda(T lambda_);
protected:
  T invM_S[Lattice<T>::q][Lattice<T>::q]; // relaxation times matrix.
  T omega; // the shear viscosity relaxatin time
  T lambda;// the bulk viscosity relaxatin time
};

/// Implementation of the entropic collision step
template<typename T, template<typename U> class Lattice>
class ForcedMRTdynamics : public MRTdynamics<T,Lattice> {
public:
  /// Constructor
  ForcedMRTdynamics(T omega_, Momenta<T,Lattice>& momenta_);
  /// Clone the object on its dynamic type.
  virtual void collide(Cell<T,Lattice>& cell,
                       LatticeStatistics<T>& statistics_);

};

/// Implementation of the entropic collision step
template<typename T, template<typename U> class Lattice>
class MRTdynamics2 : public MRTdynamics<T,Lattice> {
public:
  /// Constructor
  MRTdynamics2(T omega_, Momenta<T,Lattice>& momenta_);
  /// Clone the object on its dynamic type.
  virtual void collide(Cell<T,Lattice>& cell,
                       LatticeStatistics<T>& statistics_);

protected:
  T invM_S_2[Lattice<T>::q][Lattice<T>::q]; // relaxation times matrix.
  T omega;
};



}

#endif
