/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Orestis Malaspinas, Jonas Latt
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

#ifndef INAMURO_NEWTON_RAPHSON_DYNAMICS_H
#define INAMURO_NEWTON_RAPHSON_DYNAMICS_H

#include "dynamics/dynamics.h"
#include "io/ostreamManager.h"

namespace olb {

/**
* This class computes the inamuro BC with general dynamics. It uses the formula from the
 * paper by Inamuro et al. but since there is no explict solution
 * for a lattice different from the D2Q9 and for a speed of sound
 * c_s=q/sqrt(3), we have to use a Newton-Raphson algorithm to
 * implement these boundary conditions.
*/
template<typename T, template<typename U> class Lattice, typename Dynamics, int direction, int orientation>
class InamuroNewtonRaphsonDynamics : public BasicDynamics<T,Lattice> {
public:
  /// Constructor
  InamuroNewtonRaphsonDynamics(T omega, Momenta<T,Lattice>& momenta);
  /// Compute equilibrium distribution function
  T computeEquilibrium(int iPop, T rho, const T u[Lattice<T>::d], T uSqr) const override;
  /// Collision step
  void collide(Cell<T,Lattice>& cell,
                       LatticeStatistics<T>& statistics) override;
  /// Collide with fixed velocity
  void staticCollide(Cell<T,Lattice>& cell,
                             const T u[Lattice<T>::d],
                             LatticeStatistics<T>& statistics) override;
  /// Get local relaxation parameter of the dynamics
  T getOmega() const override;
  /// Set local relaxation parameter of the dynamics
  void setOmega(T omega) override;
private:
  void computeApproxMomentum(T approxMomentum[Lattice<T>::d],
                             const Cell<T,Lattice> &cell,
                             const T &rho, const T u[Lattice<T>::d], const T xi[Lattice<T>::d],
                             const std::vector<int> knownIndexes,const std::vector<int> missingIndexes);

  /// compute the error (L^2 norm of (u-uApprox))
  T computeError(const T &rho,const T u[Lattice<T>::d], const T approxMomentum[Lattice<T>::d]);

  void computeGradGradError(T gradGradError[Lattice<T>::d][Lattice<T>::d],
                            T gradError[Lattice<T>::d],
                            const T &rho, const T u[Lattice<T>::d],const T xi[Lattice<T>::d],
                            const T approxMomentum[Lattice<T>::d],
                            const std::vector<int> missingIndexes);

  /// compute the new xi with the newton raphson algorithm
  bool newtonRaphson(T xi[Lattice<T>::d],
                     const T gradError[Lattice<T>::d],
                     const T gradGradError[Lattice<T>::d][Lattice<T>::d]);

  bool invert(const T a[2][2],T b[2][2]);

  bool invert(const T a[3][3],T b[3][3]);

  Dynamics _boundaryDynamics;
  T _xi[Lattice<T>::d];
  mutable OstreamManager clout;
};

}

#endif
