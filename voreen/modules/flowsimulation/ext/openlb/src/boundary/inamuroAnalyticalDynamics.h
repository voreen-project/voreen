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

#ifndef INAMURO_ANALYTICAL_DYNAMICS_H
#define INAMURO_ANALYTICAL_DYNAMICS_H

#include "dynamics/dynamics.h"

namespace olb {

/**
* Implementation of Inamuro boundary condition following
 * the paper
 * "A non-slip boundary condition for lattice Boltzmann simulations",
 * Inamuro, Takaji; Yoshino, Masato; Ogino, Fumimaru, (1995).
 * This implementation works for the D2Q9 DESCRIPTOR only.
*/
template<typename T, typename DESCRIPTOR, typename Dynamics, typename MOMENTA, int direction, int orientation>
class InamuroAnalyticalDynamics : public legacy::BasicDynamics<T,DESCRIPTOR,MOMENTA> {
public:

  template<typename M>
  using exchange_momenta = InamuroAnalyticalDynamics<T,DESCRIPTOR,Dynamics,M,direction,orientation>;
  
  /// Constructor
  InamuroAnalyticalDynamics(T omega);
  /// Compute equilibrium distribution function
  T computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d], T uSqr) const override;
  /// Collision step
  CellStatistic<T> collide(Cell<T,DESCRIPTOR>& cell, LatticeStatistics<T>& statistics) override;
  /// Get local relaxation parameter of the dynamics
  T getOmega() const;
  /// Set local relaxation parameter of the dynamics
  void setOmega(T omega_);
private:
  Dynamics boundaryDynamics;
};

}

#endif
