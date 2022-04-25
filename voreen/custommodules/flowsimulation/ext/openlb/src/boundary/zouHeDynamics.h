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

#ifndef ZOU_HE_DYNAMICS_H
#define ZOU_HE_DYNAMICS_H

#include "dynamics/dynamics.h"

namespace olb {

/**
* Implementation of Zou-He boundary condition following
* the paper from Zou and He. This implementation is lattice independent.
* The implementation follow the idea proposed int the paper
* Qisu Zou, Xiaoyi He,
* "On pressure and velocity boundary conditions for the lattice Boltzmann BGK model",
* Phys. Fluids , (1997), Volume 9, Issue 6, pp. 1591-1598
*/
template<typename T, typename DESCRIPTOR, typename Dynamics, typename MOMENTA, int direction, int orientation>
class ZouHeDynamics : public legacy::BasicDynamics<T,DESCRIPTOR,MOMENTA> {
public:
  template<typename M>
  using exchange_momenta = ZouHeDynamics<T,DESCRIPTOR,Dynamics,M,direction,orientation>;

  /// Constructor
  ZouHeDynamics(T omega_);
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
