/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Asher Zarth, Lukas Baron, Mathias J. Krause, Jonas Latt, Jan E. Marquardt
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
 * BGK Dynamics for porous media -- header file.
 */
#ifndef POROUS_FORCED_BGK_DYNAMICS_H
#define POROUS_FORCED_BGK_DYNAMICS_H

#include "dynamics/dynamics.h"
#include "core/cell.h"

namespace olb {

/// Implementation of the BGK collision step for a porosity model
template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
class PorousForcedBGKdynamics : public legacy::BasicDynamics<T,DESCRIPTOR,MOMENTA> {
public:
  template<typename M>
  using exchange_momenta = PorousForcedBGKdynamics<T,DESCRIPTOR,M>;

  /// Constructor
  PorousForcedBGKdynamics(T omega_);
  void computeU (
    ConstCell<T,DESCRIPTOR>& cell,
    T u[DESCRIPTOR::d] ) const override;
  /// Compute fluid velocity and particle density on the cell.
  void computeRhoU (
    ConstCell<T,DESCRIPTOR>& cell,
    T& rho, T u[DESCRIPTOR::d]) const override;

  /// Collision step
  CellStatistic<T> collide(Cell<T,DESCRIPTOR>& cell);

  /// get relaxation parameter
  T    getOmega() const;
  /// set relaxation parameter
  void setOmega(T omega_);


private:
  using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;

  T _omega;      ///< relaxation parameter
};


/* Implementation of the BGK collision for moving porous media (HLBM approach).
 * As this scheme requires additionla data stored in an external field,
 * it is meant to be used along with a PorousParticle descriptor.
 * \param omega Lattice relaxation frequency
 * \param momenta A standard object for the momenta computation
 */
template<typename T, typename DESCRIPTOR, typename MOMENTA, bool isStatic=false>
class PorousParticleForcedBGKdynamics : public BGKdynamics<T,DESCRIPTOR,MOMENTA>, public PorousParticleDynamics<T,DESCRIPTOR,isStatic> {
public:
  template<typename M>
  using exchange_momenta = PorousParticleForcedBGKdynamics<T,DESCRIPTOR,M,isStatic>;

  /// Constructor
  PorousParticleForcedBGKdynamics(T omega_);
  /// Compute fluid velocity on the cell.
  void computeU ( ConstCell<T,DESCRIPTOR>& cell, T u[DESCRIPTOR::d] ) const override;
  /// Compute fluid velocity and particle density on the cell.
  void computeRhoU ( ConstCell<T,DESCRIPTOR>& cell,
                     T& rho, T u[DESCRIPTOR::d]) const override;
  /// extended Collision step, computes local drag in a given direction
  CellStatistic<T> collide(Cell<T,DESCRIPTOR>& cell,
               LatticeStatistics<T>& statistics_) override;
protected:
  T porousParticleBgkCollision(Cell<T,DESCRIPTOR>& cell, T rho, T u[DESCRIPTOR::d], T omega);

private:
  T KupershtokhForcing(Cell<T,DESCRIPTOR>& cell, const T rho, const T u[DESCRIPTOR::d], const T omega);
  T GuoForcing(Cell<T,DESCRIPTOR>& cell, const T rho, T u[DESCRIPTOR::d], const T omega);
  T ShanChenForcing(Cell<T,DESCRIPTOR>& cell, const T rho, T u[DESCRIPTOR::d], const T omega);
  T computeUSqr(Cell<T,DESCRIPTOR>& cell);
};

} // olb

#endif
