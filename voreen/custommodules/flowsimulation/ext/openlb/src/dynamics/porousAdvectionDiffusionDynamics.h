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
#ifndef POROUS_ADVECTION_DIFFUSION_DYNAMICS_H
#define POROUS_ADVECTION_DIFFUSION_DYNAMICS_H

#include "dynamics/dynamics.h"
#include "core/cell.h"

namespace olb {


// ===== the porous BGK advection diffusion dynamics =====//
/// This approach contains a slight error in the diffusion term.
template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::AdvectionDiffusionBulkTuple>
class PorousAdvectionDiffusionBGKdynamics : public legacy::BasicDynamics<T, DESCRIPTOR, MOMENTA> {
public:
  template<typename M>
  using exchange_momenta = PorousAdvectionDiffusionBGKdynamics<T,DESCRIPTOR,M>;

  /// Constructor
  PorousAdvectionDiffusionBGKdynamics( T omega, T tSolid );
  PorousAdvectionDiffusionBGKdynamics( const UnitConverter<T,DESCRIPTOR>& converter, T tSolid );
  /// Compute equilibrium distribution function
  T computeEquilibrium( int iPop, T rho, const T u[DESCRIPTOR::d] ) const override;
  /// Collision step
  CellStatistic<T> collide( Cell<T,DESCRIPTOR>& cell ) override;
  /// Get local relaxation parameter of the dynamics
  T getOmega() const;
  /// Set local relaxation parameter of the dynamics
  void setOmega( T omega );
private:
  using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;

  T scaleTemp(const T rho, const T porosity) const; // scales temperature relative to porosity
  T _omega;  ///< relaxation parameter
  T _tSolid; // temperature in lattice units of material with porosity 0;
};

} // namespace olb

#endif
