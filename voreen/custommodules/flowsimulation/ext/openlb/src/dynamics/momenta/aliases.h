/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006 Jonas Latt, 2021 Julius Jessberger
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
 * Instantiation of Momenta tuples which define the computation and definition
 * of momenta in a specific situation, e.g. the standard computation from
 * the populations in the bulk.
 */

#ifndef DYNAMICS_MOMENTA_ALIASES_H
#define DYNAMICS_MOMENTA_ALIASES_H

#include "interface.h"
#include "definitionRule.h"

namespace olb {

namespace momenta {

/// Standard computation of momenta from the populations in the bulk
using BulkTuple = Tuple<
  BulkDensity,
  BulkMomentum,
  BulkStress,
  DefineToNEq
>;

/// The Velocity is stored in descriptors::VELOCITY (and computed e.g. in a
/// postprocessor)
using ExternalVelocityTuple = Tuple<
  BulkDensity,
  FixedVelocityMomentum,
  BulkStress,
  DefineUSeparately
>;

/** Collection of momenta for advection diffusion situations
 * rho describes the temperature/ particle concentration, while
 * j describes the convective transport and u is the transport rate.
 * stress does not occur in this model and must not be computed.
 */
using AdvectionDiffusionBulkTuple = Tuple<
  BulkDensity,
  BulkMomentum,
  NoStress,
  DefineToEq
>;

/** Momenta collection for advection diffusion boundary with fixed temperature
 * rho describes the temperature, while j describes the (heat) transport
 * (= flux + convective transport) and u is the conduction/ flux.
 * stress does not occur in this model and must not be computed.
 */
// TODO: computing transport and conduction via compute and computeU, resp. is
// inconsequent.
template <int direction, int orientation>
using RegularizedTemperatureBoundaryTuple = Tuple<
  FixedDensity,
  FixedTemperatureMomentum<direction,orientation>,
  NoStress,
  DefineSeparately
>;

/** Momenta collection for advection diffusion boundary with fixed heat flux
 * rho describes the temperature, while j describes the (heat) transport
 * (= conduction + convective transport) and u is the conduction/ flux.
 * stress does not occur in this model and must not be computed.
 */
// TODO: computing transport and conduction via compute and computeU, resp. is
// inconsequent.
template <int direction, int orientation>
using RegularizedHeatFluxBoundaryTuple = Tuple<
  HeatFluxBoundaryDensity<direction,orientation>,
  FixedVelocityMomentumAD,
  NoStress,
  DefineToEq
>;

/// Velocity and density are stored in external fields
using EquilibriumBoundaryTuple = Tuple<
  FixedDensity,
  FixedVelocityMomentumGeneric,
  ZeroStress,
  DefineSeparately
>;

/// Regularized velocity boundary node
template <int direction, int orientation>
using RegularizedVelocityBoundaryTuple = Tuple<
  VelocityBoundaryDensity<direction,orientation>,
  FixedVelocityMomentumGeneric,
  RegularizedBoundaryStress<direction,orientation>,
  DefineSeparately
>;

/// Regularized pressure boundary node
template <int direction, int orientation>
using RegularizedPressureBoundaryTuple = Tuple<
  FixedDensity,
  FixedPressureMomentum<direction,orientation>,
  RegularizedBoundaryStress<direction,orientation>,
  DefineSeparately
>;

/// Velocity boundary node. Density and velocity are computed via boundary
/// condition, whereas stress is computed as in the bulk.
template <int direction, int orientation>
using BasicDirichletVelocityBoundaryTuple = Tuple<
  VelocityBoundaryDensity<direction,orientation>,
  FixedVelocityMomentumGeneric,
  BulkStress,
  DefineSeparately
>;

/// Pressure boundary node. Density and velocity are computed via boundary
/// condition, whereas stress is computed as in the bulk.
template <int direction, int orientation>
using BasicDirichletPressureBoundaryTuple = Tuple<
  FixedDensity,
  FixedPressureMomentum<direction,orientation>,
  BulkStress,
  DefineSeparately
>;

/** In this class, the velocity is fixed
 * As opposed to VelocityBM, the pressure is however not
 * computed from a special trick on the boundary, but the
 * same way it would be in the bulk.
 */
using FixedVelocityBoundaryTuple = Tuple<
  BulkDensity,
  FixedVelocityMomentumGeneric,
  BulkStress,
  DefineUSeparatelyTrace
>;

template <int... Normal>
using InnerCornerVelocityTuple2D = Tuple<
  InnerCornerDensity2D<Normal...>,
  FixedVelocityMomentumGeneric,
  InnerCornerStress2D<Normal...>,
  DefineSeparately
>;

template <int... PlaneAndNormal>
using InnerEdgeVelocityTuple3D = Tuple<
  InnerEdgeDensity3D<PlaneAndNormal...>,
  FixedVelocityMomentumGeneric,
  InnerEdgeStress3D<PlaneAndNormal...>,
  DefineSeparately
>;

template<int... Normal>
using InnerCornerVelocityTuple3D = Tuple<
  InnerCornerDensity3D<Normal...>,
  FixedVelocityMomentumGeneric,
  InnerCornerStress3D<Normal...>,
  DefineSeparately
>;

/// Wrapper alias for solid and wet nodes.
/*
Prototype for the usage of this class is given for NoDynamics.

Use FixedDensity & ZeroMomentum for NoDynamics, BounceBack(rho), BounceBackAnti(rho)
Use BulkDensity & ZeroMomentum for BounceBack(), BounceBackAnti(), ZeroDistributionDynamics
Use FixedDensity & FixedVelocityMomentumGeneric for BounceBackVelocity(rho, u)
Use BulkDensity & FixedVelocityMomentumGeneric for BounceBackVelocity(u)
Use FixedDensity & OffBoundaryMomentum for OffDynamics (not yet ready)

For now, the BounceBack* classes may therefore be separated into the cases whether
rho is fixed or not. They can later be reunified/ given a logical structure, when it
is clear how the dynamics will be constructed in the future.
*/
template <typename DENSITY, typename MOMENTUM>
using GenericSolidMaterialTuple = Tuple<
  DENSITY,
  MOMENTUM,
  ZeroStress,  // TODO: should we rather use NoStress and thus throw an error?
  DefineSeparately
>;

using PoissonTuple = Tuple<
  BulkDensity,
  PoissonMomentum,
  ZeroStress,
  DefineToNEq
>;

using P1Tuple = Tuple<
  BulkDensity,
  P1Momentum,
  ZeroStress,
  DefineToNEq
>;

using None = Tuple<
  ZeroDensity,
  ZeroMomentum,
  ZeroStress,
  DefineSeparately
>;

using FreeEnergyBulkTuple = Tuple<
  BulkDensity,
  FreeEnergyMomentum,
  BulkStress,
  DefineToNEq
>;

}  // namespace momenta

}  // namespace olb

#endif
