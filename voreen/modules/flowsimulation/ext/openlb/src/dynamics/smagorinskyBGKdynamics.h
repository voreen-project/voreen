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

#include "collisionLES.h"

namespace olb {

/// Smagorinsky BGK collision step
template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using SmagorinskyBGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  collision::SmagorinskyEffectiveOmega<collision::BGK>
>;

/// Smagorinsky BGK collision step with Guo forcing
template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using SmagorinskyForcedBGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  collision::SmagorinskyEffectiveOmega<collision::BGK>,
  forcing::Guo<momenta::ForcedWithStress>
>;

/// LES BGK collision using non-local TAU_EFF per-cell field
template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using ExternalTauEffLESBGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  collision::OmegaFromCellTauEff<collision::BGK>
>;

/// LES BGK collision with Guo forcing using non-local TAU_EFF per-cell field
template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using ExternalTauEffLESForcedBGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  collision::OmegaFromCellTauEff<collision::BGK>,
  forcing::Guo<momenta::Forced>
>;

/// Shear Smarorinsky BGK collision step
/**
 * Shown good results for wall-bounded flows
 * Leveque et al.: Shear-Improved Smagorinsky Model for Large-Eddy Simulation
 * of Wall-Bounded Turbulent Flows
 * DOI: http://dx.doi.org/10.1017/S0022112006003429
 **/
template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using ShearSmagorinskyBGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  collision::ShearSmagorinskyEffectiveOmega<collision::BGK>
>;

/// Shear Smarorinsky BGK collision step with Guo forcing
template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using ShearSmagorinskyForcedBGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  collision::ShearSmagorinskyEffectiveOmega<collision::BGK>,
  forcing::Guo<momenta::Forced>
>;

/// LES BGK collision for advection diffusion using non-local TAU_EFF per-cell field
template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using ExternalTauEffLESBGKadvectionDiffusionDynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::FirstOrder,
  collision::OmegaFromCellTauEff<collision::BGK>,
  AdvectionDiffusionExternalVelocityCollision
>;

/// Consistent Strain Smagorinsky BGK collision step
/**
 * Consistent subgrid scale modelling for lattice Boltzmann methods
 * Orestis Malaspinas and Pierre Sagaut
 * Journal of Fluid Mechanics / Volume / June 2012, pp 514-542
 * DOI: http://dx.doi.org/10.1017/jfm.2012.155
 **/
template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using ConStrainSmagorinskyBGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  collision::ConStrainSmagorinskyEffectiveOmega<collision::BGK>
>;

/// Consistent Smagorinsky BGK collision step
/**
 * Consistent subgrid scale modelling for lattice Boltzmann methods
 * Orestis Malaspinas and Pierre Sagaut
 * Journal of Fluid Mechanics / Volume / June 2012, pp 514-542
 * DOI: http://dx.doi.org/10.1017/jfm.2012.155
 **/
template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using ConSmagorinskyBGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  collision::ConSmagorinskyEffectiveOmega<collision::BGK>
>;

/// Smagorinsky BGK collision step with per-cell Smagorinsky constant
template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using ExternalSmagorinskyBGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  collision::ParameterFromCell<collision::LES::Smagorinsky,
                               collision::SmagorinskyEffectiveOmega<collision::BGK>>
>;

/// WALE LES BGK collision step
template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using WALEBGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  collision::WaleEffectiveOmega<collision::BGK>
>;

/// Krause BGK collision step
template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using KrauseBGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  collision::KrauseEffectiveOmega<collision::PerPopulationBGK>
>;

}

#endif
