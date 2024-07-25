/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012, 2015 Mathias J. Krause, Vojtech Cvrcek, Davide Dapelo
 *                2022 Nando Suntoyo, Adrian Kummerlaender
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
 * Strain rate similar to "J.Boyd, J. Buick and S.Green: A second-order accurate lattice Boltzmann non-Newtonian flow model"
 * Power Law similar to "Huidan Yu, Sharath S. Girimaji, Li-Shi Luo - DNS and LES of decaying isotropic turbulence with and without frame rotation using lattice Boltzmann method"
 *
 */
#ifndef POWER_LAW_BGK_DYNAMICS_H
#define POWER_LAW_BGK_DYNAMICS_H

#include "dynamics/dynamics.h"
#include "core/cell.h"

#include "collisionLES.h"
#include "porousBGKdynamics.h"

namespace olb {

// *INDENT-OFF*

namespace powerlaw {

struct OMEGA_MIN : public descriptors::FIELD_BASE<1> { };
struct OMEGA_MAX : public descriptors::FIELD_BASE<1> { };
struct M : public descriptors::FIELD_BASE<1> { };
struct N : public descriptors::FIELD_BASE<1> { };

/// Compute and update cell-wise OMEGA using Oswald-de-waele model
template <typename COLLISION>
struct OmegaFromCell {
  using parameters = typename COLLISION::parameters::template include<
    descriptors::OMEGA, OMEGA_MIN, OMEGA_MAX, M, N
  >;

  static std::string getName()
  {
    return "powerlaw::OmegaFromCell<" + COLLISION::getName() + ">";
  }

  template <typename CELL, typename PARAMETERS, typename OMEGA, typename RHO, typename PI, typename V=typename CELL::value_t>
  static V computeOmega(CELL& cell, PARAMETERS& parameters, OMEGA& omega0, RHO& rho, PI& pi) any_platform
  {
    using DESCRIPTOR = typename CELL::descriptor_t;
    V pre2 = V{0.5} * descriptors::invCs2<V,DESCRIPTOR>() * omega0 / rho; // strain rate tensor prefactor
    pre2 *= pre2;
    V gamma{};
    if constexpr (DESCRIPTOR::template provides<descriptors::FORCE>()) {
      // Cannot be done in just one line, it gives error - I don't know why. Davide Dapelo
      const auto force = cell.template getField<descriptors::FORCE>();
      gamma = util::sqrt(V{2}*pre2*lbm<DESCRIPTOR>::computePiNeqNormSqr(cell, force));
    }
    else {
      gamma = util::sqrt(V{2}*pre2*lbm<DESCRIPTOR>::computePiNeqNormSqr(cell));
    }
    V m = parameters.template get<M>();
    V n = parameters.template get<N>();
    V nuNew = m * util::pow(gamma, n-V{1}); // nu for non-Newtonian fluid
    V newOmega = V{1} / (nuNew*descriptors::invCs2<V,DESCRIPTOR>() + V{0.5});
    V omegaMax = parameters.template get<OMEGA_MAX>();
    newOmega = util::min(newOmega, omegaMax);
    V omegaMin = parameters.template get<OMEGA_MIN>();
    newOmega = util::max(newOmega, omegaMin);
    return newOmega;
  }

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  struct type {
    using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;
    using CollisionO = typename COLLISION::template type<DESCRIPTOR, MOMENTA, EQUILIBRIUM>;

    template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform
    {
      V rho, u[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR>::n] { };
      MomentaF().computeAllMomenta(cell, rho, u, pi);
      const V oldOmega = cell.template getField<descriptors::OMEGA>();
      const V newOmega = computeOmega(cell, parameters, oldOmega, rho, pi);
      cell.template setField<descriptors::OMEGA>(newOmega);
      parameters.template set<descriptors::OMEGA>(newOmega);
      return CollisionO().apply(cell, parameters);
    }
  };
};

/// Combination rule to realize a pressure drop at a periodic boundary
template <int... NORMAL>
struct PeriodicPressureOffset {
  struct OFFSET : public descriptors::FIELD_BASE<1> { };

  static std::string getName()
  {
    return "PeriodicPressureOffset";
  }

  template <typename DESCRIPTOR, typename MOMENTA>
  using combined_momenta = typename MOMENTA::template type<DESCRIPTOR>;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  using combined_equilibrium = typename EQUILIBRIUM::template type<DESCRIPTOR,MOMENTA>;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM, typename COLLISION>
  struct combined_collision {
    using CollisionO = typename COLLISION::template type<DESCRIPTOR, MOMENTA, EQUILIBRIUM>;

    template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform
    {
      static constexpr auto populations = util::populationsContributingToDirection<DESCRIPTOR, NORMAL...>();

      auto statistic = CollisionO().apply(cell, parameters);
      const V densityOffset = parameters.template get<OFFSET>();
      for (unsigned iPop : populations) {
        cell[iPop] += (cell[iPop] + descriptors::t<V,DESCRIPTOR>(iPop)) * densityOffset;
      }
      return statistic;
    }
  };

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM, typename COLLISION>
  using combined_parameters = typename COLLISION::parameters::template include<OFFSET>;
};

}

/// BGK collision using Power Law collision frequency
template <typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using PowerLawBGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  powerlaw::OmegaFromCell<collision::BGK>
>;

/// BGK collision using Power Law collision frequency with Guo forcing
template <typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using PowerLawForcedBGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  powerlaw::OmegaFromCell<collision::BGK>,
  forcing::Guo<momenta::ForcedWithStress>
>;

/// Smagorinsky BGK collision using Power Law collision frequency
template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using SmagorinskyPowerLawBGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  powerlaw::OmegaFromCell<collision::SmagorinskyEffectiveOmega<collision::BGK>>
>;

/// Smagorinsky BGK collision using Power Law collision frequency and Guo forcing
template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using SmagorinskyPowerLawForcedBGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  powerlaw::OmegaFromCell<collision::SmagorinskyEffectiveOmega<collision::BGK>>,
  forcing::Guo<momenta::Forced>
>;

/// Smagorinsky BGK collision using Power Law collision frequency for porous particles
template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using SmagorinskyPowerLawPorousParticleBGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  powerlaw::OmegaFromCell<collision::SmagorinskyEffectiveOmega<collision::PorousParticle<collision::BGK>>>
>;

// *INDENT-ON*

}

#endif
