/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2018 Robin Trunk, Sam Avis
 *                2021 Adrian Kummerlaender
 *  OpenLB e-mail contact: info@openlb.net
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

#ifndef FREE_ENERGY_DYNAMICS_H
#define FREE_ENERGY_DYNAMICS_H

#include "lbm.h"
#include "interface.h"
#include "momenta/aliases.h"

/** \file
 * In this file the dynamic calls for the free energy model is implemented. It
 * is used for the second (and third) lattices, as for the first one a BGK collision with
 * Guo forcing is applied (see ForcedBGKdynamcs).
 */
namespace olb {

namespace equilibria {

struct FreeEnergy {
  using parameters = meta::list<>;

  static std::string getName() {
    return "FreeEnergy";
  }

  template <typename DESCRIPTOR, typename MOMENTA>
  struct type {
    using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;

    template <typename RHO, typename U, typename V=RHO>
    auto compute(int iPop, const RHO& rho, const U& u) {
      V eq = equilibrium<DESCRIPTOR>::secondOrder(iPop, rho, u);
      if (iPop == 0) {
        eq += (V{1} - descriptors::t<V,DESCRIPTOR>(0)) * rho;
      } else {
        eq -= descriptors::t<V,DESCRIPTOR>(iPop) * rho;
      }
      return eq;
    }

    template <typename CELL, typename PARAMETERS, typename FEQ, typename V=typename CELL::value_t>
    CellStatistic<V> compute(CELL& cell, PARAMETERS& parameters, FEQ& fEq) {
      auto u = cell.template getField<descriptors::FORCE>();
      V rho = MomentaF().computeRho(cell);
      const V uSqr = util::normSqr<V,DESCRIPTOR::d>(u);
      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        fEq[iPop] = equilibrium<DESCRIPTOR>::secondOrder(iPop, rho, u, uSqr);
      }
      return {rho, uSqr};
    };
  };
};

}

namespace collision {

struct FreeEnergy {
  struct GAMMA : public descriptors::FIELD_BASE<1> { };

  using parameters = typename meta::list<descriptors::OMEGA,GAMMA>;

  static std::string getName() {
    return "FreeEnergy";
  }

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  struct type {
    using EquilibriumF = typename EQUILIBRIUM::template type<DESCRIPTOR,MOMENTA>;

    template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) {
      V fEq[DESCRIPTOR::q] { };
      const auto statistic = EquilibriumF().compute(cell, parameters, fEq);
      const V omega = parameters.template get<descriptors::OMEGA>();
      const V gamma = parameters.template get<GAMMA>();
      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        cell[iPop] *= V{1} - omega;
        cell[iPop] += omega * fEq[iPop];
      }
      const V tmp = gamma * descriptors::invCs2<V,DESCRIPTOR>()
                  * cell.template getField<descriptors::CHEM_POTENTIAL>();
      for (int iPop=1; iPop < DESCRIPTOR::q; ++iPop) {
        cell[iPop] -= omega * descriptors::t<V,DESCRIPTOR>(iPop) * (statistic.rho - tmp);
      }
      cell[0] += omega * (V{1} - descriptors::t<V,DESCRIPTOR>(0)) * (statistic.rho - tmp);
      return statistic;
    };
  };
};

template <int direction, int orientation>
struct FreeEnergyInletOutlet {
  using parameters = typename meta::list<descriptors::OMEGA>;

  static std::string getName() {
    return "FreeEnergyInletOutlet<"
      + std::to_string(direction) + "," + std::to_string(orientation) +
    ">";
  }

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  struct type {
    using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;

    template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) {
      // Do a standard collision neglecting the chemical potential term.
      V rho, u[DESCRIPTOR::d];
      MomentaF().computeRhoU(cell, rho, u);
      const V omega = parameters.template get<descriptors::OMEGA>();
      V uSqr = lbm<DESCRIPTOR>::bgkCollision(cell, rho, u, omega);
      for (int iPop=1; iPop < DESCRIPTOR::q; ++iPop) {
        cell[iPop] -= omega * descriptors::t<V,DESCRIPTOR>(iPop) * rho;
      }
      cell[0] += omega * (V{1} - descriptors::t<V,DESCRIPTOR>(0)) * rho;
      // Distribute the missing density to the unknown distribution functions.
      constexpr auto missingIndices = util::subIndexOutgoing<DESCRIPTOR,direction,orientation>();
      V missingRho = rho - V{1};
      V missingWeightSum = 0;
      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        if (std::find(missingIndices.begin(), missingIndices.end(), iPop) != missingIndices.end()) {
          missingWeightSum += descriptors::t<V,DESCRIPTOR>(iPop);
        } else {
          missingRho -= cell[iPop];
        }
      }
      for (unsigned iPop=0; iPop < missingIndices.size(); ++iPop) {
        cell[missingIndices[iPop]] = missingRho * descriptors::t<V,DESCRIPTOR>(missingIndices[iPop]) / missingWeightSum;
      }
      return {rho, uSqr};
    };
  };
};

}

template <typename T, typename DESCRIPTOR, typename MOMENTA=momenta::FreeEnergyBulkTuple>
using FreeEnergyBGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::FreeEnergy,
  collision::FreeEnergy
>;

template <typename T, typename DESCRIPTOR>
using FreeEnergyWallDynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  momenta::Tuple<
    momenta::BulkDensity,
    momenta::ZeroMomentum,
    momenta::ZeroStress,
    momenta::DefineSeparately
  >,
  equilibria::FreeEnergy,
  collision::Revert
>;

template <typename T, typename DESCRIPTOR, int direction, int orientation>
using FreeEnergyInletOutletDynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  momenta::Tuple<
    momenta::FreeEnergyInletOutletDensity,
    momenta::FreeEnergyInletOutletMomentum<direction,orientation>,
    momenta::RegularizedBoundaryStress<direction,orientation>,
    momenta::DefineSeparately
  >,
  equilibria::FreeEnergy,
  collision::FreeEnergyInletOutlet<direction,orientation>
>;

}

#endif

#include "freeEnergyDynamics.cse.h"
