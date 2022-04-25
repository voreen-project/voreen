/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Jonas Latt, Mathias J. Krause
 *                2021 Adrian Kummerlaender
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

#ifndef LB_DYNAMICS_H
#define LB_DYNAMICS_H

#include "interface.h"

#include "core/util.h"
#include "core/postProcessing.h"
#include "core/latticeStatistics.h"
#include "latticeDescriptors.h"

#include "momenta/interface.h"
#include "momenta/aliases.h"

#include "collision.h"
#include "equilibrium.h"
#include "forcing.h"

namespace olb {

/// Dynamics for "dead cells" doing nothing
template <typename T, typename DESCRIPTOR>
using NoDynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  // Return 1 for the 0th moment, 0 for all others
  momenta::Tuple<
    momenta::OneDensity,
    momenta::ZeroMomentum,
    momenta::ZeroStress,
    momenta::DefineSeparately
  >,
  equilibria::None,
  collision::None
>;

/// Common BGK collision step
template <typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using BGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  std::conditional_t<   (DESCRIPTOR::d == 3 && DESCRIPTOR::q == 7)
                     || (DESCRIPTOR::d == 2 && DESCRIPTOR::q == 5),
                        equilibria::FirstOrder,
                        equilibria::SecondOrder>,
  collision::BGK
>;

/// Pressure-corrected BGK collision step
template <typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using ConstRhoBGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  collision::ConstRhoBGK
>;

/// BGK collision step with external force (Guo)
/**
 * Guo Z, Zheng C, Shi B (2002) Discrete lattice effects on the
 * forcing term in the lattice Boltzmann method. Phys Rev E 65(4)
 * DOI: 10.1103/PhysRevE.65.046308
 **/
template <typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using ForcedBGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  collision::BGK,
  forcing::Guo<momenta::Forced>
>;

/// BGK collision step with external force (Kupershtokh)
/**
 * Kupershtokh A, Medvedev D, Karpov D (2009) On equations
 * of state in a lattice Boltzmann method. Comput Math Appl 58(5)
 * DOI: 10.1016/j.camwa.2009.02.024
 **/
template <typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using ForcedKupershtokhBGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  collision::BGK,
  forcing::Kupershtokh
>;

/// BGK collision step with external force (Shan and Chen)
/**
 * Shan X, Chen H (1993) Lattice Boltzmann model for simulating
 * flows with multiple phases and components. Phys Rev E. 47 (3)
 * DOI: 10.1103/PhysRevE.47.1815
 **/
template <typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using ForcedShanChenBGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  collision::BGK,
  forcing::ShanChen
>;

/// Incompressible BGK collision step
template <typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using IncBGKdynamics = dynamics::Tuple<T,DESCRIPTOR,MOMENTA,equilibria::Incompressible,collision::BGK>;

/// Incompressible BGK collision step with external force
/**
 * Using Guo forcing on incompressible BGK
 **/
template <typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using ForcedIncBGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::Incompressible,
  collision::BGK,
  forcing::Guo<momenta::Forced>
>;

/// Incompressible BGK collision step with relaxation frequency 1 / TAU_EFF and external force
template <typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using ExternalTauForcedIncBGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::Incompressible,
  collision::OmegaFromCellTauEff<collision::BGK>,
  forcing::Guo<momenta::Forced>
>;

/// Regularized BGK collision step
/**
 * This model is substantially more stable than plain BGK, and has roughly
 * the same efficiency. However, it cuts out the modes at higher Knudsen
 * numbers and can not be used in the regime of rarefied gases.
 **/
template <typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using RLBdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  collision::RLB
>;

/// Regularized BGK collision followed by any other Dynamics
/**
 * Commonly used to model boundaries consisting of a specific boundary
 * momenta and standard dynamics::Tuple-declared dynamics.
 **/
template <typename T, typename DESCRIPTOR, typename DYNAMICS, typename MOMENTA>
class CombinedRLBdynamics final : public dynamics::CustomCollision<T,DESCRIPTOR,MOMENTA> {
private:
  // Use same MOMENTA in combined and nested (boundary) dynamics
  using CORRECTED_DYNAMICS = typename DYNAMICS::template exchange_momenta<MOMENTA>;

public:
  using MomentaF     = typename MOMENTA::template type<DESCRIPTOR>;
  using ParametersD  = typename CORRECTED_DYNAMICS::ParametersD;

  template <typename M>
  using exchange_momenta = CombinedRLBdynamics<T,DESCRIPTOR,DYNAMICS,M>;

  std::type_index id() override {
    return typeid(CombinedRLBdynamics);
  }

  AbstractParameters<T,DESCRIPTOR>& getParameters(BlockLattice<T,DESCRIPTOR>& block) override {
    return block.template getData<DynamicsParameters<CombinedRLBdynamics>>();
  }

  template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
  CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
    V rho, u[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR>::n];
    MomentaF().computeAllMomenta(cell,rho,u,pi);

    for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
      cell[iPop] = typename CORRECTED_DYNAMICS::EquilibriumF().compute(iPop, rho, u)
                 + equilibrium<DESCRIPTOR>::template fromPiToFneq<V>(iPop, pi);
    }

    return typename CORRECTED_DYNAMICS::CollisionO().apply(cell, parameters);
  };

  T computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d]) const override any_platform {
    return typename CORRECTED_DYNAMICS::EquilibriumF().compute(iPop, rho, u);
  };

  std::string getName() const override {
    return "CombinedRLBdynamics<" + CORRECTED_DYNAMICS().getName() + ">";
  };

};

/// TRT collision step
template <typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using TRTdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  collision::TRT
>;

/// TRT collision step with external force
/**
 * Using Guo forcing on incompressible BGK
 **/
template <typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using ForcedTRTdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  collision::TRT,
  forcing::Guo<momenta::Forced>
>;

/// Bounce Back boundary dynamics
/**
 * This is a very popular way to implement no-slip boundary conditions,
 * due to the simplicity and due to the independence of the boundary's
 * orientation.
 */
template <typename T, typename DESCRIPTOR>
using BounceBack = dynamics::Tuple<
  T, DESCRIPTOR,
  momenta::Tuple<
    momenta::FixedDensity,
    momenta::ZeroMomentum,
    momenta::ZeroStress,
    momenta::DefineSeparately
  >,
  equilibria::None,
  collision::Revert
>;

/// Bounce Back boundary dynamics with bulk density
/**
 * Different from BounceBack these dynamics compute the
 * 0th moment from the cell's populations instead of
 * applying a fixed value.
 **/
template <typename T, typename DESCRIPTOR>
using BounceBackBulkDensity = dynamics::Tuple<
  T, DESCRIPTOR,
  momenta::Tuple<
    momenta::BulkDensity,
    momenta::ZeroMomentum,
    momenta::ZeroStress,
    momenta::DefineSeparately
  >,
  equilibria::None,
  collision::Revert
>;

/// Bounce Back boundary dynamics with Nguyen-Ladd velocity correction
/**
 * This is a very popular way to implement no-slip boundary conditions,
 * due to the simplicity and due to the independence of the boundary's
 * orientation.
 */
template <typename T, typename DESCRIPTOR>
using BounceBackVelocity = dynamics::Tuple<
  T, DESCRIPTOR,
  momenta::ExternalVelocityTuple,
  equilibria::None,
  collision::NguyenLaddCorrection<collision::Revert>
>;

/**
 * Corresponds to macro Robin boundary, micro Fresnel surface
 * Motivated by Hiorth et al. 2008; doi 10.1002/fld.1822
 */
template <typename T, typename DESCRIPTOR>
using PartialBounceBack = dynamics::Tuple<
  T, DESCRIPTOR,
  momenta::Tuple<momenta::FixedDensity,momenta::ZeroMomentum,momenta::ZeroStress,momenta::DefineSeparately>,
  equilibria::FirstOrder,
  collision::PartialBounceBack
>;

/// Poisson dynamics
template <typename T, typename DESCRIPTOR, typename MOMENTA=momenta::PoissonTuple>
using PoissonDynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::ZerothOrder,
  collision::Poisson
>;

/// P1 dynamics
template <typename T, typename DESCRIPTOR, typename MOMENTA=momenta::P1Tuple>
using P1Dynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::P1,
  collision::P1
>;

/// Models a density sink by enforcing a zero distribution on the cell
template <typename T, typename DESCRIPTOR>
struct ZeroDistributionDynamics final : public dynamics::CustomCollision<
  T, DESCRIPTOR,
  momenta::Tuple<
    momenta::BulkDensity,
    momenta::ZeroMomentum,
    momenta::ZeroStress,
    momenta::DefineSeparately
  >
> {
  using ParametersD = olb::ParametersD<T,DESCRIPTOR>;

  std::type_index id() override {
    return typeid(ZeroDistributionDynamics);
  };

  AbstractParameters<T,DESCRIPTOR>& getParameters(BlockLattice<T,DESCRIPTOR>& block) override {
    return block.template getData<DynamicsParameters<ZeroDistributionDynamics>>();
  }

  template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
  CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) {
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      cell[iPop] = -descriptors::t<T,DESCRIPTOR>(iPop);
    }
    return {-1, -1};
  };

  T computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d]) const override {
    return 0;
  };

  std::string getName() const override {
    return "ZeroDistributionDynamics";
  };
};

template <typename T, typename DESCRIPTOR, typename MOMENTA>
using ChopardDynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::Chopard,
  collision::BGK
>;

/// First order equilibrium boundary dynamics
/**
 * Applies the first order equilibrium distribution on every
 * time step using fixed density and velocity moments.
 */
template <typename T, typename DESCRIPTOR>
using EquilibriumBoundaryFirstOrder = dynamics::Tuple<
  T, DESCRIPTOR,
  momenta::EquilibriumBoundaryTuple,
  equilibria::FirstOrder,
  collision::FixedEquilibrium
>;

/// Second order equilibrium boundary dynamics
/**
 * Applies the second order equilibrium distribution on every
 * time step using fixed density and velocity moments.
 */
template <typename T, typename DESCRIPTOR>
using EquilibriumBoundarySecondOrder = dynamics::Tuple<
  T, DESCRIPTOR,
  momenta::EquilibriumBoundaryTuple,
  equilibria::SecondOrder,
  collision::FixedEquilibrium
>;

/// VANS BGK collision step with external force
/**
 * WIP: computeVANSRhoU doesn't exist (and should not)
 * Forcing to be realized using standard combinations rules (Likely Guo in this case)
 **/
template <typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
struct ForcedVANSBGKdynamics final : public dynamics::CustomCollision<T,DESCRIPTOR,MOMENTA> {
  using ParametersD = olb::ParametersD<T,DESCRIPTOR,descriptors::OMEGA>;
  using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;

  std::type_index id() override {
    return typeid(ForcedVANSBGKdynamics);
  }

  AbstractParameters<T,DESCRIPTOR>& getParameters(BlockLattice<T,DESCRIPTOR>& block) override {
    return block.template getData<DynamicsParameters<ForcedVANSBGKdynamics>>();
  }

  void computeU(ConstCell<T,DESCRIPTOR>& cell, T u[DESCRIPTOR::d]) const override {
    T rho;
    T porosity = cell.template getField<descriptors::POROSITY>();
    MomentaF().computeVANSRhoU(cell, rho, u, &porosity);
    auto force = cell.template getFieldPointer<descriptors::FORCE>();
    for (int iVel=0; iVel < DESCRIPTOR::d; ++iVel) {
      u[iVel] += force[iVel] * T{0.5};
    }
  }

  void computeRhoU(ConstCell<T,DESCRIPTOR>& cell, T& rho, T u[DESCRIPTOR::d]) const override {
    T porosity = cell.template getField<descriptors::POROSITY>();
    MomentaF().computeVANSRhoU(cell, rho, u, &porosity);
    auto force = cell.template getFieldPointer<descriptors::FORCE>();
    for (int iVel=0; iVel<DESCRIPTOR::d; ++iVel) {
      u[iVel] += force[iVel] * T{0.5};
    }
  }

  template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
  CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) {
    V rho, u[DESCRIPTOR::d];
    V porosity = cell.template getField<descriptors::POROSITY>();
    MomentaF().computeVANSRhoU(cell, rho, u, &porosity);
    auto force = cell.template getFieldPointer<descriptors::FORCE>();
    for (int iVel=0; iVel < DESCRIPTOR::d; ++iVel) {
      u[iVel] += force[iVel] * V{0.5};
      u[iVel] /= porosity;
    }
    rho *= porosity;
    const auto omega = parameters.template get<descriptors::OMEGA>();
    V uSqr = lbm<DESCRIPTOR>::bgkCollision(cell, rho, u, omega);
    lbm<DESCRIPTOR>::addExternalForce(cell, u, omega, rho);
    return {rho / porosity, uSqr * porosity * porosity};
  }

};

}

#endif

#include "legacy/dynamics.h"
