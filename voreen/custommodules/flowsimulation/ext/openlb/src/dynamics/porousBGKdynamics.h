/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Adrian Kummerlaender, Nando Suntoyo
 *
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

#ifndef POROUS_BGK_DYNAMICS_H
#define POROUS_BGK_DYNAMICS_H

#include "interface.h"
#include "collision.h"
#include "dynamics/collisionLES.h"

namespace olb {

/// Porous BGK collision step
template <typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using PorousBGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  momenta::Porous<MOMENTA>,
  equilibria::SecondOrder,
  collision::BGK
>;


namespace collision {

/* Implementation of the BGK collision for moving porous media (HLBM approach).
 * As this scheme requires additionla data stored in an external field,
 * it is meant to be used along with a PorousParticle descriptor.
 * \param omega Lattice relaxation frequency
 * \param momenta A standard object for the momenta computation
 */
template <typename COLLISION, bool isStatic=false>
struct PorousParticle {
  using parameters = typename meta::list<descriptors::OMEGA>;

  static constexpr bool is_vectorizable = false;

  static std::string getName() {
    return "PorousParticle<" + COLLISION::getName() + ">" ;
  }

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  struct type {
    using MomentaF     = typename MOMENTA::template type<DESCRIPTOR>;
    using EquilibriumF = typename EQUILIBRIUM::template type<DESCRIPTOR,MOMENTA>;
    using CollisionO   = typename COLLISION::template type<DESCRIPTOR,MOMENTA,EQUILIBRIUM>;

    constexpr static bool is_vectorizable = false;

    template <typename CELL, typename pVELOCITY, typename V=typename CELL::value_t>
    void calculate (CELL& cell, pVELOCITY& pVelocity) {
      if constexpr (isStatic) {
        for (int i=0; i<DESCRIPTOR::d; i++)  {
          pVelocity[i] -= (1.-(cell.template getField<descriptors::POROSITY>())) * pVelocity[i];
        }
      } else {
        for (int i=0; i<DESCRIPTOR::d; i++)  {
          pVelocity[i] += (1.-cell.template getField<descriptors::POROSITY>())
                          * (cell.template getFieldComponent<descriptors::VELOCITY_NUMERATOR>(i)
                             / cell.template getField<descriptors::VELOCITY_DENOMINATOR>() - pVelocity[i]);
        }
        // reset external field for next timestep
        cell.template setField<descriptors::POROSITY>(1.);
        cell.template setField<descriptors::VELOCITY_DENOMINATOR>(0.);
        cell.template setField<descriptors::VELOCITY_NUMERATOR>({0.,0.,0.});
      }
    }

    template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) {
      V rho, u[DESCRIPTOR::d];
      MomentaF().computeRhoU(cell, rho, u);

      auto statistic = CollisionO().apply(cell, parameters);

      V velDenominator = cell.template getFieldComponent<descriptors::VELOCITY_DENOMINATOR>(0);

      // use Kuperstokh forcing by default
      V uPlus[DESCRIPTOR::d]{ };
      V diff[DESCRIPTOR::q]{ };

      if (velDenominator > std::numeric_limits<V>::epsilon()) {
        for (int iDim=0; iDim<DESCRIPTOR::d; iDim++) {
          uPlus[iDim] = u[iDim];
        }
        calculate(cell, uPlus);
        for (int tmp_iPop=0; tmp_iPop < DESCRIPTOR::q; tmp_iPop++) {
          diff[tmp_iPop] +=   EquilibriumF().compute(tmp_iPop, rho, uPlus)
                            - EquilibriumF().compute(tmp_iPop, rho, u);
          cell[tmp_iPop] += diff[tmp_iPop];
        }
      }

      return statistic;
    }
  };
};

/// Implementation of the Partially Saturated Method (PSM),
/// see Krüger, Timm, et al. The Lattice Boltzmann Method. Springer, 2017. (p.447-451)
template <typename COLLISION>
struct PSM {
  using parameters = typename meta::list<descriptors::OMEGA>;

  static std::string getName() {
    return "PSM<" + COLLISION::getName() + ">" ;
  }

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  struct type {
    using MomentaF     = typename MOMENTA::template type<DESCRIPTOR>;
    using EquilibriumF = typename EQUILIBRIUM::template type<DESCRIPTOR,MOMENTA>;
    using CollisionO   = typename COLLISION::template type<DESCRIPTOR,MOMENTA,EQUILIBRIUM>;

    template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) {
      V rho, u[DESCRIPTOR::d];
      const V epsilon = 1. - cell.template getField<descriptors::POROSITY>();
      const V omega = parameters.template get<descriptors::OMEGA>();
      const V paramA = V{1.} / omega - V{0.5};
      MomentaF().computeRhoU(cell, rho, u);
      // velocity at the boundary
      auto u_s = cell.template getField<descriptors::VELOCITY_SOLID>();
      //if (epsilon < 1e-5) {
      //  return CollisionO().apply(cell, parameters);
      //}
      //else {
        // speed up paramB
        V paramB = (epsilon * paramA) / ((1. - epsilon) + paramA);
        // speed up paramC
        V paramC = (1. - paramB);
        V omega_s[DESCRIPTOR::q];
        V cell_tmp[DESCRIPTOR::q];
        for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
          cell_tmp[iPop] = cell[iPop];
          //// To be reimplemented sensibly in case one needs it
          //if constexpr (MODE == 2) {
            omega_s[iPop] = (EquilibriumF().compute(iPop, rho, u_s) - cell[iPop])
                            + (V{1} - omega) * (cell[iPop] - EquilibriumF().compute(iPop, rho, u));
          //} else if constexpr (MODE == 3) {
          //  omega_s[iPop] = (cell[descriptors::opposite<DESCRIPTOR>(iPop)] - EquilibriumF().compute(descriptors::opposite<DESCRIPTOR>(iPop), rho, u_s))
          //                  - (cell[iPop] - EquilibriumF().compute(iPop, rho, u_s));
          //}
        }
        CollisionO().apply(cell, parameters);
        for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
          cell[iPop] = cell_tmp[iPop] + paramC * (cell[iPop] - cell_tmp[iPop]);
          cell[iPop] += paramB * omega_s[iPop];
        }
        for (int iVel=0; iVel<DESCRIPTOR::d; ++iVel) {
          u[iVel] = paramC * u[iVel] + paramB * u_s[iVel];
        }
      //}
      return {rho, util::normSqr<V,DESCRIPTOR::d>(u)};
    }
  };
};

}

namespace dynamics {

struct ExposePorousParticleMomenta {
  static std::string getName() {
    return "ExposePorousParticleMomenta";
  }

  template <typename DESCRIPTOR, typename MOMENTA>
  using combined_momenta = typename momenta::PorousParticle<MOMENTA>::template type<DESCRIPTOR>;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  using combined_equilibrium = typename EQUILIBRIUM::template type<DESCRIPTOR,MOMENTA>;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM, typename COLLISION>
  using combined_collision = typename COLLISION::template type<DESCRIPTOR,MOMENTA,EQUILIBRIUM>;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM, typename COLLISION>
  using combined_parameters = typename COLLISION::parameters;
};

}


/// Porous particle BGK collision for moving particles
template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using PorousParticleBGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  collision::PorousParticle<collision::BGK,false>,
  dynamics::ExposePorousParticleMomenta
>;

/// Porous particle BGK collision for static particles
template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using StaticPorousParticleBGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  collision::PorousParticle<collision::BGK,true>,
  dynamics::ExposePorousParticleMomenta
>;

/// BGK collision step for a porosity model
template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using SmagorinskyPorousParticleBGKdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  collision::SmagorinskyEffectiveOmega<collision::PorousParticle<collision::BGK>>
>;

/// Implementation of the Partially Saturated Method (PSM),
/// see Krüger, Timm, et al. The Lattice Boltzmann Method. Springer, 2017. (p.447-451)
template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using PSMBGKdynamics = dynamics::Tuple<
  T,DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  collision::PSM<collision::BGK>
>;

/// Implementation of the Partially Saturated Method (PSM),
/// see Krüger, Timm, et al. The Lattice Boltzmann Method. Springer, 2017. (p.447-451)
template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
class ForcedPSMBGKdynamics final : public dynamics::CustomCollision<T,DESCRIPTOR,MOMENTA> {
private:
  olb::ParametersD<T,DESCRIPTOR,descriptors::OMEGA>* _parameters;

public:
  using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;
  using EquilibriumF = typename equilibria::SecondOrder::template type<DESCRIPTOR,MOMENTA>;
  using ParametersD = olb::ParametersD<T,DESCRIPTOR,descriptors::OMEGA>;

  constexpr static bool is_vectorizable = false;
  constexpr static bool has_parametrized_momenta = true;

  void setMomentaParameters(ParametersD* parameters) {
    _parameters = parameters;
  }

  std::type_index id() override {
    return typeid(ForcedPSMBGKdynamics);
  }

  AbstractParameters<T,DESCRIPTOR>& getParameters(BlockLattice<T,DESCRIPTOR>& block) override {
    return block.template getData<DynamicsParameters<ForcedPSMBGKdynamics>>();
  }

  void computeU(ConstCell<T,DESCRIPTOR>& cell, T u[DESCRIPTOR::d]) const override
  {
    MomentaF().computeU(cell, u);
    T epsilon = 1. - cell.template getField<descriptors::POROSITY>();
    T omega = _parameters->template get<descriptors::OMEGA>();
    T paramA = 1 / omega - 0.5;
    T paramB = (epsilon * paramA) / ((1. - epsilon) + paramA);
    T paramC = (1. - paramB);
    for (int iVel=0; iVel < DESCRIPTOR::d; ++iVel) {
      u[iVel] = paramC * (u[iVel] + cell.template getFieldComponent<descriptors::FORCE>(iVel) * 0.5)
              + paramB * cell.template getFieldComponent<descriptors::VELOCITY_SOLID>(iVel);
    }
  }

  void computeRhoU(ConstCell<T,DESCRIPTOR>& cell, T& rho, T u[DESCRIPTOR::d]) const override
  {
    MomentaF().computeRhoU(cell, rho, u);
    T epsilon = 1. - cell.template getField<descriptors::POROSITY>();
    T omega = _parameters->template get<descriptors::OMEGA>();
    T paramA = 1 / omega - 0.5;
    T paramB = (epsilon * paramA) / ((1. - epsilon) + paramA);
    T paramC = (1. - paramB);
    for (int iVel=0; iVel < DESCRIPTOR::d; ++iVel) {
      u[iVel] = paramC * (u[iVel] + cell.template getFieldComponent<descriptors::FORCE>(iVel) * 0.5)
              + paramB * cell.template getFieldComponent<descriptors::VELOCITY_SOLID>(iVel);
    }
  }

  template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
  CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters)
  {
    V omega = parameters.template get<descriptors::OMEGA>();
    V rho, u[DESCRIPTOR::d], uSqr;
    MomentaF().computeRhoU(cell, rho, u);
    auto force = cell.template getField<descriptors::FORCE>();
    for (int iVel=0; iVel < DESCRIPTOR::d; ++iVel) {
      u[iVel] += force[iVel] * V{0.5};
    }
    auto u_s = cell.template getField<descriptors::VELOCITY_SOLID>();

    V epsilon = V{1} - cell.template getField<descriptors::POROSITY>();
    if (epsilon < 1e-5) {
      uSqr = lbm<DESCRIPTOR>::bgkCollision(cell, rho, u, omega);
      lbm<DESCRIPTOR>::addExternalForce(cell, rho, u, omega, force);
    }
    else {
      V paramA = V{1} / omega - V{0.5};
      V paramB = (epsilon * paramA) / ((V{1} - epsilon) + paramA);
      V paramC = (V{1} - paramB);

      V omega_s[DESCRIPTOR::q];
      V cell_tmp[DESCRIPTOR::q];

      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        cell_tmp[iPop] = cell[iPop];
        omega_s[iPop] = (equilibrium<DESCRIPTOR>::secondOrder(iPop, rho, u_s) - cell[iPop])
                      + (V{1} - omega) * (cell[iPop] - equilibrium<DESCRIPTOR>::secondOrder(iPop, rho, u));
      }

      uSqr = lbm<DESCRIPTOR>::bgkCollision(cell, rho, u, omega);
      lbm<DESCRIPTOR>::addExternalForce(cell, rho, u, omega, force);

      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        cell[iPop] = cell_tmp[iPop] + paramC * (cell[iPop] - cell_tmp[iPop]);
        cell[iPop] += paramB * omega_s[iPop];
      }
      for (int iVel=0; iVel < DESCRIPTOR::d; ++iVel) {
        u[iVel] = paramC * u[iVel] + paramB * u_s[iVel];
      }
      uSqr = util::normSqr<V,DESCRIPTOR::d>(u);
    }
    return {rho, uSqr};
  }

  T computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d]) const override {
    return EquilibriumF().compute(iPop, rho, u);
  };

  std::string getName() const override {
    return "ForcedPSMBGKdynamics<" + MomentaF().getName() + ">";
  };

};

}

#endif

#include "legacy/porousBGKdynamics.h"
