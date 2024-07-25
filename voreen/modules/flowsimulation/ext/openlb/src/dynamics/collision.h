/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Adrian Kummerlaender
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

#ifndef DYNAMICS_COLLISION_H
#define DYNAMICS_COLLISION_H

#include "lbm.h"
#include "descriptorField.h"
#include "rtlbmDescriptors.h"

#include "core/latticeStatistics.h"

namespace olb {

namespace collision {

struct None {
  using parameters = typename meta::list<>;

  static std::string getName() {
    return "None";
  }

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  struct type {
    template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
      return {-1, -1};
    };
  };
};

struct BGK {
  using parameters = typename meta::list<descriptors::OMEGA>;

  static std::string getName() {
    return "BGK";
  }

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  struct type {
    using EquilibriumF = typename EQUILIBRIUM::template type<DESCRIPTOR,MOMENTA>;

    template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
      V fEq[DESCRIPTOR::q] { };
      const auto statistic = EquilibriumF().compute(cell, parameters, fEq);
      const V omega = parameters.template get<descriptors::OMEGA>();
      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        cell[iPop] *= V{1} - omega;
        cell[iPop] += omega * fEq[iPop];
      }
      return statistic;
    };
  };
};

struct RLB {
  using parameters = typename meta::list<descriptors::OMEGA>;

  static std::string getName() {
    return "RLB";
  }

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  struct type {
    using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;
    using EquilibriumF = typename EQUILIBRIUM::template type<DESCRIPTOR,MOMENTA>;

    template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
      V pi[util::TensorVal<DESCRIPTOR>::n] { };
      MomentaF().computeStress(cell, pi);
      V fEq[DESCRIPTOR::q] { };
      const auto statistic = EquilibriumF().compute(cell, parameters, fEq);
      const V omega = parameters.template get<descriptors::OMEGA>();
      cell[0] = fEq[0] + (V{1} - omega) * equilibrium<DESCRIPTOR>::template fromPiToFneq<V>(0, pi);
      for (int iPop=1; iPop <= DESCRIPTOR::q/2; ++iPop) {
        cell[iPop] = fEq[iPop];
        cell[iPop+DESCRIPTOR::q/2] = fEq[iPop+DESCRIPTOR::q/2];
        V fNeq = (V{1} - omega) * equilibrium<DESCRIPTOR>::template fromPiToFneq<V>(iPop, pi);
        cell[iPop] += fNeq;
        cell[iPop+DESCRIPTOR::q/2] += fNeq;
      }
      return statistic;
    };
  };
};

struct AdvectionDiffusionRLB {
  using parameters = typename meta::list<descriptors::OMEGA>;

  static std::string getName() {
    return "RLB";
  }

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  struct type {
    using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;
    using EquilibriumF = typename EQUILIBRIUM::template type<DESCRIPTOR,MOMENTA>;

    template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) {
      V fEq[DESCRIPTOR::q] { };
      const auto statistic = EquilibriumF().compute(cell, parameters, fEq);
      V j1[DESCRIPTOR::d] { };
      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop ) {
        for (int iD=0; iD < DESCRIPTOR::d; ++iD ) {
          j1[iD] += descriptors::c<DESCRIPTOR>(iPop,iD) * (cell[iPop] - fEq[iPop]);
        }
      }
      const V omega = parameters.template get<descriptors::OMEGA>();
      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop ) {
        V fNeq{};
        for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
          fNeq += descriptors::c<DESCRIPTOR>(iPop,iD) * j1[iD];
        }
        fNeq *= descriptors::t<V,DESCRIPTOR>(iPop) * descriptors::invCs2<V,DESCRIPTOR>();
        cell[iPop] = fEq[iPop] + (V{1} - omega) * fNeq;
      }
      return statistic;
    };
  };
};

struct ConstRhoBGK {
  using parameters = typename meta::list<descriptors::OMEGA, statistics::AVERAGE_RHO>;

  static std::string getName() {
    return "ConstRhoBGK";
  }

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  struct type {
    using MomentaF     = typename MOMENTA::template type<DESCRIPTOR>;
    using EquilibriumF = typename EQUILIBRIUM::template type<DESCRIPTOR,MOMENTA>;

    template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
      V fEq[DESCRIPTOR::q];
      const auto statistic = EquilibriumF().compute(cell, parameters, fEq);

      V deltaRho = V{1} - parameters.template get<statistics::AVERAGE_RHO>();
      V ratioRho = V{1} + deltaRho / statistic.rho;

      const V omega = parameters.template get<descriptors::OMEGA>();

      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        cell[iPop] = ratioRho * (fEq[iPop] + descriptors::t<V,DESCRIPTOR>(iPop))
                   - descriptors::t<V,DESCRIPTOR>(iPop)
                   + (V{1} - omega) * (cell[iPop] - fEq[iPop]);
      }
      return {statistic.rho+deltaRho, statistic.uSqr};
    };
  };
};

struct PerPopulationBGK {
  struct OMEGA : public descriptors::FIELD_BASE<0,0,1> { };

  using parameters = typename meta::list<OMEGA>;

  static std::string getName() {
    return "PerPopulationBGK";
  }

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  struct type {
    using EquilibriumF = typename EQUILIBRIUM::template type<DESCRIPTOR,MOMENTA>;

    template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
      V fEq[DESCRIPTOR::q] { };
      const auto statistic = EquilibriumF().compute(cell, parameters, fEq);
      const auto omega = parameters.template get<OMEGA>();
      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        cell[iPop] *= V{1} - omega[iPop];
        cell[iPop] += omega[iPop] * fEq[iPop];
      }
      return statistic;
    };
  };
};

struct TRT {
  struct MAGIC : public descriptors::FIELD_BASE<1> { };

  using parameters = typename meta::list<descriptors::OMEGA,MAGIC>;

  static std::string getName() {
    return "TRT";
  }

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  struct type {
    using MomentaF     = typename MOMENTA::template type<DESCRIPTOR>;
    using EquilibriumF = typename EQUILIBRIUM::template type<DESCRIPTOR,MOMENTA>;

    template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
      V fPlus[DESCRIPTOR::q], fMinus[DESCRIPTOR::q];
      V fEq[DESCRIPTOR::q], fEqPlus[DESCRIPTOR::q], fEqMinus[DESCRIPTOR::q];
      const auto statistic = EquilibriumF().compute(cell, parameters, fEq);
      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        fPlus[iPop]  = 0.5 * (cell[iPop] + cell[descriptors::opposite<DESCRIPTOR>(iPop)]);
        fMinus[iPop] = 0.5 * (cell[iPop] - cell[descriptors::opposite<DESCRIPTOR>(iPop)]);
      }
      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        fEqPlus[iPop]  = 0.5 * (fEq[iPop] + fEq[descriptors::opposite<DESCRIPTOR>(iPop)]);
        fEqMinus[iPop] = 0.5 * (fEq[iPop] - fEq[descriptors::opposite<DESCRIPTOR>(iPop)]);
      }

      const V omega = parameters.template get<descriptors::OMEGA>();
      const V magic = parameters.template get<MAGIC>();
      const V omega2 = V{1} / (magic / (V{1} / omega - 0.5) + V{0.5});

      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        cell[iPop] -= omega * (fPlus[iPop] - fEqPlus[iPop]) + omega2 * (fMinus[iPop] - fEqMinus[iPop]);
      }
      return statistic;
    };
  };
};

struct Revert {
  using parameters = typename meta::list<>;

  static std::string getName() {
    return "Revert";
  }

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  struct type {
    template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
      for (int iPop=1; iPop <= DESCRIPTOR::q/2; ++iPop) {
        V cell_iPop = cell[iPop];
        cell[iPop] = cell[descriptors::opposite<DESCRIPTOR>(iPop)];
        cell[descriptors::opposite<DESCRIPTOR>(iPop)] = cell_iPop;
      }
      return {-1, -1};
    };
  };
};

/// Nguyen-Ladd Velocity Correction using momenta-defined velocity
/**
 * Use in conjunction with e.g. momenta::ExternalVelocityTuple
 *
 * N.-Q. Nguyen, A. J. C. Ladd, Lubrication corrections for lattice-
 * Boltzmann simulations of particle suspensions, Physical Review
 * E 66 (4) (Oct. 2002). doi:10.1103/PhysRevE.66.046708.
**/
template <typename COLLISION>
struct NguyenLaddCorrection {
  using parameters = typename COLLISION::parameters;

  static std::string getName() {
    return "NguyenLaddCorrection<" + COLLISION::getName() + ">";
  }

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  struct type {
    using MomentaF   = typename MOMENTA::template type<DESCRIPTOR>;
    using CollisionO = typename COLLISION::template type<DESCRIPTOR,MOMENTA,EQUILIBRIUM>;

    template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) any_platform {
      CollisionO().apply(cell, parameters);
      V rho, u[DESCRIPTOR::d];
      MomentaF().computeRhoU(cell, rho, u);
      for (int iPop=1; iPop < DESCRIPTOR::q; ++iPop) {
        for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
          cell[iPop] += 2 * rho * u[iD]*descriptors::c<DESCRIPTOR>(iPop,iD) * descriptors::t<V,DESCRIPTOR>(iPop) * descriptors::invCs2<V,DESCRIPTOR>();
        }
      }
      return {rho, util::normSqr<V,DESCRIPTOR::d>(u)};
    };
  };
};

struct PartialBounceBack {
  struct RF : public descriptors::FIELD_BASE<1> { };
  using parameters = typename meta::list<RF>;

  static std::string getName() {
    return "PartialBounceBack";
  }

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  struct type {
    template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) {
      for (int iPop=1; iPop <= DESCRIPTOR::q/2; ++iPop) {
        std::swap(cell[iPop], cell[iPop+DESCRIPTOR::q/2]);
      }
      const V rf = parameters.template get<RF>();
      for (int iPop=1; iPop < DESCRIPTOR::q; ++iPop) {
        cell[iPop] = (rf - V{1}) * (cell[iPop] + descriptors::t<V,DESCRIPTOR>(iPop))
                   - descriptors::t<V,DESCRIPTOR>(iPop);
      }
      return {-1, -1};
    };
  };
};

struct Poisson {
  struct SINK : public descriptors::FIELD_BASE<1> { };
  using parameters = typename meta::list<descriptors::OMEGA,SINK>;

  static std::string getName() {
    return "Poisson";
  }

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  struct type {
    using EquilibriumF = typename EQUILIBRIUM::template type<DESCRIPTOR,MOMENTA>;

    template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) {
      V fEq[DESCRIPTOR::q] { };
      const auto statistic = EquilibriumF().compute(cell, parameters, fEq);
      const V omega = parameters.template get<descriptors::OMEGA>();
      const V sink = parameters.template get<SINK>();
      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        cell[iPop] = (cell[iPop] + descriptors::t<V,DESCRIPTOR>(iPop))
                   - omega * ((cell[iPop] + descriptors::t<V,DESCRIPTOR>(iPop)) - descriptors::t<V,DESCRIPTOR>(iPop) * statistic.rho)
                   - sink * (cell[iPop] + descriptors::t<V,DESCRIPTOR>(iPop))
                   - descriptors::t<V,DESCRIPTOR>(iPop);
      }
      return statistic;
    };
  };
};

struct P1 {
  struct SCATTERING : public descriptors::FIELD_BASE<1> { };
  struct ABSORPTION : public descriptors::FIELD_BASE<1> { };
  using parameters = typename meta::list<SCATTERING,ABSORPTION>;

  static std::string getName() {
    return "P1";
  }

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  struct type {
    using EquilibriumF = typename EQUILIBRIUM::template type<DESCRIPTOR,MOMENTA>;

    template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) {
      V fEq[DESCRIPTOR::q] { };
      const auto statistic = EquilibriumF().compute(cell, parameters, fEq);
      const V scattering = parameters.template get<SCATTERING>();
      const V absorption = parameters.template get<ABSORPTION>();
      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        cell[iPop] -= scattering * descriptors::norm_c<V,DESCRIPTOR>(iPop)*(cell[iPop] - fEq[iPop])
                    + absorption * descriptors::norm_c<V,DESCRIPTOR>(iPop)*(cell[iPop] + descriptors::t<V,DESCRIPTOR>(iPop));
      }
      return statistic;
    };
  };
};

struct FixedEquilibrium {
  using parameters = typename meta::list<>;

  static std::string getName() {
    return "FixedEquilibrium";
  }

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  struct type {
    using EquilibriumF = typename EQUILIBRIUM::template type<DESCRIPTOR,MOMENTA>;

    template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) {
      V fEq[DESCRIPTOR::q] { };
      const auto statistic = EquilibriumF().compute(cell, parameters, fEq);
      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        cell[iPop] = fEq[iPop];
      }
      return statistic;
    };
  };
};

/// Override dynamics parameter PARAMETER with cell field PARAMETER
template <typename PARAMETER, typename COLLISION>
struct ParameterFromCell {
  using parameters = typename COLLISION::parameters::template include<PARAMETER>;

  static std::string getName() {
    return "ParameterFromCell<" + COLLISION::getName() + ">";
  }

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  struct type {
    using CollisionO = typename COLLISION::template type<DESCRIPTOR, MOMENTA, EQUILIBRIUM>;

    template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) {
      parameters.template set<PARAMETER>(
        cell.template getField<PARAMETER>());
      return CollisionO().apply(cell, parameters);
    }
  };
};

/// Override dynamics parameter OMEGA with cell field OMEGA
template <typename COLLISION>
using OmegaFromCell = ParameterFromCell<descriptors::OMEGA, COLLISION>;

/// Override dynamics parameter OMEGA with inverse of cell field TAU_EFF
template <typename COLLISION>
struct OmegaFromCellTauEff {
  using parameters = typename COLLISION::parameters::template include<
    descriptors::OMEGA
  >;

  static std::string getName() {
    return "OmegaFromCellTauEff<" + COLLISION::getName() + ">";
  }

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  struct type {
    using CollisionO = typename COLLISION::template type<DESCRIPTOR, MOMENTA, EQUILIBRIUM>;

    template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) {
      parameters.template set<descriptors::OMEGA>(
        V{1} / cell.template getField<descriptors::TAU_EFF>());
      return CollisionO().apply(cell, parameters);
    }
  };
};

}

}

#endif

#include "collision.cse.h"
