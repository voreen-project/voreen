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

#ifndef DYNAMICS_COLLISION_MRT_H
#define DYNAMICS_COLLISION_MRT_H

#include "mrt.h"
#include "descriptorField.h"
#include "mrtLatticeDescriptors.h"

#include "core/latticeStatistics.h"

namespace olb {

namespace collision {

struct MRT {
  using parameters = typename meta::list<descriptors::OMEGA>;

  static std::string getName() {
    return "MRT";
  }

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  struct type {
    using MomentaF     = typename MOMENTA::template type<DESCRIPTOR>;
    using EquilibriumF = typename EQUILIBRIUM::template type<DESCRIPTOR,MOMENTA>;

    template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) {
      const V omega = parameters.template get<descriptors::OMEGA>();

      V rt[DESCRIPTOR::q] { }; // relaxation times vector.
      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        rt[iPop] = descriptors::s<V,DESCRIPTOR>(iPop);
      }
      for (int iPop=0; iPop < descriptors::shearIndexes<DESCRIPTOR>(); ++iPop) {
        rt[descriptors::shearViscIndexes<DESCRIPTOR>(iPop)] = omega;
      }
      V invM_S[DESCRIPTOR::q][DESCRIPTOR::q] { }; // relaxation times matrix
      for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
        for (int jPop = 0; jPop < DESCRIPTOR::q; ++jPop) {
          invM_S[iPop][jPop] = V{};
          for (int kPop = 0; kPop < DESCRIPTOR::q; ++kPop) {
            if (kPop == jPop) {
              invM_S[iPop][jPop] += descriptors::invM<V,DESCRIPTOR>(iPop,kPop) * rt[kPop];
            }
          }
        }
      }

      V rho, u[DESCRIPTOR::d];
      MomentaF().computeRhoU(cell, rho, u);
      V uSqr = mrt<DESCRIPTOR>::mrtCollision(cell, rho, u, invM_S);

      return {rho, uSqr};
    };
  };
};

}

namespace forcing {

/// Dynamics combination rule implementing the forcing scheme by Ladd and Verberg
/**
 * A.Ladd, R. Verberg, DESCRIPTOR-Boltzmann simulations of particle-fluid suspensions, Journal of Statistical Physics 104(2001)
 **/
struct LaddVerberg {
  static std::string getName() {
    return "LaddVerbergForcing";
  }

  template <typename DESCRIPTOR, typename MOMENTA>
  using combined_momenta = typename MOMENTA::template type<DESCRIPTOR>;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  using combined_equilibrium = typename EQUILIBRIUM::template type<DESCRIPTOR,MOMENTA>;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM, typename COLLISION>
  struct combined_collision {
    using MomentaF   = typename MOMENTA::template type<DESCRIPTOR>;
    using CollisionO = typename COLLISION::template type<DESCRIPTOR,MOMENTA,EQUILIBRIUM>;

    template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
    CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) {
      V rho, u[DESCRIPTOR::d];
      MomentaF().computeRhoU(cell, rho, u);
      const auto statistic = CollisionO().apply(cell, parameters);
      const V omega = parameters.template get<descriptors::OMEGA>();
      // While this duplication can be resolved using CSE it should be extracted into a helper
      V rt[DESCRIPTOR::q] { }; // relaxation times vector.
      for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
        rt[iPop] = descriptors::s<V,DESCRIPTOR>(iPop);
      }
      for (int iPop=0; iPop < descriptors::shearIndexes<DESCRIPTOR>(); ++iPop) {
        rt[descriptors::shearViscIndexes<DESCRIPTOR>(iPop)] = omega;
      }
      V invM_S[DESCRIPTOR::q][DESCRIPTOR::q]; // relaxation times matrix
      for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
        for (int jPop = 0; jPop < DESCRIPTOR::q; ++jPop) {
          invM_S[iPop][jPop] = V{};
          for (int kPop = 0; kPop < DESCRIPTOR::q; ++kPop) {
            if (kPop == jPop) {
              invM_S[iPop][jPop] += descriptors::invM<V,DESCRIPTOR>(iPop,kPop) * rt[kPop];
            }
          }
        }
      }
      const auto force = cell.template getField<descriptors::FORCE>();
      mrt<DESCRIPTOR>::addExternalForce(cell, rho, u, invM_S, force);
      return statistic;
    };
  };

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM, typename COLLISION>
  using combined_parameters = typename COLLISION::parameters;
};

}

}

#endif

#include "collisionMRT.cse.h"
