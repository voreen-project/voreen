/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2008 Orestis Malaspinas, Andrea Parmigiani
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

#ifndef ADVECTION_DIFFUSION_BOUNDARIES_H
#define ADVECTION_DIFFUSION_BOUNDARIES_H

#include "dynamics/latticeDescriptors.h"
#include "dynamics/advectionDiffusionDynamics.h"
#include "dynamics/dynamics.h"

namespace olb {

//===================================================================================
//================= AdvectionDiffusionDynamcison Flat Boundaries =========
//===================================================================================
template<typename T, typename DESCRIPTOR, typename DYNAMICS, typename MOMENTA, int direction, int orientation>
struct AdvectionDiffusionBoundariesDynamics final : public dynamics::CustomCollision<T,DESCRIPTOR,MOMENTA> {
  using MomentaF    = typename MOMENTA::template type<DESCRIPTOR>;
  using ParametersD = typename DYNAMICS::ParametersD;

  template <typename M>
  using exchange_momenta = AdvectionDiffusionBoundariesDynamics<T,DESCRIPTOR,DYNAMICS,M,direction,orientation>;

  std::type_index id() override {
    return typeid(AdvectionDiffusionBoundariesDynamics);
  };

  AbstractParameters<T,DESCRIPTOR>& getParameters(BlockLattice<T,DESCRIPTOR>& block) override {
    return block.template getData<DynamicsParameters<AdvectionDiffusionBoundariesDynamics>>();
  }

  template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
  CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) {
    typedef DESCRIPTOR L;

    V dirichletTemperature = MomentaF().computeRho(cell);

    // Placeholder - the index calculation of this section can and should be done at compile time
    constexpr auto tmpUnknownIndexes = util::subIndexOutgoing<L, direction, orientation>();
    std::vector<int> unknownIndexes(tmpUnknownIndexes.cbegin(), tmpUnknownIndexes.cend());
    std::vector<int> knownIndexes = util::remainingIndexes<L>(unknownIndexes);

    int missingNormal = 0;

    if constexpr ((L::d == 3 && L::q == 7)||(L::d == 2 && L::q == 5)) {
      V sum = V{0};
      for (unsigned i = 0; i < knownIndexes.size(); ++i) {
        sum += cell[knownIndexes[i]];
      }

      V difference = dirichletTemperature - V{1} - sum; // on cell there are non-shiftet values -> temperature has to be changed

      // here I know all missing and non missing f_i
      for (unsigned i = 0; i < unknownIndexes.size(); ++i) {
        int numOfNonNullComp = 0;
        for (int iDim = 0; iDim < L::d; ++iDim) {
          numOfNonNullComp += util::abs(descriptors::c<L>(unknownIndexes[i],iDim));
        }
        if (numOfNonNullComp == 1) {
          missingNormal = unknownIndexes[i];
          // here missing diagonal directions are erased
          // just the normal direction stays (D3Q7)
          unknownIndexes.erase(unknownIndexes.begin() + i);
          break;

        }
      }
      cell[missingNormal] = difference; // on cell there are non-shiftet values -> temperature has to be changed
      return typename DYNAMICS::template exchange_momenta<MOMENTA>::CollisionO().apply(cell, parameters);
    }
    else {
      auto u = cell.template getField<descriptors::VELOCITY>();
      // part for q=19 copied from AdvectionDiffusionEdgesDynamics.collide()
      // but here just all missing directions, even at border of inlet area
      // has to be checked!
      for (unsigned iteratePop = 0; iteratePop < unknownIndexes.size();
           ++iteratePop) {
        cell[unknownIndexes[iteratePop]] =
          equilibrium<DESCRIPTOR>::template firstOrder(unknownIndexes[iteratePop], dirichletTemperature, u)
          - (cell[util::opposite<L>(unknownIndexes[iteratePop])]
             - equilibrium<DESCRIPTOR>::template firstOrder(
               util::opposite<L>(unknownIndexes[iteratePop]),
               dirichletTemperature, u));
      }
      return {-1,-1};
    }
  };

  T computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d]) const override {
    return equilibrium<DESCRIPTOR>::template firstOrder(iPop, rho, u);
  };

  std::string getName() const override {
    return "AdvectionDiffusionBoundariesDynamics";
  };

};

//===================================================================================
//================= AdvectionDiffusionDynamcis On Edges =========
//===================================================================================
template<typename T, typename DESCRIPTOR, typename DYNAMICS, typename MOMENTA, int plane, int normal1, int normal2>
class AdvectionDiffusionEdgesDynamics final : public dynamics::CustomCollision<T,DESCRIPTOR,MOMENTA> {
public:
  using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;
  using ParametersD = typename DYNAMICS::ParametersD;

  template <typename M>
  using exchange_momenta = AdvectionDiffusionEdgesDynamics<T,DESCRIPTOR,DYNAMICS,M,plane,normal1,normal2>;

  std::type_index id() override {
    return typeid(AdvectionDiffusionEdgesDynamics);
  };

  AbstractParameters<T,DESCRIPTOR>& getParameters(BlockLattice<T,DESCRIPTOR>& block) override {
    return block.template getData<DynamicsParameters<AdvectionDiffusionEdgesDynamics>>();
  }

  template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
  CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) {
    typedef DESCRIPTOR L;

    V temperature = MomentaF().computeRho(cell);
    auto u = cell.template getField<descriptors::VELOCITY>();
    // I need to get Missing information on the corners !!!!
    std::vector<int> unknownIndexes = util::subIndexOutgoing3DonEdges<L,plane,normal1,normal2>();
    // here I know all missing and non missing f_i

    // The collision procedure for D2Q5 and D3Q7 lattice is the same ...
    // Given the rule f_i_neq = -f_opposite(i)_neq
    // I have the right number of equations for the number of unknowns using these lattices

    for (unsigned iPop = 0; iPop < unknownIndexes.size(); ++iPop) {
      cell[unknownIndexes[iPop]] = equilibrium<DESCRIPTOR>::template firstOrder(unknownIndexes[iPop], temperature, u)
                                   -(cell[util::opposite<L>(unknownIndexes[iPop])]
                                     - equilibrium<DESCRIPTOR>::template firstOrder(util::opposite<L>(unknownIndexes[iPop]), temperature, u) );
    }

    // Once all the f_i are known, I can call the collision for the Regularized Model.
    return typename DYNAMICS::CollisionO().apply(cell, parameters);
  };

  T computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d]) const override {
    return equilibrium<DESCRIPTOR>::template firstOrder(iPop, rho, u);
  };

  std::string getName() const override {
    return "AdvectionDiffusionEdgesDynamics";
  };
};


//===================================================================================
//================= AdvectionDiffusionDynamics on  Corners for 2D Boundaries =========
//===================================================================================
template<typename T, typename DESCRIPTOR, typename DYNAMICS, typename MOMENTA, int xNormal, int yNormal>
struct AdvectionDiffusionCornerDynamics2D final : public dynamics::CustomCollision<T,DESCRIPTOR,MOMENTA> {
  using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;
  using ParametersD = typename DYNAMICS::ParametersD;

  template <typename M>
  using exchange_momenta = AdvectionDiffusionCornerDynamics2D<T,DESCRIPTOR,DYNAMICS,M,xNormal,yNormal>;

  std::type_index id() override {
    return typeid(AdvectionDiffusionCornerDynamics2D);
  };

  AbstractParameters<T,DESCRIPTOR>& getParameters(BlockLattice<T,DESCRIPTOR>& block) override {
    return block.template getData<DynamicsParameters<AdvectionDiffusionCornerDynamics2D>>();
  }

  template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
  CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) {
    typedef DESCRIPTOR L;

    V temperature = MomentaF().computeRho(cell);
    auto u = cell.template getField<descriptors::VELOCITY>();
    // I need to get Missing information on the corners !!!!
    std::vector<int> unknownIndexes = util::subIndexOutgoing2DonCorners<L,xNormal,yNormal>();
    // here I know all missing and non missing f_i


    // The collision procedure for D2Q5 and D3Q7 lattice is the same ...
    // Given the rule f_i_neq = -f_opposite(i)_neq
    // I have the right number of equations for the number of unknowns using these lattices

    for (unsigned iPop = 0; iPop < unknownIndexes.size(); ++iPop) {
      cell[unknownIndexes[iPop]] = equilibrium<DESCRIPTOR>::template firstOrder(unknownIndexes[iPop], temperature, u)
                                   -(cell[util::opposite<L>(unknownIndexes[iPop])]
                                     - equilibrium<DESCRIPTOR>::template firstOrder(util::opposite<L>(unknownIndexes[iPop]), temperature, u) ) ;
    }

    // Once all the f_i are known, I can call the collision for the Regularized Model.
    return typename DYNAMICS::CollisionO().apply(cell, parameters);
  };

  T computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d]) const override {
    return equilibrium<DESCRIPTOR>::template firstOrder(iPop, rho, u);
  };

  std::string getName() const override {
    return "AdvectionDiffusionCornerDynamics2D";
  };

};

//===================================================================================
//================= AdvectionDiffusionDynamics on  Corners for 3D Boundaries =========
//===================================================================================
template<typename T, typename DESCRIPTOR, typename DYNAMICS, typename MOMENTA, int xNormal, int yNormal, int zNormal>
struct AdvectionDiffusionCornerDynamics3D final : public dynamics::CustomCollision<T,DESCRIPTOR,MOMENTA> {

  using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;
  using ParametersD = typename DYNAMICS::ParametersD;

  template <typename M>
  using exchange_momenta = AdvectionDiffusionCornerDynamics3D<T,DESCRIPTOR,DYNAMICS,M,xNormal,yNormal,zNormal>;

  std::type_index id() override {
    return typeid(AdvectionDiffusionCornerDynamics3D);
  };

  AbstractParameters<T,DESCRIPTOR>& getParameters(BlockLattice<T,DESCRIPTOR>& block) override {
    return block.template getData<DynamicsParameters<AdvectionDiffusionCornerDynamics3D>>();
  }

  template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
  CellStatistic<V> apply(CELL& cell, PARAMETERS& parameters) {
    typedef DESCRIPTOR L;

    V temperature = MomentaF().computeRho(cell);
    auto u = cell.template getField<descriptors::VELOCITY>();
    // I need to get Missing information on the corners !!!!
    std::vector<int> unknownIndexes = util::subIndexOutgoing3DonCorners<L,xNormal,yNormal,zNormal>();
    // here I know all missing and non missing f_i

    // The collision procedure for D2Q5 and D3Q7 lattice is the same ...
    // Given the rule f_i_neq = -f_opposite(i)_neq
    // I have the right number of equations for the number of unknowns using these lattices

    for (unsigned iPop = 0; iPop < unknownIndexes.size(); ++iPop) {
      cell[unknownIndexes[iPop]] = equilibrium<DESCRIPTOR>::template firstOrder(unknownIndexes[iPop], temperature, u)
                                   -(cell[util::opposite<L>(unknownIndexes[iPop])]
                                     - equilibrium<DESCRIPTOR>::template firstOrder(util::opposite<L>(unknownIndexes[iPop]), temperature, u) ) ;
    }

    // Once all the f_i are known, I can call the collision for the Regularized Model.
    return typename DYNAMICS::template exchange_momenta<MOMENTA>::CollisionO().apply(cell, parameters);
  };

  T computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d]) const override {
    return equilibrium<DESCRIPTOR>::template firstOrder(iPop, rho, u);
  };

  std::string getName() const override {
    return "AdvectionDiffusionCornerDynamics3D";
  };

};


}  // namespace olb

#endif
