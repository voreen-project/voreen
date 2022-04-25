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

#ifndef DYNAMICS_INTERFACE_H
#define DYNAMICS_INTERFACE_H

#include "lbm.h"
#include "core/fieldParametersD.h"
#include "momenta/interface.h"

#include <type_traits>

namespace olb {

template <typename T, typename DESCRIPTOR> class Cell;
template <typename T, typename DESCRIPTOR> class ConstCell;

/// Return value of any collision
/**
 * Used for lattice statistics if enabled
 **/
template <typename T>
struct CellStatistic {
  T rho;
  T uSqr;

  operator bool() const any_platform {
    return rho != T{-1} && uSqr != T{-1};
  }
};

/// Interface for per-cell dynamics
template <typename T, typename DESCRIPTOR>
struct Dynamics {
  using value_t = T;
  using descriptor_t = DESCRIPTOR;

  virtual ~Dynamics() any_platform { }

  /// Expose unique type-identifier for RTTI
  virtual std::type_index id() = 0;
  /// Return human-readable name
  virtual std::string getName() const {
    return "Dynamics";
  };

  /// Initialize dynamics-specific data for cell
  virtual void initialize(Cell<T,DESCRIPTOR>& cell) { };

  /// Parameters access for legacy post processors
  virtual AbstractParameters<T,DESCRIPTOR>& getParameters(BlockLattice<T,DESCRIPTOR>& block) = 0;

  /// Perform purely-local collision step on Cell interface (legacy, to be deprecated)
  virtual CellStatistic<T> collide(Cell<T,DESCRIPTOR>& cell) {
    throw std::bad_function_call();
  };

  /// Compute particle density
  virtual T computeRho(ConstCell<T,DESCRIPTOR>& cell) const = 0;
  /// Compute fluid velocity
  virtual void computeU(ConstCell<T,DESCRIPTOR>& cell, T u[DESCRIPTOR::d]) const = 0;
  /// Compute fluid momentum
  virtual void computeJ(ConstCell<T,DESCRIPTOR>& cell, T j[DESCRIPTOR::d]) const = 0;
  /// Compute stress tensor
  virtual void computeStress(ConstCell<T,DESCRIPTOR>& cell,
                             T rho, const T u[DESCRIPTOR::d], T pi[util::TensorVal<DESCRIPTOR >::n]) const = 0;
  /// Compute fluid velocity and particle density
  virtual void computeRhoU(ConstCell<T,DESCRIPTOR>& cell, T& rho, T u[DESCRIPTOR::d]) const = 0;
  /// Compute all momenta  up to second order
  virtual void computeAllMomenta(ConstCell<T,DESCRIPTOR>& cell,
                                 T& rho, T u[DESCRIPTOR::d], T pi[util::TensorVal<DESCRIPTOR >::n]) const = 0;

  /// Set particle density
  virtual void defineRho(Cell<T,DESCRIPTOR>& cell, T rho) = 0;
  /// Set fluid velocity
  virtual void defineU(Cell<T,DESCRIPTOR>& cell, const T u[DESCRIPTOR::d]) = 0;
  /// Define fluid velocity and particle density
  virtual void defineRhoU(Cell<T,DESCRIPTOR>& cell, T rho, const T u[DESCRIPTOR::d]) = 0;
  /// Define all momenta up to second order
  virtual void defineAllMomenta(Cell<T,DESCRIPTOR>& cell,
                                T rho, const T u[DESCRIPTOR::d], const T pi[util::TensorVal<DESCRIPTOR >::n]) = 0;

  /// Return iPop equilibrium for given first and second momenta
  virtual T computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d]) const any_platform = 0;
  /// Initialize to equilibrium distribution
  void iniEquilibrium(Cell<T,DESCRIPTOR>& cell, T rho, const T u[DESCRIPTOR::d]) {
    T rhoShift(rho);
    T uShift[DESCRIPTOR::d];
    util::copyN(uShift, u, DESCRIPTOR::d);
    inverseShiftRhoU(cell, rhoShift, uShift);
    for (unsigned iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      cell[iPop] = computeEquilibrium(iPop, rhoShift, uShift);
    }
  };
  /// Initialize cell to equilibrium and non-equilibrum part
  void iniRegularized(Cell<T,DESCRIPTOR>& cell,
                      T rho, const T u[DESCRIPTOR::d], const T pi[util::TensorVal<DESCRIPTOR>::n]) {
    iniEquilibrium(cell, rho, u);
    for (unsigned iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      cell[iPop] += equilibrium<DESCRIPTOR>::template fromPiToFneq<T>(iPop, pi);
    }
  };
  /// Calculate population momenta s.t. the physical momenta are reproduced by the computeRhoU
  // This is relevant for correct initialization for some forcing/ sourcing schemes
  virtual void inverseShiftRhoU(ConstCell<T,DESCRIPTOR>& cell, T& rho, T u[DESCRIPTOR::d]) const {}
};

namespace dynamics {

/// Default combination rule used by dynamics::Tuple
/**
 * Use MOMENTA, EQUILIBRIUM and COLLISION elements as they are given
 **/
struct DefaultCombination {
  static std::string getName() {
    return "Default";
  }

  template <typename DESCRIPTOR, typename MOMENTA>
  using combined_momenta = typename MOMENTA::template type<DESCRIPTOR>;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM>
  using combined_equilibrium = typename EQUILIBRIUM::template type<DESCRIPTOR,MOMENTA>;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM, typename COLLISION>
  using combined_collision = typename COLLISION::template type<DESCRIPTOR,MOMENTA,EQUILIBRIUM>;

  template <typename DESCRIPTOR, typename MOMENTA, typename EQUILIBRIUM, typename COLLISION>
  using combined_parameters = typename COLLISION::parameters;
};


/// DYNAMICS is not explicitly marked as unvectorizable
template <typename DYNAMICS, typename = void>
struct is_vectorizable : std::true_type { };

/// DYNAMICS is explicitly marked as unvectorizable
/**
 * Required as a workaround to enable non-SIMD fallback in ConcreteBlockCollisionO
 **/
template <typename DYNAMICS>
struct is_vectorizable<
  DYNAMICS,
  std::enable_if_t<!DYNAMICS::is_vectorizable>
> : std::false_type { };

/// Evaluates to true iff DYNAMICS is not explicitly marked as unvectorizable
template <typename DYNAMICS>
static constexpr bool is_vectorizable_v = is_vectorizable<DYNAMICS>::value;


/// Dynamics constructed as a tuple of momenta, equilibrium and collision
/**
 * Optionally also a combination rule thereof (e.g. a forcing scheme)
 **/
template <
  typename T, typename DESCRIPTOR,
  typename MOMENTA, typename EQUILIBRIUM, typename COLLISION,
  typename COMBINATION_RULE = DefaultCombination
>
class Tuple final : public Dynamics<T,DESCRIPTOR> {
private:
  using Parameters = typename COMBINATION_RULE::template combined_parameters<DESCRIPTOR,MOMENTA,EQUILIBRIUM,COLLISION>;

public:
  using MomentaF     = typename COMBINATION_RULE::template combined_momenta<DESCRIPTOR,MOMENTA>;
  using EquilibriumF = typename COMBINATION_RULE::template combined_equilibrium<DESCRIPTOR,MOMENTA,EQUILIBRIUM>;
  using CollisionO   = typename COMBINATION_RULE::template combined_collision<DESCRIPTOR,MOMENTA,EQUILIBRIUM,COLLISION>;
  using ParametersD  = typename Parameters::template decompose_into<
    ParametersD<T,DESCRIPTOR>::template include_fields
  >;

  constexpr static bool is_vectorizable = is_vectorizable_v<CollisionO>;

  template <typename NEW_MOMENTA>
  using exchange_momenta = Tuple<
    T, DESCRIPTOR,
    NEW_MOMENTA, EQUILIBRIUM, COLLISION,
    COMBINATION_RULE
  >;

  template <typename NEW_RULE>
  using exchange_combination_rule = Tuple<
    T, DESCRIPTOR,
    MOMENTA, EQUILIBRIUM, COLLISION,
    NEW_RULE
  >;

  std::type_index id() override {
    return typeid(Tuple);
  }

  std::string getName() const override {
    return "dynamics::Tuple<"
      + MomentaF().getName() + ","
      + EQUILIBRIUM::getName() + ","
      + COLLISION::getName() + ","
      + COMBINATION_RULE::getName() +
    ">";
  };

  /// Return true iff FIELD is a parameter
  template <typename FIELD>
  constexpr bool hasParameter() const {
    return Parameters::template contains<FIELD>();
  };

  /// Interim workaround for accessing dynamics parameters in legacy post processors
  AbstractParameters<T,DESCRIPTOR>& getParameters(BlockLattice<T,DESCRIPTOR>& block) override {
    return block.template getData<DynamicsParameters<Tuple>>();
  }

  /// Initialize MOMENTA-specific data for cell
  void initialize(Cell<T,DESCRIPTOR>& cell) override {
    MomentaF().initialize(cell);
  };

  /// Apply purely-local collision step to a generic CELL
  template <typename CELL, typename PARAMETERS, typename V=typename CELL::value_t>
  CellStatistic<typename CELL::value_t> apply(CELL& cell, PARAMETERS& parameters) any_platform {
    auto params = static_cast<ParametersD&>(parameters).template copyAs<V>();
    return CollisionO().apply(cell, params);
  };

  T computeEquilibrium(int iPop, T rho, const T u[DESCRIPTOR::d]) const override any_platform {
    return EquilibriumF().compute(iPop, rho, u);
  };

  T computeRho(ConstCell<T,DESCRIPTOR>& cell) const override {
    return MomentaF().computeRho(cell);
  };
  void computeU(ConstCell<T,DESCRIPTOR>& cell, T u[DESCRIPTOR::d]) const override {
    MomentaF().computeU(cell, u);
  };
  void computeJ(ConstCell<T,DESCRIPTOR>& cell, T j[DESCRIPTOR::d]) const override {
    MomentaF().computeJ(cell, j);
  };
  void computeStress(ConstCell<T,DESCRIPTOR>& cell,
                     T rho, const T u[DESCRIPTOR::d], T pi[util::TensorVal<DESCRIPTOR >::n]) const override{
    MomentaF().computeStress(cell, rho, u, pi);
  };
  void computeRhoU(ConstCell<T,DESCRIPTOR>& cell, T& rho, T u[DESCRIPTOR::d]) const override {
    MomentaF().computeRhoU(cell, rho, u);
  };
  void computeAllMomenta(ConstCell<T,DESCRIPTOR>& cell,
                         T& rho, T u[DESCRIPTOR::d], T pi[util::TensorVal<DESCRIPTOR >::n]) const override {
    MomentaF().computeAllMomenta(cell, rho, u, pi);
  };

  void defineRho(Cell<T,DESCRIPTOR>& cell, T rho) override {
    MomentaF().defineRho(cell, rho);
  };
  void defineU(Cell<T,DESCRIPTOR>& cell, const T u[DESCRIPTOR::d]) override {
    MomentaF().defineU(cell, u);
  };
  void defineRhoU(Cell<T,DESCRIPTOR>& cell, T rho, const T u[DESCRIPTOR::d]) override {
    MomentaF().defineRhoU(cell, rho, u);
  };
  void defineAllMomenta(Cell<T,DESCRIPTOR>& cell,
                        T rho, const T u[DESCRIPTOR::d], const T pi[util::TensorVal<DESCRIPTOR >::n]) override {
    MomentaF().defineAllMomenta(cell, rho, u, pi);
  };
  void inverseShiftRhoU(ConstCell<T,DESCRIPTOR>& cell, T& rho, T u[DESCRIPTOR::d]) const override {
    MomentaF().inverseShiftRhoU(cell, rho, u);
  }
};

template <typename T, typename DESCRIPTOR, typename MOMENTA>
struct CustomCollision : public Dynamics<T,DESCRIPTOR> {
  using MomentaF = typename MOMENTA::template type<DESCRIPTOR>;

  void initialize(Cell<T,DESCRIPTOR>& cell) override {
    MomentaF().initialize(cell);
  };

  T computeRho(ConstCell<T,DESCRIPTOR>& cell) const override {
    return MomentaF().computeRho(cell);
  };
  void computeU(ConstCell<T,DESCRIPTOR>& cell, T u[DESCRIPTOR::d]) const override {
    MomentaF().computeU(cell, u);
  };
  void computeJ(ConstCell<T,DESCRIPTOR>& cell, T j[DESCRIPTOR::d]) const override {
    MomentaF().computeJ(cell, j);
  };
  void computeStress(ConstCell<T,DESCRIPTOR>& cell,
                     T rho, const T u[DESCRIPTOR::d], T pi[util::TensorVal<DESCRIPTOR >::n]) const override {
    MomentaF().computeStress(cell, rho, u, pi);
  };
  void computeRhoU(ConstCell<T,DESCRIPTOR>& cell, T& rho, T u[DESCRIPTOR::d]) const override {
    MomentaF().computeRhoU(cell, rho, u);
  };
  void computeAllMomenta(ConstCell<T,DESCRIPTOR>& cell,
                         T& rho, T u[DESCRIPTOR::d], T pi[util::TensorVal<DESCRIPTOR >::n]) const override {
    MomentaF().computeAllMomenta(cell, rho, u, pi);
  };

  void defineRho(Cell<T,DESCRIPTOR>& cell, T rho) override {
    MomentaF().defineRho(cell, rho);
  };
  void defineU(Cell<T,DESCRIPTOR>& cell, const T u[DESCRIPTOR::d]) override {
    MomentaF().defineU(cell, u);
  };
  void defineRhoU(Cell<T,DESCRIPTOR>& cell, T rho, const T u[DESCRIPTOR::d]) override {
    MomentaF().defineRhoU(cell, rho, u);
  };
  void defineAllMomenta(Cell<T,DESCRIPTOR>& cell,
                        T rho, const T u[DESCRIPTOR::d], const T pi[util::TensorVal<DESCRIPTOR >::n]) override {
    MomentaF().defineAllMomenta(cell, rho, u, pi);
  };
  void inverseShiftRhoU(ConstCell<T,DESCRIPTOR>& cell, T& rho, T u[DESCRIPTOR::d]) const override {
    MomentaF().inverseShiftRhoU(cell, rho, u);
  }
};


/// DYNAMICS doesn't provide apply method template
template <typename DYNAMICS, typename CELL, typename PARAMETERS, typename = void>
struct is_generic : std::false_type { };

/// DYNAMICS provides apply method template
template <typename DYNAMICS, typename CELL, typename PARAMETERS>
struct is_generic<
  DYNAMICS, CELL, PARAMETERS,
  std::enable_if_t<std::is_member_function_pointer_v<decltype(&DYNAMICS::template apply<CELL,PARAMETERS>)>>
> : std::true_type { };

/// Return true iff DYNAMICS provides apply method template
/**
 * i.e. it is not legacy
 **/
template <typename T, typename DESCRIPTOR, typename DYNAMICS>
static constexpr bool is_generic_v = is_generic<DYNAMICS, Cell<T,DESCRIPTOR>, AbstractParameters<T,DESCRIPTOR>>::value;

/// DYNAMICS is not explicitly marked as requiring parameters outside DYNAMICS::apply
template <typename DYNAMICS, typename = void>
struct has_parametrized_momenta : std::false_type { };

/// DYNAMICS is explicitly marked as requiring parameters outside DYNAMICS::apply
/**
 * Required to support ForcedPSMBGKdynamics
 **/
template <typename DYNAMICS>
struct has_parametrized_momenta<
  DYNAMICS,
  std::enable_if_t<DYNAMICS::has_parametrized_momenta>
> : std::true_type { };

/// Evaluates to true iff DYNAMICS doesn't require parameters outside DYNAMICS::apply
template <typename DYNAMICS>
static constexpr bool has_parametrized_momenta_v = has_parametrized_momenta<DYNAMICS>::value;

}

}

#endif
