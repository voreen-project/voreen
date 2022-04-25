/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006 Jonas Latt, 2021 Julius Jessberger, Adrian Kummerlaender
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
 * File provides a generic interface for the computation and definition of
 * momenta (density, velocity, stress).
 * Struct Momenta provides a general (virtual) interface.
 * Struct Tuple reduces the computation and definition on the computation
 * and definition of the single density, velocity and stress "submomenta".
 * The definition of e.g. the density also depends on a rule which says how to
 * modify the cell population in that case, e.g. a new non-equilibrium state
 * could be settled. This behavior is defined in DefinitionRule structures.
 */

#ifndef DYNAMICS_MOMENTA_INTERFACE_H
#define DYNAMICS_MOMENTA_INTERFACE_H

#include <type_traits>
#include <functional>

#include "elements.h"

namespace olb {

namespace momenta {

/// Partially-specializable rho and u computation
/**
 * Currently required in order to perform basic CSE of this critical computation
 **/
template <typename BASE, typename DENSITY, typename MOMENTUM>
struct ComputeRhoU {
  template <typename CELL, typename RHO, typename U>
  void operator()(CELL& cell, RHO& rho, U& u) const any_platform
  {
    DENSITY().template compute<BASE>(cell, rho);
    MOMENTUM().template computeU<BASE>(cell, u);
  };
};

template <typename BASE>
struct ComputeRhoU<BASE,BulkDensity,BulkMomentum> {
  template <typename CELL, typename RHO, typename U, typename DESCRIPTOR=typename CELL::descriptor_t>
  void operator()(CELL& cell, RHO& rho, U& u) const any_platform
  {
    lbm<DESCRIPTOR>::computeRhoU(cell, rho, u);
  };
};

/// Tuple of momenta components forming a moment system
/**
 * Momenta are reduced on their components DENSITY, MOMETUM, STRESS
 * A generic momenta tuple consists of single "submomenta" (0th moment =
 * density, 1st velocity, 2nd moment = stress).
 * These give a rule how to compute a certain moment and how to set the
 * corresponding data, e.g. a velocity field, in case that this is necessary.
 * DefinitionRule describes how to modify the momenta, e.g. the populations
 * could be transformed to a new equilibrium or non-equilibrium state.
 */
template <
  typename DESCRIPTOR,
  typename DENSITY,
  typename MOMENTUM,
  typename STRESS,
  typename DefinitionRule
>
struct ConcreteTuple {
  using descriptor = DESCRIPTOR;
  using type = ConcreteTuple;
  using density = DENSITY;
  using momentum = MOMENTUM;
  using stress = STRESS;
  using definition = DefinitionRule;

  template <typename CELL, typename V=typename CELL::value_t>
  V computeRho(CELL& cell) const any_platform
  {
    V rho{};
    DENSITY().template compute<type>(cell, rho);
    return rho;
  }

  template <typename CELL, typename U>
  void computeU(CELL& cell, U& u) const any_platform
  {
    MOMENTUM().template computeU<type>(cell, u);
  }

  template <typename CELL, typename J>
  void computeJ(CELL& cell, J& j) const any_platform
  {
    MOMENTUM().template compute<type>(cell, j);
  }

  // TODO: Drop the rho and u args here, this is just another opportunity for
  // injecting false values. Replaces calls with computeAllMomenta
  template <typename CELL, typename RHO, typename U, typename PI>
  void computeStress(CELL& cell,
                     const RHO& rho, const U& u, PI& pi) const any_platform
  {
    STRESS().template compute<type>(cell, rho, u, pi);
  }

  //computeStress without rho and u args
  template <typename CELL, typename PI, typename V=typename CELL::value_t>
  void computeStress(CELL& cell, PI& pi) const any_platform
  {
    V rho, u[DESCRIPTOR::d];
    computeRhoU(cell, rho, u);
    STRESS().template compute<type>(cell, rho, u, pi);
  }

  template <typename CELL, typename RHO, typename J>
  void computeRhoJ(CELL& cell, RHO& rho, J& j) const any_platform
  {
    DENSITY().template compute<type>(cell, rho);
    MOMENTUM().template compute<type>(cell, j);
  }

  template <typename CELL, typename RHO, typename U>
  void computeRhoU(CELL& cell, RHO& rho, U& u) const any_platform
  {
    ComputeRhoU<type,DENSITY,MOMENTUM>()(cell, rho, u);
  }

  template <typename CELL, typename RHO, typename U, typename PI>
  void computeAllMomenta(CELL& cell, RHO& rho, U& u, PI& pi) const any_platform
  {
    computeRhoU(cell, rho, u);
    STRESS().template compute<type>(cell, rho, u, pi);
  }

  template <typename CELL, typename PINEQNORMSQR, typename V=typename CELL::value_t>
  void computePiNeqNormSqr(CELL& cell, PINEQNORMSQR& piNeqNormSqr) const any_platform
  {
    V rho, u[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR>::n] { };
    computeRhoU(cell, rho, u);
    STRESS().template compute<type>(cell, rho, u, pi);
    piNeqNormSqr = pi[0]*pi[0] + 2.*pi[1]*pi[1] + pi[2]*pi[2];
    if constexpr (util::TensorVal<DESCRIPTOR>::n == 6) {
      piNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2.*pi[4]*pi[4] +pi[5]*pi[5];
    }
  }

  template <typename CELL, typename RHO>
  void defineRho(CELL& cell, const RHO& rho)
  {
    DENSITY().template define<type>(cell, rho);

    RHO rhoShift = rho;
    inverseShiftRho(cell, rhoShift);
    DefinitionRule().template defineRho<type>(cell, rho);
  }

  template <typename CELL, typename U>
  void defineU(CELL& cell, const U& u)
  {
    MOMENTUM().template define<type>(cell, u);

    using U_ARITH = std::remove_const_t<std::remove_pointer_t<std::remove_extent_t<U>>>;
    U_ARITH uShift[DESCRIPTOR::d];
    util::copyN(uShift, u, DESCRIPTOR::d);
    inverseShiftU(cell, uShift);
    DefinitionRule().template defineU<type>(cell, u);
  }

  template <typename CELL, typename RHO, typename U>
  void defineRhoU(CELL& cell, const RHO& rho, const U& u)
  {
    DENSITY().template define<type>(cell, rho);
    MOMENTUM().template define<type>(cell, u);

    RHO rhoShift = rho;
    using U_ARITH = std::remove_const_t<std::remove_pointer_t<std::remove_extent_t<U>>>;
    U_ARITH uShift[DESCRIPTOR::d];
    util::copyN(uShift, u, DESCRIPTOR::d);
    inverseShiftRhoU(cell, rhoShift, uShift);
    DefinitionRule().template defineRhoU<type>(cell, rhoShift, uShift);
  }

  template <typename CELL, typename RHO, typename U, typename PI>
  void defineAllMomenta(CELL& cell, const RHO& rho, const U& u, const PI& pi)
  {
    DENSITY().template define<type>(cell, rho);
    MOMENTUM().template define<type>(cell, u);

    RHO rhoShift = rho;
    using U_ARITH = std::remove_const_t<std::remove_pointer_t<std::remove_extent_t<U>>>;
    U_ARITH uShift[DESCRIPTOR::d];
    util::copyN(uShift, u, DESCRIPTOR::d);
    inverseShiftRhoU(cell, rhoShift, uShift);
    DefinitionRule().template defineAllMomenta<type>(cell, rhoShift, uShift, pi);
  }

  template <typename CELL>
  void initialize(CELL& cell)
  {
    DENSITY().template initialize<type>(cell);
    MOMENTUM().template initialize<type>(cell);
  }

  template <typename CELL, typename RHO>
  void inverseShiftRho(CELL& cell, RHO& rho)
  {
    DENSITY().template inverseShift<type>(cell, rho);
  }

  template <typename CELL, typename U>
  void inverseShiftU(CELL& cell, U& u)
  {
    MOMENTUM().template inverseShift<type>(cell, u);
  }

  template <typename CELL, typename RHO, typename U>
  void inverseShiftRhoU(CELL& cell, RHO& rho, U& u)
  {
    inverseShiftRho(cell, rho);
    inverseShiftU(cell, u);
  }

  std::string getName() const
  {
    return "Momenta<"
      + density().getName() + ","
      + momentum().getName() + ","
      + stress().getName() + ","
      + definition().getName() +
    ">";
  }
};

template <
  typename DENSITY,
  typename MOMENTUM,
  typename STRESS,
  typename DefinitionRule
>
struct Tuple {
  using density    = DENSITY;
  using momentum   = MOMENTUM;
  using stress     = STRESS;
  using definition = DefinitionRule;

  template <typename DESCRIPTOR>
  using type = ConcreteTuple<DESCRIPTOR,DENSITY,MOMENTUM,STRESS,DefinitionRule>;
};

template <typename MOMENTA>
struct Forced {
  template <typename DESCRIPTOR>
  using type = ConcreteTuple<
    DESCRIPTOR,
    typename MOMENTA::template type<DESCRIPTOR>::density,
    ForcedMomentum<typename MOMENTA::template type<DESCRIPTOR>::momentum>,
    typename MOMENTA::template type<DESCRIPTOR>::stress,
    typename MOMENTA::template type<DESCRIPTOR>::definition
  >;
};

template <typename MOMENTA>
struct ForcedWithStress {
  template <typename DESCRIPTOR>
  using type = ConcreteTuple<
    DESCRIPTOR,
    typename MOMENTA::template type<DESCRIPTOR>::density,
    ForcedMomentum<typename MOMENTA::template type<DESCRIPTOR>::momentum>,
    ForcedStress<typename MOMENTA::template type<DESCRIPTOR>::stress>,
    typename MOMENTA::template type<DESCRIPTOR>::definition
  >;
};

template <typename MOMENTA>
struct Porous {
  template <typename DESCRIPTOR>
  using type = ConcreteTuple<
    DESCRIPTOR,
    typename MOMENTA::template type<DESCRIPTOR>::density,
    PorousMomentum<typename MOMENTA::template type<DESCRIPTOR>::momentum>,
    typename MOMENTA::template type<DESCRIPTOR>::stress,
    typename MOMENTA::template type<DESCRIPTOR>::definition
  >;
};

template <typename MOMENTA>
struct PorousParticle {
  template <typename DESCRIPTOR>
  using type = ConcreteTuple<
    DESCRIPTOR,
    typename MOMENTA::template type<DESCRIPTOR>::density,
    PorousParticleMomentum<typename MOMENTA::template type<DESCRIPTOR>::momentum>,
    typename MOMENTA::template type<DESCRIPTOR>::stress,
    typename MOMENTA::template type<DESCRIPTOR>::definition
  >;
};

}  // namespace momenta

}  // namespace olb

#endif
