/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Julius Jessberger
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

#ifndef DYNAMICS_MOMENTA_DEFINITION_RULE_H
#define DYNAMICS_MOMENTA_DEFINITION_RULE_H

namespace olb {

namespace momenta {

/// The momenta are defined one after the other.
struct DefineSeparately {
  template <typename TYPE, typename CELL, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void defineRho(CELL& cell, V rho)
  { }

  template <typename TYPE, typename CELL, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void defineU(CELL& cell, const V u[DESCRIPTOR::d])
  { }

  template <typename TYPE, typename CELL, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void defineRhoU(CELL& cell,
               V rho, const V u[DESCRIPTOR::d])
  { }

  template <typename TYPE, typename CELL, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void defineAllMomenta(CELL& cell,
               V rho, const V u[DESCRIPTOR::d],
               const V pi[util::TensorVal<DESCRIPTOR >::n])
  { }

  static std::string getName(){
    return "DefineSeparately";
  }
};

/// When momenta are changed, a new equilibrium state is set.
/// This struct is specific for advection-diffusion applications since it
/// applies the external field descriptors::VELOCITY directly.
struct DefineToEq {
  template <typename TYPE, typename CELL, typename RHO, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void defineRho(CELL& cell, const RHO& rho)
  {
    // get fluid velocity u
    const auto u = cell.template getField<descriptors::VELOCITY>();

    // set new equilibrium state with rho and u
    lbm<DESCRIPTOR>::defineEqFirstOrder(cell, rho, u);
  }

  template <typename TYPE, typename CELL, typename U, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void defineU(CELL& cell, const U& u)
  { }

  template <typename TYPE, typename CELL, typename RHO, typename U, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void defineRhoU(CELL& cell,
                  const RHO& rho, const U& u)
  {
    defineRho<TYPE>(cell, rho);
  }

  template <typename TYPE, typename CELL, typename RHO, typename U, typename PI, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void defineAllMomenta(CELL& cell,
                        const RHO& rho, const U& u, const PI& pi)
  {
    defineRho<TYPE>(cell, rho);
  }

  static std::string getName(){
    return "DefineToEq";
  }
};

/// When momenta are changed, the equilibrium part of the population is
/// modified while the non-equilibrium part is kept.
struct DefineToNEq {
  template <typename TYPE, typename CELL, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void defineRho(CELL& cell, V rho)
  {
    // get old equilibrium data (oldRho, u)
    V oldRho, u[DESCRIPTOR::d];
    TYPE().computeRhoU(cell, oldRho, u);
    TYPE().inverseShiftRhoU(cell, oldRho, u);

    // modify the equilibrium part of the population from (oldRho, u) to (rho, u)
    lbm<DESCRIPTOR>::defineNEq(cell, oldRho, u, rho, u);
  }

  template <typename TYPE, typename CELL, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void defineU(CELL& cell, const V u[DESCRIPTOR::d])
  {
    // get old equilibrium data (rho, oldU)
    V rho, oldU[DESCRIPTOR::d];
    TYPE().computeRhoU(cell, rho, oldU);
    TYPE().inverseShiftRhoU(cell, rho, oldU);

    // modify the equilibrium part of the population from (rho, oldU) to (rho, u)
    lbm<DESCRIPTOR>::defineNEq(cell, rho, oldU, rho, u);
  }

  template <typename TYPE, typename CELL, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void defineRhoU(CELL& cell,
                  V rho, const V u[DESCRIPTOR::d])
  {
    V oldRho, oldU[DESCRIPTOR::d];
    TYPE().computeRhoU(cell, oldRho, oldU);
    TYPE().inverseShiftRhoU(cell, oldRho, oldU);

    lbm<DESCRIPTOR>::defineNEq(cell, oldRho, oldU, rho, u);
  }

  template <typename TYPE, typename CELL, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void defineAllMomenta(CELL& cell,
               V rho, const V u[DESCRIPTOR::d],
               const V pi[util::TensorVal<DESCRIPTOR >::n])
  {
    lbm<DESCRIPTOR>::defineNEqFromPi(cell, rho, u, pi);
  }

  static std::string getName(){
    return "DefineToNEq";
  }
};

/// defineRho leads to a new non-equilibrium population, defineU
/// only sets the velocity data.
struct DefineUSeparately {
  template <typename TYPE, typename CELL, typename RHO, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void defineRho(CELL& cell, RHO& rho)
  {
    DefineToNEq().defineRho<TYPE>(cell, rho);
  }

  template <typename TYPE, typename CELL, typename U, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void defineU(CELL& cell, U& u)
  { }

  template <typename TYPE, typename CELL, typename RHO, typename U, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void defineRhoU(CELL& cell, RHO& rho, U& u)
  {
    defineRho<TYPE>(cell, rho);
  }

  template <typename TYPE, typename CELL, typename RHO, typename U, typename PI, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void defineAllMomenta(CELL& cell,
                        RHO& rho, U& u, PI& pi)
  {
    DefineToNEq().defineAllMomenta<TYPE>(cell, rho, u, pi);
  }

  static std::string getName(){
    return "DefineUSeparately";
  }
};

/// defineRho leads to a new non-equilibrium population, defineU
/// only sets the velocity data.
// In defineRho, the computation of the old momenta does not use the external
// field, but uses lbHelpers::computeRhoU as in the bulk.
struct DefineUSeparatelyTrace {
  template <typename TYPE, typename CELL, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void defineRho(CELL& cell, V rho)
  {
    // get old equilibrium data (oldRho, u)
    V oldRho, u[DESCRIPTOR::d];
    // only this line is different to DefineUSeparately
    lbm<DESCRIPTOR>::computeRhoU(cell, oldRho, u);
    TYPE().inverseShiftRhoU(cell, oldRho, u);

    // modify the equilibrium part of the population from (oldRho, u) to (rho, u)
    lbm<DESCRIPTOR>::defineNEq(cell, oldRho, u, rho, u);
  }

  template <typename TYPE, typename CELL, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void defineU(CELL& cell, const V u[DESCRIPTOR::d])
  { }

  template <typename TYPE, typename CELL, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void defineRhoU(CELL& cell,
               V rho, const V u[DESCRIPTOR::d])
  {
    defineRho<TYPE>(cell, rho);
  }

  template <typename TYPE, typename CELL, typename V=typename CELL::value_t, typename DESCRIPTOR=typename CELL::descriptor_t>
  void defineAllMomenta(CELL& cell,
               V rho, const V u[DESCRIPTOR::d],
               const V pi[util::TensorVal<DESCRIPTOR >::n])
  {
    DefineUSeparately().defineAllMomenta<TYPE>(cell, rho, u, pi);
  }

  static std::string getName(){
    return "DefineUSeparatelyTrace";
  }
};

}

}

#endif
