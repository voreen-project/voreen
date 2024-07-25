/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007, 2013 Jonas Latt, Mathias J. Krause, Geng Liu
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
 * Helper functions for the implementation of LB dynamics. This file is all
 * about efficiency. The generic template code is specialized for commonly
 * used Lattices, so that a maximum performance can be taken out of each
 * case.
 */
#ifndef MRT_H
#define MRT_H

#include "lbm.h"
#include "mrtLatticeDescriptors.h"

namespace olb {

template<typename DESCRIPTOR>
struct mrt {
  static_assert(
    std::is_same<typename DESCRIPTOR::category_tag, descriptors::tag::MRT>::value,
    "DESCRIPTOR is tagged as MRT");

  /// Computation of equilibrium distribution (in momenta space)
  template <typename RHO, typename U, typename V=RHO>
  static V equilibrium(int iPop, const RHO& rho, const U& u)
  {
    V equ{};
    const V uSqr = util::normSqr<U,DESCRIPTOR::d>(u);
    for (int jPop = 0; jPop < DESCRIPTOR::q; ++jPop) {
      equ += descriptors::m<V,DESCRIPTOR>(iPop,jPop) *
             (olb::equilibrium<DESCRIPTOR>::secondOrder(jPop,rho,u,uSqr) +
              descriptors::t<V,DESCRIPTOR>(jPop));
    }
    return equ;
  }

  /// Computation of all equilibrium distribution (in momenta space)
  template <typename MOMENTAEQ, typename RHO, typename U, typename V=RHO>
  static void computeEquilibrium(MOMENTAEQ& momentaEq,
                                 const RHO& rho, const U& u)
  {
    const V uSqr = util::normSqr<U,DESCRIPTOR::d>(u);
    for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
      momentaEq[iPop] = V{};
      for (int jPop = 0; jPop < DESCRIPTOR::q; ++jPop) {
        momentaEq[iPop] += descriptors::m<V,DESCRIPTOR>(iPop,jPop) *
                           (olb::equilibrium<DESCRIPTOR>::secondOrder(jPop,rho,u,uSqr) +
                            descriptors::t<V,DESCRIPTOR>(jPop));
      }
    }
  }

  template <typename MOMENTA, typename CELL, typename V=typename CELL::value_t>
  static void computeMomenta(MOMENTA& momenta, CELL& cell)
  {
    for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
      momenta[iPop] = V{};
      for (int jPop = 0; jPop < DESCRIPTOR::q; ++jPop) {
        momenta[iPop] += descriptors::m<V,DESCRIPTOR>(iPop,jPop) *
                         (cell[jPop] + descriptors::t<V,DESCRIPTOR>(jPop));
      }
    }
  }

  /// MRT collision step
  template <typename CELL, typename RHO, typename U, typename INVM_S, typename V=typename CELL::value_t>
  static V mrtCollision(CELL& cell,
                        const RHO& rho, const U& u,
                        const INVM_S& invM_S)
  {
    V uSqr = util::normSqr<V,DESCRIPTOR::d>(u);
    V momenta[DESCRIPTOR::q] { };
    V momentaEq[DESCRIPTOR::q] { };

    computeMomenta(momenta, cell);
    computeEquilibrium(momentaEq, rho, u);

    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      V collisionTerm{};
      for (int jPop = 0; jPop < DESCRIPTOR::q; ++jPop) {
        collisionTerm += invM_S[iPop][jPop] * (momenta[jPop] - momentaEq[jPop]);
      }
      cell[iPop] -= collisionTerm;
    }
    return uSqr;
  }

  /// MRT SGS collision step
  template <typename CELL, typename RHO, typename U, typename OMEGA, typename INVM_S_SGS, typename V=typename CELL::value_t>
  static V mrtSGSCollision(CELL& cell,
                           const RHO& rho, const U& u,
                           const OMEGA& omega,
                           const INVM_S_SGS& invM_S_SGS)
  {
    V uSqr = util::normSqr<V,DESCRIPTOR::d>(u);
    V momenta[DESCRIPTOR::q];
    V momentaEq[DESCRIPTOR::q];

    computeMomenta(momenta, cell);
    computeEquilibrium(momentaEq, rho, u);

    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      V collisionTerm{};
      for (int jPop = 0; jPop < DESCRIPTOR::q; ++jPop) {
        collisionTerm += invM_S_SGS[iPop][jPop] * (momenta[jPop] - momentaEq[jPop]);
      }
      cell[iPop] -= collisionTerm;
    }

    return uSqr;
  }


  /// Ladd-Verberg-I body force model for MRT
  /// A.Ladd, R. Verberg, DESCRIPTOR-Boltzmann simulations of particle-fluid suspensions, Journal of Statistical Physics 104(2001)
  template <typename CELL, typename RHO, typename U, typename INVM_S, typename FORCE, typename V=typename CELL::value_t>
  static void addExternalForce(CELL& cell,
                               const RHO& rho, const U& u,
                               const INVM_S& invM_S,
                               const FORCE& force)
  {
    V f_u{};
    for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
      f_u += force[iD] * u[iD];
    }

    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      V c_u{};
      V c_f{};
      for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
        c_u += descriptors::c<DESCRIPTOR>(iPop,iD)*u[iD];
        c_f += descriptors::c<DESCRIPTOR>(iPop,iD)*force[iD];
      }
      V f1 = descriptors::t<V,DESCRIPTOR>(iPop)*rho*c_f*descriptors::invCs2<V,DESCRIPTOR>();
      V f2 = descriptors::t<V,DESCRIPTOR>(iPop)*rho*(c_u*c_f*descriptors::invCs2<V,DESCRIPTOR>()-f_u)*descriptors::invCs2<V,DESCRIPTOR>();

      V invMsM{};
      for (int jPop=0; jPop < DESCRIPTOR::q; ++jPop) {
        invMsM += invM_S[iPop][jPop]*descriptors::m<V,DESCRIPTOR>(jPop,iPop);
      }
      cell[iPop] += f1 + f2 - invMsM*f2*V{0.5};
    }
  }

};

}

#endif
