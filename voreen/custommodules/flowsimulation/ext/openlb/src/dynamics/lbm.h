/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Jonas Latt
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

#ifndef DYNAMICS_LBM_H
#define DYNAMICS_LBM_H

#include "core/util.h"
#include "descriptorFunction.h"

namespace olb {

template <typename DESCRIPTOR>
struct equilibrium {
  /// Computation of equilibrium distribution, first order in u
  template <typename RHO, typename U, typename V=RHO>
  static V firstOrder(int iPop, const RHO& rho, const U& u) any_platform
  {
    V c_u{};
    for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
      c_u += descriptors::c<DESCRIPTOR>(iPop,iD)*u[iD];
    }
    return rho
           * descriptors::t<V,DESCRIPTOR>(iPop)
           * ( V{1} + c_u * descriptors::invCs2<V,DESCRIPTOR>() )
           - descriptors::t<V,DESCRIPTOR>(iPop);
  }

  /// Computation of equilibrium distribution, second order in u
  template <typename RHO, typename U, typename USQR, typename V=RHO>
  static V secondOrder(int iPop, const RHO& rho, const U& u, const USQR& uSqr) any_platform
  {
    V c_u{};
    for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
      c_u += descriptors::c<DESCRIPTOR>(iPop,iD)*u[iD];
    }
    return rho
           * descriptors::t<V,DESCRIPTOR>(iPop)
           * ( V{1}
               + descriptors::invCs2<V,DESCRIPTOR>() * c_u
               + descriptors::invCs2<V,DESCRIPTOR>() * descriptors::invCs2<V,DESCRIPTOR>() * V{0.5} * c_u * c_u
               - descriptors::invCs2<V,DESCRIPTOR>() * V{0.5} * uSqr)
           - descriptors::t<V,DESCRIPTOR>(iPop);
  }
  /// Computation of equilibrium distribution, second order in u
  template <typename RHO, typename U, typename V=RHO>
  static V secondOrder(int iPop, const RHO& rho, const U& u) any_platform
  {
    const V uSqr = util::normSqr<U,DESCRIPTOR::d>(u);
    return secondOrder(iPop, rho, u, uSqr);
  }

  // compute equilibrium f^eq_i eq. (5.32) from DOI:10.1002/9780470177013
  template <typename RHO, typename U, typename V=RHO>
  static V P1(int iPop, const RHO& rho, const U& u) any_platform
  {
    V c_u{};
    // compute scalar product of c[iPop]*u
    for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
      c_u += descriptors::c<DESCRIPTOR>(iPop,iD)*u[iD];
    }
    return descriptors::t<V,DESCRIPTOR>(iPop) * (rho + c_u)
         - descriptors::t<V,DESCRIPTOR>(iPop);
  }

  template <typename J, typename JSQR, typename PRESSURE, typename V=PRESSURE>
  static V incompressible(int iPop, const J& j, const JSQR& jSqr, const PRESSURE& pressure) any_platform
  {
    V c_j{};
    for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
      c_j += descriptors::c<DESCRIPTOR>(iPop,iD)*j[iD];
    }
    return descriptors::t<V,DESCRIPTOR>(iPop)
           * ( descriptors::invCs2<V,DESCRIPTOR>() * pressure
               + descriptors::invCs2<V,DESCRIPTOR>() * c_j
               + descriptors::invCs2<V,DESCRIPTOR>() * descriptors::invCs2<V,DESCRIPTOR>()/V{2} * c_j * c_j
               - descriptors::invCs2<V,DESCRIPTOR>()/V{2} * jSqr )
           - descriptors::t<V,DESCRIPTOR>(iPop);
  }
  template <typename J, typename PRESSURE, typename V=PRESSURE>
  static V incompressible(int iPop, const J& j, const PRESSURE& pressure) any_platform
  {
    const V jSqr = util::normSqr<J,DESCRIPTOR::d>(j);
    return incompressible(iPop, j, jSqr, pressure);
  }

  /// compute off-equilibrium part of the populations from gradient of the flux
  /// for asymmetric regularization init to circumvent pi computation
  template <typename V>
  static V fromJgradToFneq(int iPop,
                           const V Jgrad[DESCRIPTOR::d * DESCRIPTOR::d],
                           V omega) any_platform
  {
    using L = DESCRIPTOR;
    V fNeq{};
    int iJgrad = 0;
    for (int iAlpha=0; iAlpha < L::d; ++iAlpha) {
      for (int iBeta=0; iBeta < L::d; ++iBeta) {
        V toAdd = descriptors::c<L>(iPop,iAlpha) * descriptors::c<L>(iPop,iBeta);
        if (iAlpha == iBeta) {
          toAdd -= 1./descriptors::invCs2<V,L>();
        }
        else {
          toAdd *= V{2};
        }
        toAdd *= Jgrad[iJgrad++];
        fNeq += toAdd;
      }
    }
    fNeq *= - descriptors::t<V,L>(iPop) * descriptors::invCs2<V,L>() / omega;
    return fNeq;
  }

  /// Compute off-equilibrium part of the f's from the stress tensor Pi.
  /** Implements the following formula (with Einstein index contraction):
   * \f[ f_i^{neq} = t_i / (2 c_s^4) *
   *                 (c_{ia} c_{ib} - c_s^2 \delta_{ab}) \Pi_{ab} \f]
   * By Pi we mean the tensor computed from the off-equilibrium functions:
   * \f[ \Pi = \sum c_i c_i f_i^{neq}
   *         = \sum c_i c_i f_i - \rho u u - c_s^2 \rho\ Id \f]
   */
  template <typename V, typename PI>
  static V fromPiToFneq(int iPop, const PI& pi) any_platform
  {
    using L = DESCRIPTOR;
    V fNeq{};
    int iPi = 0;
    // Iterate only over superior triangle + diagonal, and add
    // the elements under the diagonal by symmetry
    for (int iAlpha=0; iAlpha < L::d; ++iAlpha) {
      for (int iBeta=iAlpha; iBeta < L::d; ++iBeta) {
        V toAdd = descriptors::c<L>(iPop,iAlpha)*descriptors::c<L>(iPop,iBeta);
        if (iAlpha == iBeta) {
          toAdd -= V{1}/descriptors::invCs2<V,L>();
        }
        else {
          toAdd *= V{2}; // multiply off-diagonal elements by 2
        }                // because the Q tensor is symmetric
        toAdd *= pi[iPi++];
        fNeq += toAdd;
      }
    }
    fNeq *= descriptors::t<V,L>(iPop) * descriptors::invCs2<V,L>() * descriptors::invCs2<V,L>() / V{2};
    return fNeq;
  }

  template <typename V>
  static V fromJneqToFneq(int iPop, const V jNeq[DESCRIPTOR::d]) any_platform
  {
    V fNeq{};
    for (int iD = 0; iD < DESCRIPTOR::d; ++iD) {
      fNeq += descriptors::c<DESCRIPTOR>(iPop,iD) * jNeq[iD];
    }
    fNeq *= descriptors::t<V,DESCRIPTOR>(iPop) * descriptors::invCs2<V,DESCRIPTOR>();
    return fNeq;
  }

};

/// Collection of common computations for LBM
template <typename DESCRIPTOR>
struct lbm {
  /// Computation of density
  template <typename CELL, typename V=typename CELL::value_t>
  static V computeRho(CELL& cell) any_platform
  {
    V rho = V();
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      rho += cell[iPop];
    }
    rho += V{1};
    return rho;
  }

  /// Computation of momentum
  template <typename CELL, typename J, typename V=typename CELL::value_t>
  static void computeJ(CELL& cell, J& j) any_platform
  {
    for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
      j[iD] = V();
    }
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
        j[iD] += cell[iPop]*descriptors::c<DESCRIPTOR>(iPop,iD);
      }
    }
  }

  /// Computation of hydrodynamic variables
  template <typename CELL, typename RHO, typename J, typename V=typename CELL::value_t>
  static void computeRhoJ(CELL& cell, RHO& rho, J& j) any_platform
  {
    rho = computeRho(cell);
    computeJ(cell, j);
  }

  /// Computation of hydrodynamic variables
  template <typename CELL, typename RHO, typename U, typename V=typename CELL::value_t>
  static void computeRhoU(CELL& cell, RHO& rho, U& u) any_platform
  {
    computeRhoJ(cell, rho, u);
    for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
      u[iD] /= rho;
    }
  }

  /// Computation of stress tensor
  template <typename CELL, typename RHO, typename U, typename PI, typename V=typename CELL::value_t>
  static void computeStress(CELL& cell, const RHO& rho, const U& u, PI& pi) any_platform
  {
    int iPi = 0;
    for (int iAlpha=0; iAlpha < DESCRIPTOR::d; ++iAlpha) {
      for (int iBeta=iAlpha; iBeta < DESCRIPTOR::d; ++iBeta) {
        pi[iPi] = V();
        for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
          pi[iPi] += descriptors::c<DESCRIPTOR>(iPop,iAlpha)*
                     descriptors::c<DESCRIPTOR>(iPop,iBeta) * cell[iPop];
        }
        // stripe off equilibrium contribution
        pi[iPi] -= rho*u[iAlpha]*u[iBeta];
        if (iAlpha==iBeta) {
          pi[iPi] -= V{1} / descriptors::invCs2<V,DESCRIPTOR>()*(rho-V{1});
        }
        ++iPi;
      }
    }
  }

  /// Computation of all hydrodynamic variables
  template <typename CELL, typename RHO, typename U, typename PI, typename V=typename CELL::value_t>
  static void computeAllMomenta(CELL& cell, RHO& rho, U& u, PI& pi) any_platform
  {
    computeRhoU(cell, rho, u);
    computeStress(cell, rho, u, pi);
  }

  template <typename CELL, typename FEQ, typename V=typename CELL::value_t>
  static void computeFeq(CELL& cell, FEQ& fEq) any_platform
  {
    V rho {};
    V u[DESCRIPTOR::d] {};
    computeRhoU(cell, rho, u);
    const V uSqr = util::normSqr<V,DESCRIPTOR::d>(u);
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      fEq[iPop] = equilibrium<DESCRIPTOR>::secondOrder(iPop, rho, u, uSqr);
    }
  }

  /// Computation of non-equilibrium distribution
  template <typename CELL, typename FNEQ, typename RHO, typename U, typename V=typename CELL::value_t>
  static void computeFneq(CELL& cell, FNEQ& fNeq, const RHO& rho, const U& u) any_platform
  {
    const V uSqr = util::normSqr<U,DESCRIPTOR::d>(u);
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      fNeq[iPop] = cell[iPop] - equilibrium<DESCRIPTOR>::secondOrder(iPop, rho, u, uSqr);
    }
  }

  template <typename CELL, typename FNEQ, typename V=typename CELL::value_t>
  static void computeFneq(CELL& cell, FNEQ& fNeq) any_platform
  {
    V rho{};
    V u[DESCRIPTOR::d] {};
    computeRhoU(cell, rho, u);
    computeFneq(cell, fNeq, rho, u);
  }

  /// BGK collision step
  template <typename CELL, typename RHO, typename VELOCITY, typename OMEGA, typename V=typename CELL::value_t>
  static V bgkCollision(CELL& cell, const RHO& rho, const VELOCITY& u, const OMEGA& omega) any_platform
  {
    const V uSqr = util::normSqr<VELOCITY,DESCRIPTOR::d>(u);
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      cell[iPop] *= V{1} - omega;
      cell[iPop] += omega * equilibrium<DESCRIPTOR>::secondOrder(iPop, rho, u, uSqr);
    }
    return uSqr;
  }

  /// Advection diffusion BGK collision step
  template <typename CELL, typename RHO, typename VELOCITY, typename OMEGA, typename V=typename CELL::value_t>
  static V adeBgkCollision(CELL& cell, const RHO& rho, const VELOCITY& u, const OMEGA& omega)
  {
    const V uSqr = util::normSqr<VELOCITY,DESCRIPTOR::d>(u);
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      cell[iPop] *= V{1} - omega;
      cell[iPop] += omega * equilibrium<DESCRIPTOR>::firstOrder(iPop, rho, u);
    }
    return uSqr;
  }

  /// Incompressible BGK collision step
  template <typename CELL, typename PRESSURE, typename J, typename OMEGA, typename V=typename CELL::value_t>
  static V incBgkCollision(CELL& cell, const PRESSURE& pressure, const J& j, const OMEGA& omega) any_platform
  {
    const V jSqr = util::normSqr<J,DESCRIPTOR::d>(j);
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      cell[iPop] *= V{1} - omega;
      cell[iPop] += omega * equilibrium<DESCRIPTOR>::template incompressible<J,V,V>(iPop, j, jSqr, pressure);
    }
    return jSqr;
  }

  /// BGK collision step with density correction
  template <typename CELL, typename RHO, typename U, typename RATIORHO, typename OMEGA, typename V=typename CELL::value_t>
  static V constRhoBgkCollision(CELL& cell, const RHO& rho, const U& u,
                                const RATIORHO& ratioRho, const OMEGA& omega) any_platform
  {
    const V uSqr = util::normSqr<U,DESCRIPTOR::d>(u);
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      V feq = equilibrium<DESCRIPTOR>::secondOrder(iPop, rho, u, uSqr);
      cell[iPop] =
        ratioRho*(feq+descriptors::t<V,DESCRIPTOR>(iPop))-descriptors::t<V,DESCRIPTOR>(iPop) +
        (V{1}-omega)*(cell[iPop]-feq);
    }
    return uSqr;
  }

  /// RLB advection diffusion collision step
  template <typename CELL, typename RHO, typename U, typename OMEGA, typename V=typename CELL::value_t>
  static V rlbCollision(CELL& cell, const RHO& rho, const U& u, const OMEGA& omega) any_platform
  {
    const V uSqr = util::normSqr<U,DESCRIPTOR::d>(u);
    // First-order moment for the regularization
    V j1[DESCRIPTOR::d];
    for ( int iD = 0; iD < DESCRIPTOR::d; ++iD ) {
      j1[iD] = V();
    }

    V fEq[DESCRIPTOR::q];
    for ( int iPop = 0; iPop < DESCRIPTOR::q; ++iPop ) {
      fEq[iPop] = equilibrium<DESCRIPTOR>::firstOrder( iPop, rho, u );
      for ( int iD = 0; iD < DESCRIPTOR::d; ++iD ) {
        j1[iD] += descriptors::c<DESCRIPTOR>(iPop,iD) * ( cell[iPop] - fEq[iPop] );
      }
    }

    // Collision step
    for ( int iPop = 0; iPop < DESCRIPTOR::q; ++iPop ) {
      V fNeq = V();
      for ( int iD = 0; iD < DESCRIPTOR::d; ++iD ) {
        fNeq += descriptors::c<DESCRIPTOR>(iPop,iD) * j1[iD];
      }
      fNeq *= descriptors::t<V,DESCRIPTOR>(iPop) * descriptors::invCs2<V,DESCRIPTOR>();
      cell[iPop] = fEq[iPop] + ( V{1} - omega ) * fNeq;
    }
    return uSqr;
  }

  /// Renormalized DESCRIPTOR Boltzmann collision operator, fIn --> fOut
  template <typename CELL, typename RHO, typename U, typename PI, typename OMEGA, typename V=typename CELL::value_t>
  static V rlbCollision(CELL& cell, const RHO& rho, const U& u, const PI& pi, const OMEGA& omega) any_platform
  {
    const V uSqr = util::normSqr<U,DESCRIPTOR::d>(u);
    cell[0] = equilibrium<DESCRIPTOR>::secondOrder(0, rho, u, uSqr)
              + (V{1}-omega) * equilibrium<DESCRIPTOR>::template fromPiToFneq<V>(0, pi);
    for (int iPop=1; iPop <= DESCRIPTOR::q/2; ++iPop) {
      cell[iPop] = equilibrium<DESCRIPTOR>::secondOrder(iPop, rho, u, uSqr);
      cell[iPop+DESCRIPTOR::q/2] = equilibrium<DESCRIPTOR>::secondOrder(iPop+DESCRIPTOR::q/2, rho, u, uSqr);

      V fNeq = (V{1}-omega) * equilibrium<DESCRIPTOR>::template fromPiToFneq<V>(iPop, pi);
      cell[iPop] += fNeq;
      cell[iPop+DESCRIPTOR::q/2] += fNeq;
    }
    return uSqr;
  }

  template <typename CELL, typename NEWRHO, typename NEWU, typename V=typename CELL::value_t>
  static void defineEqFirstOrder(CELL& cell, const NEWRHO& newRho, const NEWU& newU)
  {
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      cell[iPop] = equilibrium<DESCRIPTOR>::firstOrder(iPop, newRho, newU);
    }
  }

  template <typename CELL, typename OLDRHO, typename OLDU, typename NEWRHO, typename NEWU, typename V=typename CELL::value_t>
  static void defineNEq(CELL& cell,
                        const OLDRHO& oldRho, const OLDU& oldU,
                        const NEWRHO& newRho, const NEWU& newU)
  {
    const V oldUSqr = util::normSqr<OLDU,DESCRIPTOR::d>(oldU);
    const V newUSqr = util::normSqr<NEWU,DESCRIPTOR::d>(newU);
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      cell[iPop] += equilibrium<DESCRIPTOR>::secondOrder(iPop, newRho, newU, newUSqr)
                  - equilibrium<DESCRIPTOR>::secondOrder(iPop, oldRho, oldU, oldUSqr);
    }
  }

  template <typename CELL, typename RHO, typename U, typename PI, typename V=typename CELL::value_t>
  static void defineNEqFromPi(CELL& cell,
                              const RHO& rho,
                              const U& u,
                              const PI& pi)
  {
    const V uSqr = util::normSqr<U,DESCRIPTOR::d>(u);
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      cell[iPop] = equilibrium<DESCRIPTOR>::secondOrder(iPop, rho, u, uSqr)
                 + equilibrium<DESCRIPTOR>::template fromPiToFneq<V>(iPop, pi);
    }
  }

  /// Computes squared norm of non-equilibrium part of 2nd momentum for forced dynamics
  template <typename CELL, typename FORCE, typename V=typename CELL::value_t>
  static V computePiNeqNormSqr(CELL& cell, const FORCE& force)
  {
    V rho, u[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR>::n];
    computeAllMomenta(cell, rho, u, pi);
    V ForceTensor[util::TensorVal<DESCRIPTOR>::n];
    // Creation of body force tensor (rho/2.)*(G_alpha*U_beta + U_alpha*G_Beta)
    int iPi = 0;
    for (int Alpha=0; Alpha<DESCRIPTOR::d; ++Alpha) {
      for (int Beta=Alpha; Beta<DESCRIPTOR::d; ++Beta) {
        ForceTensor[iPi] = rho/2.*(force[Alpha]*u[Beta] + u[Alpha]*force[Beta]);
        ++iPi;
      }
    }
    // Creation of second-order moment off-equilibrium tensor
    for (int iPi=0; iPi < util::TensorVal<DESCRIPTOR >::n; ++iPi) {
      pi[iPi] += ForceTensor[iPi];
    }
    V PiNeqNormSqr = pi[0]*pi[0] + 2.*pi[1]*pi[1] + pi[2]*pi[2];
    if constexpr (util::TensorVal<DESCRIPTOR>::n == 6) {
      PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2.*pi[4]*pi[4] +pi[5]*pi[5];
    }
    return PiNeqNormSqr;
  }

  /// Computes squared norm of non-equilibrium part of 2nd momentum for standard (non-forced) dynamics
  template <typename CELL, typename V=typename CELL::value_t>
  static V computePiNeqNormSqr(CELL& cell)
  {
    V rho, u[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR>::n];
    computeAllMomenta(cell, rho, u, pi);
    V PiNeqNormSqr = pi[0]*pi[0] + 2.*pi[1]*pi[1] + pi[2]*pi[2];
    if constexpr (util::TensorVal<DESCRIPTOR >::n == 6) {
      PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2.*pi[4]*pi[4] +pi[5]*pi[5];
    }
    return PiNeqNormSqr;
  }

  /// Add a force term after BGK collision
  template <typename CELL, typename RHO, typename U, typename OMEGA, typename FORCE, typename V=typename CELL::value_t>
  static void addExternalForce(CELL& cell, const RHO& rho, const U& u, const OMEGA& omega, const FORCE& force)
  {
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      V c_u{};
      for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
        c_u += descriptors::c<DESCRIPTOR>(iPop,iD)*u[iD];
      }
      c_u *= descriptors::invCs2<V,DESCRIPTOR>() * descriptors::invCs2<V,DESCRIPTOR>();
      V forceTerm{};
      for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
        forceTerm +=
          (   (V{descriptors::c<DESCRIPTOR>(iPop,iD)} - u[iD]) * descriptors::invCs2<V,DESCRIPTOR>()
              + c_u * descriptors::c<DESCRIPTOR>(iPop,iD)
          )
          * force[iD];
      }
      forceTerm *= descriptors::t<V,DESCRIPTOR>(iPop);
      forceTerm *= V{1} - omega * V{0.5};
      forceTerm *= rho;
      cell[iPop] += forceTerm;
    }
  }

};

}

#endif

#include "lbm.cse.h"
