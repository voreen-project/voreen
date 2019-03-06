/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Jonas Latt
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
 * Specialized helper functions for advanced techniques around LB
 * implementations. They implement the physics of the first-order terms
 * of the Chapman-Enskog expansion and are useful whenever a transition
 * from hydrodynamical variables (rho, u) to kinetic variables (f) si to
 * be implemented. Additionally, they are used for the implementation of
 * the stable RLB dynamics.
 *
 * This file is all about efficiency. The generic
 * template code is specialized for commonly used Lattices, so that a
 * maximum performance can be taken out of each case.
 */
#ifndef FIRST_ORDER_LB_HELPERS_H
#define FIRST_ORDER_LB_HELPERS_H

#include "latticeDescriptors.h"
#include "core/cell.h"
#include "core/util.h"
#include "lbHelpers.h"

namespace olb {

/// General first-order functions
template<typename T, template<typename U> class Lattice>
struct firstOrderLbHelpers {

  /// Compute off-equilibrium part of the f's from the stress tensor Pi.
  /** Implements the following formula (with Einstein index contraction):
   * \f[ f_i^{neq} = t_i / (2 c_s^4) *
   *                 (c_{ia} c_{ib} - c_s^2 \delta_{ab}) \Pi_{ab} \f]
   * By Pi we mean the tensor computed from the off-equilibrium functions:
   * \f[ \Pi = \sum c_i c_i f_i^{neq}
   *         = \sum c_i c_i f_i - \rho u u - c_s^2 \rho\ Id \f]
   */
  static T fromPiToFneq (
    int iPop, const T pi[util::TensorVal<Lattice<T> >::n] )
  {
    typedef Lattice<T> L;
    T fNeq = T();
    int iPi = 0;
    // Iterate only over superior triangle + diagonal, and add
    // the elements under the diagonal by symmetry
    for (int iAlpha=0; iAlpha<L::d; ++iAlpha) {
      for (int iBeta=iAlpha; iBeta<L::d; ++iBeta) {
        T toAdd = L::c[iPop][iAlpha]*L::c[iPop][iBeta];
        if (iAlpha==iBeta) {
          toAdd -= 1./L::invCs2;
        } else {
          toAdd *= (T)2; // multiply off-diagonal elements by 2
        }                  // because the Q tensor is symmetric
        toAdd *= pi[iPi++];
        fNeq += toAdd;
      }
    }
    fNeq *= L::t[iPop] * L::invCs2 * L::invCs2 / (T)2;
    return fNeq;
  }

  /// Compute off-equilibrium part of the f's from the strain rate tensor S.
  /** Implements the following formula:
   * \f[ f_i^{neq} = - t_i / (c_s^2\omega) *
   *                 (c_{ia} c_{ib} - c_s^2 \delta_{ab}) S_{ab} \f]
   * By S we mean the tensor computed from the velocity gradients:
   * \f[ S_{\alpha\beta} = 1/2 (
   *     \partial_\alpha(\rho u_\beta) + \partial_\beta(\rho u_\alpha) ) \f]
   */
  static T fromStrainToFneq (
    int iPop, const T S[util::TensorVal<Lattice<T> >::n], T omega )
  {
    typedef Lattice<T> L;
    T fNeq = fromPiToFneq(iPop,S) * (-(T)2 / L::invCs2 / omega);
    return fNeq;
  }

};  // struct firstOrderLbHelpers

/// Specific helper functions for RLB dynamics
template<typename T, template<typename U> class Lattice>
struct rlbHelpers {
  /// Renormalized Lattice Boltzmann collision operator, fIn --> fOut
  static T rlbCollision (
    Cell<T,Lattice>& cell, T rho, const T u[Lattice<T>::d],
    const T pi[util::TensorVal<Lattice<T> >::n], T omega )
  {
    typedef Lattice<T> L;
    const T uSqr = util::normSqr<T,L::d>(u);
    cell[0] = lbHelpers<T,Lattice>::equilibrium(0, rho, u, uSqr)
              + ((T)1-omega) *
              firstOrderLbHelpers<T,Lattice>::fromPiToFneq(0, pi);
    for (int iPop=1; iPop<=L::q/2; ++iPop) {
      cell[iPop] =
        lbHelpers<T,Lattice>::equilibrium(iPop, rho, u, uSqr);
      cell[iPop+L::q/2] =
        lbHelpers<T,Lattice>::equilibrium(iPop+L::q/2, rho, u, uSqr);

      T fNeq = ((T)1-omega) *
               firstOrderLbHelpers<T,Lattice>::fromPiToFneq(iPop, pi);
      cell[iPop] += fNeq;
      cell[iPop+L::q/2] += fNeq;
    }
    return uSqr;
  }

};  // struct rlbHelpers

}  // namespace olb

// The specialized code is directly included. That is because we never want
// it to be precompiled so that in both the precompiled and the
// "include-everything" version, the compiler can apply all the
// optimizations it wants.
#include "firstOrderLbHelpers2D.h"
#include "firstOrderLbHelpers3D.h"

#endif
