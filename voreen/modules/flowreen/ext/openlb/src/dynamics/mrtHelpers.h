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
#ifndef MRT_HELPERS_H
#define MRT_HELPERS_H

#include "lbHelpers.h"
#include "mrtLatticeDescriptors.h"

namespace olb {

using namespace descriptors;

/// All helper functions are inside this structure
template<typename T, template<typename U> class Lattice>
struct mrtHelpers {
  /// Computation of equilibrium distribution (in momenta space)
  static T equilibrium( int iPop, T rho, const T u[Lattice<T>::d],
                        const T uSqr )
  {
    T equ = T();
    for (int jPop = 0; jPop < Lattice<T>::q; ++jPop) {
      equ += Lattice<T>::M[iPop][jPop] *
             (lbHelpers<T,Lattice>::equilibrium(jPop,rho,u,uSqr) +
              Lattice<T>::t[jPop]);
    }

    return equ;
  }

  /// Computation of all equilibrium distribution (in momenta space)
  static void computeEquilibrium( T momentaEq[Lattice<T>::q],
                                  T rho, const T u[Lattice<T>::d],
                                  const T uSqr )
  {
    for (int iPop = 0; iPop < Lattice<T>::q; ++iPop) {
      momentaEq[iPop] = T();
      for (int jPop = 0; jPop < Lattice<T>::q; ++jPop) {
        momentaEq[iPop] += Lattice<T>::M[iPop][jPop] *
                           (lbHelpers<T,Lattice>::equilibrium(jPop,rho,u,uSqr) +
                            Lattice<T>::t[jPop]);
      }
    }
  }

  static void computeMomenta(T momenta[Lattice<T>::q], Cell<T,Lattice> &cell)
  {
    for (int iPop = 0; iPop < Lattice<T>::q; ++iPop) {
      momenta[iPop] = T();
      for (int jPop = 0; jPop < Lattice<T>::q; ++jPop) {
        momenta[iPop] += Lattice<T>::M[iPop][jPop] *
                         (cell[jPop] + Lattice<T>::t[jPop]);
      }
    }
  }

  /// MRT collision step
  static T mrtCollision( Cell<T,Lattice>& cell,
                         T rho, const T u[Lattice<T>::d],
                         T invM_S[Lattice<T>::q][Lattice<T>::q])
  {
    T uSqr = util::normSqr<T,Lattice<T>::d>(u);
    T momenta[Lattice<T>::q];
    T momentaEq[Lattice<T>::q];

    computeMomenta(momenta,cell);
    computeEquilibrium(momentaEq,rho,u,uSqr);

    for (int iPop=0; iPop < Lattice<T>::q; ++iPop) {
      T collisionTerm = T();
      for (int jPop = 0; jPop < Lattice<T>::q; ++jPop) {
        collisionTerm += invM_S[iPop][jPop] *
                         (momenta[jPop] - momentaEq[jPop]);
      }
      cell[iPop] -= collisionTerm;
    }

    return uSqr;
  }

  /// MRT SGS collision step
  static T mrtSGSCollision( Cell<T,Lattice>& cell,
                            T rho, const T u[Lattice<T>::d],
                            T omega,
                            T invM_S_SGS[Lattice<T>::q][Lattice<T>::q])
  {
    T uSqr = util::normSqr<T,Lattice<T>::d>(u);
    T momenta[Lattice<T>::q];
    T momentaEq[Lattice<T>::q];


    computeMomenta(momenta,cell);
    computeEquilibrium(momentaEq,rho,u,uSqr);

    for (int iPop=0; iPop < Lattice<T>::q; ++iPop) {
      T collisionTerm = T();
      for (int jPop = 0; jPop < Lattice<T>::q; ++jPop) {

        // cout << "wert_in helpers"<<iPop <<jPop << "= "<<  invM_S_SGS[iPop][jPop]<< endl;
        collisionTerm += invM_S_SGS[iPop][jPop] *
                         (momenta[jPop] - momentaEq[jPop]);



        //      if (iPop==jPop && (iPop==Lattice<T>::shearViscIndexes[iPop]))
        //      {
        //        collisionTerm += invM_S_SGS[iPop][jPop] /(invM_S_SGS[iPop][jPop])*omega*(Lattice<T>::invM[iPop][jPop])*
        //                                 (momenta[jPop] - momentaEq[jPop]);
        //        cout << "omege: " << omega<< endl;
        //      }
        //      else
        //      {
        //        collisionTerm += invM_S_SGS[iPop][jPop] *
        //                                 (momenta[jPop] - momentaEq[jPop]);
        //      }

      }
      cell[iPop] -= collisionTerm;
    }

    return uSqr;
  }


  /// Ladd-Verberg-I body force model for MRT
  /// A.Ladd, R. Verberg, Lattice-Boltzmann simulations of particle-fluid suspensions, Journal of Statistical Physics 104(2001)
  static void addExternalForce( Cell<T,Lattice>& cell,
                                T rho,
                                const T u[Lattice<T>::d],
                                T invM_S[Lattice<T>::q][Lattice<T>::q])
  {
    static const int forceBeginsAt = Lattice<T>::ExternalField::forceBeginsAt;
    T* force = cell.getExternal(forceBeginsAt);
    T f_u = T();
    for (int iD=0; iD < Lattice<T>::d; ++iD) {
      f_u += force[iD]*u[iD];
    }

    for (int iPop=0; iPop < Lattice<T>::q; ++iPop) {
      T c_u = T();
      T c_f = T();
      for (int iD=0; iD < Lattice<T>::d; ++iD) {
        c_u += Lattice<T>::c[iPop][iD]*u[iD];
        c_f += Lattice<T>::c[iPop][iD]*force[iD];
      }
      T f1 = Lattice<T>::t[iPop]*rho*c_f*Lattice<T>::invCs2;
      T f2 = Lattice<T>::t[iPop]*rho*(c_u*c_f*Lattice<T>::invCs2-f_u)*Lattice<T>::invCs2;

      T invMsM = T();
      for (int jPop=0; jPop < Lattice<T>::q; ++jPop) {
        invMsM += invM_S[iPop][jPop]*Lattice<T>::M[jPop][iPop];
      }
      cell[iPop] += f1 + f2 - invMsM*f2/2.;
    }

    return;
  }

};  // struct mrtHelpers

}  // namespace olb

// The specialized code is directly included. That is because we never want
// it to be precompiled so that in both the precompiled and the
// "include-everything" version, the compiler can apply all the
// optimizations it wants.
#include "mrtHelpers2D.h"
#include "mrtHelpers3D.h"

#endif
