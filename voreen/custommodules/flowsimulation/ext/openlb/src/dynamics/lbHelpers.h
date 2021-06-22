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
 * Helper functions for the implementation of LB dynamics. This file is all
 * about efficiency. The generic template code is specialized for commonly
 * used Lattices, so that a maximum performance can be taken out of each
 * case.
 */
#ifndef LB_HELPERS_H
#define LB_HELPERS_H

#include "latticeDescriptors.h"
#include "core/cell.h"
#include "core/util.h"


namespace olb {


// Forward declarations
template<typename T, typename DESCRIPTOR> struct lbDynamicsHelpers;
template<typename T, typename DESCRIPTOR> struct lbExternalHelpers;
template<typename T, typename DESCRIPTOR> struct lbLatticeHelpers;

/// This structure forwards the calls to the appropriate helper class
template<typename T, typename DESCRIPTOR>
struct lbHelpers {

  static T equilibrium(int iPop, T rho, const T u[DESCRIPTOR::d], const T uSqr)
  {
    return lbDynamicsHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho, u, uSqr);
  }

  static T equilibriumFirstOrder(int iPop, T rho, const T u[DESCRIPTOR::d])
  {
    return lbDynamicsHelpers<T,DESCRIPTOR>
           ::equilibriumFirstOrder(iPop, rho, u);
  }

  static T incEquilibrium(int iPop, const T j[DESCRIPTOR::d], const T jSqr, const T pressure)
  {
    return lbDynamicsHelpers<T,DESCRIPTOR>
           ::incEquilibrium(iPop, j, jSqr, pressure);
  }

  static T equilibriumP1(int iPop, T rho, std::array<T,DESCRIPTOR::d>& u)
  {
    return lbDynamicsHelpers<T,DESCRIPTOR>::equilibriumP1(iPop, rho, u);
  }

  static void computeFneq ( ConstCell<T,DESCRIPTOR>& cell,
                            T fNeq[DESCRIPTOR::q], T rho, const T u[DESCRIPTOR::d] )
  {
    lbDynamicsHelpers<T,DESCRIPTOR>::computeFneq(cell, fNeq, rho, u);
  }

  static T bgkCollision(Cell<T,DESCRIPTOR>& cell, T const& rho, const T u[DESCRIPTOR::d], T const& omega)
  {
    return lbDynamicsHelpers<T,DESCRIPTOR>
           ::bgkCollision(cell, rho, u, omega);
  }

  static T bgkCollision(Cell<T,DESCRIPTOR>& cell, T const& rho, const T u[DESCRIPTOR::d], const T omega[DESCRIPTOR::q])
  {
    const T uSqr = util::normSqr<T,DESCRIPTOR::d>(u);
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      cell[iPop] *= (T)1-omega[iPop];
      cell[iPop] += omega[iPop] * lbDynamicsHelpers<T,DESCRIPTOR>::equilibrium (iPop, rho, u, uSqr );
    }
    return uSqr;
  }

  static T rlbCollision( Cell<T,DESCRIPTOR>& cell, T rho, const T u[DESCRIPTOR::d], T omega )
  {
    return lbDynamicsHelpers<T, DESCRIPTOR>
           ::rlbCollision( cell, rho, u, omega );
  }

  static T mrtCollision( Cell<T,DESCRIPTOR>& cell, T const& rho, const T u[DESCRIPTOR::d], T invM_S[DESCRIPTOR::q][DESCRIPTOR::q] )
  {

    return lbDynamicsHelpers<T, DESCRIPTOR>
           ::mrtCollision(cell, rho, u, invM_S );
  }

  static T incBgkCollision(Cell<T,DESCRIPTOR>& cell, T pressure, const T j[DESCRIPTOR::d], T omega)
  {
    return lbDynamicsHelpers<T,DESCRIPTOR>
           ::incBgkCollision(cell, pressure, j, omega);
  }

  static T constRhoBgkCollision(Cell<T,DESCRIPTOR>& cell,
                                T rho, const T u[DESCRIPTOR::d], T ratioRho, T omega)
  {
    return lbDynamicsHelpers<T,DESCRIPTOR>
           ::constRhoBgkCollision(cell, rho, u, ratioRho, omega);
  }

  static T computeRho(ConstCell<T,DESCRIPTOR>& cell)
  {
    return lbDynamicsHelpers<T,DESCRIPTOR>
           ::computeRho(cell);
  }

  static void computeJ(ConstCell<T,DESCRIPTOR>& cell, T j[DESCRIPTOR::d] )
  {
    lbDynamicsHelpers<T,DESCRIPTOR>
    ::computeJ(cell, j);
  }

  static void computeRhoU(ConstCell<T,DESCRIPTOR>& cell, T& rho, T u[DESCRIPTOR::d])
  {
    lbDynamicsHelpers<T,DESCRIPTOR>
    ::computeRhoU(cell, rho, u);
  }

  static void computeFeq(ConstCell<T,DESCRIPTOR>& cell, T fEq[DESCRIPTOR::q])
  {
    T rho{};
    T u[2] {};
    computeRhoU(cell, rho, u);
    const T uSqr = u[0]*u[0] + u[1]*u[1];
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      fEq[iPop] = equilibrium(iPop, rho, u, uSqr);
    }
  }

  static void computeFneq(ConstCell<T,DESCRIPTOR>& cell, T fNeq[DESCRIPTOR::q])
  {
    T rho{};
    T u[2] {};
    computeRhoU(cell, rho, u);
    computeFneq(cell, fNeq, rho, u);
  }

  static void computeRhoJ(ConstCell<T,DESCRIPTOR>& cell, T& rho, T j[DESCRIPTOR::d])
  {
    lbDynamicsHelpers<T,DESCRIPTOR>::computeRhoJ(cell, rho, j);
  }

  static void computeStress(ConstCell<T,DESCRIPTOR>& cell, T rho, const T u[DESCRIPTOR::d],
                            T pi[util::TensorVal<DESCRIPTOR >::n] )
  {
    lbDynamicsHelpers<T,DESCRIPTOR>::computeStress(cell, rho, u, pi);
  }

  static void computeAllMomenta(ConstCell<T,DESCRIPTOR>& cell, T& rho, T u[DESCRIPTOR::d],
                                T pi[util::TensorVal<DESCRIPTOR >::n] )
  {
    lbDynamicsHelpers<T,DESCRIPTOR>::computeAllMomenta(cell, rho, u, pi);
  }

  /// Computes squared norm of non-equilibrium part of 2nd momentum for standard (non-forced) dynamics
  template <typename X = DESCRIPTOR>
  static std::enable_if_t<!X::template provides<descriptors::FORCE>(), T>
  computePiNeqNormSqr(ConstCell<T,DESCRIPTOR>& cell)
  {
    //return lbDynamicsHelpers<T,DESCRIPTOR>
    //::computePiNeqNormSqr(cell);
    T rho, u[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR >::n];
    computeAllMomenta(cell, rho, u, pi);
    T PiNeqNormSqr = pi[0]*pi[0] + 2.*pi[1]*pi[1] + pi[2]*pi[2];
    if (util::TensorVal<DESCRIPTOR >::n == 6) {
      PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2.*pi[4]*pi[4] +pi[5]*pi[5];
    }
    return PiNeqNormSqr;
  }

  /// Computes squared norm of non-equilibrium part of 2nd momentum for forced dynamics
  template <typename X = DESCRIPTOR>
  static std::enable_if_t<X::template provides<descriptors::FORCE>(), T>
  computePiNeqNormSqr(ConstCell<T,DESCRIPTOR>& cell)
  {
    //return lbDynamicsHelpers<T,DESCRIPTOR>
    //::computeForcedPiNeqNormSqr(cell);
    T rho, u[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR >::n];
    computeAllMomenta(cell, rho, u, pi);
    //Creation of body force tensor (rho/2.)*(G_alpha*U_beta + U_alpha*G_Beta)
    T ForceTensor[util::TensorVal<DESCRIPTOR >::n];
    int iPi = 0;
    auto force = cell.template getField<descriptors::FORCE>();
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
    T PiNeqNormSqr = pi[0]*pi[0] + 2.*pi[1]*pi[1] + pi[2]*pi[2];
    if (util::TensorVal<DESCRIPTOR >::n == 6) {
      PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2.*pi[4]*pi[4] +pi[5]*pi[5];
    }
    return PiNeqNormSqr;
  }

  static void modifyVelocity(Cell<T,DESCRIPTOR>& cell, const T newU[DESCRIPTOR::d])
  {
    lbDynamicsHelpers<T,DESCRIPTOR>
    ::modifyVelocity(cell, newU);
  }

  static void addExternalForce(Cell<T,DESCRIPTOR>& cell, const T u[DESCRIPTOR::d], T omega, T amplitude=(T)1)
  {
    lbExternalHelpers<T,DESCRIPTOR>::addExternalForce(cell, u, omega, amplitude);
  }

  static void swapAndStream2D(BlockLattice2D<T,DESCRIPTOR>& lattice, int iX, int iY)
  {
    lbLatticeHelpers<T,DESCRIPTOR>::swapAndStream2D(lattice, iX, iY);
  }

  static void swapAndStream3D(BlockLattice3D<T,DESCRIPTOR>& lattice, int iX, int iY, int iZ)
  {
    lbLatticeHelpers<T,DESCRIPTOR>::swapAndStream3D(lattice, iX, iY, iZ);
  }

};  // struct lbHelpers


/// All helper functions are inside this structure
template<typename T, typename DESCRIPTOR>
struct lbDynamicsHelpers {
  /// Computation of equilibrium distribution, second order in u
  static T equilibrium(int iPop, T rho, const T u[DESCRIPTOR::d], const T uSqr)
  {
    T c_u = T();
    for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
      c_u += descriptors::c<DESCRIPTOR>(iPop,iD)*u[iD];
    }
    return rho
           * descriptors::t<T,DESCRIPTOR>(iPop)
           * ( T{1}
               + descriptors::invCs2<T,DESCRIPTOR>() * c_u
               + descriptors::invCs2<T,DESCRIPTOR>() * descriptors::invCs2<T,DESCRIPTOR>() * T{0.5} * c_u *c_u
               - descriptors::invCs2<T,DESCRIPTOR>() * T{0.5} * uSqr )
           - descriptors::t<T,DESCRIPTOR>(iPop);
  }

  /// Computation of equilibrium distribution, first order in u
  static T equilibriumFirstOrder(int iPop, T rho, const T u[DESCRIPTOR::d])
  {
    T c_u = T();
    for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
      c_u += descriptors::c<DESCRIPTOR>(iPop,iD)*u[iD];
    }
    return rho
           * descriptors::t<T,DESCRIPTOR>(iPop)
           * ( T{1} + c_u * descriptors::invCs2<T,DESCRIPTOR>() )
           - descriptors::t<T,DESCRIPTOR>(iPop);
  }

  // compute equilibrium f^eq_i eq. (5.32) from DOI:10.1002/9780470177013
  static T equilibriumP1(int iPop, T rho, std::array<T,DESCRIPTOR::d>& u)
  {
    T c_u = T();
    // compute scalar product of c[iPop]*u
    for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
      c_u += descriptors::c<DESCRIPTOR>(iPop,iD)*u[iD];
    }
    return descriptors::t<T,DESCRIPTOR>(iPop) * (rho + c_u)
           - descriptors::t<T,DESCRIPTOR>(iPop);
  }

  static T incEquilibrium( int iPop, const T j[DESCRIPTOR::d],
                           const T jSqr, const T pressure )
  {
    T c_j = T();
    for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
      c_j += descriptors::c<DESCRIPTOR>(iPop,iD)*j[iD];
    }

    return descriptors::t<T,DESCRIPTOR>(iPop)
           * ( descriptors::invCs2<T,DESCRIPTOR>() * pressure
               + descriptors::invCs2<T,DESCRIPTOR>() * c_j
               + descriptors::invCs2<T,DESCRIPTOR>() * descriptors::invCs2<T,DESCRIPTOR>()/T{2} * c_j*c_j
               - descriptors::invCs2<T,DESCRIPTOR>()/T{2} * jSqr )
           - descriptors::t<T,DESCRIPTOR>(iPop);
  }

  /// Computation of non-equilibrium distribution
  static void computeFneq(ConstCell<T,DESCRIPTOR>& cell, T fNeq[DESCRIPTOR::q], T rho, const T u[DESCRIPTOR::d])
  {
    const T uSqr = util::normSqr<T,DESCRIPTOR::d>(u);
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      fNeq[iPop] = cell[iPop] - equilibrium(iPop, rho, u, uSqr);
    }
  }

  /// BGK collision step
  static T bgkCollision(Cell<T,DESCRIPTOR>& cell, T const& rho, const T u[DESCRIPTOR::d], T const& omega)
  {
    const T uSqr = util::normSqr<T,DESCRIPTOR::d>(u);
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      cell[iPop] *= (T)1-omega;
      cell[iPop] += omega * lbDynamicsHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho, u, uSqr );
    }
    return uSqr;
  }

  /// Incompressible BGK collision step
  static T incBgkCollision(Cell<T,DESCRIPTOR>& cell, T pressure, const T j[DESCRIPTOR::d], T omega)
  {
    const T jSqr = util::normSqr<T,DESCRIPTOR::d>(j);
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      cell[iPop] *= (T)1-omega;
      cell[iPop] += omega * lbDynamicsHelpers<T,DESCRIPTOR>::incEquilibrium(iPop, j, jSqr, pressure);
    }
    return jSqr;
  }

  /// BGK collision step with density correction
  static T constRhoBgkCollision(Cell<T,DESCRIPTOR>& cell, T rho, const T u[DESCRIPTOR::d], T ratioRho, T omega)
  {
    const T uSqr = util::normSqr<T,DESCRIPTOR::d>(u);
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      T feq = lbDynamicsHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho, u, uSqr );
      cell[iPop] =
        ratioRho*(feq+descriptors::t<T,DESCRIPTOR>(iPop))-descriptors::t<T,DESCRIPTOR>(iPop) +
        ((T)1-omega)*(cell[iPop]-feq);
    }
    return uSqr;
  }

  /// RLB advection diffusion collision step
  static T rlbCollision(Cell<T,DESCRIPTOR>& cell, T rho, const T u[DESCRIPTOR::d], T omega )
  {
    const T uSqr = util::normSqr<T, DESCRIPTOR::d>( u );
    // First-order moment for the regularization
    T j1[DESCRIPTOR::d];
    for ( int iD = 0; iD < DESCRIPTOR::d; ++iD ) {
      j1[iD] = T();
    }

    T fEq[DESCRIPTOR::q];
    for ( int iPop = 0; iPop < DESCRIPTOR::q; ++iPop ) {
      fEq[iPop] = lbDynamicsHelpers<T, DESCRIPTOR>::equilibriumFirstOrder( iPop, rho, u );
      for ( int iD = 0; iD < DESCRIPTOR::d; ++iD ) {
        j1[iD] += descriptors::c<DESCRIPTOR>(iPop,iD) * ( cell[iPop] - fEq[iPop] );
      }
    }

    // Collision step
    for ( int iPop = 0; iPop < DESCRIPTOR::q; ++iPop ) {
      T fNeq = T();
      for ( int iD = 0; iD < DESCRIPTOR::d; ++iD ) {
        fNeq += descriptors::c<DESCRIPTOR>(iPop,iD) * j1[iD];
      }
      fNeq *= descriptors::t<T,DESCRIPTOR>(iPop) * descriptors::invCs2<T,DESCRIPTOR>();
      cell[iPop] = fEq[iPop] + ( (T)1 - omega ) * fNeq;
    }
    return uSqr;
  }


///// Computation of all equilibrium distribution (in momenta space)
  static void computeMomentaEquilibrium( T momentaEq[DESCRIPTOR::q], T rho, const T u[DESCRIPTOR::d], T uSqr )
  {
    for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
      momentaEq[iPop] = T();
      for (int jPop = 0; jPop < DESCRIPTOR::q; ++jPop) {
        momentaEq[iPop] += DESCRIPTOR::M[iPop][jPop] *
                           (lbDynamicsHelpers<T, DESCRIPTOR>::equilibrium(jPop,rho,u,uSqr) + DESCRIPTOR::t[jPop]);
      }
    }
  }

  static void computeMomenta(T momenta[DESCRIPTOR::q], Cell<T,DESCRIPTOR>& cell)
  {
    for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
      momenta[iPop] = T();
      for (int jPop = 0; jPop < DESCRIPTOR::q; ++jPop) {
        momenta[iPop] += DESCRIPTOR::M[iPop][jPop] *
                         (cell[jPop] + DESCRIPTOR::t[jPop]);
      }
    }
  }

  static T mrtCollision( Cell<T,DESCRIPTOR>& cell, T const& rho, const T u[DESCRIPTOR::d], T invM_S[DESCRIPTOR::q][DESCRIPTOR::q] )
  {
    //// Implemented in advectionDiffusionMRTlbHelpers2D.h and advectionDiffusionMRTlbHelpers3D.h
    T uSqr = util::normSqr<T, DESCRIPTOR::d>(u);
    T momenta[DESCRIPTOR::q];
    T momentaEq[DESCRIPTOR::q];

    computeMomenta(momenta, cell);
    computeMomentaEquilibrium(momentaEq, rho, u, uSqr);

//    std::cout << "momenta = ";
//    for (int i=0; i < DESCRIPTOR::q; ++i) {
//        std::cout << momenta[i] << ", ";
//    }
//    std::cout << std::endl;

//    std::cout << "momentaEq = ";
//    for (int i=0; i < DESCRIPTOR::q; ++i) {
//        std::cout << momentaEq[i] << ", ";
//    }
//    std::cout << std::endl;

    for (int iPop = 0; iPop < DESCRIPTOR::q; ++iPop) {
      T collisionTerm = T();
      for (int jPop = 0; jPop < DESCRIPTOR::q; ++jPop) {
        collisionTerm += invM_S[iPop][jPop] * (momenta[jPop] - momentaEq[jPop]);
      }
      cell[iPop] -= collisionTerm;
    }
    return uSqr;
  }

  /// Computation of density
  static T computeRho(ConstCell<T,DESCRIPTOR>& cell)
  {
    T rho = T();
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      rho += cell[iPop];
    }
    rho += (T)1;
    return rho;
  }

  /// Computation of momentum
  static void computeJ(ConstCell<T,DESCRIPTOR>& cell, T j[DESCRIPTOR::d])
  {
    for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
      j[iD] = T();
    }
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
        j[iD] += cell[iPop]*descriptors::c<DESCRIPTOR>(iPop,iD);
      }
    }
  }

  /// Computation of hydrodynamic variables
  static void computeRhoU(ConstCell<T,DESCRIPTOR>& cell, T& rho, T u[DESCRIPTOR::d])
  {
    rho = T();
    for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
      u[iD] = T();
    }
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      rho += cell[iPop];
      for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
        u[iD] += cell[iPop]*descriptors::c<DESCRIPTOR>(iPop,iD);
      }
    }
    rho += (T)1;
    for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
      u[iD] /= rho;
    }
  }

  /// Computation of hydrodynamic variables
  static void computeRhoJ(ConstCell<T,DESCRIPTOR>& cell, T& rho, T j[DESCRIPTOR::d])
  {
    rho = T();
    for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
      j[iD] = T();
    }
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      rho += cell[iPop];
      for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
        j[iD] += cell[iPop]*descriptors::c<DESCRIPTOR>(iPop,iD);
      }
    }
    rho += (T)1;
  }

  /// Computation of stress tensor
  static void computeStress(ConstCell<T,DESCRIPTOR>& cell,
                            T rho, const T u[DESCRIPTOR::d],
                            T pi[util::TensorVal<DESCRIPTOR>::n] )
  {
    int iPi = 0;
    for (int iAlpha=0; iAlpha < DESCRIPTOR::d; ++iAlpha) {
      for (int iBeta=iAlpha; iBeta < DESCRIPTOR::d; ++iBeta) {
        pi[iPi] = T();
        for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
          pi[iPi] += descriptors::c<DESCRIPTOR>(iPop,iAlpha)*
                     descriptors::c<DESCRIPTOR>(iPop,iBeta) * cell[iPop];
        }
        // stripe off equilibrium contribution
        pi[iPi] -= rho*u[iAlpha]*u[iBeta];
        if (iAlpha==iBeta) {
          pi[iPi] -= 1./descriptors::invCs2<T,DESCRIPTOR>()*(rho-(T)1);
        }
        ++iPi;
      }
    }
  }

  /// Computation of all hydrodynamic variables
  static void computeAllMomenta(ConstCell<T,DESCRIPTOR>& cell, T& rho, T u[DESCRIPTOR::d],
                                T pi[util::TensorVal<DESCRIPTOR>::n] )
  {
    computeRhoU(cell, rho, u);
    computeStress(cell, rho, u, pi);
  }

  /*
  /// Computes squared norm of non-equilibrium part of 2nd momentum
  static T computePiNeqNormSqr(ConstCell<T,DESCRIPTOR>& cell)
  {
    T rho, u[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR >::n];
    computeAllMomenta(cell, rho, u, pi);
    T PiNeqNormSqr = pi[0]*pi[0] + 2.*pi[1]*pi[1] + pi[2]*pi[2];
    if (util::TensorVal<DESCRIPTOR >::n == 6) {
      PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2.*pi[4]*pi[4] +pi[5]*pi[5];
    }
    return PiNeqNormSqr;
  }

  /// Computes squared norm of forced non-equilibrium part of 2nd momentum
  static T computeForcedPiNeqNormSqr(ConstCell<T,DESCRIPTOR>& cell)
  {
    T rho, u[DESCRIPTOR::d], pi[util::TensorVal<DESCRIPTOR >::n];
    computeAllMomenta(cell, rho, u, pi);
    //Creation of body force tensor (rho/2.)*(G_alpha*U_beta + U_alpha*G_Beta)
    T ForceTensor[util::TensorVal<DESCRIPTOR >::n];
    int iPi = 0;
    for (int Alpha=0; Alpha<DESCRIPTOR::d; ++Alpha) {
      for (int Beta=Alpha; Beta<DESCRIPTOR::d; ++Beta) {
        ForceTensor[iPi] = rho/2.*(cell[DESCRIPTOR::template index<descriptors::FORCE>()][Alpha]*u[Beta] + u[Alpha]*cell[DESCRIPTOR::template index<descriptors::FORCE>()][Beta]);
        ++iPi;
      }
    }
    // Creation of second-order moment off-equilibrium tensor
    for (int iPi=0; iPi < util::TensorVal<DESCRIPTOR >::n; ++iPi){
      pi[iPi] += ForceTensor[iPi];
    }
    T PiNeqNormSqr = pi[0]*pi[0] + 2.*pi[1]*pi[1] + pi[2]*pi[2];
    if (util::TensorVal<DESCRIPTOR >::n == 6) {
      PiNeqNormSqr += pi[2]*pi[2] + pi[3]*pi[3] + 2.*pi[4]*pi[4] +pi[5]*pi[5];
    }
    return PiNeqNormSqr;
  }
  */

  static void modifyVelocity(Cell<T,DESCRIPTOR>& cell, const T newU[DESCRIPTOR::d])
  {
    T rho, oldU[DESCRIPTOR::d];
    computeRhoU(cell, rho, oldU);
    const T oldUSqr = util::normSqr<T,DESCRIPTOR::d>(oldU);
    const T newUSqr = util::normSqr<T,DESCRIPTOR::d>(newU);
    for (int iPop=0; iPop<DESCRIPTOR::q; ++iPop) {
      cell[iPop] = cell[iPop]
                   - equilibrium(iPop, rho, oldU, oldUSqr)
                   + equilibrium(iPop, rho, newU, newUSqr);
    }
  }

};  // struct lbDynamicsHelpers

/// Helper functions for dynamics that access external field
template<typename T, typename DESCRIPTOR>
struct lbExternalHelpers {
  /// Add a force term after BGK collision
  static void addExternalForce(Cell<T,DESCRIPTOR>& cell, const T u[DESCRIPTOR::d], T omega, T amplitude)
  {
    const auto force = cell.template getField<descriptors::FORCE>();
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      T c_u = T();
      for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
        c_u += descriptors::c<DESCRIPTOR>(iPop,iD)*u[iD];
      }
      c_u *= descriptors::invCs2<T,DESCRIPTOR>()*descriptors::invCs2<T,DESCRIPTOR>();
      T forceTerm = T();
      for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
        forceTerm +=
          (   ((T)descriptors::c<DESCRIPTOR>(iPop,iD)-u[iD]) * descriptors::invCs2<T,DESCRIPTOR>()
              + c_u * descriptors::c<DESCRIPTOR>(iPop,iD)
          )
          * force[iD];
      }
      forceTerm *= descriptors::t<T,DESCRIPTOR>(iPop);
      forceTerm *= T(1) - omega/T(2);
      forceTerm *= amplitude;
      cell[iPop] += forceTerm;
    }
  }
};  // struct externalFieldHelpers

/// Helper functions with full-lattice access
template<typename T, typename DESCRIPTOR>
struct lbLatticeHelpers {
  /// Swap ("bounce-back") values of a cell (2D), and apply streaming step
  static void swapAndStream2D(BlockLattice2D<T,DESCRIPTOR>& lattice, int iX, int iY)
  {
    const int half = DESCRIPTOR::q/2;
    for (int iPop=1; iPop<=half; ++iPop) {
      const int jX = iX + descriptors::c<DESCRIPTOR>(iPop,0);
      const int jY = iY + descriptors::c<DESCRIPTOR>(iPop,1);
      auto iCell = lattice.get(iX,iY);
      auto jCell = lattice.get(jX,jY);
      T fTmp            = iCell[iPop];
      iCell[iPop]       = iCell[iPop+half];
      iCell[iPop+half]  = jCell[iPop];
      jCell[iPop]       = fTmp;
    }
  }

  /// Swap ("bounce-back") values of a cell (3D), and apply streaming step
  static void swapAndStream3D(BlockLattice3D<T,DESCRIPTOR>& lattice,
                              int iX, int iY, int iZ)
  {
    const int half = DESCRIPTOR::q/2;
    for (int iPop=1; iPop<=half; ++iPop) {
      const int jX = iX + descriptors::c<DESCRIPTOR>(iPop,0);
      const int jY = iY + descriptors::c<DESCRIPTOR>(iPop,1);
      const int jZ = iZ + descriptors::c<DESCRIPTOR>(iPop,2);
      auto iCell = lattice.get(iX,iY,iZ);
      auto jCell = lattice.get(jX,jY,jZ);
      T fTmp           = iCell[iPop];
      iCell[iPop]      = iCell[iPop+half];
      iCell[iPop+half] = jCell[iPop];
      jCell[iPop]      = fTmp;
    }
  }
};

/// All boundary helper functions are inside this structure
template<typename T, typename DESCRIPTOR, int direction, int orientation>
struct BoundaryHelpers {
  static void computeStress (
    ConstCell<T,DESCRIPTOR>& cell, T rho, const T u[DESCRIPTOR::d],
    T pi[util::TensorVal<DESCRIPTOR >::n] )
  {
    typedef DESCRIPTOR L;
    const T uSqr = util::normSqr<T,L::d>(u);

    std::vector<int> const& onWallIndices = util::subIndex<L, direction, 0>();
    std::vector<int> const& normalIndices = util::subIndex<L, direction, orientation>();

    T fNeq[DESCRIPTOR::q];
    for (unsigned fIndex=0; fIndex<onWallIndices.size(); ++fIndex) {
      int iPop = onWallIndices[fIndex];
      if (iPop == 0) {
        fNeq[0] = T();  // fNeq[0] will not be used anyway
      }
      else {
        fNeq[iPop] =
          cell[iPop] -
          lbHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho, u, uSqr);
      }
    }
    for (unsigned fIndex=0; fIndex<normalIndices.size(); ++fIndex) {
      int iPop = normalIndices[fIndex];
      fNeq[iPop] =
        cell[iPop] -
        lbHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho, u, uSqr);
    }

    int iPi = 0;
    for (int iAlpha=0; iAlpha<L::d; ++iAlpha) {
      for (int iBeta=iAlpha; iBeta<L::d; ++iBeta) {
        pi[iPi] = T();
        for (unsigned fIndex=0; fIndex<onWallIndices.size(); ++fIndex) {
          const int iPop = onWallIndices[fIndex];
          pi[iPi] +=
            descriptors::c<L>(iPop,iAlpha)*descriptors::c<L>(iPop,iBeta)*fNeq[iPop];
        }
        for (unsigned fIndex=0; fIndex<normalIndices.size(); ++fIndex) {
          const int iPop = normalIndices[fIndex];
          pi[iPi] += (T)2 * descriptors::c<L>(iPop,iAlpha)*descriptors::c<L>(iPop,iBeta)*
                     fNeq[iPop];
        }
        ++iPi;
      }
    }
  }

};  // struct boundaryHelpers

}  // namespace olb

// The specialized code is directly included. That is because we never want
// it to be precompiled so that in both the precompiled and the
// "include-everything" version, the compiler can apply all the
// optimizations it wants.
#include "lbHelpersD2Q5.h"
#include "lbHelpersD2Q9.h"
#include "lbHelpersD3Q7.h"
#include "lbHelpersD3Q15.h"
#include "lbHelpersD3Q19.h"
#include "lbHelpersD3Q27.h"

#endif
