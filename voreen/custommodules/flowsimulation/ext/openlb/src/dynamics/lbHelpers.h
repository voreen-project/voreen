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
template<typename T, class Descriptor> struct lbDynamicsHelpers;
template<typename T, template<typename U> class Lattice> struct lbExternalHelpers;
template<typename T, template<typename U> class Lattice> struct lbLatticeHelpers;

/// This structure forwards the calls to the appropriate helper class
template<typename T, template<typename U> class Lattice>
struct lbHelpers {

  static T equilibrium(int iPop, T rho, const T u[Lattice<T>::d], const T uSqr)
  {
    return lbDynamicsHelpers<T,typename Lattice<T>::BaseDescriptor>
           ::equilibrium(iPop, rho, u, uSqr);
  }

  static T equilibriumFirstOrder(int iPop, T rho, const T u[Lattice<T>::d])
  {
    return lbDynamicsHelpers<T,typename Lattice<T>::BaseDescriptor>
           ::equilibriumFirstOrder(iPop, rho, u);
  }

  static T incEquilibrium(int iPop, const T j[Lattice<T>::d], const T jSqr, const T pressure)
  {
    return lbDynamicsHelpers<T,typename Lattice<T>::BaseDescriptor>
           ::incEquilibrium(iPop, j, jSqr, pressure);
  }

  static void computeFneq ( Cell<T,Lattice> const& cell,
                            T fNeq[Lattice<T>::q], T rho, const T u[Lattice<T>::d] )
  {
    lbDynamicsHelpers<T,typename Lattice<T>::BaseDescriptor>::computeFneq(cell, fNeq, rho, u);
  }

  static T bgkCollision(Cell<T,Lattice>& cell, T const& rho, const T u[Lattice<T>::d], T const& omega)
  {
    return lbDynamicsHelpers<T,typename Lattice<T>::BaseDescriptor>
           ::bgkCollision(cell, rho, u, omega);
  }

  static T bgkCollision(Cell<T,Lattice>& cell, T const& rho, const T u[Lattice<T>::d], const T omega[Lattice<T>::q])
  {
    const T uSqr = util::normSqr<T,Lattice<T>::d>(u);
    for (int iPop=0; iPop < Lattice<T>::q; ++iPop) {
      cell[iPop] *= (T)1-omega[iPop];
      cell[iPop] += omega[iPop] * lbDynamicsHelpers<T,typename Lattice<T>::BaseDescriptor>::equilibrium (iPop, rho, u, uSqr );
    }
    return uSqr;
  }

  //====================================
  /*
   * Helpers included from the adLbDynamicsHelpers
   */
  static T rlbCollision( Cell<T, Lattice>& cell, T rho, const T u[Lattice<T>::d], T omega )
  {
    return lbDynamicsHelpers<T, typename Lattice<T>::BaseDescriptor>
           ::rlbCollision( cell, rho, u, omega );
  }

  static T sinkCollision( Cell<T, Lattice>& cell, T intensity, T omega, T sink )
  {
    return lbDynamicsHelpers<T, typename Lattice<T>::BaseDescriptor>
           ::sinkCollision( cell, intensity, omega, sink );
  }

  static T anisoCollision( Cell<T, Lattice>& cell, T intensity, T absorption, T scattering )
  {
    return lbDynamicsHelpers<T, typename Lattice<T>::BaseDescriptor>
           ::anisoCollision( cell, intensity, absorption, scattering );
  }

  //====================================
  /*
   * Helpers included from the advectionDiffusionMRTlbHelpers
   */
  static T mrtCollision( Cell<T, Lattice>& cell, T const& rho, const T u[Lattice<T>::d], T invM_S[Lattice<T>::q][Lattice<T>::q] ) {

    return lbDynamicsHelpers<T, typename Lattice<T>::BaseDescriptor>
           ::mrtCollision(cell, rho, u, invM_S );
  }
  //====================================

  static T incBgkCollision(Cell<T,Lattice>& cell, T pressure, const T j[Lattice<T>::d], T omega)
  {
    return lbDynamicsHelpers<T,typename Lattice<T>::BaseDescriptor>
           ::incBgkCollision(cell, pressure, j, omega);
  }

  static T constRhoBgkCollision(Cell<T,Lattice>& cell,
                                T rho, const T u[Lattice<T>::d], T ratioRho, T omega)
  {
    return lbDynamicsHelpers<T,typename Lattice<T>::BaseDescriptor>
           ::constRhoBgkCollision(cell, rho, u, ratioRho, omega);
  }

  static T computeRho(Cell<T,Lattice> const& cell)
  {
    return lbDynamicsHelpers<T,typename Lattice<T>::BaseDescriptor>
           ::computeRho(cell);
  }

  static void computeJ(Cell<T,Lattice> const& cell, T j[Lattice<T>::d] )
  {
    lbDynamicsHelpers<T,typename Lattice<T>::BaseDescriptor>
    ::computeJ(cell, j);
  }

  static void computeRhoU(Cell<T,Lattice> const& cell, T& rho, T u[Lattice<T>::d])
  {
    lbDynamicsHelpers<T,typename Lattice<T>::BaseDescriptor>
    ::computeRhoU(cell, rho, u);
  }

  static void computeRhoJ(Cell<T,Lattice> const& cell, T& rho, T j[Lattice<T>::d])
  {
    lbDynamicsHelpers<T,typename Lattice<T>::BaseDescriptor>
    ::computeRhoJ(cell, rho, j);
  }

  static void computeStress(Cell<T,Lattice> const& cell, T rho, const T u[Lattice<T>::d],
                            T pi[util::TensorVal<Lattice<T> >::n] )
  {
    lbDynamicsHelpers<T,typename Lattice<T>::BaseDescriptor>
    ::computeStress(cell, rho, u, pi);
  }

  static void computeAllMomenta(Cell<T,Lattice> const& cell, T& rho, T u[Lattice<T>::d],
                                T pi[util::TensorVal<Lattice<T> >::n] )
  {
    lbDynamicsHelpers<T,typename Lattice<T>::BaseDescriptor>
    ::computeAllMomenta(cell, rho, u, pi);
  }

  static void modifyVelocity(Cell<T,Lattice>& cell, const T newU[Lattice<T>::d])
  {
    lbDynamicsHelpers<T,typename Lattice<T>::BaseDescriptor>
    ::modifyVelocity(cell, newU);
  }

  static void addExternalForce(Cell<T,Lattice>& cell, const T u[Lattice<T>::d], T omega, T amplitude=(T)1)
  {
    lbExternalHelpers<T,Lattice>::addExternalForce(cell, u, omega, amplitude);
  }

  static void swapAndStream2D(Cell<T,Lattice> **grid, int iX, int iY)
  {
    lbLatticeHelpers<T,Lattice>::swapAndStream2D(grid, iX, iY);
  }

  static void swapAndStream3D(Cell<T,Lattice> ***grid, int iX, int iY, int iZ)
  {
    lbLatticeHelpers<T,Lattice>::swapAndStream3D(grid, iX, iY, iZ);
  }

};  // struct lbHelpers


/// All helper functions are inside this structure
template<typename T, class Descriptor>
struct lbDynamicsHelpers {
  /// Computation of equilibrium distribution
  static T equilibrium(int iPop, T rho, const T u[Descriptor::d], const T uSqr)
  {
    T c_u = T();
    for (int iD=0; iD < Descriptor::d; ++iD) {
      c_u += Descriptor::c[iPop][iD]*u[iD];
    }
    return rho * Descriptor::t[iPop] * ( (T)1 + Descriptor::invCs2 * c_u
                                              + Descriptor::invCs2 * Descriptor::invCs2 * (T)0.5 * c_u *c_u
                                              - Descriptor::invCs2 * (T)0.5 * uSqr )
                                              - Descriptor::t[iPop];
  }

  /// First order computation of equilibrium distribution
  // TODO AP: check if this equilibrium function is the needed one
  // -> 12.09.2017: the implementation is taken from advectionDiffusionLbHelpers.h
  static T equilibriumFirstOrder(int iPop, T rho, const T u[Descriptor::d])
  {
    T c_u = T();
    for (int iD=0; iD < Descriptor::d; ++iD) {
      c_u += Descriptor::c[iPop][iD]*u[iD];
    }
    return rho * Descriptor::t[iPop] * ( (T)1 + c_u * Descriptor::invCs2 ) - Descriptor::t[iPop];
  }

  static T incEquilibrium( int iPop, const T j[Descriptor::d],
                           const T jSqr, const T pressure )
  {
    T c_j = T();
    for (int iD=0; iD < Descriptor::d; ++iD) {
      c_j += Descriptor::c[iPop][iD]*j[iD];
    }
    T rho = (T)1 + pressure*Descriptor::invCs2;

    return Descriptor::t[iPop] * ( Descriptor::invCs2 * pressure +
                                   Descriptor::invCs2 * c_j +
                                   Descriptor::invCs2 * Descriptor::invCs2/(T)2 * c_j*c_j -
                                   Descriptor::invCs2/(T)2 * jSqr
                                 ) - Descriptor::t[iPop];
  }

  static void computeFneq(CellBase<T,Descriptor> const& cell, T fNeq[Descriptor::q], T rho, const T u[Descriptor::d])
  {
    const T uSqr = util::normSqr<T,Descriptor::d>(u);
    for (int iPop=0; iPop < Descriptor::q; ++iPop) {
      fNeq[iPop] = cell[iPop] - equilibrium(iPop, rho, u, uSqr);
    }
  }

  /// BGK collision step
  static T bgkCollision(CellBase<T,Descriptor>& cell, T const& rho, const T u[Descriptor::d], T const& omega)
  {
    const T uSqr = util::normSqr<T,Descriptor::d>(u);
    for (int iPop=0; iPop < Descriptor::q; ++iPop) {
      cell[iPop] *= (T)1-omega;
      cell[iPop] += omega * lbDynamicsHelpers<T,Descriptor>::equilibrium(iPop, rho, u, uSqr );
    }
    return uSqr;
  }

  /// Incompressible BGK collision step
  static T incBgkCollision(CellBase<T,Descriptor>& cell, T pressure, const T j[Descriptor::d], T omega)
  {
    const T jSqr = util::normSqr<T,Descriptor::d>(j);
    for (int iPop=0; iPop < Descriptor::q; ++iPop) {
      cell[iPop] *= (T)1-omega;
      cell[iPop] += omega * lbDynamicsHelpers<T,Descriptor>::incEquilibrium (
                      iPop, j, jSqr, pressure );
    }
    return jSqr;
  }

  /// BGK collision step with density correction
  static T constRhoBgkCollision(CellBase<T,Descriptor>& cell, T rho, const T u[Descriptor::d], T ratioRho, T omega)
  {
    const T uSqr = util::normSqr<T,Descriptor::d>(u);
    for (int iPop=0; iPop < Descriptor::q; ++iPop) {
      T feq = lbDynamicsHelpers<T,Descriptor>::equilibrium(iPop, rho, u, uSqr );
      cell[iPop] =
        ratioRho*(feq+Descriptor::t[iPop])-Descriptor::t[iPop] +
        ((T)1-omega)*(cell[iPop]-feq);
    }
    return uSqr;
  }

  //========================================

  /// RLB advection diffusion collision step
  static T rlbCollision(CellBase<T, Descriptor>& cell, T rho, const T u[Descriptor::d], T omega )
  {
    const T uSqr = util::normSqr<T, Descriptor::d>( u );
    // First-order moment for the regularization
    T j1[Descriptor::d];
    for ( int iD = 0; iD < Descriptor::d; ++iD ) {
      j1[iD] = T();
    }

    T fEq[Descriptor::q];
    for ( int iPop = 0; iPop < Descriptor::q; ++iPop ) {
      fEq[iPop] = lbDynamicsHelpers<T, Descriptor>::equilibriumFirstOrder( iPop, rho, u );
      for ( int iD = 0; iD < Descriptor::d; ++iD ) {
        j1[iD] += Descriptor::c[iPop][iD] * ( cell[iPop] - fEq[iPop] );
      }
    }

    // Collision step
    for ( int iPop = 0; iPop < Descriptor::q; ++iPop ) {
      T fNeq = T();
      for ( int iD = 0; iD < Descriptor::d; ++iD ) {
        fNeq += Descriptor::c[iPop][iD] * j1[iD];
      }
      fNeq *= Descriptor::t[iPop] * Descriptor::invCs2;
      cell[iPop] = fEq[iPop] + ( (T)1 - omega ) * fNeq;
    }
    return uSqr;
  }

  /// Paper: Mink et al. 2016 DOI: 10.1016/j.jocs.2016.03.014
  static T sinkCollision(CellBase<T, Descriptor>& cell, T intensity, T omega, T sink){
    // collision step only valid for D3Q7 lattices
    // spezialisation implemented in lbHelpersD3Q7.h
    for ( int iPop = 0; iPop < Descriptor::q; ++iPop ) {
      cell[iPop] = (1 - omega) * (cell[iPop] + Descriptor::t[iPop])
                    + omega * Descriptor::t[iPop] * intensity
                    - sink*(cell[iPop] + Descriptor::t[iPop]) - Descriptor::t[iPop];
    }
    return T(0);
  }

  static T anisoCollision(CellBase<T, Descriptor>& cell, T intensity, T absorption, T scattering)
  {
    //    T rho[Descriptor::q];
    //    for (int iPop = 0; iPop < Descriptor::q; ++iPop) {
    //      rho[iPop] = T(0);
    //    }
    //
    //    for (int iPop = 0; iPop < Descriptor::q; iPop++) {
    //      for ( int jPop = 0; jPop < Descriptor::q; ++jPop ) {
    //        rho[iPop] += (cell[jPop] + Descriptor::t[jPop]); // * Descriptor::henyeyPhaseFunction[jPop][iPop];
    //      }
    //    }

    //    double feq = adLbDynamicsHelpers<T, Descriptor>::equilibriumFirstOrder(iPop, rhoBGK, 0.0);
    //    for (int iPop = 0; iPop < Descriptor::q; ++iPop) {
    //      cell[iPop] = (cell[iPop] + Descriptor::t[iPop])
    //                   - ( cell[iPop] + Descriptor::t[iPop] - intensity * Descriptor::t[iPop] )
    //                   - absorption/scattering * scaling * ( cell[iPop] + Descriptor::t[iPop] )
    //                   - Descriptor::t[iPop];
    //    }

    for (int iPop = 0; iPop < Descriptor::q; ++iPop) {
      cell[iPop] = (cell[iPop] + Descriptor::t[iPop])
                  - scattering * ( cell[iPop] + Descriptor::t[iPop] - intensity * Descriptor::t[iPop] )
                  - absorption * ( cell[iPop] + Descriptor::t[iPop] )
                  - Descriptor::t[iPop];
    }
    return T(0);
  }

  //====================================
  /*
   * Helpers included from the advectionDiffusionMRTlbHelpers
   */
///// Computation of all equilibrium distribution (in momenta space)
  static void computeMomentaEquilibrium( T momentaEq[Descriptor::q], T rho, const T u[Descriptor::d], T uSqr )
  {
    for (int iPop = 0; iPop < Descriptor::q; ++iPop) {
      momentaEq[iPop] = T();
      for (int jPop = 0; jPop < Descriptor::q; ++jPop) {
        momentaEq[iPop] += Descriptor::M[iPop][jPop] *
                (lbDynamicsHelpers<T, Descriptor>::equilibrium(jPop,rho,u,uSqr) + Descriptor::t[jPop]);
      }
    }
  }

  static void computeMomenta(T momenta[Descriptor::q], CellBase<T, Descriptor>& cell)
  {
    for (int iPop = 0; iPop < Descriptor::q; ++iPop) {
      momenta[iPop] = T();
      for (int jPop = 0; jPop < Descriptor::q; ++jPop) {
        momenta[iPop] += Descriptor::M[iPop][jPop] *
                         (cell[jPop] + Descriptor::t[jPop]);
      }
    }
  }

  static T mrtCollision( CellBase<T,Descriptor>& cell, T const& rho, const T u[Descriptor::d], T invM_S[Descriptor::q][Descriptor::q] )
  {
    //// Implemented in advectionDiffusionMRTlbHelpers2D.h and advectionDiffusionMRTlbHelpers3D.h
    T uSqr = util::normSqr<T, Descriptor::d>(u);
    T momenta[Descriptor::q];
    T momentaEq[Descriptor::q];

    computeMomenta(momenta, cell);
    computeMomentaEquilibrium(momentaEq, rho, u, uSqr);

//    std::cout << "momenta = ";
//    for (int i=0; i < Descriptor::q; ++i) {
//        std::cout << momenta[i] << ", ";
//    }
//    std::cout << std::endl;

//    std::cout << "momentaEq = ";
//    for (int i=0; i < Descriptor::q; ++i) {
//        std::cout << momentaEq[i] << ", ";
//    }
//    std::cout << std::endl;

    for (int iPop = 0; iPop < Descriptor::q; ++iPop) {
      T collisionTerm = T();
      for (int jPop = 0; jPop < Descriptor::q; ++jPop) {
        collisionTerm += invM_S[iPop][jPop] * (momenta[jPop] - momentaEq[jPop]);
      }
      cell[iPop] -= collisionTerm;
    }
    return uSqr;
  }
  //========================================

  /// Computation of density
  static T computeRho(CellBase<T,Descriptor> const& cell)
  {
    T rho = T();
    for (int iPop=0; iPop < Descriptor::q; ++iPop) {
      rho += cell[iPop];
    }
    rho += (T)1;
    return rho;
  }

  /// Computation of momentum
  static void computeJ(CellBase<T,Descriptor> const& cell, T j[Descriptor::d])
  {
    for (int iD=0; iD < Descriptor::d; ++iD) {
      j[iD] = T();
    }
    for (int iPop=0; iPop < Descriptor::q; ++iPop) {
      for (int iD=0; iD < Descriptor::d; ++iD) {
        j[iD] += cell[iPop]*Descriptor::c[iPop][iD];
      }
    }
  }

  /// Computation of hydrodynamic variables
  static void computeRhoU(CellBase<T,Descriptor> const& cell, T& rho, T u[Descriptor::d])
  {
    rho = T();
    for (int iD=0; iD < Descriptor::d; ++iD) {
      u[iD] = T();
    }
    for (int iPop=0; iPop < Descriptor::q; ++iPop) {
      rho += cell[iPop];
      for (int iD=0; iD < Descriptor::d; ++iD) {
        u[iD] += cell[iPop]*Descriptor::c[iPop][iD];
      }
    }
    rho += (T)1;
    for (int iD=0; iD < Descriptor::d; ++iD) {
      u[iD] /= rho;
    }
  }

  /// Computation of hydrodynamic variables
  static void computeRhoJ(CellBase<T,Descriptor> const& cell, T& rho, T j[Descriptor::d])
  {
    rho = T();
    for (int iD=0; iD < Descriptor::d; ++iD) {
      j[iD] = T();
    }
    for (int iPop=0; iPop < Descriptor::q; ++iPop) {
      rho += cell[iPop];
      for (int iD=0; iD < Descriptor::d; ++iD) {
        j[iD] += cell[iPop]*Descriptor::c[iPop][iD];
      }
    }
    rho += (T)1;
  }

  /// Computation of stress tensor
  static void computeStress(CellBase<T,Descriptor> const& cell, T rho, const T u[Descriptor::d],
                            T pi[util::TensorVal<Descriptor>::n] )
  {
    int iPi = 0;
    for (int iAlpha=0; iAlpha < Descriptor::d; ++iAlpha) {
      for (int iBeta=iAlpha; iBeta < Descriptor::d; ++iBeta) {
        pi[iPi] = T();
        for (int iPop=0; iPop < Descriptor::q; ++iPop) {
          pi[iPi] += Descriptor::c[iPop][iAlpha]*
                     Descriptor::c[iPop][iBeta] * cell[iPop];
        }
        // stripe off equilibrium contribution
        pi[iPi] -= rho*u[iAlpha]*u[iBeta];
        if (iAlpha==iBeta) {
          pi[iPi] -= 1./Descriptor::invCs2*(rho-(T)1);
        }
        ++iPi;
      }
    }
  }

  /// Computation of all hydrodynamic variables
  static void computeAllMomenta(CellBase<T,Descriptor> const& cell, T& rho, T u[Descriptor::d],
                                T pi[util::TensorVal<Descriptor>::n] )
  {
    computeRhoU(cell, rho, u);
    computeStress(cell, rho, u, pi);
  }

  static void modifyVelocity(CellBase<T,Descriptor>& cell, const T newU[Descriptor::d])
  {
    T rho, oldU[Descriptor::d];
    computeRhoU(cell, rho, oldU);
    const T oldUSqr = util::normSqr<T,Descriptor::d>(oldU);
    const T newUSqr = util::normSqr<T,Descriptor::d>(newU);
    for (int iPop=0; iPop<Descriptor::q; ++iPop) {
      cell[iPop] = cell[iPop]
                   - equilibrium(iPop, rho, oldU, oldUSqr)
                   + equilibrium(iPop, rho, newU, newUSqr);
    }
  }

};  // struct lbDynamicsHelpers

/// Helper functions for dynamics that access external field
template<typename T, template<typename U> class Lattice>
struct lbExternalHelpers {
  /// Add a force term after BGK collision
  static void addExternalForce(Cell<T,Lattice>& cell, const T u[Lattice<T>::d], T omega, T amplitude)
  {
    static const int forceBeginsAt = Lattice<T>::ExternalField::forceBeginsAt;
    T* force = cell.getExternal(forceBeginsAt);
    for (int iPop=0; iPop < Lattice<T>::q; ++iPop) {
      T c_u = T();
      for (int iD=0; iD < Lattice<T>::d; ++iD) {
        c_u += Lattice<T>::c[iPop][iD]*u[iD];
      }
      c_u *= Lattice<T>::invCs2*Lattice<T>::invCs2;
      T forceTerm = T();
      for (int iD=0; iD < Lattice<T>::d; ++iD) {
        forceTerm +=
          (   ((T)Lattice<T>::c[iPop][iD]-u[iD]) * Lattice<T>::invCs2
              + c_u * Lattice<T>::c[iPop][iD]
          )
          * force[iD];
      }
      forceTerm *= Lattice<T>::t[iPop];
      forceTerm *= T(1) - omega/T(2);
      forceTerm *= amplitude;
      cell[iPop] += forceTerm;
    }
  }
};  // struct externalFieldHelpers

/// Helper functions with full-lattice access
template<typename T, template<typename U> class Lattice>
struct lbLatticeHelpers {
  /// Swap ("bounce-back") values of a cell (2D), and apply streaming step
  static void swapAndStream2D(Cell<T,Lattice> **grid, int iX, int iY)
  {
    const int half = Lattice<T>::q/2;
    for (int iPop=1; iPop<=half; ++iPop) {
      int nextX = iX + Lattice<T>::c[iPop][0];
      int nextY = iY + Lattice<T>::c[iPop][1];
      T fTmp                   = grid[iX][iY][iPop];
      grid[iX][iY][iPop]       = grid[iX][iY][iPop+half];
      grid[iX][iY][iPop+half]  = grid[nextX][nextY][iPop];
      grid[nextX][nextY][iPop] = fTmp;
    }
  }

  /// Swap ("bounce-back") values of a cell (3D), and apply streaming step
  static void swapAndStream3D(Cell<T,Lattice> ***grid,
                              int iX, int iY, int iZ)
  {
    const int half = Lattice<T>::q/2;
    for (int iPop=1; iPop<=half; ++iPop) {
      int nextX = iX + Lattice<T>::c[iPop][0];
      int nextY = iY + Lattice<T>::c[iPop][1];
      int nextZ = iZ + Lattice<T>::c[iPop][2];
      T fTmp                          = grid[iX][iY][iZ][iPop];
      grid[iX][iY][iZ][iPop]          = grid[iX][iY][iZ][iPop+half];
      grid[iX][iY][iZ][iPop+half]     = grid[nextX][nextY][nextZ][iPop];
      grid[nextX][nextY][nextZ][iPop] = fTmp;
    }
  }
};

/// All boundary helper functions are inside this structure
template<typename T, template<typename U> class Lattice, int direction, int orientation>
struct BoundaryHelpers {
  static void computeStress (
    Cell<T,Lattice> const& cell, T rho, const T u[Lattice<T>::d],
    T pi[util::TensorVal<Lattice<T> >::n] )
  {
    typedef Lattice<T> L;
    const T uSqr = util::normSqr<T,L::d>(u);

    std::vector<int> const& onWallIndices = util::subIndex<L, direction, 0>();
    std::vector<int> const& normalIndices = util::subIndex<L, direction, orientation>();

    T fNeq[Lattice<T>::q];
    for (unsigned fIndex=0; fIndex<onWallIndices.size(); ++fIndex) {
      int iPop = onWallIndices[fIndex];
      fNeq[iPop] =
        cell[iPop] -
        lbHelpers<T,Lattice>::equilibrium(iPop, rho, u, uSqr);
    }
    for (unsigned fIndex=0; fIndex<normalIndices.size(); ++fIndex) {
      int iPop = normalIndices[fIndex];
      if (iPop == 0) {
        fNeq[iPop] = T();  // fNeq[0] will not be used anyway
      } else {
        fNeq[iPop] =
          cell[iPop] -
          lbHelpers<T,Lattice>::equilibrium(iPop, rho, u, uSqr);
      }
    }

    int iPi = 0;
    for (int iAlpha=0; iAlpha<L::d; ++iAlpha) {
      for (int iBeta=iAlpha; iBeta<L::d; ++iBeta) {
        pi[iPi] = T();
        for (unsigned fIndex=0; fIndex<onWallIndices.size(); ++fIndex) {
          const int iPop = onWallIndices[fIndex];
          pi[iPi] +=
            L::c[iPop][iAlpha]*L::c[iPop][iBeta]*fNeq[iPop];
        }
        for (unsigned fIndex=0; fIndex<normalIndices.size(); ++fIndex) {
          const int iPop = normalIndices[fIndex];
          pi[iPi] += (T)2 * L::c[iPop][iAlpha]*L::c[iPop][iBeta]*
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
