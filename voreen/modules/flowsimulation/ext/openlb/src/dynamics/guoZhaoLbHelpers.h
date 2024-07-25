/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2016-2017 Davide Dapelo, Mathias J. Krause
 *  OpenLB e-mail contact: info@openlb.net
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
 * Helper functions for the implementation of the
 * Guo-ZhapLB dynamics.
 */

#ifndef LB_GUOZHAO_HELPERS_H
#define LB_GUOZHAO_HELPERS_H

#include "dynamics/descriptorAlias.h"
#include "core/cell.h"
#include "core/util.h"


namespace olb {


// Forward declarations
template<typename T, typename DESCRIPTOR> struct GuoZhaoLbDynamicsHelpers;
template<typename T, typename DESCRIPTOR> struct GuoZhaoLbExternalHelpers;

/// This structure forwards the calls to the appropriate Guo Zhao helper class
template<typename T, typename DESCRIPTOR>
struct GuoZhaoLbHelpers {

  static T equilibrium(int iPop, T epsilon, T rho, const T u[DESCRIPTOR::d], const T uSqr)
  {
    return GuoZhaoLbDynamicsHelpers<T,DESCRIPTOR>
           ::equilibrium(iPop, epsilon, rho, u, uSqr);
  }

  static T bgkCollision(Cell<T,DESCRIPTOR>& cell, T const& epsilon, T const& rho, const T u[DESCRIPTOR::d], T const& omega)
  {
    return GuoZhaoLbDynamicsHelpers<T,DESCRIPTOR>
           ::bgkCollision(cell, epsilon, rho, u, omega);
  }

  static T bgkCollisionAndForcing(Cell<T,DESCRIPTOR>& cell, T const& rho, const T u[DESCRIPTOR::d], T const& omega) {
    return GuoZhaoLbDynamicsHelpers<T,DESCRIPTOR>
           ::bgkCollisionAndForcing(cell, rho, u, omega);
  }

  static void updateGuoZhaoForce(Cell<T,DESCRIPTOR>& cell, const T u[DESCRIPTOR::d]) 
  {
    GuoZhaoLbExternalHelpers<T,DESCRIPTOR>::updateGuoZhaoForce(cell, u);
  }

  static void updateGuoForceSimple(Cell<T,DESCRIPTOR>& cell, const T u[DESCRIPTOR::d]) {
    GuoZhaoLbExternalHelpers<T,DESCRIPTOR>::updateGuoForceSimple(cell, u);
  }

  static void addExternalForce(Cell<T,DESCRIPTOR>& cell, const T u[DESCRIPTOR::d], T omega, T rho, T epsilon=(T)1)
  {
    GuoZhaoLbExternalHelpers<T,DESCRIPTOR>::addExternalForce(cell, u, omega, rho, epsilon);
  }

};  // struct GuoZhaoLbHelpers


/// All Guo Zhao helper functions are inside this structure
template<typename T, typename DESCRIPTOR>
struct GuoZhaoLbDynamicsHelpers {

  /// Computation of Guo Zhao equilibrium distribution - original (compressible) formulation following Guo and Zhao (2002).
  static T forceEquilibrium(int iPop, T epsilon, T rho, const T u[DESCRIPTOR::d], const T force[DESCRIPTOR::d], T nu)
  {
  }

  /// Computation of Guo Zhao equilibrium distribution - original (compressible) formulation following Guo and Zhao (2002).
  static T equilibrium(int iPop, T epsilon, T rho, const T u[DESCRIPTOR::d], const T uSqr)
  {
    T c_u = T();
    for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
      c_u += descriptors::c<DESCRIPTOR>(iPop,iD)*u[iD];
    }
    return rho * descriptors::t<T,DESCRIPTOR>(iPop) * (
             (T)1 + descriptors::invCs2<T,DESCRIPTOR>() * c_u +
             descriptors::invCs2<T,DESCRIPTOR>() * descriptors::invCs2<T,DESCRIPTOR>()/((T)2*epsilon) * c_u*c_u -
             descriptors::invCs2<T,DESCRIPTOR>()/((T)2*epsilon) * uSqr
           ) - descriptors::t<T,DESCRIPTOR>(iPop);
  }

  /// Computation of Guo Zhao equilibrium distribution with increased efficiency (Hafen2020)
  static T equilibriumReturnCu(int iPop, T epsilon, T rho, const T u[DESCRIPTOR::d], const T uSqr, T& c_u) {
    for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
      c_u += descriptors::c<DESCRIPTOR>(iPop,iD)*u[iD];
    }
    return rho * descriptors::t<T,DESCRIPTOR>(iPop) * (
             (T)1 + descriptors::invCs2<T,DESCRIPTOR>() * c_u +
             descriptors::invCs2<T,DESCRIPTOR>() * descriptors::invCs2<T,DESCRIPTOR>()/((T)2*epsilon) * c_u*c_u -
             descriptors::invCs2<T,DESCRIPTOR>()/((T)2*epsilon) * uSqr
           ) - descriptors::t<T,DESCRIPTOR>(iPop);
  }

  /// Guo Zhao BGK collision step
  static T bgkCollision(Cell<T,DESCRIPTOR>& cell, T const& epsilon, T const& rho, const T u[DESCRIPTOR::d], T const& omega)
  {
    const T uSqr = util::normSqr<T,DESCRIPTOR::d>(u);
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      cell[iPop] *= (T)1-omega;
      cell[iPop] += omega * GuoZhaoLbDynamicsHelpers<T,DESCRIPTOR>::equilibrium (
                      iPop, epsilon, rho, u, uSqr );
    }
    return uSqr;
  }

  // Guo Zhao BGK collision and forcing step for simplified Guo dynamics (Hafen2020)
  static T bgkCollisionAndForcing(Cell<T,DESCRIPTOR>& cell, T const& rho, const T u[DESCRIPTOR::d], T const& omega) {
    auto force = cell.template getFieldPointer<descriptors::FORCE>();
    const T epsilon = cell.template getField<descriptors::EPSILON>();
    const T uSqr = util::normSqr<T,DESCRIPTOR::d>(u);
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      T c_u = T();
      //Do BGK collision
      cell[iPop] *= (T)1-omega;
      cell[iPop] += omega * GuoZhaoLbDynamicsHelpers<T,DESCRIPTOR>::equilibriumReturnCu (
                      iPop, epsilon, rho, u, uSqr, c_u );
      //Add Guo source term
      c_u *= descriptors::invCs2<T,DESCRIPTOR>()*descriptors::invCs2<T,DESCRIPTOR>()/epsilon;
      T forceTerm = T();
      for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
        forceTerm +=
          (   (epsilon*(T)descriptors::c<DESCRIPTOR>(iPop,iD)-u[iD]) * descriptors::invCs2<T,DESCRIPTOR>()/epsilon
              + c_u * descriptors::c<DESCRIPTOR>(iPop,iD)
          )
          * force[iD];
      }
      forceTerm *= descriptors::t<T,DESCRIPTOR>(iPop);
      forceTerm *= T(1) - omega/T(2);
      forceTerm *= rho;
      cell[iPop] += forceTerm;
    }
    return uSqr;
  }

};  // struct GuoZhaoLbDynamicsHelpers

/// Helper functions for dynamics that access external field
template<typename T, typename DESCRIPTOR>
/// Updates Guo Zhao porous force
struct GuoZhaoLbExternalHelpers {
  static void updateGuoZhaoForce(Cell<T,DESCRIPTOR>& cell, const T u[DESCRIPTOR::d])
  {
    const T epsilon = cell.template getField<descriptors::EPSILON>();
    const T k       = cell.template getField<descriptors::K>();
    const T nu      = cell.template getField<descriptors::NU>();
    auto bodyF      = cell.template getField<descriptors::BODY_FORCE>();

    auto force = cell.template getFieldPointer<descriptors::FORCE>();

    const T uMag = util::sqrt( util::normSqr<T,DESCRIPTOR::d>(u) );
    const T Fe = 0;//1.75/util::sqrt(150.*util::pow(epsilon,3));

    // Linear Darcy term, nonlinear Forchheimer term and body force
    for (int iDim=0; iDim <DESCRIPTOR::d; iDim++) {
      force[iDim] = -u[iDim]*epsilon*nu/k - epsilon*Fe/util::sqrt(k)*uMag*u[iDim] + bodyF[iDim]*epsilon;
    }
  }

  /// Simplified version
  static void updateGuoForceSimple(Cell<T,DESCRIPTOR>& cell, const T u[DESCRIPTOR::d]) {
    const T porosity = cell.template getField<descriptors::POROSITY>();
    auto force = cell.template getFieldPointer<descriptors::FORCE>();
    for (int iDim=0; iDim <DESCRIPTOR::d; iDim++) {
      force[iDim] = -u[iDim]*2.*(1./porosity-1.);
    }
  }


  /// Add a force term scaled by physical porosity epsilon after BGK collision
  static void addExternalForce(Cell<T,DESCRIPTOR>& cell, const T u[DESCRIPTOR::d], T omega, T rho, T epsilon)
  {
    auto force = cell.template getFieldPointer<descriptors::FORCE>();
    for (int iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      T c_u = T();
      for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
        c_u += descriptors::c<DESCRIPTOR>(iPop,iD)*u[iD];
      }
      c_u *= descriptors::invCs2<T,DESCRIPTOR>()*descriptors::invCs2<T,DESCRIPTOR>()/epsilon;
      T forceTerm = T();
      for (int iD=0; iD < DESCRIPTOR::d; ++iD) {
        forceTerm +=
          (   (epsilon*(T)descriptors::c<DESCRIPTOR>(iPop,iD)-u[iD]) * descriptors::invCs2<T,DESCRIPTOR>()/epsilon
              + c_u * descriptors::c<DESCRIPTOR>(iPop,iD)
          )
          * force[iD];
      }
      forceTerm *= descriptors::t<T,DESCRIPTOR>(iPop);
      forceTerm *= T(1) - omega/T(2);
      forceTerm *= rho;
      cell[iPop] += forceTerm;
    }
  }
};  // struct GuoZhaoLbExternalHelpers

}  // namespace olb


#endif
