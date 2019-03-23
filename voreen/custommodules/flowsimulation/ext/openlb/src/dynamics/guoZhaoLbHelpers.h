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

#include "dynamics/guoZhaoLatticeDescriptors.h"
#include "core/cell.h"
#include "core/util.h"


namespace olb {


// Forward declarations
template<typename T, class Descriptor> struct GuoZhaoLbDynamicsHelpers;
template<typename T, template<typename U> class Lattice> struct GuoZhaoLbExternalHelpers;

/// This structure forwards the calls to the appropriate Guo Zhao helper class
template<typename T, template<typename U> class Lattice>
struct GuoZhaoLbHelpers {

  static T equilibrium(int iPop, T epsilon, T rho, const T u[Lattice<T>::d], const T uSqr) {
    return GuoZhaoLbDynamicsHelpers<T,typename Lattice<T>::BaseDescriptor>
           ::equilibrium(iPop, epsilon, rho, u, uSqr);
  }

  static T bgkCollision(Cell<T,Lattice>& cell, T const& epsilon, T const& rho, const T u[Lattice<T>::d], T const& omega) {
    return GuoZhaoLbDynamicsHelpers<T,typename Lattice<T>::BaseDescriptor>
           ::bgkCollision(cell, epsilon, rho, u, omega);
  }

  static void updateGuoZhaoForce(Cell<T,Lattice>& cell, const T u[Lattice<T>::d]) {
    GuoZhaoLbExternalHelpers<T,Lattice>::updateGuoZhaoForce(cell, u);
  }

};  // struct GuoZhaoLbHelpers


/// All Guo Zhao helper functions are inside this structure
template<typename T, class Descriptor>
struct GuoZhaoLbDynamicsHelpers {

  /// Computation of Guo Zhao equilibrium distribution - original (compressible) formulation following Guo and Zhao (2002).
  static T equilibrium(int iPop, T epsilon, T rho, const T u[Descriptor::d], const T uSqr) {
    T c_u = T();
    for (int iD=0; iD < Descriptor::d; ++iD) {
      c_u += Descriptor::c[iPop][iD]*u[iD];
    }
    return rho * Descriptor::t[iPop] * (
             (T)1 + Descriptor::invCs2 * c_u +
             Descriptor::invCs2 * Descriptor::invCs2/((T)2*epsilon) * c_u*c_u -
             Descriptor::invCs2/((T)2*epsilon) * uSqr
           ) - Descriptor::t[iPop];
  }

  /// Guo Zhao BGK collision step
  static T bgkCollision(CellBase<T,Descriptor>& cell, T const& epsilon, T const& rho, const T u[Descriptor::d], T const& omega) {
    const T uSqr = util::normSqr<T,Descriptor::d>(u);
    for (int iPop=0; iPop < Descriptor::q; ++iPop) {
      cell[iPop] *= (T)1-omega;
      cell[iPop] += omega * GuoZhaoLbDynamicsHelpers<T,Descriptor>::equilibrium (
                      iPop, epsilon, rho, u, uSqr );
    }
    return uSqr;
  }

};  // struct GuoZhaoLbDynamicsHelpers

/// Helper functions for dynamics that access external field
template<typename T, template<typename U> class Lattice>
/// Updates Guo Zhao porous force
struct GuoZhaoLbExternalHelpers {
  static void updateGuoZhaoForce(Cell<T,Lattice>& cell, const T u[Lattice<T>::d]) {
    T epsilon = *cell.getExternal(Lattice<T>::ExternalField::epsilonAt);
    T k       = *cell.getExternal(Lattice<T>::ExternalField::KAt);
    T nu      = *cell.getExternal(Lattice<T>::ExternalField::nuAt);
    T bodyF0  = *cell.getExternal(Lattice<T>::ExternalField::bodyForceBeginsAt);
    T bodyF1  = *cell.getExternal(Lattice<T>::ExternalField::bodyForceBeginsAt+1);

    T* force0 = cell.getExternal(Lattice<T>::ExternalField::forceBeginsAt);
    T* force1 = cell.getExternal(Lattice<T>::ExternalField::forceBeginsAt+1);

    *force0 = -u[0]*epsilon*nu/k + bodyF0*epsilon;
    *force1 = -u[1]*epsilon*nu/k + bodyF1*epsilon;
  }
};  // struct GuoZhaoLbExternalHelpers


}  // namespace olb


#endif
