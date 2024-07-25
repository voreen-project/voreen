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

#ifndef BOUNDARY_HELPER_H
#define BOUNDARY_HELPER_H

#include "core/cell.h"
#include "core/util.h"

namespace olb {

/// All boundary helper functions are inside this structure
template<typename T, typename DESCRIPTOR, int direction, int orientation>
struct BoundaryHelpers {
  template <typename CELL>
  static void computeStress(
    CELL& cell, T rho, const T u[DESCRIPTOR::d],
    T pi[util::TensorVal<DESCRIPTOR>::n] ) any_platform
  {
    const T uSqr = util::normSqr<T,DESCRIPTOR::d>(u);

    constexpr auto onWallIndices = util::populationsContributingToVelocity<DESCRIPTOR,direction,0>();
    constexpr auto normalIndices = util::populationsContributingToVelocity<DESCRIPTOR,direction,orientation>();

    T fNeq[DESCRIPTOR::q];
    for (unsigned fIndex=0; fIndex<onWallIndices.size(); ++fIndex) {
      int iPop = onWallIndices[fIndex];
      if (iPop == 0) {
        fNeq[0] = T();  // fNeq[0] will not be used anyway
      }
      else {
        fNeq[iPop] =
          cell[iPop] -
          equilibrium<DESCRIPTOR>::secondOrder(iPop, rho, u, uSqr);
      }
    }
    for (unsigned fIndex=0; fIndex<normalIndices.size(); ++fIndex) {
      int iPop = normalIndices[fIndex];
      fNeq[iPop] =
        cell[iPop] -
        equilibrium<DESCRIPTOR>::secondOrder(iPop, rho, u, uSqr);
    }

    int iPi = 0;
    for (int iAlpha=0; iAlpha < DESCRIPTOR::d; ++iAlpha) {
      for (int iBeta=iAlpha; iBeta < DESCRIPTOR::d; ++iBeta) {
        pi[iPi] = T();
        for (unsigned fIndex=0; fIndex < onWallIndices.size(); ++fIndex) {
          const int iPop = onWallIndices[fIndex];
          pi[iPi] += descriptors::c<DESCRIPTOR>(iPop,iAlpha)
                   * descriptors::c<DESCRIPTOR>(iPop,iBeta)
                   * fNeq[iPop];
        }
        for (unsigned fIndex=0; fIndex < normalIndices.size(); ++fIndex) {
          const int iPop = normalIndices[fIndex];
          pi[iPi] += T{2}
                   * descriptors::c<DESCRIPTOR>(iPop,iAlpha)
                   * descriptors::c<DESCRIPTOR>(iPop,iBeta)
                   * fNeq[iPop];
        }
        ++iPi;
      }
    }
  }

};

}

#endif
