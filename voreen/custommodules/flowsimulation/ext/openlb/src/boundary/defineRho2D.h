/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007 the OpenLB project
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

#ifndef DEFINE_RHO_2D_H
#define DEFINE_RHO_2D_H

#include <vector>
#include "io/ostreamManager.h"
#include "utilities/functorPtr.h"
#include "geometry/superGeometry2D.h"
#include "geometry/blockGeometryStatistics2D.h"
#include "core/superLattice2D.h"
#include "core/blockLatticeStructure2D.h"
#include "functors/lattice/indicator/superIndicatorF2D.h"
#include "functors/lattice/indicator/blockIndicatorF2D.h"
#include "dynamics/dynamics.h"
#include "geometry/blockGeometry2D.h"
#include "defineU2D.h"  ///for getBondaryIntersection
//defineU for offLatticeBoundaryConditions
namespace olb {
////////// SuperLattice Domain  /////////////////////////////////////////
template<typename T, typename DESCRIPTOR>
void defineRho(SuperLattice2D<T, DESCRIPTOR>& sLattice, SuperGeometry2D<T>& superGeometry, int material,
               AnalyticalF2D<T,T>& rho,
               std::vector<int> bulkMaterials = std::vector<int>(1,1));

template<typename T, typename DESCRIPTOR>
void defineRho(SuperLattice2D<T, DESCRIPTOR>& sLattice, FunctorPtr<SuperIndicatorF2D<T>>&& indicator,
               FunctorPtr<SuperIndicatorF2D<T>>&& bulkIndicator,
               AnalyticalF2D<T,T>&                rho);
////////// BlockLattice Domain  /////////////////////////////////////////

template<typename T, typename DESCRIPTOR>
void defineRho(BlockLatticeStructure2D<T,DESCRIPTOR>& block, BlockIndicatorF2D<T>& indicator, BlockIndicatorF2D<T>& bulkIndicator, AnalyticalF2D<T,T>& rho);

template<typename T, typename DESCRIPTOR>
void defineRho(BlockLatticeStructure2D<T,DESCRIPTOR>& block, int iX, int iY, int iPop, const T rho);


}//namespace olb


#endif
