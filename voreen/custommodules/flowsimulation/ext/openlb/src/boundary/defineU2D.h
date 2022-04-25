/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2020 Alexander Schulz
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

#ifndef DEFINE_U_2D_H
#define DEFINE_U_2D_H

#include <vector>
#include "io/ostreamManager.h"
#include "utilities/functorPtr.h"
#include "geometry/superGeometry.h"
#include "geometry/blockGeometryStatistics2D.h"
#include "core/superLattice2D.h"
#include "functors/lattice/indicator/superIndicatorF2D.h"
#include "functors/lattice/indicator/blockIndicatorF2D.h"
#include "dynamics/dynamics.h"
#include "geometry/blockGeometry.h"

//defineU for offLatticeBoundaryConditions
namespace olb {
////////// SuperLattice Domain  /////////////////////////////////////////
template<typename T, typename DESCRIPTOR>
void defineUBouzidi(SuperLattice<T, DESCRIPTOR>& sLattice, SuperGeometry<T,2>& superGeometry, int material,
                    AnalyticalF2D<T,T>& u,
                    olb::Vector<int, 2> bulkMaterials = olb::Vector<int, 2>(1,1));

template<typename T, typename DESCRIPTOR>
void defineUBouzidi(SuperLattice<T, DESCRIPTOR>& sLattice, FunctorPtr<SuperIndicatorF2D<T>>&& indicator,
                    FunctorPtr<SuperIndicatorF2D<T>>&& bulkIndicator,
                    AnalyticalF2D<T,T>& u);

////////// BlockLattice Domain  /////////////////////////////////////////

template<typename T, typename DESCRIPTOR>
void defineUBouzidi(BlockLattice<T,DESCRIPTOR>& block, BlockIndicatorF2D<T>& indicator, BlockIndicatorF2D<T>& bulkIndicator, AnalyticalF2D<T,T>& u);



template<typename T, typename DESCRIPTOR>
void defineUBouzidi(BlockLattice<T,DESCRIPTOR>& block, int iX, int iY, int iPop, const T u[DESCRIPTOR::d]);


template<typename T, typename DESCRIPTOR>
bool getBoundaryIntersection(BlockLattice<T,DESCRIPTOR>& block, int iX, int iY, int iPop, T point[DESCRIPTOR::d]);

template<typename T, typename DESCRIPTOR>
void setBoundaryIntersection(BlockLattice<T,DESCRIPTOR>& block, int iX, int iY, int iPop, T distance);

}//namespace olb


#endif
