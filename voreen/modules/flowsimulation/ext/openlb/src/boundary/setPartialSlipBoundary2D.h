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

//This file contains the Partial Slip Boundary
//This is a new version of the Boundary, which only contains free floating functions
#ifndef SET_PARTIAL_SLIP_BOUNDARY_2D_H
#define SET_PARTIAL_SLIP_BOUNDARY_2D_H


#include <vector>
#include "io/ostreamManager.h"
#include "utilities/functorPtr.h"
#include "geometry/superGeometry.h"
#include "geometry/blockGeometryStatistics2D.h"
#include "geometry/blockGeometry.h"
#include "core/superLattice2D.h"
#include "functors/lattice/indicator/superIndicatorF2D.h"
#include "functors/lattice/indicator/blockIndicatorF2D.h"
#include "dynamics/dynamics.h"
#include "boundaryPostProcessors2D.h"
#include "setBoundary2D.h"

namespace olb {

///Initialising the Partial slip boundary on the superLattice domain
template<typename T, typename DESCRIPTOR>
void setPartialSlipBoundary(SuperLattice<T, DESCRIPTOR>& sLattice, T tuner, SuperGeometry<T,2>& superGeometry, int material);

///Initialising the Partial slip boundary on the superLattice domain
template<typename T, typename DESCRIPTOR>
void setPartialSlipBoundary(SuperLattice<T, DESCRIPTOR>& sLattice, T tuner, FunctorPtr<SuperIndicatorF2D<T>>&& indicator);


///Set Partial Slip boundary for any indicated cells inside the block domain
template<typename T, typename DESCRIPTOR>
void setPartialSlipBoundary(BlockLattice<T,DESCRIPTOR>& block, T tuner, BlockIndicatorF2D<T>& indicator, bool includeOuterCells=false);

}//namespace olb


#endif
