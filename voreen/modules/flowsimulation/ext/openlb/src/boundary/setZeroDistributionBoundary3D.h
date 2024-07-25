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

//This file contains the Zero Distribution Boundary.
//This is a new version of the Boundary, which only contains free floating functions

#ifndef SET_ZERO_DISTRIBUTION_BOUNDARY_3D_H
#define SET_ZERO_DISTRIBUTION_BOUNDARY_3D_H

#include <vector>
#include "io/ostreamManager.h"
#include "utilities/functorPtr.h"
#include "geometry/superGeometry.h"
#include "core/superLattice.h"
#include "functors/lattice/indicator/superIndicatorBaseF3D.h"
#include "functors/lattice/indicator/blockIndicatorF3D.h"
#include "dynamics/dynamics.h"
#include "dynamics/freeEnergyDynamics.h"
#include "advectionDiffusionBoundaryPostProcessor3D.h"
#include "advectionDiffusionBoundaries.h"
#include "extendedFiniteDifferenceBoundary3D.h"
#include "setBoundary3D.h"

namespace olb {

//Initialising the setZeroDistributionBoundary function on the superLattice domain
//advection diffusion boundary condition therefore mostly --> MixinDynamics = AdvectionDiffusionRLBdynamics<T,DESCRIPTOR>>
template<typename T, typename DESCRIPTOR>
void setZeroDistributionBoundary(SuperLattice<T, DESCRIPTOR>& sLattice,SuperGeometry<T,3>& superGeometry, int material);

//setZeroDistributionBoundary function on the superLattice domain
//depending on the application, the first function can be skipped and this function can be called directly in the app
//the addPoints2CommBC is currently located inside setBoundary3D.h
template<typename T, typename DESCRIPTOR>
void setZeroDistributionBoundary(SuperLattice<T, DESCRIPTOR>& sLattice, FunctorPtr<SuperIndicatorF3D<T>>&& indicator);


//sets the ZeroDistributionBoundary for any indicated cell inside the block domain
template<typename T, typename DESCRIPTOR>
void setZeroDistributionBoundary(BlockLattice<T,DESCRIPTOR>& _block, BlockIndicatorF3D<T>& indicator, bool includeOuterCells=false);

}//namespace olb


#endif

