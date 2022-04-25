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

//This file contains the Free Energy Wall Boundary
//This is a new version of the Boundary, which only contains free floating functions
#ifndef SET_FREE_ENERGY_WALL_BOUNDARY_2D_H
#define SET_FREE_ENERGY_WALL_BOUNDARY_2D_H


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
#include "dynamics/freeEnergyDynamics.h"
#include "boundaryPostProcessors2D.h"
#include "setBoundary2D.h"


namespace olb {
    
/// Implementation of a wetting boundary condition for the ternary free energy model, consisting of a BounceBack
/// dynamics and an FreeEnergyWall PostProcessor.
/// \param[in] alpha_ - Parameter related to the interface width. [lattice units]
/// \param[in] kappa1_ - Parameter related to surface tension. [lattice units]
/// \param[in] kappa2_ - Parameter related to surface tension. [lattice units]
/// \param[in] h1_ - Parameter related to resulting contact angle of the boundary. [lattice units]
/// \param[in] h2_ - Parameter related to resulting contact angle of the boundary. [lattice units]
/// \param[in] latticeNumber - determines the number of the free energy lattice to set the boundary accordingly


///Initialising the Free Energy Wall Boundary on the superLattice domain
template<typename T, typename DESCRIPTOR>
void setFreeEnergyWallBoundary(SuperLattice<T, DESCRIPTOR>& sLattice, SuperGeometry<T,2>& superGeometry,
                               int material, T alpha, T kappa1, T kappa2, T h1, T h2, int latticeNumber);

///Initialising the Free Energy Wall Boundary on the superLattice domain
template<typename T, typename DESCRIPTOR>
void setFreeEnergyWallBoundary(SuperLattice<T, DESCRIPTOR>& sLattice, FunctorPtr<SuperIndicatorF2D<T>>&& indicator,
                               T alpha, T kappa1, T kappa2, T h1, T h2, int latticeNumber);


/// Implementation of a wetting boundary condition for the ternary free energy model, consisting of a BounceBack
/// dynamics and an FreeEnergyWall PostProcessor.
/// \param[in] alpha_ - Parameter related to the interface width. [lattice units]
/// \param[in] kappa1_ - Parameter related to surface tension. [lattice units]
/// \param[in] kappa2_ - Parameter related to surface tension. [lattice units]
/// \param[in] kappa3_ - Parameter related to surface tension. [lattice units]
/// \param[in] h1_ - Parameter related to resulting contact angle of the boundary. [lattice units]
/// \param[in] h2_ - Parameter related to resulting contact angle of the boundary. [lattice units]
/// \param[in] h3_ - Parameter related to resulting contact angle of the boundary. [lattice units]
/// \param[in] latticeNumber - determines the number of the free energy lattice to set the boundary accordingly

///Initialising the Free Energy Wall Boundary on the superLattice domain
template<typename T, typename DESCRIPTOR>
void setFreeEnergyWallBoundary(SuperLattice<T, DESCRIPTOR>& sLattice, SuperGeometry<T,2>& superGeometry,
                               int material, T alpha, T kappa1, T kappa2, T kappa3, T h1, T h2, T h3, int latticeNumber);

///Initialising the Free Energy Wall Boundary on the superLattice domain
template<typename T, typename DESCRIPTOR>
void setFreeEnergyWallBoundary(SuperLattice<T, DESCRIPTOR>& sLattice, FunctorPtr<SuperIndicatorF2D<T>>&& indicator,
                               T alpha, T kappa1, T kappa2, T kappa3, T h1, T h2, T h3, int latticeNumber);



//set FreeEnergyWallBoundary on block domain.
//This function works for the setFreeEnergyWallBoundary with h1,h2,h3 Parameters and h1,h2 Parameters
template<typename T, typename DESCRIPTOR>
void setFreeEnergyWallBoundary(BlockLattice<T,DESCRIPTOR>& block, BlockIndicatorF2D<T>& indicator,
                               T addend, int latticeNumber, bool includeOuterCells=false);

}//namespace olb

#endif

