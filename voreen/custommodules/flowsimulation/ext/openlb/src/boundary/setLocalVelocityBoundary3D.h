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

///This file contains the Local Velocity Boundary
///This is a new version of the Boundary, which only contains free floating functions
#ifndef SET_LOCAL_VELOCITY_BOUNDARY_H
#define SET_LOCAL_VELOCITY_BOUNDARY_H

#include <vector>
#include "utilities/functorPtr.h"
#include "extendedFiniteDifferenceBoundary3D.h"
#include "geometry/superGeometry3D.h"
#include "extendedFiniteDifferenceBoundary3D.h"
#include "core/superLattice3D.h"
#include "functors/lattice/indicator/superIndicatorBaseF3D.h"
#include "dynamics/dynamics.h"
#include "functors/lattice/indicator/blockIndicatorF3D.h"
#include "momentaOnBoundaries3D.h"
#include "io/ostreamManager.h"
#include "functors/lattice/indicator/blockIndicatorF3D.h"
#include "dynamics/freeEnergyDynamics.h"
#include "boundaryPostProcessors3D.h"

namespace olb {
////////// SuperLattice Domain  /////////////////////////////////////////

///Initialising the setLocalVelocityBoundary function on the superLattice domain
//Local Boundaries use the RLBdynamics collision operator
template<typename T, typename DESCRIPTOR, typename MixinDynamics=RLBdynamics<T,DESCRIPTOR>>
void setLocalVelocityBoundary(SuperLattice3D<T, DESCRIPTOR>& sLattice, T omega, SuperGeometry3D<T>& superGeometry, int material);

///Initialising the setLocalVelocityBoundary function on the superLattice domain
template<typename T, typename DESCRIPTOR, typename MixinDynamics=RLBdynamics<T,DESCRIPTOR>>
void setLocalVelocityBoundary(SuperLattice3D<T, DESCRIPTOR>& sLattice, T omega, FunctorPtr<SuperIndicatorF3D<T>>&& indicator);

/// Adds needed Cells to the Communicator _commBC in SuperLattice
template<typename T, typename DESCRIPTOR>
void addPoints2CommBC(SuperLattice3D<T, DESCRIPTOR>& sLattice, FunctorPtr<SuperIndicatorF3D<T>>&& indicator, int _overlap);

////////// BlockLattice Domain  /////////////////////////////////////////

/// Set local velocity boundary for any indicated cells inside the block domain
template<typename T, typename DESCRIPTOR, typename MixinDynamics=RLBdynamics<T,DESCRIPTOR>>
void setLocalVelocityBoundary(BlockLatticeStructure3D<T,DESCRIPTOR>& _block, T omega, BlockIndicatorF3D<T>& indicator, bool includeOuterCells);

///sets boundary on indicated cells. This is a function, which can be used on many boundaries. It is not specific to the local velocity boundary
template<int direction, int orientation,typename T, typename DESCRIPTOR, typename MixinDynamics>
void setBoundary(BlockLatticeStructure3D<T,DESCRIPTOR>& _block, T omega, int iX,int iY,int iZ,
                 Momenta<T,DESCRIPTOR>* momenta,Dynamics<T,DESCRIPTOR>* dynamics,PostProcessorGenerator3D<T,DESCRIPTOR>* postProcessor);

}//namespace olb


#endif

