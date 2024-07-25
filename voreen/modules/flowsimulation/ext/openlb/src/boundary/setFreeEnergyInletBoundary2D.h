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

//This file contains the Free Energy Inlet Boundary
//This is a new version of the Boundary, which only contains free floating functions
#ifndef SET_FREE_ENERGY_INLET_BOUNDARY_2D_H
#define SET_FREE_ENERGY_INLET_BOUNDARY_2D_H

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

///Implementation of a inlet boundary condition for the partner lattices of the binary or ternary free energy model.
///FreeEnergyInletBoundary is local boundary Condition --> MixinDynamics = RLBdynamics
///Initialising the Free Energy Inlet Boundary on the superLattice domain
template<typename T, typename DESCRIPTOR, typename MixinDynamics = RLBdynamics<T, DESCRIPTOR>>
void setFreeEnergyInletBoundary(SuperLattice<T, DESCRIPTOR>& sLattice, T omega, SuperGeometry<T,2>& superGeometry, int material, 
                                std::string type, int latticeNumber);

///Initialising the Free Energy Inlet Boundary on the superLattice domain
template<typename T, typename DESCRIPTOR, typename MixinDynamics = RLBdynamics<T,DESCRIPTOR>>
void setFreeEnergyInletBoundary(SuperLattice<T, DESCRIPTOR>& sLattice, T omega, FunctorPtr<SuperIndicatorF2D<T>>&& indicator, 
                                std::string type, int latticeNumber);

///Set FreeEnergyInlet boundary for any indicated cells inside the block domain
template<typename T, typename DESCRIPTOR, typename MixinDynamics>
void setFreeEnergyInletBoundary(BlockLattice<T,DESCRIPTOR>& block,T omega, BlockIndicatorF2D<T>& indicator, std::string type,
                                int latticeNumber, bool includeOuterCells=false);

}//namespace olb
#endif
