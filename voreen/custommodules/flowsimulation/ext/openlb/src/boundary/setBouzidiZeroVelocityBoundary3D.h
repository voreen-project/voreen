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

//This file contains the BouzidiZeroVelocityBoundary
//This is an offLattice Boundary
//This is a new version of the Boundary, which only contains free floating functions
#ifndef SET_BOUZIDI_ZERO_VELOCITY_BOUNDARY_H
#define SET_BOUZIDI_ZERO_VELOCITY_BOUNDARY_H

#include <vector>
#include "utilities/functorPtr.h"
#include "extendedFiniteDifferenceBoundary3D.h"
#include "geometry/superGeometry.h"
#include "extendedFiniteDifferenceBoundary3D.h"

#include "functors/lattice/indicator/superIndicatorBaseF3D.h"
#include "dynamics/dynamics.h"
#include "functors/lattice/indicator/blockIndicatorF3D.h"
#include "io/ostreamManager.h"
#include "functors/lattice/indicator/blockIndicatorF3D.h"
#include "dynamics/freeEnergyDynamics.h"
#include "offBoundaryPostProcessors3D.h"


namespace olb {
/// Set offDynamics with boundary links and post processors using indicators
/**
 * Add offDynamics with initialisation of boundary links and the corresponding
 * post processors
 * Note: Uses information of the second neighbours of the cell (x,y,z)
 * Add post processors. Ensure that offDynamics are defined!
 *
 * \param boundaryIndicator Indicator describing boundary cells
 * \param bulkIndicator     Indicator describing bulk cells
 * \param geometryIndicator Indicator describing the geometry to be bounded
 **/

////////// SuperLattice Domain  /////////////////////////////////////////

///Initialising the setLocalVelocityBoundary function on the superLattice domain
template<typename T, typename DESCRIPTOR>
void setBouzidiZeroVelocityBoundary( SuperLattice<T, DESCRIPTOR>& sLattice,SuperGeometry<T,3>& superGeometry, int material,IndicatorF3D<T>& geometryIndicator,
                                     std::vector<int> bulkMaterials = std::vector<int>(1,1));

template<typename T, typename DESCRIPTOR>
void setBouzidiZeroVelocityBoundary(SuperLattice<T, DESCRIPTOR>& sLattice, FunctorPtr<SuperIndicatorF3D<T>>&& boundaryIndicator, IndicatorF3D<T>& geometryIndicator,
                                    std::vector<int> bulkMaterials);
template<typename T, typename DESCRIPTOR>
void setBouzidiZeroVelocityBoundary(SuperLattice<T, DESCRIPTOR>& sLattice, FunctorPtr<SuperIndicatorF3D<T>>&& boundaryIndicator,
                                    FunctorPtr<SuperIndicatorF3D<T>>&& bulkIndicator,
                                    IndicatorF3D<T>& geometryIndicator);


////////// BlockLattice Domain  /////////////////////////////////////////


template<typename T, typename DESCRIPTOR>
void setBouzidiZeroVelocityBoundary(BlockLattice<T, DESCRIPTOR>& block, BlockIndicatorF3D<T>& boundaryIndicator, BlockIndicatorF3D<T>& bulkIndicator, IndicatorF3D<T>& geometryIndicator, T _epsFraction);

template<typename T, typename DESCRIPTOR>
void setBouzidiZeroVelocityBoundary1(BlockLattice<T, DESCRIPTOR>& block, BlockGeometry<T,3>& blockGeometryStructure, int iX, int iY, int iZ, IndicatorF3D<T>& geometryIndicator, BlockIndicatorF3D<T>& bulkIndicator, T _epsFraction);

template<typename T, typename DESCRIPTOR>
void setBouzidiZeroVelocityBoundary(BlockLattice<T, DESCRIPTOR>& block, BlockGeometry<T,3>& blockGeometryStructure, int x, int y, int z, T distances[DESCRIPTOR::q]);

template<typename T, typename DESCRIPTOR>
void setBouzidiZeroVelocityBoundary(BlockLattice<T, DESCRIPTOR>& block, BlockGeometry<T,3>& blockGeometryStructure, int x, int y, int z, int iPop, T dist);
} //namespace olb


#endif
