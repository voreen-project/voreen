/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2020 Alexander Schulz, Davide Dapelo
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

//This file contains the Interpolated Velocity Boundary
//This is a new version of the Boundary, which only contains free floating functions

#ifndef SET_NO_PENETRATION_BOUNDARY_2D_H
#define SET_NO_PENETRATION_BOUNDARY_2D_H

#include "fdBoundaryPostProcessors2D.h"


namespace olb {
///Initialising the setFdNeumannZeroBoundary function on the superLattice domain
///Interpolated Boundaries use the BGKdynamics collision-operator
template<typename T, typename DESCRIPTOR, typename SCHEME_ADV, typename FIELD=descriptors::AD_FIELD, typename SOURCE=void>
void setFdNeumannZeroBoundary(SuperLattice<T,DESCRIPTOR>& sLattice, std::size_t& iT, std::shared_ptr<FdModel<T,DESCRIPTOR>> model,
                                SuperGeometry<T,2>& superGeometry, int material);

///Initialising the setFdNeumannZeroBoundary function on the superLattice domain
template<typename T, typename DESCRIPTOR, typename SCHEME_ADV, typename FIELD=descriptors::AD_FIELD, typename SOURCE=void>
void setFdNeumannZeroBoundary(SuperLattice<T,DESCRIPTOR>& sLattice, std::size_t& iT, std::shared_ptr<FdModel<T,DESCRIPTOR>> model,
                                FunctorPtr<SuperIndicatorF2D<T>>&& indicator);

/// Set interpolated velocity boundary for any indicated cells inside the block domain
template<typename T, typename DESCRIPTOR, typename SCHEME_ADV, typename FIELD=descriptors::AD_FIELD, typename SOURCE=void>
void setFdNeumannZeroBoundary(BlockLattice<T,DESCRIPTOR>& block, std::size_t& iT, std::shared_ptr<FdModel<T,DESCRIPTOR>> model,
                                BlockIndicatorF2D<T>& indicator, bool includeOuterCells=false);

}//namespace olb
#endif
