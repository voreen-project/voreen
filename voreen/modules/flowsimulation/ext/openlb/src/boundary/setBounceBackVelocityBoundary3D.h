/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2020 Max Gaedtke
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
#ifndef SET_BOUNCE_BACK_VELOCITY_BOUNDARY_3D_H
#define SET_BOUNCE_BACK_VELOCITY_BOUNDARY_3D_H

#include <vector>
#include "utilities/functorPtr.h"
#include "extendedFiniteDifferenceBoundary3D.h"
#include "geometry/superGeometry.h"
#include "extendedFiniteDifferenceBoundary3D.h"
#include "core/superLattice.h"
#include "functors/lattice/indicator/superIndicatorBaseF3D.h"
#include "dynamics/dynamics.h"
#include "functors/lattice/indicator/blockIndicatorF3D.h"
#include "io/ostreamManager.h"
#include "functors/lattice/indicator/blockIndicatorF3D.h"
#include "dynamics/freeEnergyDynamics.h"
#include "setBoundary3D.h"

namespace olb {

///Initialising the setLocalVelocityBoundary function on the superLattice domain
template<typename T, typename DESCRIPTOR>
void setBounceBackVelocityBoundary(SuperGeometry<T,3>& superGeometry, int
                                   material, T omega, SuperLattice<T, DESCRIPTOR>& sLattice);

template<typename T, typename DESCRIPTOR>
void setBounceBackVelocityBoundary(FunctorPtr<SuperIndicatorF3D<T>>&& indicator, T omega,SuperLattice<T, DESCRIPTOR>& sLattice);


///Set BounceBackVelocityBoundary for any indicated cells inside the block domain
template<typename T, typename DESCRIPTOR>
void setBounceBackVelocityBoundary(BlockIndicatorF3D<T>& indicator, T omega, bool includeOuterCells,BlockLattice<T,DESCRIPTOR>& _block);

}//namespace olb

#endif
