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
#ifndef DEFINE_U_3D_H
#define DEFINE_U_3D_H

#include <vector>
#include "utilities/functorPtr.h"
#include "geometry/superGeometry.h"

#include "functors/lattice/indicator/superIndicatorBaseF3D.h"
#include "functors/lattice/indicator/blockIndicatorF3D.h"
#include "dynamics/dynamics.h"
#include "io/ostreamManager.h"


//defineU for offLatticeBoundaryConditions

namespace olb {
////////// SuperLattice Domain  /////////////////////////////////////////
template<typename T, typename DESCRIPTOR>
void defineUBouzidi(SuperLattice<T,DESCRIPTOR>& sLattice, SuperGeometry<T,3>& superGeometry, int material,
                    AnalyticalF3D<T,T>& u, std::vector<int> bulkMaterials = std::vector<int>(1,1) );


template<typename T, typename DESCRIPTOR>
void defineUBouzidi(SuperLattice<T,DESCRIPTOR>& sLattice, FunctorPtr<SuperIndicatorF3D<T>>&& boundaryIndicator,
                    AnalyticalF3D<T,T>& u, std::vector<int> bulkMaterials = std::vector<int>(1,1) );

template<typename T, typename DESCRIPTOR>
void defineUBouzidi(SuperLattice<T,DESCRIPTOR>& sLattice, FunctorPtr<SuperIndicatorF3D<T>>&& boundaryIndicator,
                    FunctorPtr<SuperIndicatorF3D<T>>&& bulkIndicator, AnalyticalF3D<T,T>& u);

////////// BlockLattice Domain  /////////////////////////////////////////
template<typename T, typename DESCRIPTOR>
void defineUBouzidi(BlockLattice<T,DESCRIPTOR>& _block, BlockIndicatorF3D<T>& indicator, BlockIndicatorF3D<T>& bulkIndicator, AnalyticalF3D<T,T>& u);


template<typename T, typename DESCRIPTOR>
void defineUBouzidi(BlockLattice<T,DESCRIPTOR>& _block, int iX, int iY, int iZ, int iPop, const T u[DESCRIPTOR::d]);


template<typename T, typename DESCRIPTOR>
bool getBoundaryIntersection(BlockLattice<T,DESCRIPTOR>& _block, int iX, int iY, int iZ, int iPop, T point[DESCRIPTOR::d]);


template<typename T, typename DESCRIPTOR>
void setBoundaryIntersection(BlockLattice<T,DESCRIPTOR>& _block, int iX, int iY, int iZ, int iPop, T distance);



}//namespace olb


#endif
