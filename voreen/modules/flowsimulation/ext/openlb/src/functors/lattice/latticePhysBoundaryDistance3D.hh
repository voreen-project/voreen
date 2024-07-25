/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Lukas Baron, Tim Dornieden, Mathias J. Krause,
 *  Albert Mink
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

#ifndef LATTICE_PHYS_BOUNDARY_DISTANCE_3D_HH
#define LATTICE_PHYS_BOUNDARY_DISTANCE_3D_HH

#include<vector>    // for generic i/o
#include<cmath>     // for lpnorm
#include<math.h>

#include "latticePhysBoundaryDistance3D.h"
#include "superBaseF3D.h"
#include "functors/analytical/indicator/indicatorBaseF3D.h"
#include "indicator/superIndicatorF3D.h"
#include "dynamics/lbm.h"  // for computation of lattice rho and velocity
#include "geometry/superGeometry.h"
#include "blockBaseF3D.h"
#include "communication/mpiManager.h"
#include "utilities/vectorHelpers.h"

namespace olb {

template <typename T, typename DESCRIPTOR>
SuperLatticePhysBoundaryDistance3D<T,DESCRIPTOR>::SuperLatticePhysBoundaryDistance3D
(SuperLattice<T,DESCRIPTOR>& sLattice, SuperGeometry<T,3>& superGeometry,
 XMLreader const& xmlReader)
  : SuperLatticeF3D<T,DESCRIPTOR>(sLattice,1),
    _superGeometry(superGeometry)
{
  this->getName() = "physBoundaryDistance";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back( new BlockLatticePhysBoundaryDistance3D<T,DESCRIPTOR>(this->_sLattice.getBlock(iC), this->_superGeometry.getBlockGeometry(iC), xmlReader));
  }
}

template<typename T, typename DESCRIPTOR>
BlockLatticePhysBoundaryDistance3D<T, DESCRIPTOR>::BlockLatticePhysBoundaryDistance3D(
  BlockLattice<T, DESCRIPTOR>& blockLattice, BlockGeometry<T,3>& blockGeometry, XMLreader const& xmlReader)
  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, 1), _blockGeometry(blockGeometry)
{
  this->getName() = "physBoundaryDistance";

  for (XMLreader* child : xmlReader) {
    // clout << "iterator to xml-child: " << child->getName() << std::endl;
    _tmpIndicator = createIndicatorSphere3D<T>(*child);
    _indicatorList.push_back(_tmpIndicator);
  }
}

template<typename T, typename DESCRIPTOR>
bool BlockLatticePhysBoundaryDistance3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  T minDistance = std::numeric_limits<T>::util::max();
  T origin[3];
  _blockGeometry.getPhysR(origin, input);
  for (auto &indicator : _indicatorList) {
    bool inside[1] = {false};
    (*indicator)(inside, origin);
    if (inside[0]) {
      minDistance = -1;
      break;
    }
    T distance = 0;
    indicator->distance(distance, origin);
    // clout << "sphere distance = " << distance << std::endl;
    if ( distance < minDistance ) {
      minDistance = distance;
    }
  }
  // clout << "min distance = " << minDistance << std::endl;

  output[0] = minDistance;
  return true;
}

}
#endif
