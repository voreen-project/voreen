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

#ifndef LATTICE_PHYS_PORE_SIZE_DISTRIBUTION_3D_HH
#define LATTICE_PHYS_PORE_SIZE_DISTRIBUTION_3D_HH

#include<vector>    // for generic i/o
#include<cmath>     // for lpnorm
#include<math.h>

#include "latticePhysPoreSizeDistribution3D.h"
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
SuperLatticePhysPoreSizeDistribution3D<T,DESCRIPTOR>::SuperLatticePhysPoreSizeDistribution3D
(SuperLattice<T,DESCRIPTOR>& sLattice, SuperGeometry<T,3>& superGeometry, int material,
 XMLreader const& xmlReader)
  : SuperLatticeF3D<T,DESCRIPTOR>(sLattice,1),
    _superGeometry(superGeometry)
{
  this->getName() = "physPoreSizeDistribution";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back( new BlockLatticePhysPoreSizeDistribution3D<T,DESCRIPTOR>(this->_sLattice.getBlock(iC), this->_superGeometry.getBlockGeometry(iC), material, xmlReader));
  }
}

template<typename T, typename DESCRIPTOR>
BlockLatticePhysPoreSizeDistribution3D<T, DESCRIPTOR>::BlockLatticePhysPoreSizeDistribution3D(
  BlockLattice<T, DESCRIPTOR>& blockLattice, BlockGeometry<T,3>& blockGeometry, int material, XMLreader const& xmlReader)
  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, 1), _blockGeometry(blockGeometry), _material(material),
    _distanceFunctor(blockLattice, blockGeometry, xmlReader), _distanceCache(_distanceFunctor)
{
  this->getName() = "physPoreSizeDistribution";

  for (XMLreader* child : xmlReader) {
    // clout << "iterator to xml-child: " << child->getName() << std::endl;
    _tmpIndicator = createIndicatorSphere3D<T>(*child);
    _indicatorList.push_back(_tmpIndicator);
  }

}

template<typename T, typename DESCRIPTOR>
bool BlockLatticePhysPoreSizeDistribution3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  if (this->_blockGeometry.get({input[0],input[1],input[2]}) == _material) {
    T localDistance[1] = {0.};
    _distanceCache(localDistance, input);
    // cout << localDistance[0] << std::endl;

    // filter by local maximum (compare to 26 neighbours)
    for (int iPop = 1; iPop < 27; iPop++) {
      int neighbourInput[3] = {0,0,0};
      for (int iDim = 0; iDim < 3; iDim++) {
        neighbourInput[iDim] = input[iDim] + descriptors::c<descriptors::D3Q27<>>(iPop,iDim);
        // cout << neighbourInput[iDim] << " ";
      }
      // cout << iPop << ": ";

      T neighbourDistance[1] = {0.};
      if ( _distanceCache(neighbourDistance, neighbourInput) ) {
        if ( neighbourDistance[0] > localDistance[0] ) {
          // cout << "neighbour larger ";
          // cout << iPop << std::endl;
          output[0] = -1;
          return true;
        }
      }
      else {
        // cout << "distance not calculated ";
        // cout << iPop << std::endl;
        output[0] = -1;
        return true;
      }
    }

    output[0] = localDistance[0];
    return true;
  }
  else {
    output[0] = -1;
    return false;
  }
}

}
#endif
