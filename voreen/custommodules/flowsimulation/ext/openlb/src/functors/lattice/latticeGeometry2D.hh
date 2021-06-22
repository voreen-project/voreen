/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2019 Albert Mink, Mathias J. Krause, Lukas Baron
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

#ifndef LATTICE_GEOMETRY_2D_HH
#define LATTICE_GEOMETRY_2D_HH

#include <vector>
#include <cmath>
#include <limits>

#include "latticeGeometry2D.h"
#include "dynamics/lbHelpers.h"  // for computation of lattice rho and velocity
#include "geometry/superGeometry2D.h"
#include "indicator/superIndicatorF2D.h"
#include "blockBaseF2D.h"
#include "functors/genericF.h"
#include "functors/analytical/analyticalF.h"
#include "functors/analytical/indicator/indicatorF2D.h"
#include "core/blockLattice2D.h"
#include "communication/mpiManager.h"
#include "core/blockLatticeStructure2D.h"


namespace olb {

template<typename T,typename DESCRIPTOR>
SuperLatticeGeometry2D<T,DESCRIPTOR>::SuperLatticeGeometry2D(
  SuperLattice2D<T,DESCRIPTOR>& sLattice, SuperGeometry2D<T>& superGeometry,
  const int material)
  : SuperLatticeF2D<T,DESCRIPTOR>(sLattice, 1), _superGeometry(superGeometry),
    _material(material)
{
  this->getName() = "geometry";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back( new  BlockLatticeGeometry2D<T,DESCRIPTOR>(
                                  this->_sLattice.getBlockLattice(iC),
                                  this->_superGeometry.getBlockGeometry(iC),
                                  _material) );
  }
}

template <typename T, typename DESCRIPTOR>
BlockLatticeGeometry2D<T,DESCRIPTOR>::BlockLatticeGeometry2D
(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice, BlockGeometryStructure2D<T>& blockGeometry, int material)
  : BlockLatticeF2D<T,DESCRIPTOR>(blockLattice,1), _blockGeometry(blockGeometry), _material(material)
{
  this->getName() = "geometry";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticeGeometry2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  int materialTmp = _blockGeometry.getMaterial( input[0], input[1] );

  if (_material != -1) {
    if (_material == materialTmp) {
      output[0] = T(1);
      return true;
    }
    else {
      output[0] = T();
      return true;
    }
  }
  output[0]=T(materialTmp);
  return false;
}

}
#endif
