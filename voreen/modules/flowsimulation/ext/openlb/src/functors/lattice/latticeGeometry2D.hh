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
#include "utilities/omath.h"
#include <limits>

#include "latticeGeometry2D.h"
#include "dynamics/lbm.h"  // for computation of lattice rho and velocity
#include "geometry/superGeometry.h"
#include "indicator/superIndicatorF2D.h"
#include "blockBaseF2D.h"
#include "functors/genericF.h"
#include "functors/analytical/analyticalF.h"
#include "functors/analytical/indicator/indicatorF2D.h"
#include "communication/mpiManager.h"


namespace olb {

template<typename T,typename DESCRIPTOR>
SuperLatticeGeometry2D<T,DESCRIPTOR>::SuperLatticeGeometry2D(
  SuperLattice<T,DESCRIPTOR>& sLattice, SuperGeometry<T,2>& superGeometry,
  const int material)
  : SuperLatticeGeometry2D<T,DESCRIPTOR>(superGeometry, material)
{ }

template<typename T,typename DESCRIPTOR>
SuperLatticeGeometry2D<T,DESCRIPTOR>::SuperLatticeGeometry2D(
  SuperGeometry<T,2>& superGeometry, const int material)
  : SuperF2D<T>(superGeometry, 1), _superGeometry(superGeometry),
    _material(material)
{
  this->getName() = "geometry";
  const int maxC = superGeometry.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back( new  BlockLatticeGeometry2D<T,DESCRIPTOR>(
                                  this->_superGeometry.getBlockGeometry(iC),
                                  _material) );
  }
}

template <typename T, typename DESCRIPTOR>
BlockLatticeGeometry2D<T,DESCRIPTOR>::BlockLatticeGeometry2D
(BlockLattice<T,DESCRIPTOR>& blockLattice, BlockGeometry<T,2>& blockGeometry, int material)
  : BlockLatticeGeometry2D<T,DESCRIPTOR>(blockGeometry, material)
{ }

template <typename T, typename DESCRIPTOR>
BlockLatticeGeometry2D<T,DESCRIPTOR>::BlockLatticeGeometry2D
(BlockGeometry<T,2>& blockGeometry, int material)
  : BlockF2D<T>(blockGeometry,1), _blockGeometry(blockGeometry), _material(material)
{
  this->getName() = "geometry";
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticeGeometry2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  const int materialTmp = _blockGeometry.getMaterial( {input[0], input[1]} );

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
