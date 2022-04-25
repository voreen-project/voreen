/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2020 Clara Schragmann
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

#ifndef LATTICE_DISCRETE_NORMAL_3D_HH
#define LATTICE_DISCRETE_NORMAL_3D_HH

#include<vector>    // for generic i/o

#include "latticeDiscreteNormal3D.h"
#include "indicator/superIndicatorF3D.h"
#include "geometry/superGeometry.h"


namespace olb {

template<typename T, typename DESCRIPTOR>
SuperLatticeDiscreteNormal3D<T, DESCRIPTOR>::SuperLatticeDiscreteNormal3D(SuperLattice<T, DESCRIPTOR>& sLattice,  SuperGeometry<T,3>& superGeometry, FunctorPtr<SuperIndicatorF3D<T>>&& indicatorF)
  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, DESCRIPTOR::d), _superGeometry(superGeometry),
    _indicatorF(std::move(indicatorF))
{
  this->getName() = "discreteNormal";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new  BlockLatticeDiscreteNormal3D<T, DESCRIPTOR>(this->_sLattice.getBlock(iC),
                               this->_superGeometry.getBlockGeometry(iC),
                               indicatorF->getBlockIndicatorF(iC) ));
  }
}

template<typename T, typename DESCRIPTOR>
BlockLatticeDiscreteNormal3D<T, DESCRIPTOR>::BlockLatticeDiscreteNormal3D(BlockLattice<T, DESCRIPTOR>& blockLattice, BlockGeometry<T,3>& blockGeometry, BlockIndicatorF3D<T>& indicatorF)
  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, DESCRIPTOR::d),
    _blockGeometry(blockGeometry),
    _indicatorF(indicatorF)
{
  this->getName() = "discreteNormal";
}

template<typename T, typename DESCRIPTOR>
bool BlockLatticeDiscreteNormal3D<T, DESCRIPTOR>::operator()(T output[], const int input[])
{
  if (_indicatorF(input)) {
    std::vector<int> normalVector = _blockGeometry.getStatistics().getType(input);
    output[0] = normalVector[1];
    output[1] = normalVector[2];
    output[2] = normalVector[3];
  }
  else {
    output[0] = 0.;
    output[1] = 0.;
    output[2] = 0.;
  }

  return true;
}

template<typename T, typename DESCRIPTOR>
SuperLatticeDiscreteNormalType3D<T, DESCRIPTOR>::SuperLatticeDiscreteNormalType3D(SuperLattice<T, DESCRIPTOR>& sLattice,  SuperGeometry<T,3>& superGeometry, FunctorPtr<SuperIndicatorF3D<T>>&& indicatorF)
  : SuperLatticeF3D<T, DESCRIPTOR>(sLattice, DESCRIPTOR::d), _superGeometry(superGeometry),
    _indicatorF(std::move(indicatorF))
{
  this->getName() = "discreteNormalType";
  int maxC = this->_sLattice.getLoadBalancer().size();
  this->_blockF.reserve(maxC);
  for (int iC = 0; iC < maxC; iC++) {
    this->_blockF.emplace_back(new  BlockLatticeDiscreteNormalType3D<T, DESCRIPTOR>(this->_sLattice.getBlock(iC),
                               this->_superGeometry.getBlockGeometry(iC),
                               indicatorF->getBlockIndicatorF(iC) ));
  }
}

template<typename T, typename DESCRIPTOR>
BlockLatticeDiscreteNormalType3D<T, DESCRIPTOR>::BlockLatticeDiscreteNormalType3D(BlockLattice<T, DESCRIPTOR>& blockLattice, BlockGeometry<T,3>& blockGeometry, BlockIndicatorF3D<T>& indicatorF)
  : BlockLatticeF3D<T, DESCRIPTOR>(blockLattice, DESCRIPTOR::d),
    _blockGeometry(blockGeometry),
    _indicatorF(indicatorF)
{
  this->getName() = "discreteNormalType";
}

template<typename T, typename DESCRIPTOR>
bool BlockLatticeDiscreteNormalType3D<T, DESCRIPTOR>::operator()(T output[1], const int input[])
{
  if (_indicatorF(input)) {
    std::vector<int> normalVector = _blockGeometry.getStatistics().getType(input);
    output[0] = normalVector[0];
  }
  else {
    output[0] = 0.;
  }

  return true;
}

}
#endif
