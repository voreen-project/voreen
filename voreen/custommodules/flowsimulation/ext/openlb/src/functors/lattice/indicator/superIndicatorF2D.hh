/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2017 Adrian Kummerländer
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

#ifndef SUPER_INDICATOR_F_2D_HH
#define SUPER_INDICATOR_F_2D_HH

#include <numeric>
#include <algorithm>

#include "core/util.h"
#include "superIndicatorF2D.h"
#include "blockIndicatorF2D.h"

namespace olb {


template <typename T>
SuperIndicatorFfromIndicatorF2D<T>::SuperIndicatorFfromIndicatorF2D(
  IndicatorF2D<T>& indicatorF, SuperGeometry2D<T>& geometry)
  : SuperIndicatorF2D<T>(geometry),
    _indicatorF(indicatorF)
{
  this->getName() = "SuperIndicator_from_" + _indicatorF.getName();

  LoadBalancer<T>&     load   = this->getSuperStructure().getLoadBalancer();
  CuboidGeometry2D<T>& cuboid = this->getSuperStructure().getCuboidGeometry();

  for (int iC = 0; iC < load.size(); ++iC) {
    this->_blockF.emplace_back(
      new BlockIndicatorFfromIndicatorF2D<T>(
        _indicatorF, geometry.getExtendedBlockGeometry(iC), cuboid.get(load.glob(iC)))
    );
  }
}

template <typename T>
bool SuperIndicatorFfromIndicatorF2D<T>::operator() (bool output[], const int input[])
{
  T physR[2];
  this->_superStructure.getCuboidGeometry().getPhysR(physR, input);
  return _indicatorF(output, physR);
}


template <typename T>
SuperIndicatorMaterial2D<T>::SuperIndicatorMaterial2D(
  SuperGeometry2D<T>& geometry, std::vector<int> materials)
  : SuperIndicatorF2D<T>(geometry),
    _superGeometry(geometry),
    _materialNumbers(materials)
{
  const std::string matString = std::accumulate(
                                  materials.begin()+1,
                                  materials.end(),
                                  std::to_string(materials[0]),
  [](const std::string& a, int b) {
    return a + '_' + std::to_string(b);
  });
  this->getName() = "SuperIndicator_on_Material_" + matString;

  for (int iC = 0; iC < _superGeometry.getLoadBalancer().size(); ++iC) {
    this->_blockF.emplace_back(
      new BlockIndicatorMaterial2D<T>(_superGeometry.getExtendedBlockGeometry(iC),
                                      _superGeometry.getOverlap(),
                                      &_materialNumbers)
    );
  }
}

template <typename T>
bool SuperIndicatorMaterial2D<T>::operator() (bool output[], const int input[])
{
  output[0] = false;

  LoadBalancer<T>& load = _superGeometry.getLoadBalancer();

  if (!this->_blockF.empty() && load.isLocal(input[0])) {
    this->getBlockF(load.loc(input[0]))(output, &input[1]);
  } else {
    // query material number locally instead of performing unnecessary synchronization
    output[0] = std::find(_materialNumbers.cbegin(),
                          _materialNumbers.cend(),
                          _superGeometry.get(input[0], input[1], input[2]))
                != _materialNumbers.cend();
  }

  return output[0];
}


} // namespace olb

#endif
