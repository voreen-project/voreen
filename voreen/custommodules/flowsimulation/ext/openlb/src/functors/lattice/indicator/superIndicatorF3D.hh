/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2016 Benjamin Förster, Adrian Kummerländer
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

#ifndef SUPER_INDICATOR_F_3D_HH
#define SUPER_INDICATOR_F_3D_HH

#include <numeric>

#include "superIndicatorF3D.h"
#include "blockIndicatorF3D.h"
#include "core/util.h"

namespace olb {

template <typename T>
SuperIndicatorFfromIndicatorF3D<T>::SuperIndicatorFfromIndicatorF3D(
  IndicatorF3D<T>& indicatorF, SuperGeometry3D<T>& geometry)
  : SuperIndicatorF3D<T>(geometry),
    _indicatorF(indicatorF)
{
  this->getName() = "SuperIndicator_from_" + _indicatorF.getName();

  LoadBalancer<T>&     load   = this->getSuperStructure().getLoadBalancer();
  CuboidGeometry3D<T>& cuboid = this->getSuperStructure().getCuboidGeometry();

  for (int iC = 0; iC < load.size(); ++iC) {
    this->_blockF.emplace_back(
      new BlockIndicatorFfromIndicatorF3D<T>(
        _indicatorF, geometry.getExtendedBlockGeometry(iC), cuboid.get(load.glob(iC)))
    );
  }
}

template <typename T>
bool SuperIndicatorFfromIndicatorF3D<T>::operator() (bool output[], const int input[])
{
  T physR[3];
  this->_superStructure.getCuboidGeometry().getPhysR(physR, input);
  _indicatorF(output, physR);
  return true;
}


template <typename T>
SuperIndicatorMaterial3D<T>::SuperIndicatorMaterial3D(
  SuperGeometry3D<T>& geometry, std::vector<int> materials)
  : SuperIndicatorF3D<T>(geometry),
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

  for (int iC = 0; iC < this->_superGeometry.getLoadBalancer().size(); ++iC) {
    this->_blockF.emplace_back(
      new BlockIndicatorMaterial3D<T>(this->_superGeometry.getExtendedBlockGeometry(iC),
                                      this->_superGeometry.getOverlap(),
                                      &_materialNumbers)
    );
  }
}

template <typename T>
SuperIndicatorMaterial3D<T>::SuperIndicatorMaterial3D(
  SuperGeometry3D<T>& geometry, std::list<int> materials)
  : SuperIndicatorMaterial3D(geometry,
                             std::vector<int>(materials.begin(), materials.end())) { }

template <typename T>
bool SuperIndicatorMaterial3D<T>::operator() (bool output[], const int input[])
{
  output[0] = false;

  LoadBalancer<T>& load = this->_superGeometry.getLoadBalancer();

  if (!this->_blockF.empty() && load.isLocal(input[0])) {
    // query material number of appropriate block indicator
    this->getBlockF(load.loc(input[0]))(output,&input[1]);
  } else {
    // Try to query material number locally as a fallback if no block indicators
    // were instantiated. This will terminate if data is unavailable.
    output[0] = std::find(_materialNumbers.cbegin(),
                          _materialNumbers.cend(),
                          this->_superGeometry.get(input))
                != _materialNumbers.cend();
  }

  return true;
}


} // namespace olb

#endif
