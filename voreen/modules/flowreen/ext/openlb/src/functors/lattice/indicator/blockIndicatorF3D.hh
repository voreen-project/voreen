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

#ifndef BLOCK_INDICATOR_F_3D_HH
#define BLOCK_INDICATOR_F_3D_HH

#include <algorithm>

#include "blockIndicatorF3D.h"
#include "core/util.h"

namespace olb {

template <typename T>
BlockIndicatorFfromIndicatorF3D<T>::BlockIndicatorFfromIndicatorF3D(
  IndicatorF3D<T>&    indicatorF,
  BlockGeometry3D<T>& blockGeometry,
  Cuboid3D<T>&        cuboid)
  : BlockIndicatorF3D<T>(blockGeometry),
    _indicatorF(indicatorF),
    _cuboid(cuboid)
{};

template <typename T>
bool BlockIndicatorFfromIndicatorF3D<T>::operator() (bool output[], const int input[])
{
  T physR[3] = {};
  _cuboid.getPhysR(physR,input);
  return _indicatorF(output,physR);
}


template <typename T>
BlockIndicatorMaterial3D<T>::BlockIndicatorMaterial3D(
  BlockGeometry3D<T>& blockGeometry, int overlap, const std::vector<int>* materials)
  : BlockIndicatorF3D<T>(blockGeometry),
    _overlap(overlap),
    _materialNumbers(materials)
{};

template <typename T>
bool BlockIndicatorMaterial3D<T>::operator() (bool output[], const int input[])
{
  const int blockInput[3] = {
    input[0] + _overlap,
    input[1] + _overlap,
    input[2] + _overlap,
  };

  OLB_PRECONDITION(blockInput[0] < this->_blockGeometry.getNx());
  OLB_PRECONDITION(blockInput[1] < this->_blockGeometry.getNy());
  OLB_PRECONDITION(blockInput[2] < this->_blockGeometry.getNz());

  output[0] = std::find(_materialNumbers->cbegin(),
                        _materialNumbers->cend(),
                        // read material number explicitly using the const version
                        // of BlockGeometry3D<T>::get to avoid resetting geometry
                        // statistics:
                        const_cast<const BlockGeometry3D<T>&>(this->_blockGeometry).get(blockInput[0], blockInput[1], blockInput[2]))
              != _materialNumbers->cend();

  return output[0];
}


} // namespace olb

#endif
