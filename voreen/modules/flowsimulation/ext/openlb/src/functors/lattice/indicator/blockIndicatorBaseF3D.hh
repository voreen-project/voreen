/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2017 Adrian Kummerlaender
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

#ifndef BLOCK_INDICATOR_BASE_F_3D_HH
#define BLOCK_INDICATOR_BASE_F_3D_HH

#include "blockIndicatorBaseF3D.h"
#include "geometry/blockGeometry.h"

namespace olb {

template <typename T>
BlockIndicatorF3D<T>::BlockIndicatorF3D(BlockGeometry<T,3>& geometry)
  : BlockF3D<bool>(geometry, 1),
    _block(geometry),
    _cachedData{nullptr}
{ }

template <typename T>
BlockGeometry<T,3>& BlockIndicatorF3D<T>::getBlockGeometry()
{
  return _block;
}

template <typename T>
bool BlockIndicatorF3D<T>::operator() (const int input[])
{
  bool output{};
  if (_cachedData == nullptr) {
    this->operator()(&output, input);
  }
  else {
    _cachedData->get(input);
  }
  return output;
}

template <typename T>
bool BlockIndicatorF3D<T>::operator() (int iX, int iY, int iZ)
{
  int latticeR[3] { iX, iY, iZ };
  return this->operator()(latticeR);
}

template <typename T>
bool BlockIndicatorF3D<T>::operator() (LatticeR<3> loc)
{
  return operator()(loc[0], loc[1], loc[2]);
}

template <typename T>
void BlockIndicatorF3D<T>::setCache(const BlockData<3,T,bool>& cache)
{
  _cachedData = &cache;
}

template <typename T>
bool BlockIndicatorF3D<T>::isEmpty()
{
  // There is no way to determine domain emptyness in a fashion that is both
  // generic and efficient.
  return false;
}


} // namespace olb

#endif
