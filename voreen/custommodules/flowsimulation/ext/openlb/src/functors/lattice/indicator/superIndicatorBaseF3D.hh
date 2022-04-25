/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2016 Benjamin Förster
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

#ifndef SUPER_INDICATOR_BASE_F_3D_HH
#define SUPER_INDICATOR_BASE_F_3D_HH

#include "superIndicatorBaseF3D.h"
#include "geometry/superGeometry.h"

namespace olb {


template <typename T>
SuperIndicatorF3D<T>::SuperIndicatorF3D(SuperGeometry<T,3>& geometry)
  : SuperF3D<T, bool>(geometry, 1),
    _superGeometry(geometry)
{ }

template <typename T>
BlockIndicatorF3D<T>& SuperIndicatorF3D<T>::getBlockIndicatorF(int iCloc)
{
  OLB_ASSERT(iCloc < int(this->_blockF.size()) && iCloc >= 0,
             "block functor index within bounds");
  // Note: The type system doesn't guarantee this operation to be valid
  //       as blockF may contain any implementation of the BlockF3D interface.
  return *static_cast<BlockIndicatorF3D<T>*>(this->_blockF[iCloc].get());
}

template <typename T>
SuperGeometry<T,3>& SuperIndicatorF3D<T>::getSuperGeometry()
{
  return this->_superGeometry;
}

template <typename T>
bool SuperIndicatorF3D<T>::operator() (const int input[])
{
  bool output{};
  if (_cachedData) {
    output = _cachedData->getBlock(input[0]).get(input+1);
  }
  else {
    this->operator()(&output, input);
  }
  return output;
}

template <typename T>
bool SuperIndicatorF3D<T>::operator() (int iC, int iX, int iY, int iZ)
{
  bool output{};
  if (_cachedData) {
    output = _cachedData->getBlock(iC).get({iX, iY, iZ});
  }
  else {
    this->operator()(&output, iC, iX, iY, iZ);
  }
  return output;
}

template <typename T>
void SuperIndicatorF3D<T>::cache()
{
  _cachedData = std::unique_ptr<SuperData<3,T,bool>>(
                  new SuperData<3,T,bool>(*this));
  for (unsigned iC = 0; iC < this->_blockF.size(); ++iC) {
    getBlockIndicatorF(iC).setCache(_cachedData->getBlock(iC));
  }
}


} // namespace olb

#endif
