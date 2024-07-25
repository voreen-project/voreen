/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2019 Tom Braun
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

#ifndef BLOCK_DISCRETIZATION_F_2D_HH
#define BLOCK_DISCRETIZATION_F_2D_HH

#include "blockDiscretizationF2D.h"
#include <algorithm>

namespace olb {

template <typename T>
BlockDiscretizationF2D<T>::BlockDiscretizationF2D(BlockF2D<T>& f,
    T bottomBoundary,
    T topBoundary,
    int n)
  :BlockF2D<T>(f.getBlockStructure(),f.getTargetDim()),_f(f),_bottomBoundary(bottomBoundary),_topBoundary(topBoundary),_n(n)
{
  this->getName() = "Discretization(" + _f.getName() +  " ) ";
}

template <typename T>
bool BlockDiscretizationF2D<T>::operator()(T output[], const int input[])
{
  const T delta = (_topBoundary - _bottomBoundary) / (util::max(_n-1,1));
  _f(output, input);

  for (int i=0; i < _f.getTargetDim(); ++i) {
    if ( output[i] <= _bottomBoundary ) {
      output[i] = _bottomBoundary;
    }
    else if ( output[i] >= _topBoundary ) {
      output[i] = _topBoundary;
    }
    else if ( _n != 1 ) {
      output[i] -= _bottomBoundary;
      output[i]  = _bottomBoundary + std::nearbyint(output[i] / delta) * delta;
    }
  }

  return true;
}

}
#endif
