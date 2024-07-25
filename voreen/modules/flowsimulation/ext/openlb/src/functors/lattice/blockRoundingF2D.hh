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

#ifndef BLOCK_ROUNDING_F_2D_HH
#define BLOCK_ROUNDING_F_2D_HH

#include "utilities/roundingMode.h"
#include "blockRoundingF2D.h"

namespace olb {

template <typename T>
BlockRoundingF2D<T>::BlockRoundingF2D(BlockF2D<T>& f,
                                      RoundingMode rounding)
  :BlockF2D<T>(f.getBlockStructure(),f.getTargetDim()),_f(f),_rounding(rounding)
{
  this->getName() = "Rounding(" + _f.getName() + ")";
}

template <typename T>
bool BlockRoundingF2D<T>::operator()(T output[], const int input[])
{

  _f(output, input);


  switch ( _rounding ) {
  case RoundingMode::NearestInteger:
    for ( int i=0; i < _f.getTargetDim(); i++ ) {
      output[i] = std::nearbyint(output[i]);
    }
    break;
  case RoundingMode::None:
    for ( int i=0; i < _f.getTargetDim(); i++ ) {
      output[i] = output[i];
    }
    break;
  case RoundingMode::Floor:
    for ( int i=0; i < _f.getTargetDim(); i++ ) {
      output[i] = util::floor(output[i]);
    }
    break;
  case RoundingMode::Ceil:
    for ( int i=0; i < _f.getTargetDim(); i++ ) {
      output[i] = util::ceil(output[i]);
    }
    break;
  }

  return true;
}

}
#endif
