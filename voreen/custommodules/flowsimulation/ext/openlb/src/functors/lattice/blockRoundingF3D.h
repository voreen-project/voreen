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

#ifndef BLOCK_ROUNDING_F_3D_H
#define BLOCK_ROUNDING_F_3D_H

#include "blockBaseF3D.h"
#include "utilities/roundingMode.h"

namespace olb {

/**
 * Block functor for rounding the value in a certain mode:
 * None :=           No rounding
 * NearestInteger := rounding to nearest integer
 * Floor:=           rounding to nearest lower integer
 * Ceil :=           rounding to nearest higher integer
**/
template <typename T>
class BlockRoundingF3D: public BlockF3D<T> {
private:
  BlockF3D<T>& _f;
  const RoundingMode _rounding;
public:

  /**
   * \param f           Block functor
   * \param rounding    Mode how the value shall be rounded
   **/

  BlockRoundingF3D(BlockF3D<T>& f, RoundingMode rounding = RoundingMode::NearestInteger );

  bool operator() (T output[], const int input []);
};

}
#endif
