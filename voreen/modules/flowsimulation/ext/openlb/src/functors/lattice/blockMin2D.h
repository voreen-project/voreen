/*  This file is part of the OpenLB library
 *
 *  Copyright (C) Adrian Kummerlaender
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

#ifndef BLOCK_MIN_2D_H
#define BLOCK_MIN_2D_H

#include "blockBaseF2D.h"
#include "geometry/cuboid2D.h"
#include "indicator/blockIndicatorBaseF2D.h"

namespace olb {


/// BlockMin2D returns the min in each component of f on a indicated subset
template <typename T, typename W = T>
class BlockMin2D final : public BlockF2D<W> {
private:
  BlockF2D<W>&          _f;
  BlockIndicatorF2D<T>& _indicatorF;
  Cuboid2D<T>&          _cuboid;
public:
  BlockMin2D(BlockF2D<W>&          f,
             BlockIndicatorF2D<T>& indicatorF,
             Cuboid2D<T>&          cuboid);
  bool operator() (W output[], const int input[]) override;
};


}

#endif
