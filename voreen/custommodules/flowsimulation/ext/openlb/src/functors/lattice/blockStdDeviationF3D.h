/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2019 Jakob Mangold, Mathias J. Krause
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

#ifndef BLOCK_STDDEVIATION_F3D_H
#define BLOCK_STDDEVIATION_F3D_H

#include "geometry/cuboid3D.h"
#include "indicator/blockIndicatorBaseF3D.h"
#include "blockAverage3D.h"
#include "blockBaseF3D.h"


namespace olb {


/// BlockStdDeviationF3D returns the Deviation in each component of f on a indicated subset
template <typename T, typename W = T>
class BlockStdDeviationF3D final : public BlockAverage3D<W> {
private:
  BlockF3D<W>&          _f;
  BlockIndicatorF3D<T>& _indicatorF;
  Cuboid3D<T>&      _cuboid;
  T _expectedValue;
public:
  BlockStdDeviationF3D(BlockF3D<W>&          f,
                       BlockIndicatorF3D<T>& indicatorF,
                       Cuboid3D<T>&      cuboid,
                       T expectedValue);
  bool operator() (W output[], const int input[]) override;
};


}

#endif
