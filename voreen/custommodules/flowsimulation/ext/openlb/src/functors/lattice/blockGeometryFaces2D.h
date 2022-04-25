/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2018 Adrian Kummerlaender
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

#ifndef BLOCK_GEOMETRY_FACES_2D_H
#define BLOCK_GEOMETRY_FACES_2D_H

#include "blockBaseF2D.h"
#include "indicator/blockIndicatorBaseF2D.h"

namespace olb {


template <typename T>
class BlockGeometryFaces2D final : public BlockF2D<T> {
private:
  BlockIndicatorF2D<T>& _indicatorF;
  T _latticeL;
public:
  BlockGeometryFaces2D(BlockIndicatorF2D<T>& indicatorF, T latticeL);
  bool operator() (T output[], const int input[]) override;
};

template <typename T, bool HLBM=false>
class BlockGeometryFacesIndicator2D final : public GenericF<T,int> {
private:
  BlockGeometry<T,2>& _blockGeometry;
  SmoothIndicatorF2D<T,T,HLBM>& _indicator;
  int _material;
  T _latticeL;
public:
  BlockGeometryFacesIndicator2D(BlockGeometry<T,2>& blockGeometry,
                                SmoothIndicatorF2D<T,T,HLBM>& indicator,
                                int material, T deltaX);
  bool operator() (T output[], const int input[]) override;
};

}

#endif
