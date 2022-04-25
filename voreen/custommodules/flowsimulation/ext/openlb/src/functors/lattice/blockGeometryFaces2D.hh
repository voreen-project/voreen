/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2013-2018 Mathias Krause, Albert Mink, Adrian Kummerlaender
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

#ifndef BLOCK_GEOMETRY_FACES_2D_HH
#define BLOCK_GEOMETRY_FACES_2D_HH

#include "blockGeometryFaces2D.h"

namespace olb {


template <typename T>
BlockGeometryFaces2D<T>::BlockGeometryFaces2D(BlockIndicatorF2D<T>& indicatorF, T latticeL)
  : BlockF2D<T>(indicatorF.getBlockStructure(), 5),
    _indicatorF(indicatorF),
    _latticeL(latticeL)
{
  this->getName() = "blockGeometryFaces";
}

template <typename T>
bool BlockGeometryFaces2D<T>::operator() (T output[], const int input[])
{
  for (int i=0; i<5; ++i) {
    output[i] = T();
  }

  std::size_t counter[5] = {0};

  if (!_indicatorF.isEmpty()) {
    auto& blockGeometry = _indicatorF.getBlockGeometry();
    const Vector<int,2> min = _indicatorF.getMin();
    const Vector<int,2> max = _indicatorF.getMax();

    // Iterate over all cells and count the cells of the face
    for (int iX = min[0]; iX <= max[0]; ++iX) {
      for (int iY = min[1]; iY <= max[1]; ++iY) {
          // Lock at solid nodes only
          if (_indicatorF(iX, iY)) {
            if (blockGeometry.getMaterial({iX-1, iY}) == 1) {
              counter[0]++;
            }
            if (blockGeometry.getMaterial({iX, iY-1}) == 1) {
              counter[1]++;
            }
            if (blockGeometry.getMaterial({iX+1, iY}) == 1) {
              counter[2]++;
            }
            if (blockGeometry.getMaterial({iX, iY+1}) == 1) {
              counter[3]++;
            }
          }
      }
    }

    const T dx2 = _latticeL*_latticeL;
    T total = T(); //added
    for (int i=0; i<4; ++i) {
      output[i]  = (T) counter[i] * dx2;
      total+= (T) counter[i] * dx2; //was output[4]
    }
    output[4]=total;
    return true;
  } else {
    for (int i=0; i<5; ++i) {
      output[i]=T();
    }
    return true;
  }
  return false;
}

template <typename T, bool HLBM>
BlockGeometryFacesIndicator2D<T,HLBM>::BlockGeometryFacesIndicator2D(
  BlockGeometry<T,2>& blockGeometry, SmoothIndicatorF2D<T,T,HLBM>& indicator,
  int material, T latticeL)
  : GenericF<T,int>(5,0), _blockGeometry(blockGeometry), _indicator(indicator),
    _material(material), _latticeL(latticeL)  // _latticeLsqr(latticeL*latticeL)
{
      this->getName() = "facesSmoothInd";
}
template <typename T, bool HLBM>
bool BlockGeometryFacesIndicator2D<T,HLBM>::operator() (T output[], const int input[])
{
  int counter[4] = {0,0,0,0};
  T inside[1];
  T physR[2];
  if (_blockGeometry.getStatistics().getNvoxel(_material)!=0) {
    const int x0 = _blockGeometry.getStatistics().getMinLatticeR(_material)[0];
    const int y0 = _blockGeometry.getStatistics().getMinLatticeR(_material)[1];
    const int x1 = _blockGeometry.getStatistics().getMaxLatticeR(_material)[0];
    const int y1 = _blockGeometry.getStatistics().getMaxLatticeR(_material)[1];

    // Iterate over all cells and count the cells of the face
    for (int iX = x0; iX <= x1; ++iX) {
      for (int iY = y0; iY <= y1; ++iY) {
          // Look at solid nodes only
        _blockGeometry.getPhysR(physR, {iX, iY});
          _indicator(inside, physR);
          if ( !util::nearZero(inside[0]) ) {
            _blockGeometry.getPhysR(physR, {iX-1, iY});
            _indicator(inside, physR);
            if ( util::nearZero(inside[0]) )
              counter[0]++;
            _blockGeometry.getPhysR(physR, {iX, iY-1});
            _indicator(inside, physR);
            if ( util::nearZero(inside[0]) )
              counter[1]++;
            _blockGeometry.getPhysR(physR, {iX+1, iY});
            _indicator(inside, physR);
            if ( util::nearZero(inside[0]) )
              counter[2]++;
            _blockGeometry.getPhysR(physR, {iX, iY+1});
            _indicator(inside, physR);
            if ( util::nearZero(inside[0]) )
              counter[3]++;
          }
      }
    }

    T total = T();
    for (int i=0; i<4; ++i) {
      output[i]= ((T) counter[i]) * _latticeL;
      total+= ((T) counter[i]) * _latticeL;
    }
    output[4]=total;
    return true;
  } else {
    for (int i=0; i<5; ++i) {
      output[i]=T();
    }
    return true;
  }
  return false;

}

}
#endif
