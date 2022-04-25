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

#ifndef BLOCK_GEOMETRY_FACES_3D_HH
#define BLOCK_GEOMETRY_FACES_3D_HH

#include "blockGeometryFaces3D.h"

namespace olb {


template <typename T>
BlockGeometryFaces3D<T>::BlockGeometryFaces3D(BlockIndicatorF3D<T>& indicatorF, T latticeL)
  : BlockF3D<T>(indicatorF.getBlockStructure(), 7),
    _indicatorF(indicatorF),
    _latticeL(latticeL)
{
  this->getName() = "blockGeometryFaces";
}

template <typename T>
bool BlockGeometryFaces3D<T>::operator() (T output[], const int input[])
{
  for (int i=0; i<7; ++i) {
    output[i] = T();
  }

  std::size_t counter[7] = {0};

  if (!_indicatorF.isEmpty()) {
    auto& blockGeometry = _indicatorF.getBlockGeometry();
    const Vector<int,3> min = _indicatorF.getMin();
    const Vector<int,3> max = _indicatorF.getMax();

    // Iterate over all cells and count the cells of the face
    for (int iX = min[0]; iX <= max[0]; ++iX) {
      for (int iY = min[1]; iY <= max[1]; ++iY) {
        for (int iZ = min[2]; iZ <= max[2]; ++iZ) {
          // Lock at solid nodes only
          if (_indicatorF(iX, iY, iZ)) {
            if (blockGeometry.getMaterial({iX-1, iY, iZ}) == 1) {
              counter[0]++;
            }
            if (blockGeometry.getMaterial({iX, iY-1, iZ}) == 1) {
              counter[1]++;
            }
            if (blockGeometry.getMaterial({iX, iY, iZ-1}) == 1) {
              counter[2]++;
            }
            if (blockGeometry.getMaterial({iX+1, iY, iZ}) == 1) {
              counter[3]++;
            }
            if (blockGeometry.getMaterial({iX, iY+1, iZ}) == 1) {
              counter[4]++;
            }
            if (blockGeometry.getMaterial({iX, iY, iZ+1}) == 1) {
              counter[5]++;
            }
          }
        }
      }
    }

    const T dx2 = _latticeL*_latticeL;
    for (int i=0; i<6; ++i) {
      output[i]  = (T) counter[i] * dx2;
      output[6] += (T) counter[i] * dx2;
    }
  }

  return true;
}

template <typename T, bool HLBM>
BlockGeometryFacesIndicator3D<T,HLBM>::BlockGeometryFacesIndicator3D(
  BlockGeometry<T,3>& blockGeometry, SmoothIndicatorF3D<T,T,HLBM>& indicator,
  int material, T latticeL)
  : GenericF<T,int>(7,0), _blockGeometry(blockGeometry), _indicator(indicator),
    _material(material), _latticeLsqr(latticeL*latticeL)
{
  this->getName() = "facesSmoothInd";
}
template <typename T, bool HLBM>
bool BlockGeometryFacesIndicator3D<T,HLBM>::operator() (T output[], const int input[])
{
  int counter[6] = {0,0,0,0,0,0};
  T inside[1];
  T physR[3];
  if (_blockGeometry.getStatistics().getNvoxel(_material)!=0) {
    const int x0 = _blockGeometry.getStatistics().getMinLatticeR(_material)[0];
    const int y0 = _blockGeometry.getStatistics().getMinLatticeR(_material)[1];
    const int z0 = _blockGeometry.getStatistics().getMinLatticeR(_material)[2];
    const int x1 = _blockGeometry.getStatistics().getMaxLatticeR(_material)[0];
    const int y1 = _blockGeometry.getStatistics().getMaxLatticeR(_material)[1];
    const int z1 = _blockGeometry.getStatistics().getMaxLatticeR(_material)[2];

    // Iterate over all cells and count the cells of the face
    for (int iX = x0; iX <= x1; ++iX) {
      for (int iY = y0; iY <= y1; ++iY) {
        for (int iZ = z0; iZ <= z1; ++iZ) {
          // Look at solid nodes only
          _blockGeometry.getPhysR(physR, {iX, iY, iZ});
          _indicator(inside, physR);
          if ( !util::nearZero(inside[0]) ) {
            _blockGeometry.getPhysR(physR, {iX-1, iY, iZ});
            _indicator(inside, physR);
            if ( util::nearZero(inside[0]) ) {
              counter[0]++;
            }
            _blockGeometry.getPhysR(physR, {iX, iY-1, iZ});
            _indicator(inside, physR);
            if ( util::nearZero(inside[0]) ) {
              counter[1]++;
            }
            _blockGeometry.getPhysR(physR, {iX, iY, iZ-1});
            _indicator(inside, physR);
            if ( util::nearZero(inside[0]) ) {
              counter[2]++;
            }
            _blockGeometry.getPhysR(physR, {iX+1, iY, iZ});
            _indicator(inside, physR);
            if ( util::nearZero(inside[0]) ) {
              counter[3]++;
            }
            _blockGeometry.getPhysR(physR, {iX, iY+1, iZ});
            _indicator(inside, physR);
            if ( util::nearZero(inside[0]) ) {
              counter[4]++;
            }
            _blockGeometry.getPhysR(physR, {iX, iY, iZ+1});
            _indicator(inside, physR);
            if ( util::nearZero(inside[0]) ) {
              counter[5]++;
            }
          }
        }
      }
    }

    T total = T();
    for (int i=0; i<6; ++i) {
      output[i]= ((T) counter[i]) * _latticeLsqr;
      total+= ((T) counter[i]) * _latticeLsqr;
    }
    output[6]=total;
    return true;
  }
  else {
    for (int i=0; i<7; ++i) {
      output[i]=T();
    }
    return true;
  }
  return false;

}

}
#endif
