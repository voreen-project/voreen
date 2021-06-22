/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014 Mathias J. Krause
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

/** \file
 * Representation of the 3D block geometry view -- generic implementation.
 */

#ifndef BLOCK_GEOMETRY_VIEW_3D_HH
#define BLOCK_GEOMETRY_VIEW_3D_HH


#include <vector>
#include "geometry/blockGeometryView3D.h"


namespace olb {

template<typename T>
BlockGeometryView3D<T>::BlockGeometryView3D(
  BlockGeometryStructure3D<T>& originalBlockGeometry,
  int x0, int x1, int y0,
  int y1, int z0, int z1)
  : BlockGeometryStructure3D<T>(x1-x0+1, y1-y0+1, z1-z0+1, originalBlockGeometry.getIcGlob()),
    _originalBlockGeometry(&originalBlockGeometry),
    _x0(x0), _y0(y0), _z0(z0)
{
  this->_statistics = BlockGeometryStatistics3D<T>(this);
  addToStatisticsList( &(this->_statistics.getStatisticsStatus()) );
}

template<typename T>
BlockGeometryView3D<T>::~BlockGeometryView3D()
{
  removeFromStatisticsList( &(this->_statistics.getStatisticsStatus()) );
}

template<typename T>
BlockGeometryStatistics3D<T>& BlockGeometryView3D<T>::getStatistics(bool verbose)
{
  return this->_statistics;
}

template<typename T>
BlockGeometryStatistics3D<T> const& BlockGeometryView3D<T>::getStatistics(bool verbose) const
{
  return this->_statistics;
}

template<typename T>
Vector<T,3> BlockGeometryView3D<T>::getOrigin() const
{
  Vector<T,3> origin(_originalBlockGeometry->getOrigin());
  origin[0] += _x0 * getDeltaR();
  origin[1] += _y0 * getDeltaR();
  origin[2] += _z0 * getDeltaR();
  return origin;
}

template<typename T>
T BlockGeometryView3D<T>::getDeltaR() const
{
  return _originalBlockGeometry->getDeltaR();
}

template<typename T>
int& BlockGeometryView3D<T>::get(int iX, int iY, int iZ)
{
  return _originalBlockGeometry->get(_x0+iX, _y0+iY, _z0+iZ);
}

template<typename T>
const int& BlockGeometryView3D<T>::get(int iX, int iY, int iZ) const
{
  return _originalBlockGeometry->get(_x0+iX, _y0+iY, _z0+iZ);
}

template<typename T>
int BlockGeometryView3D<T>::getMaterial(int iX, int iY, int iZ) const
{
  return _originalBlockGeometry->getMaterial(_x0+iX, _y0+iY, _z0+iZ);
}


template<typename T>
void BlockGeometryView3D<T>::getPhysR(T physR[3], const int& iX, const int& iY, const int& iZ) const
{
  _originalBlockGeometry->getPhysR(physR, _x0 + iX, _y0 + iY, _z0 + iZ);
  return;
}


template<typename T>
void BlockGeometryView3D<T>::addToStatisticsList(bool* statisticStatus)
{
  _originalBlockGeometry->addToStatisticsList(statisticStatus);
}

template<typename T>
void BlockGeometryView3D<T>::removeFromStatisticsList(bool* statisticStatus)
{
  _originalBlockGeometry->removeFromStatisticsList(statisticStatus);
}

} // namespace olb

#endif
