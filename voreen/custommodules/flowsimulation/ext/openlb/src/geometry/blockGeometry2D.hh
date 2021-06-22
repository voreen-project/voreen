/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2011, 2014 Mathias J. Krause, Simon Zimny
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
 * Representation of the 2D block geometry -- generic implementation.
 */

#ifndef BLOCK_GEOMETRY_2D_HH
#define BLOCK_GEOMETRY_2D_HH


#include "geometry/blockGeometry2D.h"

namespace olb {

template<typename T>
BlockGeometry2D<T>::BlockGeometry2D(T x0, T y0, T h, int nX, int nY, int iCglob)
  : BlockGeometryStructure2D<T>(nX, nY, iCglob),
    _data(this->getNcells()),
    _cuboid(x0, y0, h, nX, nY)
{
  addToStatisticsList( &(this->_statistics.getStatisticsStatus()) );
}

template<typename T>
BlockGeometry2D<T>::BlockGeometry2D(Cuboid2D<T>& cuboid, int iCglob)
  : BlockGeometryStructure2D<T>(cuboid.getNx(), cuboid.getNy(), iCglob),
    _data(this->getNcells()),
    _cuboid(cuboid)
{
  addToStatisticsList( &(this->_statistics.getStatisticsStatus()) );
}

template<typename T>
BlockGeometryStatistics2D<T>& BlockGeometry2D<T>::getStatistics(bool verbose)
{
  return this->_statistics;
}

template<typename T>
BlockGeometryStatistics2D<T> const& BlockGeometry2D<T>::getStatistics(bool verbose) const
{
  return this->_statistics;
}

template<typename T>
Vector<T,2> BlockGeometry2D<T>::getOrigin() const
{
  return _cuboid.getOrigin();
}

template<typename T>
T BlockGeometry2D<T>::getDeltaR() const
{
  return _cuboid.getDeltaR();
}

template<typename T>
int& BlockGeometry2D<T>::get(int iX, int iY)
{
  resetStatistics();
  return _data[0][this->getCellId(iX,iY)];
}

template<typename T>
int const& BlockGeometry2D<T>::get(int iX, int iY) const
{
  return _data[0][this->getCellId(iX,iY)];
}

template<typename T>
int BlockGeometry2D<T>::getMaterial(int iX, int iY) const
{
  int material;
  if (iX < 0 || iX + 1 > this->getNx() || iY < 0 || iY + 1 > this->getNy() ) {
    material = 0;
  }
  else {
    material = _data[0][this->getCellId(iX,iY)];
  }
  return material;
}

template<typename T>
int& BlockGeometry2D<T>::get(std::size_t iCell)
{
  resetStatistics();
  return _data[0][iCell];
}

template<typename T>
void BlockGeometry2D<T>::getPhysR(T physR[2], const int& iX, const int& iY) const
{
  _cuboid.getPhysR(physR, iX, iY);
}

template<typename T>
void BlockGeometry2D<T>::addToStatisticsList(bool* statisticStatus)
{
  _statisticsUpdateNeeded.push_back(statisticStatus);
}

template<typename T>
void BlockGeometry2D<T>::removeFromStatisticsList(bool* statisticStatus)
{
  _statisticsUpdateNeeded.remove(statisticStatus);
}

template<typename T>
void BlockGeometry2D<T>::printLayer(int x0, int x1, int y0, int y1, bool linenumber)
{
  for (int x = x0; x <= x1; x++) {
    if (linenumber) {
      this->clout << x << ": ";
    }
    for (int y = y0; y <= y1; y++) {
      this->clout << getMaterial(x, y) << " ";
    }
    if (x1 - x0 != 0) {
      this->clout << std::endl;
    }
  }
  this->clout << std::endl;
}

template<typename T>
void BlockGeometry2D<T>::printLayer(int direction, int layer, bool linenumber)
{
  assert(direction >= 0 && direction <= 2);
  switch (direction) {
  case 0:
    printLayer(layer, layer, 0, this->getNy() - 1, linenumber);
    break;
  case 1:
    printLayer(0, this->getNx() - 1, layer, layer, linenumber);
    break;
  }
}

template<typename T>
void BlockGeometry2D<T>::printNode(int x0, int y0)
{
  for (int x = x0 - 1; x <= x0 + 1; x++) {
    this->clout << "x=" << x << std::endl;
    for (int y = y0 - 1; y <= y0 + 1; y++) {
      this->clout << getMaterial(x, y) << " ";
    }
    this->clout << std::endl;
  }
  this->clout << std::endl;
}

template<typename T>
void BlockGeometry2D<T>::resetStatistics()
{
  for (bool* update : _statisticsUpdateNeeded) {
    *update = true;
  }
}

} // namespace olb

#endif
