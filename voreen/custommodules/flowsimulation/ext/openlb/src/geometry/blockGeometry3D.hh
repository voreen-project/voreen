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
 * Representation of the 3D block geometry -- generic implementation.
 */

#ifndef BLOCK_GEOMETRY_3D_HH
#define BLOCK_GEOMETRY_3D_HH


#include "geometry/blockGeometry3D.h"


namespace olb {

template<typename T>
BlockGeometry3D<T>::BlockGeometry3D(T x0, T y0, T z0, T h, int nX, int nY, int nZ, int iCglob)
  : BlockGeometryStructure3D<T>(nX, nY, nZ, iCglob),
    _data(this->getNcells()),
    _cuboid(x0, y0, z0, h, nX, nY, nZ)
{
  this->_statistics = BlockGeometryStatistics3D<T>(this);
  addToStatisticsList( &(this->_statistics.getStatisticsStatus()) );
}

template<typename T>
BlockGeometry3D<T>::BlockGeometry3D(Cuboid3D<T>& cuboid, int iCglob)
  : BlockGeometryStructure3D<T>(cuboid.getNx(), cuboid.getNy(), cuboid.getNz(), iCglob),
    _data(this->getNcells()),
    _cuboid(cuboid)
{
  this->_statistics = BlockGeometryStatistics3D<T>(this);
  addToStatisticsList( &(this->_statistics.getStatisticsStatus()) );
}

template<typename T>
BlockGeometryStatistics3D<T>& BlockGeometry3D<T>::getStatistics(bool verbose)
{
  return this->_statistics;
}

template<typename T>
BlockGeometryStatistics3D<T> const& BlockGeometry3D<T>::getStatistics(bool verbose) const
{
  return this->_statistics;
}

template<typename T>
Vector<T,3> BlockGeometry3D<T>::getOrigin() const
{
  return _cuboid.getOrigin();
}

template<typename T>
T BlockGeometry3D<T>::getDeltaR() const
{
  return _cuboid.getDeltaR();
}

template<typename T>
int& BlockGeometry3D<T>::get(int iX, int iY, int iZ)
{
  resetStatistics();
  return _data[0][this->getCellId(iX,iY,iZ)];
}

template<typename T>
int& BlockGeometry3D<T>::get(std::size_t iCell)
{
  resetStatistics();
  return _data[0][iCell];
}

template<typename T>
int const& BlockGeometry3D<T>::get(int iX, int iY, int iZ) const
{
  return _data[0][this->getCellId(iX,iY,iZ)];
}

template<typename T>
int const& BlockGeometry3D<T>::get(std::size_t iCell) const
{
  return _data[0][iCell];
}

template<typename T>
int BlockGeometry3D<T>::getMaterial(int iX, int iY, int iZ) const
{
  int material;
  if (iX < 0 || iX + 1 > this->getNx() || iY < 0 || iY
      + 1 > this->getNy() || iZ < 0 || iZ + 1 > this->getNz()) {
    material = 0;
  }
  else {
    material = _data[0][this->getCellId(iX,iY,iZ)];
  }
  return material;
}

template<typename T>
void BlockGeometry3D<T>::getPhysR(T physR[3], const int& iX, const int& iY,
                                  const int& iZ) const
{
  _cuboid.getPhysR(physR, iX, iY, iZ);
  return;
}

template<typename T>
void BlockGeometry3D<T>::addToStatisticsList(bool* statisticStatus)
{
  _statisticsUpdateNeeded.push_back(statisticStatus);
}

template<typename T>
void BlockGeometry3D<T>::removeFromStatisticsList(bool* statisticStatus)
{
  _statisticsUpdateNeeded.remove(statisticStatus);
}

template<typename T>
void BlockGeometry3D<T>::printLayer(int x0, int x1, int y0, int y1, int z0,
                                    int z1, bool linenumber)
{
  for (int x = x0; x <= x1; x++) {
    if (linenumber) {
      this->clout << x << ": ";
    }
    for (int y = y0; y <= y1; y++) {
      for (int z = z0; z <= z1; z++) {
        this->clout << getMaterial(x, y, z) << " ";
      }
      if (y1 - y0 != 0 && z1 - z0 != 0) {
        this->clout << std::endl;
      }
    }
    if (x1 - x0 != 0) {
      this->clout << std::endl;
    }
  }
  this->clout << std::endl;
}

template<typename T>
void BlockGeometry3D<T>::printLayer(int direction, int layer, bool linenumber)
{
  assert(direction >= 0 && direction <= 2);
  switch (direction) {
  case 0:
    printLayer(layer, layer, 0, this->getNy() - 1, 0, this->getNz() - 1, linenumber);
    break;
  case 1:
    printLayer(0, this->getNx() - 1, layer, layer, 0, this->getNz() - 1, linenumber);
    break;
  case 2:
    printLayer(0, this->getNx() - 1, 0, this->getNy() - 1, layer, layer, linenumber);
    break;
  }
}

template<typename T>
void BlockGeometry3D<T>::printNode(int x0, int y0, int z0)
{
  for (int x = x0 - 1; x <= x0 + 1; x++) {
    this->clout << "x=" << x << std::endl;
    for (int y = y0 - 1; y <= y0 + 1; y++) {
      for (int z = z0 - 1; z <= z0 + 1; z++) {
        this->clout << getMaterial(x, y, z) << " ";
      }
      this->clout << std::endl;
    }
    this->clout << std::endl;
  }
  this->clout << std::endl;
}

template<typename T>
void BlockGeometry3D<T>::resetStatistics()
{
  for (std::list<bool*>::iterator it = _statisticsUpdateNeeded.begin(); it != _statisticsUpdateNeeded.end(); ++it) {
    **it = true;
  }
}

} // namespace olb

#endif
