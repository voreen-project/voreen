/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2015 Mathias J. Krause
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
 * Dynamics for a generic 3D block -- header file.
 */
#ifndef BLOCK_STRUCTURE_3D_H
#define BLOCK_STRUCTURE_3D_H

#include <cstdint>

#include "core/olbDebug.h"
#include "core/vector.h"

namespace olb {

/** An empty hull with left bottom corner at (0,0,0).
 *
 * \param _nx extension in x direction
 * \param _ny extension in y direction
 * \param _nz extension in z direction
 *
 */
class BlockStructure3D {
protected:
  /// Block width
  int _nx;
  /// Block height
  int _ny;
  /// Block lenght
  int _nz;
public:
  BlockStructure3D(int nx, int ny, int nz):
    _nx(nx), _ny(ny), _nz(nz) {};
  BlockStructure3D(int nx, int ny, int nz, int overlap):
    _nx(nx+2*overlap), _ny(ny+2*overlap), _nz(nz+2*overlap) {};

  /// Read only access to block width
  int getNx() const
  {
    return _nx;
  };
  /// Read only access to block height
  int getNy() const
  {
    return _ny;
  };
  /// Read only access to block height
  int getNz() const
  {
    return _nz;
  };
  /// Get number of cells
  std::size_t getNcells() const
  {
    // The conversions to std::size_t ensure 64-bit compatibility. Note that
    // Nx, Ny and Nz are of type int, which might be 32-bit types, even on
    // 64-bit platforms. Therefore, Nx*Ny*Nz may lead to a type overflow.
    return static_cast<std::size_t>(getNx())
         * static_cast<std::size_t>(getNy())
         * static_cast<std::size_t>(getNz());
  }
  /// Get 1D cell ID
  std::size_t getCellId(int iX, int iY, int iZ) const
  {
    OLB_PRECONDITION(iX >= 0 && iX < this->_nx);
    OLB_PRECONDITION(iY >= 0 && iY < this->_ny);
    OLB_PRECONDITION(iZ >= 0 && iZ < this->_nz);
    return iX*_ny*_nz + iY*_nz + iZ;
  }
  /// Get 1D neighbor distance
  std::ptrdiff_t getNeighborDistance(int iX, int iY, int iZ) const
  {
    return iX*_ny*_nz + iY*_nz + iZ;
  }
  /// Get 1D neighbor distance
  std::ptrdiff_t getNeighborDistance(Vector<int,3> c) const
  {
    return getNeighborDistance(c[0], c[1], c[2]);
  }
  /// Return whether location is valid
  bool isInside(int iX, int iY, int iZ) const
  {
    return 0 <= iX && iX < getNx() &&
           0 <= iY && iY < getNy() &&
           0 <= iZ && iZ < getNz();
  };
};

}  // namespace olb

#endif
