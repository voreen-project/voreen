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
 * Dynamics for a generic 2D block -- header file.
 */
#ifndef BLOCK_STRUCTURE_2D_H
#define BLOCK_STRUCTURE_2D_H

#include <cstdint>
#include <array>

#include "core/olbDebug.h"
#include "core/vector.h"

namespace olb {

/** An empty hull with left bottom corner at (0,0).
 *
 * \param _nx extension in x direction
 * \param _ny extension in y direction
 *
 */

  static int foo_x = 0;
  static int foo_y = 0;
class BlockStructure2D {
protected:
  /// Block width
  int _nx;
  /// Block height
  int _ny;
public:
  BlockStructure2D(int nx, int ny) : _nx(nx), _ny(ny) {};
  BlockStructure2D(int nx, int ny, int overlap):
    _nx(nx+2*overlap), _ny(ny+2*overlap) { };
  BlockStructure2D(): BlockStructure2D(1,1) {};
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
  /// Get number of cells
  std::size_t getNcells() const
  {
    // The conversions to std::size_t ensure 64-bit compatibility. Note that
    // Nx and Ny are of type int, which might be 32-bit types, even on
    // 64-bit platforms. Therefore, Nx*Ny may lead to a type overflow.
    return static_cast<std::size_t>(getNx())
           * static_cast<std::size_t>(getNy());
  }

  /// Get 1D cell ID
  std::size_t getCellId(int iX, int iY) const
  {


    if(!(iX >= 0 && iX < this->_nx && iY >= 0 && iY < this->_ny)){
      throw std::runtime_error{"Bad access at: " + std::to_string(iX) + " " + std::to_string(iY) + " " + std::to_string(foo_x) + " " + std::to_string(foo_y) + " "+ std::to_string(_nx) + " " + std::to_string(_ny)};

    }

    foo_x = iX;
    foo_y = iY;

    OLB_PRECONDITION(iX >= 0 && iX < this->_nx);
    OLB_PRECONDITION(iY >= 0 && iY < this->_ny);
    return iX*_ny + iY;
  }
  /// Get 1D neighbor distance
  std::ptrdiff_t getNeighborDistance(int iX, int iY) const
  {
    return iX*_ny + iY;
  }
  /// Get 1D neighbor distance
  std::ptrdiff_t getNeighborDistance(Vector<int,2> c) const
  {
    return getNeighborDistance(c[0], c[1]);
  }
  /// Return whether location is valid
  bool isInside(int iX, int iY) const
  {
    return 0 <= iX && iX < getNx() &&
           0 <= iY && iY < getNy();
  };
};

}  // namespace olb

#endif
