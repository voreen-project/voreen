/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Mathias Krause, Adrian Kummerlaender
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

#ifndef BLOCK_STRUCTURE_H
#define BLOCK_STRUCTURE_H

#include <cstdint>
#include <type_traits>

#include "vector.h"
#include "olbDebug.h"

namespace olb {

/// Type for sequential block-local cell indices
using CellID = std::uint32_t;

/// Type for in-memory distance of block-local cell indices
using CellDistance = std::int64_t;

/// Type for spatial block-local lattice coordinates
template <unsigned D>
using LatticeR = Vector<std::int32_t,D>;

/// Type for spatial (physical) coordinates
template <typename T, unsigned D>
using PhysR = Vector<T,D>;

/// Base of a regular block
/**
 * With extent, optional padding and memory bijection for spatial locations
 **/
template <unsigned D>
class BlockStructureD {
protected:
  LatticeR<D> _core;
  LatticeR<D> _size;
  LatticeR<D> _projection;

  int _padding;

public:
  static_assert(D == 2 || D == 3, "Only D=2 and D=3 are supported");

  BlockStructureD(Vector<int,D> size, int padding=0):
    _core(size),
    _size(size + 2*padding),
    _padding(padding)
  {
    if constexpr (D == 3) {
      _projection = {_size[1]*_size[2], _size[2], 1};
    } else {
      _projection = {_size[1], 1};
    }

    if (getNcells()-1 > std::numeric_limits<CellID>::max()) {
      throw std::invalid_argument("Cell count must not exceed cell index space");
    }
  };

  BlockStructureD():
    BlockStructureD(1, 0)
  { };

  /// Read only access to block width
  int getNx() const
  {
    return _core[0];
  };
  /// Read only access to block height
  int getNy() const
  {
    return _core[1];
  };
  /// Read only access to block height
  int getNz() const
  {
    static_assert(D == 3, "z-component only available in 3D");
    return _core[2];
  };

  LatticeR<D> getExtent() const
  {
    return _core;
  }

  /// Read only access to padding
  int getPadding() const
  {
    return _padding;
  };

  /// Get number of cells
  std::size_t getNcells() const
  {
    if constexpr (D == 3) {
      return _size[0] * _size[1] * _size[2];
    } else {
      return _size[0] * _size[1];
    }
	__builtin_unreachable();
  }

  /// Get 1D cell ID
  CellID getCellId(LatticeR<D> latticeR) const
  {
    OLB_PRECONDITION(isInside(latticeR));
    return (latticeR+_padding) * _projection;
  }

  template <typename... L>
  std::enable_if_t<sizeof...(L) == D, CellID>
  getCellId(L... latticeR) const
  {
    return LatticeR<D>{latticeR+_padding...} * _projection;
  }

  /// Get 1D neighbor distance
  CellDistance getNeighborDistance(LatticeR<D> dir) const
  {
    return dir * _projection;
  }

  /// Return whether location is valid
  bool isInside(LatticeR<D> latticeR) const
  {
    return latticeR >= -_padding && latticeR < _core + _padding;
  };

  /// Return whether location is valid
  bool isPadding(LatticeR<D> latticeR) const
  {
    return isInside(latticeR) && !(latticeR >= 0 && latticeR < _core);
  };

  template <typename... L>
  std::enable_if_t<sizeof...(L) == D, bool>
  isInside(L... latticeR) const
  {
    return isInside({latticeR...});
  };

  /// Return maximum valid neighborhood sphere radius w.r.t. latticeR
  CellDistance getNeighborhoodRadius(LatticeR<D> latticeR) const
  {
    auto lower = latticeR + _padding;
    auto upper = LatticeR<D>([&](unsigned iDim) -> int {
      int x = lower[iDim] - _size[iDim] + 1;
      return x < 0 ? -x : 0;
    });
    auto y = std::min(*std::min_element(lower.data(), lower.data() + D),
                      *std::min_element(upper.data(), upper.data() + D));
    return y;
  };

  template <typename F>
  void forSpatialLocations(F f) const
  {
    using loc = typename LatticeR<D>::value_t;
    for (loc iX=-_padding; iX < _core[0] + _padding; ++iX) {
      for (loc iY=-_padding; iY < _core[1] + _padding; ++iY) {
        if constexpr (D == 3) {
          for (loc iZ=-_padding; iZ < _core[2] + _padding; ++iZ) {
            if constexpr (std::is_invocable_v<F, LatticeR<D>>) {
              f({iX,iY,iZ});
            } else {
              f(iX,iY,iZ);
            }
          }
        } else {
          if constexpr (std::is_invocable_v<F, LatticeR<D>>) {
            f({iX,iY});
          } else {
            f(iX,iY);
          }
        }
      }
    }
  };

  template <typename F>
  void forSpatialLocations(LatticeR<D> min, LatticeR<D> max, F f) const
  {
    using loc = typename LatticeR<D>::value_t;
    for (loc iX=std::max(-_padding, min[0]); iX < std::min(_core[0] + _padding, max[0]+1); ++iX) {
      for (loc iY=std::max(-_padding, min[1]); iY < std::min(_core[1] + _padding, max[1]+1); ++iY) {
        if constexpr (D == 3) {
          for (loc iZ=std::max(-_padding, min[2]); iZ < std::min(_core[2] + _padding, max[2]+1); ++iZ) {
            if constexpr (std::is_invocable_v<F, LatticeR<D>>) {
              f({iX,iY,iZ});
            } else {
              f(iX,iY,iZ);
            }
          }
        } else {
          if constexpr (std::is_invocable_v<F, LatticeR<D>>) {
            f({iX,iY});
          } else {
            f(iX,iY);
          }
        }
      }
    }
  };

  template <typename F>
  void forCoreSpatialLocations(F f) const
  {
    using loc = typename LatticeR<D>::value_t;
    for (loc iX=0; iX < _core[0]; ++iX) {
      for (loc iY=0; iY < _core[1]; ++iY) {
        if constexpr (D == 3) {
          for (loc iZ=0; iZ < _core[2]; ++iZ) {
            if constexpr (std::is_invocable_v<F, LatticeR<D>>) {
              f({iX,iY,iZ});
            } else {
              f(iX,iY,iZ);
            }
          }
        } else {
          if constexpr (std::is_invocable_v<F, LatticeR<D>>) {
            f({iX,iY});
          } else {
            f(iX,iY);
          }
        }
      }
    }
  };

  template <typename F>
  void forCellIndices(F f) const
  {
    for (CellID iCell=0; iCell < getNcells(); ++iCell) {
      f(iCell);
    }
  }

};

template <typename DESCRIPTOR>
using BlockStructure = BlockStructureD<DESCRIPTOR::d>;


}

#endif
