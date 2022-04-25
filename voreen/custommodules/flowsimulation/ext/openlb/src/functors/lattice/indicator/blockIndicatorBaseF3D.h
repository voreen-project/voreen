/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2017 Adrian Kummerlaender
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

#ifndef BLOCK_INDICATOR_BASE_F_3D_H
#define BLOCK_INDICATOR_BASE_F_3D_H

#include "functors/lattice/superBaseF3D.h"
#include "functors/lattice/blockBaseF3D.h"
#include "core/blockData.h"
#include "core/blockStructure.h"
#include "geometry/blockGeometry.h"

namespace olb {

/// Base block indicator functor
/**
 * Derived functors implement a indicator function on the full domain given by BlockGeometry3D.
 **/
template <typename T>
class BlockIndicatorF3D : public BlockF3D<bool> {
protected:
  BlockGeometry<T,3>& _block;
  const BlockData<3,T,bool>* _cachedData;
public:
  using BlockF3D<bool>::operator();

  BlockIndicatorF3D(BlockGeometry<T,3>& geometry);

  /// Get underlying block geometry structure
  /**
   * \returns _block
   **/
  BlockGeometry<T,3>& getBlockGeometry();

  /// Block indicator specific function operator overload
  /**
   * The boolean return value of `operator()(T output[], S input[])` describes
   * the call's success and by convention must not describe the indicated domain.
   *
   * \return Domain indicator i.e. `true` iff the input lies within the described subset.
   **/
  bool operator() (const int input[]);
  bool operator() (int iX, int iY, int iZ);
  bool operator() (LatticeR<3> loc);

  /// Set bool-mask cache to be used by indicator operator overloads
  void setCache(const BlockData<3,T,bool>& cache);

  /// Returns true only if the indicated domain subset is empty
  /**
   * May return false even if the indicated domain subset is in fact empty.
   * Primarily implemented to minimize block accesses if an empty domain can
   * be inferred by e.g. BlockGeometryStatistics3D data.
   *
   * i.e. only override this method if the domain can be checked for emptyness
   *      in an efficient fashion.
   **/
  virtual bool isEmpty();
  /// Returns min lattice position of the indicated subset's bounding box
  virtual Vector<int,3> getMin() = 0;
  /// Returns max lattice position of the indicated subset's bounding box
  virtual Vector<int,3> getMax() = 0;

};

} // namespace olb

#endif
