/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2017 Adrian Kummerländer
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

#include "functors/lattice/blockBaseF3D.h"
#include "core/blockStructure3D.h"

namespace olb {

template<typename T> class BlockF3D;
template<typename T> class BlockGeometry3D;

/// Base block indicator functor (discrete)
template <typename T>
class BlockIndicatorF3D : public BlockF3D<bool> {
protected:
  BlockGeometry3D<T>& _blockGeometry;
public:
  using BlockF3D<bool>::operator();

  BlockIndicatorF3D(BlockGeometry3D<T>& blockGeometry);
  /**
   * Get underlying block geometry
   *
   * \returns _blockGeometry
   **/
  BlockGeometry3D<T>& getBlockGeometry();
  /**
   * Block indicator specific function operator overload
   *
   * The boolean return value of `operator()(T output[], S input[])` describes
   * the call's success and by convention must not describe the indicated domain.
   *
   * \return Domain indicator i.e. `true` iff the input lies within the described domain.
   **/
  virtual bool operator() (const int input[]);
};

} // namespace olb

#endif
