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

#ifndef BLOCK_INDICATOR_F_2D_H
#define BLOCK_INDICATOR_F_2D_H

#include "blockIndicatorBaseF2D.h"
#include "geometry/blockGeometryView2D.h"

namespace olb {

/// BlockIndicatorF2D from IndicatorF2D
template <typename T>
class BlockIndicatorFfromIndicatorF2D : public BlockIndicatorF2D<T> {
protected:
  IndicatorF2D<T>& _indicatorF;
  Cuboid2D<T>&     _cuboid;
public:
  BlockIndicatorFfromIndicatorF2D(IndicatorF2D<T>&  indicatorF,
                                  BlockStructure2D& blockStructure,
                                  Cuboid2D<T>&      cuboid);
  bool operator() (bool output[], const int input[]) override;
};


/// Block indicator functor from material numbers
template <typename T>
class BlockIndicatorMaterial2D : public BlockIndicatorF2D<T> {
protected:
  BlockGeometry2D<T>&           _blockGeometry;
  const int                     _overlap;
  const std::vector<int>* const _materialNumbers;
public:
  /**
   * \param blockGeometry Block geometry to be queried, accessible via SuperGeometry2D::getExtendedBlockGeometry
   * \param overlap       Overlap of given block geometry
   * \param materials     Material number vector is accepted via pointer to the actual vector in SuperIndicatorMaterial2D
   **/
  BlockIndicatorMaterial2D(BlockGeometry2D<T>&           blockGeometry,
                           int                           overlap,
                           const std::vector<int>* const materials);
  bool operator() (bool output[], const int input[]) override;
};

} // namespace olb

#endif
