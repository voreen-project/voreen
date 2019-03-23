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

#ifndef BLOCK_INDICATOR_F_3D_H
#define BLOCK_INDICATOR_F_3D_H

#include "blockIndicatorBaseF3D.h"
#include "geometry/blockGeometryView3D.h"

namespace olb {

/// BlockIndicatorF3D from IndicatorF3D
template <typename T>
class BlockIndicatorFfromIndicatorF3D : public BlockIndicatorF3D<T> {
protected:
  IndicatorF3D<T>& _indicatorF;
  Cuboid3D<T>&     _cuboid;
public:
  BlockIndicatorFfromIndicatorF3D(IndicatorF3D<T>&    indicatorF,
                                  BlockGeometry3D<T>& blockGeometry,
                                  Cuboid3D<T>&        cuboid);
  bool operator() (bool output[], const int input[]) override;
};


/// Block indicator functor from material numbers
template <typename T>
class BlockIndicatorMaterial3D : public BlockIndicatorF3D<T> {
protected:
  const int                     _overlap;
  const std::vector<int>* const _materialNumbers;
public:
  /**
   * \param blockGeometry Block geometry to be queried, accessible via SuperGeometry3D::getExtendedBlockGeometry
   * \param overlap       Overlap of given block geometry
   * \param materials     Material number vector is accepted via pointer to the actual vector in SuperIndicatorMaterial3D
   **/
  BlockIndicatorMaterial3D(BlockGeometry3D<T>&           blockGeometry,
                           int                           overlap,
                           const std::vector<int>* const materials);
  bool operator() (bool output[], const int input[]) override;
};

} // namespace olb

#endif
