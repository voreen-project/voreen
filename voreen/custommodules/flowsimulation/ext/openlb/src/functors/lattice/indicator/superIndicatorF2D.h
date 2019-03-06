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

#ifndef SUPER_INDICATOR_F_2D_H
#define SUPER_INDICATOR_F_2D_H

#include "superIndicatorBaseF2D.h"
#include "geometry/superGeometry2D.h"

namespace olb {


/// SuperIndicatorF2D from IndicatorF2D
/**
 * Maintains block level BlockIndicatorFfromIndicatorF2D instances in SuperF2D::_blockF.
 **/
template <typename T>
class SuperIndicatorFfromIndicatorF2D : public SuperIndicatorF2D<T> {
protected:
  IndicatorF2D<T>& _indicatorF;
public:
  /**
   * \param indicatorF Indicator to be converted into a super indicator
   * \param geometry   Super geometry required for block indicator construction
   **/
  SuperIndicatorFfromIndicatorF2D(IndicatorF2D<T>& indicatorF, SuperGeometry2D<T>& geometry);
  bool operator() (bool output[], const int input[]) override;
};

/// Indicator functor from material numbers
/**
 * Maintains block level BlockIndicatorMaterial2D instances in SuperF2D::_blockF.
 **/
template <typename T>
class SuperIndicatorMaterial2D : public SuperIndicatorF2D<T> {
protected:
  SuperGeometry2D<T>&    _superGeometry;
  const std::vector<int> _materialNumbers;
public:
  /**
   * \param geometry  Super geometry required for block indicator construction
   *                  and global material number queries
   * \param materials Vector of material numbers to be indicated
   **/
  SuperIndicatorMaterial2D(SuperGeometry2D<T>& geometry, std::vector<int> materials);
  bool operator() (bool output[], const int input[]) override;
};


} // namespace olb

#endif
