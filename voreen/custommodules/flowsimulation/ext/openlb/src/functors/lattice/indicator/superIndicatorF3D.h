/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2016 Benjamin Förster, Adrian Kummerländer
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

#ifndef SUPER_INDICATOR_F_3D_H
#define SUPER_INDICATOR_F_3D_H

#include "geometry/superGeometry3D.h"
#include "indicatorBaseF3D.h"
#include "superIndicatorBaseF3D.h"

namespace olb {

template<typename T> class SuperGeometry3D;

/// SuperIndicatorF3D from IndicatorF3D
/**
 * Maintains block level BlockIndicatorFfromIndicatorF3D instances in SuperF3D::_blockF.
 **/
template <typename T>
class SuperIndicatorFfromIndicatorF3D : public SuperIndicatorF3D<T> {
protected:
  IndicatorF3D<T>& _indicatorF;
public:
  /**
   * \param indicatorF Indicator to be converted into a super indicator
   * \param geometry   Super geometry required for block indicator construction
   **/
  SuperIndicatorFfromIndicatorF3D(IndicatorF3D<T>& indicatorF, SuperGeometry3D<T>& geometry);
  bool operator() (bool output[], const int input[]) override;
};


/// Indicator functor from material numbers
/**
 * Maintains block level BlockIndicatorMaterial3D instances in SuperF3D::_blockF.
 **/
template <typename T>
class SuperIndicatorMaterial3D : public SuperIndicatorF3D<T> {
protected:
  const std::vector<int> _materialNumbers;
public:
  /**
   * \param geometry  Super geometry required for block indicator construction
   *                  and global material number queries
   * \param materials Vector of material numbers to be indicated
   **/
  SuperIndicatorMaterial3D(SuperGeometry3D<T>& geometry, std::vector<int> materials);
  /**
   * \param geometry  Super geometry required for block indicator construction
   *                  and global material number queries
   * \param materials List of material numbers to be indicated
   **/
  SuperIndicatorMaterial3D(SuperGeometry3D<T>& geometry, std::list<int> materials);
  bool operator() (bool output[], const int input[]) override;
};


} // namespace olb

#endif
