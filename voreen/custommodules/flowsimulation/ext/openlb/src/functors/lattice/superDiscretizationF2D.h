/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2019 Tom Braun
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

#ifndef SUPER_DISCRETIZATION_F_2D_H
#define SUPER_DISCRETIZATION_F_2D_H

#include "utilities/functorPtr.h"
#include "superBaseF2D.h"

namespace olb {

/**
 * Super functor for discretizing values by an interval (bottomBoundary,topBoundary),
 * as well as restricting the value by setting n equal-distributed points and rounding the value to the nearest point
 * If n = 1, there won't be restricting, and for n>=1 there will be n-1 restricting points.
 **/
template <typename T, typename W=T>
class SuperDiscretizationF2D: public SuperF2D<T, W> {
private:
  FunctorPtr<SuperF2D<T,W>>& _f;
  const T _bottomBoundary;
  const T _topBoundary;
  const int _n;
public:

  /**
    * \param f               Block functor
    * \param bottomBoundary     bottom border of the interval the value will be restricted to
    * \param topBoundary     top border of the interval the value will be restricted to
    * \param n               number of points the interval will be split by, and the value discretizied to
    **/

  SuperDiscretizationF2D(FunctorPtr<SuperF2D<T, W>>&& f,
                         T bottomBoundary,
                         T topBoundary,
                         int n = 1 );
};

}
#endif
