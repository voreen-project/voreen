/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2018 Adrian Kummerlaender
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

#ifndef SUPER_AVERAGE_2D_H
#define SUPER_AVERAGE_2D_H

#include "superBaseF2D.h"
#include "blockAverage2D.h"
#include "indicator/superIndicatorBaseF2D.h"
#include "utilities/functorPtr.h"

namespace olb {


/// SuperAverage2D returns the average in each component of f on a indicated subset
template <typename T, typename W = T>
class SuperAverage2D final : public SuperF2D<T,W> {
private:
  FunctorPtr<SuperF2D<T,W>>        _f;
  FunctorPtr<SuperIndicatorF2D<T>> _indicatorF;
public:
  /// Constructor for determining the average of f on a indicated subset
  /**
   * \param f          functor of which the average is to be determined
   * \param indicatorF indicator describing the subset on which to evaluate f
   **/
  SuperAverage2D(FunctorPtr<SuperF2D<T,W>>&&        f,
                 FunctorPtr<SuperIndicatorF2D<T>>&& indicatorF);
  /// Constructor for determining the average of f on a given material
  /**
   * \param f             functor of which the average is to be determined
   * \param superGeometry super geometry for constructing material indicator
   * \param material      number of the relevant material
   **/
  SuperAverage2D(FunctorPtr<SuperF2D<T,W>>&& f,
                 SuperGeometry<T,2>& superGeometry,
                 const int material);

  /// Global average operator
  /**
   * Note: While this functor exposes BlockAverage3D functors if possible, a call to
   * this function will not use them but calculate the global average by summing all
   * components and voxel counts.
   * Calling BlockAverage2D in this situation would unnecessarily complicate this as
   * we would have to weight the aggregated averages according to their share in the
   * global average.
   **/
  bool operator() (W output[], const int input[]) override;
};


}

#endif
