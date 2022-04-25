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

#ifndef SUPER_INTEGRAL_F_2D_H
#define SUPER_INTEGRAL_F_2D_H

#include "functors/lattice/superBaseF2D.h"
#include "functors/lattice/indicator/superIndicatorBaseF2D.h"
#include "functors/analytical/indicator/indicatorBaseF2D.h"
#include "functors/lattice/indicator/superIndicatorF2D.h"
#include "geometry/superGeometry.h"
#include "utilities/functorPtr.h"

namespace olb {


/// SuperSum2D sums all components of f over a indicated subset
template <typename T, typename W = T>
class SuperSum2D final : public SuperF2D<T,W> {
private:
  FunctorPtr<SuperF2D<T,W>>        _f;
  FunctorPtr<SuperIndicatorF2D<T>> _indicatorF;
public:
  /// Constructor for summing f on a indicated subset
  /**
   * \param f          functor to be summed
   * \param indicatorF indicator describing the subset on which to evaluate f
   **/
  SuperSum2D(FunctorPtr<SuperF2D<T,W>>&&        f,
             FunctorPtr<SuperIndicatorF2D<T>>&& indicatorF);
  /// Constructor for summing f on a given material
  /**
   * \param f             functor to be summed
   * \param superGeometry super geometry for constructing material indicator
   * \param material      number of the relevant material
   **/
  SuperSum2D(FunctorPtr<SuperF2D<T,W>>&& f,
             SuperGeometry<T,2>& superGeometry,
             const int material);

  bool operator() (W output[], const int input[]) override;
};


/// SuperIntegral2D integrates f on a indicated subset
template <typename T, typename W = T>
class SuperIntegral2D final : public SuperF2D<T,W> {
private:
  FunctorPtr<SuperF2D<T,W>>        _f;
  FunctorPtr<SuperIndicatorF2D<T>> _indicatorF;
public:
  /// Constructor for integrating f on a indicated subset
  /**
   * \param f          functor to be integrated
   * \param indicatorF indicator describing the subset on which to evaluate f
   **/
  SuperIntegral2D(FunctorPtr<SuperF2D<T,W>>&&        f,
                  FunctorPtr<SuperIndicatorF2D<T>>&& indicatorF);
  /// Constructor for integrating f on a given material
  /**
   * \param f             functor to be integrated
   * \param superGeometry super geometry for constructing material indicator
   * \param material      number of the relevant material
   **/
  SuperIntegral2D(FunctorPtr<SuperF2D<T,W>>&& f,
                  SuperGeometry<T,2>& superGeometry,
                  const int material);

  bool operator() (W output[], const int input[]) override;
};


}

#endif
