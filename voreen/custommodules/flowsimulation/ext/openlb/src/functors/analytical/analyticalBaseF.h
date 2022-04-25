/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Lukas Baron, Tim Dornieden, Mathias J. Krause, Albert Mink
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

#ifndef ANALYTICAL_BASE_F_H
#define ANALYTICAL_BASE_F_H

#include "functors/genericF.h"
#include "functors/analytical/indicator/indicatorBaseF3D.h"

/**
 *  The functor dimensions are given by F: S^m -> T^n  (S=source, T=target)
 *  and are implemented via GenericF(n,m).
 *  Don't get confused by the flipped order of source and target.
 */

namespace olb {

////////////////////////////////////////////////////////////////////////////////
// 2nd level classes
// note: for LatticeFunctions the number indicates the SOURCE dimension,
//       target dim depends on return variable type, so std::vector<T> is used

template<unsigned D, typename T, typename S> class AnalyticalIdentity;

/// AnalyticalF are applications from DD to XD, where X is set by the constructor.
template<unsigned D, typename T, typename S> class AnalyticalF : public GenericF<T,S> {
protected:
  // n denotes the target dimension
  AnalyticalF(int n): GenericF<T,S>(n,D)
  {
    static_assert(D==1 || D==2 || D==3, "Only D=1,2,3 allowed.");
  }
public:
  using identity_functor_type = AnalyticalIdentity<D,T,S>;
  static constexpr unsigned dim = D;

  AnalyticalF<D,T,S>& operator-(AnalyticalF<D,T,S>& rhs);
  AnalyticalF<D,T,S>& operator+(AnalyticalF<D,T,S>& rhs);
  AnalyticalF<D,T,S>& operator*(AnalyticalF<D,T,S>& rhs);
  AnalyticalF<D,T,S>& operator/(AnalyticalF<D,T,S>& rhs);
};

/// AnalyticalIdentity stores vectors, result of addition,multiplication, ...
template <unsigned D, typename T, typename S>
class AnalyticalIdentity : public AnalyticalF<D,T,S> {
protected:
  AnalyticalF<D,T,S>& _f;
public:
  AnalyticalIdentity(AnalyticalF<D,T,S>& f);
  bool operator() (T output[], const S input[]) override;
};

////////////// CONVERSION FROM NEW TO OLD IMPLEMENTATION //////////////////////////

template <typename T, typename S>
using AnalyticalF1D = AnalyticalF<1,T,S>;
template <typename T, typename S>
using AnalyticalF2D = AnalyticalF<2,T,S>;
template <typename T, typename S>
using AnalyticalF3D = AnalyticalF<3,T,S>;

template <typename T, typename S>
using AnalyticalIdentity1D = AnalyticalIdentity<1,T,S>;
template <typename T, typename S>
using AnalyticalIdentity2D = AnalyticalIdentity<2,T,S>;
template <typename T, typename S>
using AnalyticalIdentity3D = AnalyticalIdentity<3,T,S>;


/// Converts IndicatorF to AnalyticalF (used for Analytical operands for Identity)
template <typename T, typename S>
class AnalyticalFfromIndicatorF3D : public AnalyticalF3D<T,S> {
protected:
  IndicatorF3D<T>& _indicatorF;
public:
  AnalyticalFfromIndicatorF3D(IndicatorF3D<T>& indicatorF);
  bool operator() (T output[], const S input[]) override;
};



} // end namespace olb

#endif
