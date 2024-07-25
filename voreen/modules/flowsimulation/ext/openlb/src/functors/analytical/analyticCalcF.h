/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012-2018 Lukas Baron, Tim Dornieden, Mathias J. Krause,
 *  Albert Mink, Adrian Kummerlaender
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

#ifndef ANALYTICAL_CALC_F_H
#define ANALYTICAL_CALC_F_H


#include "analyticalBaseF.h"
#include "utilities/functorPtr.h"
#include "utilities/arithmetic.h"

namespace olb {


/// arithmetic helper class for analytical functors
template <unsigned D, typename T, typename S, template<typename> class F>
class AnalyticCalcF : public AnalyticalF<D,T,S> {
protected:
  FunctorPtr<AnalyticalF<D,T,S>> _f;
  FunctorPtr<AnalyticalF<D,T,S>> _g;
public:
  AnalyticCalcF(FunctorPtr<AnalyticalF<D,T,S>>&& f,
                FunctorPtr<AnalyticalF<D,T,S>>&& g);

  AnalyticCalcF(T scalar, FunctorPtr<AnalyticalF<D,T,S>>&& g);
  AnalyticCalcF(FunctorPtr<AnalyticalF<D,T,S>>&& f, T scalar);

  bool operator() (T output[], const S input[]) override;
};

/// addition functor
template <unsigned D, typename T, typename S>
using AnalyticCalcPlus = AnalyticCalcF<D,T,S,util::plus>;

template <typename T, typename S>
using AnalyticCalcPlus1D = AnalyticCalcPlus<1,T,S>;
template <typename T, typename S>
using AnalyticCalcPlus2D = AnalyticCalcPlus<2,T,S>;
template <typename T, typename S>
using AnalyticCalcPlus3D = AnalyticCalcPlus<3,T,S>;

/// subtraction functor
template <unsigned D, typename T, typename S>
using AnalyticCalcMinus = AnalyticCalcF<D,T,S,util::minus>;

template <typename T, typename S>
using AnalyticCalcMinus1D = AnalyticCalcMinus<1,T,S>;
template <typename T, typename S>
using AnalyticCalcMinus2D = AnalyticCalcMinus<2,T,S>;
template <typename T, typename S>
using AnalyticCalcMinus3D = AnalyticCalcMinus<3,T,S>;

/// multiplication functor
template <unsigned D, typename T, typename S>
using AnalyticCalcMultiplication = AnalyticCalcF<D,T,S,util::multiplies>;

template <typename T, typename S>
using AnalyticCalcMultiplication1D = AnalyticCalcMultiplication<1,T,S>;
template <typename T, typename S>
using AnalyticCalcMultiplication2D = AnalyticCalcMultiplication<2,T,S>;
template <typename T, typename S>
using AnalyticCalcMultiplication3D = AnalyticCalcMultiplication<3,T,S>;

/// division functor
template <unsigned D, typename T, typename S>
using AnalyticCalcDivision = AnalyticCalcF<D,T,S,util::divides>;

template <typename T, typename S>
using AnalyticCalcDivision1D = AnalyticCalcDivision<1,T,S>;
template <typename T, typename S>
using AnalyticCalcDivision2D = AnalyticCalcDivision<2,T,S>;
template <typename T, typename S>
using AnalyticCalcDivision3D = AnalyticCalcDivision<3,T,S>;

template <unsigned D, typename T, typename S>
std::shared_ptr<AnalyticalF<D,T,S>> operator+(std::shared_ptr<AnalyticalF<D,T,S>> lhs, std::shared_ptr<AnalyticalF<D,T,S>> rhs);
template <unsigned D, typename T, typename S>
std::shared_ptr<AnalyticalF<D,T,S>> operator+(std::shared_ptr<AnalyticalF<D,T,S>> lhs, T rhs);
template <unsigned D, typename T, typename S>
std::shared_ptr<AnalyticalF<D,T,S>> operator+(T lhs, std::shared_ptr<AnalyticalF<D,T,S>> rhs);

template <unsigned D, typename T, typename S>
std::shared_ptr<AnalyticalF<D,T,S>> operator-(std::shared_ptr<AnalyticalF<D,T,S>> lhs, std::shared_ptr<AnalyticalF<D,T,S>> rhs);
template <unsigned D, typename T, typename S>
std::shared_ptr<AnalyticalF<D,T,S>> operator-(std::shared_ptr<AnalyticalF<D,T,S>> lhs, T rhs);
template <unsigned D, typename T, typename S>
std::shared_ptr<AnalyticalF<D,T,S>> operator-(T lhs, std::shared_ptr<AnalyticalF<D,T,S>> rhs);

template <unsigned D, typename T, typename S>
std::shared_ptr<AnalyticalF<D,T,S>> operator*(std::shared_ptr<AnalyticalF<D,T,S>> lhs, std::shared_ptr<AnalyticalF<D,T,S>> rhs);
template <unsigned D, typename T, typename S>
std::shared_ptr<AnalyticalF<D,T,S>> operator*(std::shared_ptr<AnalyticalF<D,T,S>> lhs, T rhs);
template <unsigned D, typename T, typename S>
std::shared_ptr<AnalyticalF<D,T,S>> operator*(T lhs, std::shared_ptr<AnalyticalF<D,T,S>> rhs);

template <unsigned D, typename T, typename S>
std::shared_ptr<AnalyticalF<D,T,S>> operator/(std::shared_ptr<AnalyticalF<D,T,S>> lhs, std::shared_ptr<AnalyticalF<D,T,S>> rhs);
template <unsigned D, typename T, typename S>
std::shared_ptr<AnalyticalF<D,T,S>> operator/(std::shared_ptr<AnalyticalF<D,T,S>> lhs, T rhs);
template <unsigned D, typename T, typename S>
std::shared_ptr<AnalyticalF<D,T,S>> operator/(T lhs, std::shared_ptr<AnalyticalF<D,T,S>> rhs);


} // end namespace olb

#endif
