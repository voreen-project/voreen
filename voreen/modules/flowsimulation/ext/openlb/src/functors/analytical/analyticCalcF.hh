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

#ifndef ANALYTICAL_CALC_F_HH
#define ANALYTICAL_CALC_F_HH


#include "analyticCalcF.h"
#include "analyticalF.h"
#include "core/olbDebug.h"

namespace olb {


template <unsigned D, typename T, typename S, template<typename> class F>
AnalyticCalcF<D,T,S,F>::AnalyticCalcF(FunctorPtr<AnalyticalF<D,T,S>>&& f,
                                      FunctorPtr<AnalyticalF<D,T,S>>&& g)
  : AnalyticalF<D,T,S>(f->getTargetDim()),
    _f(std::move(f)),
    _g(std::move(g))
{
  OLB_ASSERT(g->getTargetDim() == f->getTargetDim(),
             "the dimensions of both functors need to be equal");
  std::swap(f->_ptrCalcC, this->_ptrCalcC);
  this->getName() = "(" + _f->getName() + F<T>::symbol + _g->getName() + ")";
}

template <unsigned D, typename T, typename S, template<typename> class F>
AnalyticCalcF<D,T,S,F>::AnalyticCalcF(T scalar, FunctorPtr<AnalyticalF<D,T,S>>&& g)
  : AnalyticCalcF(
      std::unique_ptr<AnalyticalF<D,T,S>>(
        new AnalyticalConst<D,T,S>(std::vector<T>(g->getTargetDim(), scalar))),
      std::forward<decltype(g)>(g))
{ }

template <unsigned D, typename T, typename S, template<typename> class F>
AnalyticCalcF<D,T,S,F>::AnalyticCalcF(FunctorPtr<AnalyticalF<D,T,S>>&& f, T scalar)
  : AnalyticCalcF(
      std::forward<decltype(f)>(f),
      std::unique_ptr<AnalyticalF<D,T,S>>(
        new AnalyticalConst<D,T,S>(std::vector<T>(f->getTargetDim(), scalar))))
{ }

template <unsigned D, typename T, typename S, template<typename> class F>
bool AnalyticCalcF<D,T,S,F>::operator()(T output[], const S input[])
{
  T outputTmp[this->_g->getTargetDim()];
  this->_g(outputTmp, input);
  this->_f(output, input);
  for (int i = 0; i < this->_f->getTargetDim(); ++i) {
    output[i] = F<T>()(output[i], outputTmp[i]);
  }
  return true;
}


template <unsigned D, typename T, typename S>
std::shared_ptr<AnalyticalF<D,T,S>> operator+(std::shared_ptr<AnalyticalF<D,T,S>> lhs, std::shared_ptr<AnalyticalF<D,T,S>> rhs)
{
  return std::shared_ptr<AnalyticalF<D,T,S>>(
           new AnalyticCalcPlus<D,T,S>(std::move(lhs), std::move(rhs)));
}

template <unsigned D, typename T, typename S>
std::shared_ptr<AnalyticalF<D,T,S>> operator+(std::shared_ptr<AnalyticalF<D,T,S>> lhs, T rhs)
{
  return std::shared_ptr<AnalyticalF<D,T,S>>(
           new AnalyticCalcPlus<D,T,S>(std::move(lhs), rhs));
}

template <unsigned D, typename T, typename S>
std::shared_ptr<AnalyticalF<D,T,S>> operator+(T lhs, std::shared_ptr<AnalyticalF<D,T,S>> rhs)
{
  return std::shared_ptr<AnalyticalF<D,T,S>>(
           new AnalyticCalcPlus<D,T,S>(lhs, std::move(rhs)));
}

template <unsigned D, typename T, typename S>
std::shared_ptr<AnalyticalF<D,T,S>> operator-(std::shared_ptr<AnalyticalF<D,T,S>> lhs, std::shared_ptr<AnalyticalF<D,T,S>> rhs)
{
  return std::shared_ptr<AnalyticalF<D,T,S>>(
           new AnalyticCalcMinus<D,T,S>(std::move(lhs), std::move(rhs)));
}

template <unsigned D, typename T, typename S>
std::shared_ptr<AnalyticalF<D,T,S>> operator-(std::shared_ptr<AnalyticalF<D,T,S>> lhs, T rhs)
{
  return std::shared_ptr<AnalyticalF<D,T,S>>(
           new AnalyticCalcMinus<D,T,S>(std::move(lhs), rhs));
}

template <unsigned D, typename T, typename S>
std::shared_ptr<AnalyticalF<D,T,S>> operator-(T lhs, std::shared_ptr<AnalyticalF<D,T,S>> rhs)
{
  return std::shared_ptr<AnalyticalF<D,T,S>>(
           new AnalyticCalcMinus<D,T,S>(lhs, std::move(rhs)));
}

template <unsigned D, typename T, typename S>
std::shared_ptr<AnalyticalF<D,T,S>> operator*(std::shared_ptr<AnalyticalF<D,T,S>> lhs, std::shared_ptr<AnalyticalF<D,T,S>> rhs)
{
  return std::shared_ptr<AnalyticalF<D,T,S>>(
           new AnalyticCalcMultiplication<D,T,S>(std::move(lhs), std::move(rhs)));
}

template <unsigned D, typename T, typename S>
std::shared_ptr<AnalyticalF<D,T,S>> operator*(std::shared_ptr<AnalyticalF<D,T,S>> lhs, T rhs)
{
  return std::shared_ptr<AnalyticalF<D,T,S>>(
           new AnalyticCalcMultiplication<D,T,S>(std::move(lhs), rhs));
}

template <unsigned D, typename T, typename S>
std::shared_ptr<AnalyticalF<D,T,S>> operator*(T lhs, std::shared_ptr<AnalyticalF<D,T,S>> rhs)
{
  return std::shared_ptr<AnalyticalF<D,T,S>>(
           new AnalyticCalcMultiplication<D,T,S>(lhs, std::move(rhs)));
}

template <unsigned D, typename T, typename S>
std::shared_ptr<AnalyticalF<D,T,S>> operator/(std::shared_ptr<AnalyticalF<D,T,S>> lhs, std::shared_ptr<AnalyticalF<D,T,S>> rhs)
{
  return std::shared_ptr<AnalyticalF<D,T,S>>(
           new AnalyticCalcDivision<D,T,S>(std::move(lhs), std::move(rhs)));
}

template <unsigned D, typename T, typename S>
std::shared_ptr<AnalyticalF<D,T,S>> operator/(std::shared_ptr<AnalyticalF<D,T,S>> lhs, T rhs)
{
  return std::shared_ptr<AnalyticalF<D,T,S>>(
           new AnalyticCalcDivision<D,T,S>(std::move(lhs), rhs));
}

template <unsigned D, typename T, typename S>
std::shared_ptr<AnalyticalF<D,T,S>> operator/(T lhs, std::shared_ptr<AnalyticalF<D,T,S>> rhs)
{
  return std::shared_ptr<AnalyticalF<D,T,S>>(
           new AnalyticCalcDivision<D,T,S>(lhs, std::move(rhs)));
}

/////////////////////////////////operator()/// ////////////////////////////////
template <unsigned D, typename T, typename S>
AnalyticalF<D,T,S>& AnalyticalF<D,T,S>::operator+(AnalyticalF<D,T,S>& rhs)
{
  auto tmp = std::make_shared< AnalyticCalcPlus<D,T,S> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}

template <unsigned D, typename T, typename S>
AnalyticalF<D,T,S>& AnalyticalF<D,T,S>::operator-(AnalyticalF<D,T,S>& rhs)
{
  auto tmp = std::make_shared< AnalyticCalcMinus<D,T,S> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}

template <unsigned D, typename T, typename S>
AnalyticalF<D,T,S>& AnalyticalF<D,T,S>::operator*(AnalyticalF<D,T,S>& rhs)
{
  auto tmp = std::make_shared< AnalyticCalcMultiplication<D,T,S> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}

template <unsigned D, typename T, typename S>
AnalyticalF<D,T,S>& AnalyticalF<D,T,S>::operator/(AnalyticalF<D,T,S>& rhs)
{
  auto tmp = std::make_shared< AnalyticCalcDivision<D,T,S> >(*this,rhs);
  this->_ptrCalcC = tmp;
  return *tmp;
}


} // end namespace olb

#endif
