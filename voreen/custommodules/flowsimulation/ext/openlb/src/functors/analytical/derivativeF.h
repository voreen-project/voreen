/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Julius Jessberger, Mathias J. Krause
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


#ifndef DERIVATIVE_F_H
#define DERIVATIVE_F_H


#include "analyticalBaseF.h"
#include "../../utilities/aDiff.h"
#include "../../utilities/adHelpers.h"

#include <stdlib.h>




namespace olb {


/// Class for computing the derivative of a given 1D functor with a finite difference
template <typename T>
class AnalyticalDerivativeFD1D : public AnalyticalF1D<T,T> {
protected:
  AnalyticalF1D<T,T>& _f;
  T _eps {1.e-10};
public:
  explicit AnalyticalDerivativeFD1D(AnalyticalF1D<T,T>& f) : AnalyticalF1D<T,T>(1), _f(f) { };
  AnalyticalDerivativeFD1D(AnalyticalF1D<T,T>& f, T eps) : AnalyticalF1D<T,T>(1), _f(f), _eps(eps) { };
  bool operator() (T output[], const T input[]) override;
};

/// Finite Difference computation
template <typename T>
bool AnalyticalDerivativeFD1D<T>::operator() (T output[], const T input[])
{
  _f(output,input);
  T x = output[0];
  T input2[1];
  input2[0] = input[0] + _eps;
  _f(output,input2);
  output[0] -= x;
  output[0] /= _eps;
  return true;
}

// Conversion to old implementation
template class AnalyticalDerivativeFD1D<double>;




/// Class for AD Differentiation of 1-dim Functor F: S -> T.
template <unsigned D, typename T, typename S>
class AnalyticalDerivativeAD1D : public AnalyticalF<D,T,S> {
public:
  using SAD = util::ADf<S, 1>;
  using TAD = util::ADf<T, 1>;

  explicit AnalyticalDerivativeAD1D(AnalyticalF<D,TAD,SAD>& f) : AnalyticalF<D,T,S>(f.getTargetDim()), _f(f) { };
  bool operator() (T output[], const S input[]) override;

protected:
  AnalyticalF<D,TAD,SAD>& _f;   // the functor which is to be differentiated
};

/// AD Differentiation operator implementation
template <unsigned D, typename T, typename S>
bool AnalyticalDerivativeAD1D<D,T,S>::operator() (T output[], const S input[]) {
  // Define corresponding AD-data types
  SAD inpAD(input[0]);
  inpAD.setDiffVariable(0);
  TAD outAD(output[0]);

  // Apply f to the AD versions of input, output.
  _f(&outAD, &inpAD);

  // Return the derivative.
  output[0] = outAD.d(0);

  return true;
}


template<
  class F,
  typename T, typename S, unsigned sourceDIM,
  typename... ARGS
>
class AnalyticalDerivativeAD : public AnalyticalF<sourceDIM,T,S> {

public:
  using SAD = util::ADf<S,sourceDIM>;
  using TAD = util::ADf<T,sourceDIM>;

protected:
  typename F::template exchange_type<TAD,SAD> _f;

public:
  explicit AnalyticalDerivativeAD(unsigned targetDIM, ARGS&&... args)
   : AnalyticalF<sourceDIM,T,S>(targetDIM*sourceDIM), _f(args...) { }

  explicit AnalyticalDerivativeAD(const F& f)
   : AnalyticalF<sourceDIM,T,S>(f.getTargetDim()*sourceDIM),
     _f(f.template copyAs<TAD,SAD>()) { }

  /** Computes the partial derivatives of _f and writes them into output.
   * \result output[j*sourceDIM+i] = d/dx_i F_j
   * --> output is an arithmetic array of length targetDIM*sourceDIM
   */
  bool operator()(T output[], const S input[]) override {
    SAD inputAD[sourceDIM];
    util::copyN(inputAD, input, sourceDIM);
    util::iniDiagonal(inputAD, sourceDIM);
    TAD outputAD[_f.getTargetDim()];

    _f(outputAD, inputAD);

    util::copyDerivatives(output, outputAD, _f.getTargetDim());

    return true;
  }
};

/* // Note: this cannot work since the mapping AFGS -> F is not unique
template <class F, typename... ARGS>
AnalyticalDerivativeAD(unsigned, ARGS&&...) -> AnalyticalDerivativeAD <F,
  typename F::targetType, typename F::sourceType, F::dim, ARGS...>;
*/

template <class F>
AnalyticalDerivativeAD(const F&) -> AnalyticalDerivativeAD <F,
  typename F::targetType, typename F::sourceType, F::dim>;



}  // end namespace olb

#endif
