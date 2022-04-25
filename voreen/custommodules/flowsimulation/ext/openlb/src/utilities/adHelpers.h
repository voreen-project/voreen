/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Julius Jeßberger
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

/** \file
 * Some helper functions for the ADf data type.
 */

#ifndef AD_HELPERS_H
#define AD_HELPERS_H

#include "aDiff.h"


namespace olb {

namespace util {

/** The variables of an array are set to be the differential variables.
 * --> d/di a_j = 1 if (i==j) else 0.
 */
template <typename SAD>
void iniDiagonal(SAD* a, int dim) {
  for (int i = 0; i < dim; ++i) {
    a[i].setDiffVariable(i);
  }
}

template <typename SAD>
void iniDiagonal(std::vector<SAD>& a) {
  for (unsigned i = 0; i < a.size(); ++i) {
    a[i].setDiffVariable(i);
  }
}


/** Copy the derivatives from an ADf array into an array.
 * \result target[j*dim+i] = source[j].d(i) = d/dx_i (s_j)
 * \param length = length of the array
 */
template <typename T, typename TAD>
void copyDerivatives(T* target, const TAD* source, int length) {
  for (int j = 0; j < length; ++j) {
    for (unsigned i = 0; i < TAD::dim; ++i) {
      target[j*(TAD::dim)+i] = source[j].d(i);
    }
  }
}


/// copy vector with specified typecast
template<typename T, typename S>
std::vector<T> copyAs(const std::vector<S>& input) {
  std::vector<T> result(input.size(), T());

  for(std::size_t it = 0; it < input.size(); ++it) {
    result[it] = T(input[it]);
  }
  return result;
}


/// copy value and initialize derivative
template<typename S>
ADf<S,1> iniAD(const S source) {
  ADf<S,1> result(source);
  result.setDiffVariable(0);
  return result;
}

/// copy array values and initialize derivatives
template<unsigned n, typename S>
ADf<S,n>* iniAD(const S* source) {
  ADf<S,n>* result = new ADf<S,n>[n];
  util::copyN (result, source, n);
  util::iniDiagonal (result, n);
  return result;
}

/// copy values and initialize derivatives
template<unsigned n, typename S>
std::vector<ADf<S,n>> iniAD(const std::vector<S>& source) {
  auto result = copyAs<ADf<S,n>,S> (source);
  util::iniDiagonal (result);
  return result;
}



/** Compute derivatives of a function f: S -> T via forward AD
 * Signature is "V f(U)" so f is expected to accept and return a single value.
 * It must be (captured by) a lambda function that handles arbitrary types,
 * so "auto f(auto x)" is to be written.
 */
template <typename S, typename F>
auto derivativeFAD (F f, const S input) {
  ADf<S,1> inputAD = input;
  inputAD.setDiffVariable(0);

  return f(inputAD).d(0);
}

/** Compute derivatives of a function f: S^sourceDIM -> T via forward AD
 * Signature is "V f(U*)" so f is expected to return a single value.
 * It must be (captured by) a lambda function that handles arbitrary types,
 * so "auto f(auto x)" is to be written.
 */
template <unsigned sourceDIM, typename S, typename F>
auto derivativeFAD (F f, const S* input) {
  using SAD = ADf<S,sourceDIM>;
  using T = decltype(f(input));
  using TAD = ADf<T,sourceDIM>;

  SAD inputAD[sourceDIM];
  util::copyN(inputAD, input, sourceDIM);
  util::iniDiagonal(inputAD, sourceDIM);

  TAD resultAD = f(inputAD);
  T* result = new T[sourceDIM];
  copyDerivatives(result, &resultAD, 1);
  return result;
}

/** Compute derivatives of a function f: S^sourceDIM -> T^targetDIM via forward AD
 * Signature of f is "void f(V*, U*)", it is expected to modify the first argument
 * in-place.
 * f must be (captured by) a lambda function that handles arbitrary types,
 * so "void f(auto y, auto x)" is to be written.
 */
template <unsigned sourceDIM, typename T, typename S, typename F>
void derivativeFAD (F f, T* output, const S* input, unsigned targetDIM) {
  using SAD = ADf<S,sourceDIM>;
  using TAD = ADf<T,sourceDIM>;

  SAD inputAD[sourceDIM];
  util::copyN(inputAD, input, sourceDIM);
  util::iniDiagonal(inputAD, sourceDIM);
  TAD outputAD[targetDIM];

  f(outputAD, inputAD);
  copyDerivatives(output, outputAD, targetDIM);
  return;
}

}

}

#endif
