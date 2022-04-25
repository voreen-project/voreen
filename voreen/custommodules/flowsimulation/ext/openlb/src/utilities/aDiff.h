/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012-2021 Mathias J. Krause, Jan E. Marquardt, Luca Heim
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
 * The description of a algoritmic differentiation data type
 * using the forward method -- header file.
 */

#ifndef A_DIFF_H
#define A_DIFF_H

#include <cstdlib>
#include <type_traits>
#include <array>
#include <iostream>
#include <cassert>

#include "core/baseType.h"
#include "core/vector.h"
#include "core/util.h"

struct AD { };

// All OpenLB code is contained in this namespace.
namespace olb {

namespace util {

/// Definition of a description of a algoritmic differentiation
/// data type using the forward method
/** An ADf is a data type which enables the calculation of
 * derivatives by means of algorithmic differentiation.
 *
 * This class is not intended to be derived from.
 */

// Used to enable ADf instances in scalar-vector operations
/**
 * See olb::meta::is_arithmetic
 **/

template <class T, unsigned DIM>
class ADf final : public AD {
private:
  inline void checkDataType();
public:
  using base_t = T;
  static constexpr unsigned dim = DIM;

  /// value
  T _v = T();
  /// derivatives
  Vector<T,DIM> _d = ( T( 0 ) );

  // class functions
  inline ADf();

  //explicit ADf(const T& v);
  inline ADf(const ADf& a);
  inline ADf(const ADf& a, bool b);
  inline ADf(const T& v);
  inline ADf(const T& v, const Vector<T,DIM>& d);

  ADf(ADf&& a);

  inline T& v();
  inline T& d(unsigned i);
  inline const T& d(unsigned i) const;
  inline Vector<T,DIM>& d();

  inline ADf& operator = (const ADf& a);
  inline ADf& operator = (ADf&& a);
  inline ADf& operator = (const T& v);

  inline ADf& operator += (const ADf& a);
  inline ADf& operator += (const T& v);

  inline ADf& operator -= (const ADf& a);
  inline ADf& operator -= (const T& v);

  inline ADf& operator *= (const ADf& a);
  inline ADf& operator *= (const T& v);

  inline ADf& operator /= (const ADf& a);
  inline ADf& operator /= (const T& v);

  inline void setDiffVariable(unsigned iD);

  inline bool hasZeroDerivative() const;
  inline operator base_t() const;
};

// class functions
template <class T, unsigned DIM>
inline ADf<T,DIM>::ADf():ADf<T,DIM>( T() )
{ }

template <class T, unsigned DIM>
inline ADf<T,DIM>::ADf(const ADf& a):ADf<T,DIM>(a._v, a._d)
{ }

template <class T, unsigned DIM>
inline ADf<T,DIM>::ADf(const T& v):ADf<T,DIM>(v,Vector<T,DIM> {})
{ }

template <class T, unsigned DIM>
inline ADf<T,DIM>::ADf(const T& v, const Vector<T,DIM>& d): _v(v), _d(d)
{
  checkDataType();
}

template <class T, unsigned DIM>
inline ADf<T,DIM>::ADf(ADf&& a):ADf<T,DIM>(a._v, a._d)
{ }

/// Check whether the base type is of floating point type
template <class T, unsigned DIM>
inline void ADf<T,DIM>::checkDataType()
{
  static_assert(std::is_floating_point<T>::value,
                "Use floating point types (float, double or long double). Do not use integer types.");
}

template <class T, unsigned DIM>
inline ADf<T,DIM>& ADf<T,DIM>::operator= (const ADf& a)
{
  _v=a._v;
  _d=a._d;
  return *this;
}

template <class T, unsigned DIM>
inline ADf<T,DIM>& ADf<T,DIM>::operator= (ADf&& a)
{
  _v=a._v;
  _d=a._d;
  return *this;
}

template <class T, unsigned DIM>
inline ADf<T,DIM>& ADf<T,DIM>::operator= (const T& v)
{
  _v=v;
  _d = ( T() );
  return *this;
}

template <class T, unsigned DIM>
inline ADf<T,DIM>& ADf<T,DIM>::operator += (const ADf<T,DIM>& a)
{
  _v+=a._v;
  _d+=a._d;
  return *this;
}

template <class T, unsigned DIM>
inline ADf<T,DIM>& ADf<T,DIM>::operator += (const T& v)
{
  _v+=v;
  return *this;
}

template <class T, unsigned DIM>
inline ADf<T,DIM>& ADf<T,DIM>::operator -= (const ADf<T,DIM>& a)
{
  _v-=a._v;
  _d-=a._d;
  return *this;
}

template <class T, unsigned DIM>
inline ADf<T,DIM>& ADf<T,DIM>::operator -= (const T& v)
{
  _v-=v;
  return *this;
}

template <class T, unsigned DIM>
inline ADf<T,DIM>& ADf<T,DIM>::operator *= (const ADf<T,DIM>& a)
{
  _d=_d*a._v+_v*a._d;
  _v*=a._v;
  return *this;
}

template <class T, unsigned DIM>
inline ADf<T,DIM>& ADf<T,DIM>::operator *= (const T& v)
{
  _v*=v;
  _d*=v;
  return *this;
}

template <class T, unsigned DIM>
inline ADf<T,DIM>& ADf<T,DIM>::operator /= (const ADf<T,DIM>& a)
{
  T tmp(T(1)/a._v);
  _v*=tmp;
  _d=tmp*(_d-_v*a._d);
  return *this;
}

template <class T, unsigned DIM>
inline ADf<T,DIM>& ADf<T,DIM>::operator /= (const T& v)
{
  T tmp(T(1)/v);
  _v*=tmp;
  _d*=tmp;
  return *this;
}

template <class T, unsigned DIM>
inline T& ADf<T,DIM>::v()
{
  return _v;
}

template <class T, unsigned DIM>
inline T& ADf<T,DIM>::d(unsigned i)
{
  return _d[i];
}

template <class T, unsigned DIM>
inline const T& ADf<T,DIM>::d(unsigned i) const
{
  return _d[i];
}

template <class T, unsigned DIM>
inline Vector<T,DIM>& ADf<T,DIM>::d()
{
  return _d;
}

template <class T, unsigned DIM>
inline void ADf<T,DIM>::setDiffVariable(unsigned iD)
{
  _d =( T(0) );
  _d[iD] = T(1);
}

template <class T, unsigned DIM>
inline std::ostream& operator << (std::ostream& os, const ADf<T,DIM>& o)
{
  os << o._v;
  if constexpr (DIM>0) {
    os << o._d;
  }
  return os;
}

// WARNING: Reads only the value, not the derivatives
template <class T, unsigned DIM>
inline std::istream& operator >> (std::istream& is, ADf<T,DIM>& in)
{
  is >> in._v;
  in._d = T{0};
  return is;
}


// indirect helper functions, providing general arithmetics
// overloading of operators for interaction with BaseType and itself
// addition, subtraction, multiplication and division

template <class T, unsigned DIM, typename U, typename std::enable_if<std::is_integral<U>::value | std::is_floating_point<U>::value,int>::type = 0>
inline ADf<T,DIM> operator+ (const U& a, const ADf<T,DIM>& b)
{
  return ADf<T,DIM>(b) += olb::BaseType<ADf<T,DIM>>(a);
}

template <class T, unsigned DIM, typename U, typename std::enable_if<std::is_integral<U>::value | std::is_floating_point<U>::value,int>::type = 0>
inline ADf<T,DIM> operator+ (const ADf<T,DIM>& a,const U& b)
{
  return ADf<T,DIM>(a) += olb::BaseType<ADf<T,DIM>>(b);
}

template <class T, unsigned DIM>
inline ADf<T,DIM> operator+ (const ADf<T,DIM>& a, const ADf<T,DIM>& b)
{
  return ADf<T,DIM>(a) += b;
}

template <class T, unsigned DIM, typename U, typename std::enable_if<std::is_integral<U>::value | std::is_floating_point<U>::value,int>::type = 0>
inline ADf<T,DIM> operator- (const U& a, const ADf<T,DIM>& b)
{
  return ADf<T,DIM>(-ADf<T,DIM>(b) + olb::BaseType<ADf<T,DIM>>(a));
}

template <class T, unsigned DIM, typename U, typename std::enable_if<std::is_integral<U>::value | std::is_floating_point<U>::value,int>::type = 0>
inline ADf<T,DIM> operator- (const ADf<T,DIM>& a,const U& b)
{
  return ADf<T,DIM>(a) -= olb::BaseType<ADf<T,DIM>>(b);
}

template <class T, unsigned DIM>
inline ADf<T,DIM> operator- (const ADf<T,DIM>& a, const ADf<T,DIM>& b)
{
  return ADf<T,DIM>(a)-=b;
}

template <class T, unsigned DIM, typename U, typename std::enable_if<std::is_integral<U>::value | std::is_floating_point<U>::value,int>::type = 0>
inline ADf<T,DIM> operator* (const U& a, const ADf<T,DIM>& b)
{
  return ADf<T,DIM>(b) *= olb::BaseType<ADf<T,DIM>>(a);
}

template <class T, unsigned DIM, typename U, typename std::enable_if<std::is_integral<U>::value | std::is_floating_point<U>::value,int>::type = 0>
inline ADf<T,DIM> operator* (const ADf<T,DIM>& a,const U& b)
{
  return ADf<T,DIM>(a) *= olb::BaseType<ADf<T,DIM>>(b);
}

template <class T, unsigned DIM>
inline ADf<T,DIM> operator* (const ADf<T,DIM>& a, const ADf<T,DIM>& b)
{
  return ADf<T,DIM>(a) *= b;
}

template <class T, unsigned DIM, typename U, typename std::enable_if<std::is_integral<U>::value | std::is_floating_point<U>::value,int>::type = 0>
inline ADf<T,DIM> operator/ (const U& a, const ADf<T,DIM>& b)
{
  return ADf<T,DIM>(a) /= b;
}

template <class T, unsigned DIM, typename U, typename std::enable_if<std::is_integral<U>::value | std::is_floating_point<U>::value,int>::type = 0>
inline ADf<T,DIM> operator/ (const ADf<T,DIM>& a,const U& b)
{
  return ADf<T,DIM>(a) /= olb::BaseType<ADf<T,DIM>>(b);
}


template <class T, unsigned DIM>
inline ADf<T,DIM> operator/ (const ADf<T,DIM>& a, const ADf<T,DIM>& b)
{
  return ADf<T,DIM>(a) /= b;
}

/* Unary operators */
template <class T, unsigned DIM>
inline ADf<T,DIM> operator+ (const ADf<T,DIM>& a)
{
  return ADf<T,DIM>(a);
}

template <class T, unsigned DIM>
inline ADf<T,DIM> operator- (const ADf<T,DIM>& a)
{
  return ADf<T,DIM>(a) *= -1;
}










// General arithmetics for interaction with objects of the first template
// parameter. This enables computations of second derivatives by pyramiding
// ADf.
// addition, subtraction, multiplication and division
template <class T, unsigned DIM>
inline ADf<T,DIM> add1(const T& a, const ADf<T,DIM>& b)
{
  ADf<T,DIM> c(a+b._v, b._d);
  return c;
}

template <class T, unsigned DIM>
inline ADf<T,DIM> add2(const ADf<T,DIM>& a,const T& b)
{
  ADf<T,DIM> c(a._v+b, a._d);
  return c;
}

template <class T, unsigned DIM>
inline ADf<T,DIM> sub1(const T& a, const ADf<T,DIM>& b)
{
  ADf<T,DIM> c(a-b._v, b._d);
  c._d *= -1;
  return c;
}

template <class T, unsigned DIM>
inline ADf<T,DIM> sub2(const ADf<T,DIM>& a,const T& b)
{
  ADf<T,DIM> c(a._v-b, a._d);
  return c;
}

template <class T, unsigned DIM>
inline ADf<T,DIM> mul1(const T& a, const ADf<T,DIM>& b)
{
  ADf<T,DIM> c(a*b._v, b._d);
  c._d *= a;
  return c;
}

template <class T, unsigned DIM>
inline ADf<T,DIM> mul2(const ADf<T,DIM>& a,const T& b)
{
  ADf<T,DIM> c(a._v*b, a._d);
  c._d*=b;
  return c;
}

template <class T, unsigned DIM>
inline ADf<T,DIM> di_v1(const T& a, const ADf<T,DIM>& b)
{
  ADf<T,DIM> c(a/b._v, b._d);
  T tmp(-c._v/b._v);
  c._d*=tmp;
  return c;
}

template <class T, unsigned DIM>
inline ADf<T,DIM> di_v2(const ADf<T,DIM>& a,const T& b)
{
  ADf<T,DIM> c(a._v/b, a._d);
  T tmp = T(1)/b;
  c._d*=tmp;
  return c;
}

template <class T, unsigned DIM>
inline ADf<T,DIM> pow1(const T& a, const ADf<T,DIM>& b)
{
  ADf<T,DIM> c(pow(a,b._v), b._d);
  T tmp(b._v*pow(a,b._v-T(1)));//,tmp1(c._v*log10(a));
  c._d*=tmp;
  return c;
}
template <class T, unsigned DIM>
inline ADf<T,DIM> pow2(const ADf<T,DIM>& a, const T& b)
{
  ADf<T,DIM> c(pow(a._v,b), a._d);
  T tmp(b*pow(a._v,b-T(1)));//,tmp1(c._v*log10(a._v));
  c._d*=tmp;
  return c;
}










// extended math functions:
// pow, sqr, exp, log, log10, sqrt, sin, cos, tan, asin, acos, atan
template <class T, unsigned DIM, typename U, typename std::enable_if<std::is_integral<U>::value | std::is_floating_point<U>::value,int>::type = 0>
inline ADf<T,DIM> pow (const U& a, const ADf<T,DIM>& b)
{
  ADf<T,DIM> c(std::pow(T(a),b._v), b._d);
  T tmp(c._v*std::log(T(a)));
  c._d*=tmp;
  return c;
}

template <class T, unsigned DIM, typename U, typename std::enable_if<std::is_integral<U>::value | std::is_floating_point<U>::value,int>::type = 0>
inline ADf<T,DIM> pow (const ADf<T,DIM>& a,const U& b)
{
  ADf<T,DIM> c(std::pow(a._v,T(b)), a._d);
  T tmp(T(b)*std::pow(a._v,b-T(1)));
  c._d*=tmp;
  return c;
}


template <class T, unsigned DIM>
inline ADf<T,DIM> pow (const ADf<T,DIM>& a, const ADf<T,DIM>& b)
{
  if (!util::nearZero(a._v)) {
    return exp(b*log(a));
  }
  else {
    ADf<T,DIM> c(a._v, Vector<T,DIM> {std::numeric_limits<T>::quiet_NaN()});
#ifdef AdWarnings
    std::cout << "ADf WARNING: pow(Adf) - pow evaluated at non-differentiable point" << std::endl;
#endif
    if (std::isfinite(c._v)) {
      for (unsigned i = 0; i < DIM; ++i) {
        if (util::nearZero(a.d(i))) {
          c._d[i] = T(0);
        }
      }
    }
    return c;
  }
}

template <class T, unsigned DIM>
inline ADf<T,DIM> pow (const ADf<T,DIM>& a, int b)
{
  ADf<T,DIM> c(std::pow(a._v,b), a._d);
  T tmp(b*std::pow(a._v, b-1));
  c._d*=tmp;
  return c;
}

template <class T, unsigned DIM>
inline ADf<T,DIM> sqr (const ADf<T,DIM>& a)
{
  ADf<T,DIM> c((a._v)*(a._v), a._d);
  T tmp(T(2)*a._v);
  c._d*=tmp;
  return c;
}

template <class T, unsigned DIM>
inline ADf<T,DIM> exp (const ADf<T,DIM>& a)
{
  ADf<T,DIM> c(std::exp(a._v), a._d);
  c._d*=c._v;
  return c;
}

template <unsigned DIM>
inline ADf<float,DIM> expf (const ADf<float,DIM>& a)
{
  return exp(a);
}

template <unsigned DIM>
inline ADf<long double,DIM> expl (const ADf<long double,DIM>& a)
{
  return exp(a);
}

template <class T, unsigned DIM>
inline ADf<T,DIM> log (const ADf<T,DIM>& a)
{
  ADf<T,DIM> c(std::log(a._v), a._d);
  T tmp (T(1)/a._v);
  c._d*=tmp;
  return c;
}

template <unsigned DIM>
inline ADf<float,DIM> logf (const ADf<float,DIM>& a)
{
  return log(a);
}

template <unsigned DIM>
inline ADf<long double,DIM> logl (const ADf<long double,DIM>& a)
{
  return log(a);
}

template <class T, unsigned DIM>
inline ADf<T,DIM> log10 (const ADf<T,DIM>& a)
{
  ADf<T,DIM> c(std::log10(a._v), a._d);
  T tmp (T(1)/(a._v*std::log(T(10))));
  c._d*=tmp;
  return c;
}

template <unsigned DIM>
inline ADf<float,DIM> log10f (const ADf<float,DIM>& a)
{
  return log10(a);
}

template <unsigned DIM>
inline ADf<long double,DIM> log10l (const ADf<long double,DIM>& a)
{
  return log10(a);
}

template <class T, unsigned DIM>
inline ADf<T,DIM> log2 (const ADf<T,DIM>& a)
{
  ADf<T,DIM> c(std::log2(a._v), a._d);
  T tmp (T(1)/(a._v*std::log(T(2))));
  c._d*=tmp;
  return c;
}

template <unsigned DIM>
inline ADf<float,DIM> log2f (const ADf<float,DIM>& a)
{
  return log2(a);
}

template <unsigned DIM>
inline ADf<long double,DIM> log2l (const ADf<long double,DIM>& a)
{
  return log2(a);
}

template <class T, unsigned DIM>
inline ADf<T,DIM> log1p (const ADf<T,DIM>& a)
{
  return log(a+T(1));
}

template <unsigned DIM>
inline ADf<float,DIM> log1pf (const ADf<float,DIM>& a)
{
  return log1p(a);
}

template <unsigned DIM>
inline ADf<long double,DIM> log1pl (const ADf<long double,DIM>& a)
{
  return log1p(a);
}

template <class T, unsigned DIM>
inline ADf<T,DIM> sqrt (const ADf<T,DIM>& a)
{
  ADf<T,DIM> c(std::sqrt(a._v), a._d);
  T tmp(T(1.)/(c._v*T(2)));
  for (unsigned i = 0; i < DIM; ++i) {
    if (!util::nearZero(a.d(i))) {
      c._d[i] *= tmp;
    }
  }
  return c;
}

template <class T, unsigned DIM>
inline ADf<T,DIM> sin (const ADf<T,DIM>& a)
{
  ADf<T,DIM> c(std::sin(a._v), a._d);
  T tmp(std::cos(a._v));
  c._d*=tmp;
  return c;
}

template <class T, unsigned DIM>
inline ADf<T,DIM> cos (const ADf<T,DIM>& a)
{
  ADf<T,DIM> c(std::cos(a._v), a._d);
  T tmp(-std::sin(a._v));
  c._d*=tmp;
  return c;
}

template <class T, unsigned DIM>
inline ADf<T,DIM> tan (const ADf<T,DIM>& a)
{
  ADf<T,DIM> c(std::tan(a._v), a._d);
  T tmp(T(1)+c._v*c._v);
  c._d*=tmp;
  return c;
}

template <class T, unsigned DIM>
inline ADf<T,DIM> asin (const ADf<T,DIM>& a)
{
  ADf<T,DIM> c(std::asin(a._v), a._d);
  T tmp(T(1)/std::sqrt(T(1)-(a._v)*(a._v)));
  c._d*=tmp;
  return c;
}

template <class T, unsigned DIM>
inline ADf<T,DIM> acos (const ADf<T,DIM>& a)
{
  ADf<T,DIM> c(std::acos(a._v), a._d);
  T tmp(-T(1)/std::sqrt(T(1)-(a._v)*(a._v)));
  c._d*=tmp;
  return c;
}

template <class T, unsigned DIM>
inline ADf<T,DIM> atan (const ADf<T,DIM>& a)
{
  ADf<T,DIM> c(std::atan(a._v), a._d);
  T tmp(T(1)/(T(1)+(a._v)*(a._v)));
  c._d*=tmp;
  return c;
}

template <class T, unsigned DIM>
inline ADf<T,DIM> atan2 (const T& y, const ADf<T,DIM>& x)
{
  ADf<T,DIM> c(std::atan2(y, x._v));
  T tmpB(-y / (x._v*x._v + y*y));
  c._d = tmpB*x._d;
  return c;
}

template <class T, unsigned DIM>
inline ADf<T,DIM> atan2 (const ADf<T,DIM>& y, const T& x)
{
  ADf<T,DIM> c(std::atan2(y._v, x));
  T tmpA(x / (x*x + y._v*y._v));
  c._d = tmpA*y._d;
  return c;
}

template <class T, unsigned DIM>
inline ADf<T,DIM> atan2 (const ADf<T,DIM>& y, const ADf<T,DIM>& x)
{
  ADf<T,DIM> c(std::atan2(y._v, x._v));
  T tmpA(x._v / (x._v*x._v + y._v*y._v));
  T tmpB(-y._v / (x._v*x._v + y._v*y._v));
  c._d = tmpA*y._d + tmpB*x._d;
  return c;
}

template <class T, unsigned DIM>
inline ADf<T,DIM> sinh (const ADf<T,DIM>& a)
{
  return 0.5*(exp(a)-exp(-a));
}

template <class T, unsigned DIM>
inline ADf<T,DIM> cosh (const ADf<T,DIM>& a)
{
  return 0.5*(exp(a)+exp(-a));
}

template <class T, unsigned DIM>
inline ADf<T,DIM> tanh (const ADf<T,DIM>& a)
{
  return 1 - 2.0 / (exp(2*a) + 1);
}

template <class T, unsigned DIM>
inline ADf<T,DIM> asinh (const ADf<T,DIM>& a)
{
  return log( a + sqrt(a*a+1) );
}

template <class T, unsigned DIM>
inline ADf<T,DIM> acosh (const ADf<T,DIM>& a)
{
  if (a._v >= 1) {
    return log( a + sqrt(a*a-1) );
  }
  else {
    ADf<T,DIM> c;
    c._v = std::acosh(a._v);
    c._d = ( std::numeric_limits<T>::quiet_NaN() );
    return c;
  }
}

template <class T, unsigned DIM>
inline ADf<T,DIM> atanh (const ADf<T,DIM>& a)
{
  if (std::abs(a._v) < 1) {
    return 0.5 * log( (1+a) / (1-a) );
  }
  else {
    ADf<T,DIM> c;
    c._v = std::atanh(a._v);
    c._d = ( std::numeric_limits<T>::quiet_NaN() );
    return c;
  }
}

template <class T, unsigned DIM>
inline ADf<T,DIM> fmod (const ADf<T,DIM>& a, const ADf<T,DIM>& b)
{
  const T k = std::floor(T(a) / T(b));
  return ADf<T,DIM>(a - k * b);
}

template <class T, unsigned DIM, typename U, typename std::enable_if<std::is_integral<U>::value | std::is_floating_point<U>::value,int>::type = 0>
inline ADf<T,DIM> fmod (const ADf<T,DIM>& a, const U& b)
{
  return fmod(a, ADf<T,DIM>(b));
}

template <class T, unsigned DIM, typename U, typename std::enable_if<std::is_integral<U>::value | std::is_floating_point<U>::value,int>::type = 0>
inline ADf<T,DIM> fmod (const U& a, const ADf<T,DIM>& b)
{
  return fmod( ADf<T,DIM>(a),b);
}









/// tests if ADf has only zero derivatives
template <class T, unsigned DIM>
inline bool ADf<T,DIM>::hasZeroDerivative() const
{
  T sum = 0;
  for (unsigned i=0; i<DIM; i++) {
    sum += std::fabs(_d[i]);
  }
  return sum==0;
}



// boolean comparison operators: ==, !=, >, >=, <, <=
template <class T, unsigned DIM>
inline bool operator == (const ADf<T,DIM> &a, const ADf<T,DIM> &b)
{
  return a._v==b._v;
}

template <class T, unsigned DIM, typename U, typename std::enable_if<std::is_integral<U>::value | std::is_floating_point<U>::value,int>::type = 0>
inline bool operator == (const ADf<T,DIM> &a, const U &b)
{
  return a._v==olb::BaseType<ADf<T,DIM>>(b);
}

template <class T, unsigned DIM, typename U, typename std::enable_if<std::is_integral<U>::value | std::is_floating_point<U>::value,int>::type = 0>
inline bool operator == (const U &a, const ADf<T,DIM> &b)
{
  return olb::BaseType<ADf<T,DIM>>(a)==b._v;
}

template <class T, unsigned DIM>
inline bool operator != (const ADf<T,DIM> &a, const ADf<T,DIM> &b)
{
  return a._v!=b._v;
}

template <class T, unsigned DIM, typename U, typename std::enable_if<std::is_integral<U>::value | std::is_floating_point<U>::value,int>::type = 0>
inline bool operator != (const ADf<T,DIM> &a, const U &b)
{
  return a._v!=olb::BaseType<ADf<T,DIM>>(b);
}

template <class T, unsigned DIM, typename U, typename std::enable_if<std::is_integral<U>::value | std::is_floating_point<U>::value,int>::type = 0>
inline bool operator != (const U &a, const ADf<T,DIM> &b)
{
  return olb::BaseType<ADf<T,DIM>>(a)!=b._v;
}


template <class T, unsigned DIM>
inline bool operator > (const ADf<T,DIM> &a, const ADf<T,DIM> &b)
{
  return a._v>b._v;
}

template <class T, unsigned DIM, typename U, typename std::enable_if<std::is_integral<U>::value | std::is_floating_point<U>::value,int>::type = 0>
inline bool operator > (const ADf<T,DIM> &a, const U &b)
{
  return a._v>olb::BaseType<ADf<T,DIM>>(b);
}

template <class T, unsigned DIM, typename U, typename std::enable_if<std::is_integral<U>::value | std::is_floating_point<U>::value,int>::type = 0>
inline bool operator > (const U &a, const ADf<T,DIM> &b)
{
  return olb::BaseType<ADf<T,DIM>>(a)>b._v;
}

template <class T, unsigned DIM>
inline bool operator >= (const ADf<T,DIM> &a, const ADf<T,DIM> &b)
{
  return a._v>=b._v;
}

template <class T, unsigned DIM, typename U, typename std::enable_if<std::is_integral<U>::value | std::is_floating_point<U>::value,int>::type = 0>
inline bool operator >= (const ADf<T,DIM> &a, const U &b)
{
  return a._v>=olb::BaseType<ADf<T,DIM>>(b);
}

template <class T, unsigned DIM, typename U, typename std::enable_if<std::is_integral<U>::value | std::is_floating_point<U>::value,int>::type = 0>
inline bool operator >= (const U &a, const ADf<T,DIM> &b)
{
  return olb::BaseType<ADf<T,DIM>>(a)>=b._v;
}


template <class T, unsigned DIM>
inline bool operator < (const ADf<T,DIM> &a, const ADf<T,DIM> &b)
{
  return a._v<b._v;
}

template <class T, unsigned DIM, typename U, typename std::enable_if<std::is_integral<U>::value | std::is_floating_point<U>::value,int>::type = 0>
inline bool operator < (const ADf<T,DIM> &a, const U &b)
{
  return a._v<olb::BaseType<ADf<T,DIM>>(b);
}

template <class T, unsigned DIM, typename U, typename std::enable_if<std::is_integral<U>::value | std::is_floating_point<U>::value,int>::type = 0>
inline bool operator < (const U &a, const ADf<T,DIM> &b)
{
  return olb::BaseType<ADf<T,DIM>>(a)<b._v;
}

template <class T, unsigned DIM>
inline bool operator <= (const ADf<T,DIM> &a, const ADf<T,DIM> &b)
{
  return a._v<=b._v;
}

template <class T, unsigned DIM, typename U, typename std::enable_if<std::is_integral<U>::value | std::is_floating_point<U>::value,int>::type = 0>
inline bool operator <= (const ADf<T,DIM> &a, const U &b)
{
  return a._v<=olb::BaseType<ADf<T,DIM>>(b);
}

template <class T, unsigned DIM, typename U, typename std::enable_if<std::is_integral<U>::value | std::is_floating_point<U>::value,int>::type = 0>
inline bool operator <= (const U &a, const ADf<T,DIM> &b)
{
  return olb::BaseType<ADf<T,DIM>>(a)<=b._v;
}

// typecast functions
/// typecast from ADf --> olb::BaseType
template <class T, unsigned DIM>
inline ADf<T,DIM>::operator ADf<T,DIM>::base_t() const
{
#ifdef AdWarnings
  if (this->hasZeroDerivative() != true) {
    std::cout << "ADf WARNING: non-zero derivative in typecast ADf -> Basetype" << std::endl;
    std::cout << *this << std::endl;
    //exit(1);
  }
#endif
  return base_t(_v);
}

template <class T, unsigned DIM>
inline ADf<T,DIM> floor (const ADf<T,DIM>& a)
{
  ADf<T,DIM> c;
  if (a._v == std::floor(a._v)) {
    c._v = a._v;
#ifdef AdWarnings
    std::cout << "ADf WARNING: floor(Adf) - floor evaluated at non-differentiable point" << std::endl;
#endif
    for (unsigned i = 0; i < DIM; ++i) {
      c._d[i] = (a.d(i) == 0) ? T(0) : std::numeric_limits<T>::quiet_NaN();
    }
    return c;
  }
  c._v = std::floor(a._v);
  c._d = ( T(0) );
  return c;
}

template <unsigned DIM>
inline ADf<float,DIM> floorf (const ADf<float,DIM>& a)
{
  return floor(a);
}

template <unsigned DIM>
inline ADf<long double,DIM> floorl (const ADf<long double,DIM>& a)
{
  return floor(a);
}

template <class T, unsigned DIM>
inline ADf<T,DIM> ceil (const ADf<T,DIM>& a)
{
  // Similar to floor
  ADf<T,DIM> c (std::ceil(a._v));
  if (a._v == c._v) {
#ifdef AdWarnings
    std::cout << "ADf WARNING: ceil(Adf) - ceil evaluated at non-differentiable point" << std::endl;
#endif
    for (unsigned i = 0; i < DIM; ++i) {
      c._d[i] = (a.d(i) == 0) ? T(0) : std::numeric_limits<T>::quiet_NaN();
    }
  }
  return c;
}

template <unsigned DIM>
inline ADf<float,DIM> ceilf (const ADf<float,DIM>& a)
{
  return ceil(a);
}

template <unsigned DIM>
inline ADf<long double,DIM> ceill (const ADf<long double,DIM>& a)
{
  return ceil(a);
}

template <class T, unsigned DIM>
inline ADf<T,DIM> round (const ADf<T,DIM>& a)
{
  ADf<T,DIM> c (std::round(a._v));
  constexpr T EPSILON (std::numeric_limits<T>::epsilon());
  if (std::round(a._v - EPSILON) != std::round(a._v + EPSILON)) {
#ifdef AdWarnings
    std::cout << "ADf WARNING: round(Adf) - round evaluated at non-differentiable point" << std::endl;
#endif
    for (unsigned i = 0; i < DIM; ++i) {
      c._d[i] = (a.d(i) == 0) ? T(0) : std::numeric_limits<T>::quiet_NaN();
    }
  }
  return c;
}

template <unsigned DIM>
inline ADf<float,DIM> roundf (const ADf<float,DIM>& a)
{
  return round(a);
}

template <unsigned DIM>
inline ADf<long double,DIM> roundl (const ADf<long double,DIM>& a)
{
  return round(a);
}

// lround, llround, etc. throw error, since we do not support integers with ADf
template <class T, unsigned DIM>
inline ADf<T,DIM> lround (const ADf<T,DIM>& a)
{
#ifdef AdWarnings
  std::cout << "ADf WARNING: ADf does not support integer types. Using round() instead." << std::endl;
#endif
  return round(a);
}

template <unsigned DIM>
inline ADf<float,DIM> lroundf (const ADf<float,DIM>& a)
{
  return lround(a);
}

template <unsigned DIM>
inline ADf<long double,DIM> lroundl (const ADf<long double,DIM>& a)
{
  return lround(a);
}

template <class T, unsigned DIM>
inline ADf<T,DIM> llround (const ADf<T,DIM>& a)
{
  return lround(a);
}

template <unsigned DIM>
inline ADf<float,DIM> llroundf (const ADf<float,DIM>& a)
{
  return llround(a);
}

template <unsigned DIM>
inline ADf<long double,DIM> llroundl (const ADf<long double,DIM>& a)
{
  return llround(a);
}

template <class T, unsigned DIM>
inline ADf<T,DIM> fabs (const ADf<T,DIM>& a)
{
#ifdef AdWarnings
  if (a._v == T(0)) {
    std::cout << "ADf WARNING: fabs(Adf) - fabs evaluated at non-differentiable point" << std::endl;
  }
#endif
  return ADf<T,DIM>( util::sign(a._v) * a );
}

template <unsigned DIM>
inline ADf<float,DIM> fabsf (const ADf<float,DIM>& a)
{
  return fabs(a);
}

template <unsigned DIM>
inline ADf<long double,DIM> fabsl (const ADf<long double,DIM>& a)
{
  return fabs(a);
}

template <class T, unsigned DIM>
inline ADf<T,DIM> abs (const ADf<T,DIM>& a)
{
  return fabs(a);
}

template <class T, unsigned DIM>
inline ADf<T,DIM> labs (const ADf<T,DIM>& a)
{
#ifdef AdWarnings
  std::cout << "ADf WARNING: ADf is not supporting integer types. Using abs() instead." << std::endl;
#endif
  return abs(a);
}

template <class T, unsigned DIM>
inline ADf<T,DIM> llabs (const ADf<T,DIM>& a)
{
  return labs(a);
}

template <class T, unsigned DIM>
inline ADf<T,DIM> max (const olb::BaseType<ADf<T,DIM>>& a, const ADf<T,DIM>& b)
{
  ADf<T,DIM> c(a);
  c = max(c, b);
  return c;
}

template <class T, unsigned DIM>
inline ADf<T,DIM> max(const ADf<T,DIM>& a,const olb::BaseType<ADf<T,DIM>>& b)
{
  ADf<T,DIM> c(b);
  c = max(a, c);
  return c;
}

template <class T, unsigned DIM>
inline ADf<T,DIM> max(const ADf<T,DIM>& a, const ADf<T,DIM>& b)
{
  ADf<T,DIM> c((a + b + fabs(a - b) ) * T(0.5));
  return c;
}

template <class T, unsigned DIM>
inline ADf<T,DIM> min (const olb::BaseType<ADf<T,DIM>>& a, const ADf<T,DIM>& b)
{
  ADf<T,DIM> c(a);
  c = min(c, b);
  return c;
}

template <class T, unsigned DIM>
inline ADf<T,DIM> min(const ADf<T,DIM>& a,const olb::BaseType<ADf<T,DIM>>& b)
{
  ADf<T,DIM> c(b);
  c = min(a, c);
  return c;
}

template <class T, unsigned DIM>
inline ADf<T,DIM> min(const ADf<T,DIM>& a, const ADf<T,DIM>& b)
{
  ADf<T,DIM> c((a + b - fabs(a - b) ) * T(0.5));
  return c;
}

// overload nearZero to prevent bad behaviour due to implicit casting
template <class T, unsigned DIM>
inline bool nearZero(const ADf<T,DIM>& a)
{
  const T EPSILON = std::numeric_limits<T>::epsilon();
  if (a._v > -EPSILON && a._v < EPSILON) {
    return true;
  }
  else {
    return false;
  }
}

template <typename T>
inline void throwADfException(T arg)
{
  static_assert(std::is_floating_point<T>::value,
                "The data type ADf has slipped in here. If you really want to use it here and lose all derivation information, then perform an explicit typecast.");
}

} // namespace util

} // namespace olb

#endif
