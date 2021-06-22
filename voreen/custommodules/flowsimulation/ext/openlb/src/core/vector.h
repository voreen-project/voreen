/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2015 Asher Zarth, Mathias J. Krause, Albert Mink
 *                2020 Adrian Kummerlaender
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
 * efficient implementation of a vector class
 */

#ifndef VECTOR_H
#define VECTOR_H

#include <cstring>
#include <algorithm>
#include <type_traits>

#include "olbDebug.h"
#include "scalarVector.h"
#include "util.h"

namespace olb {


/// Plain old scalar vector
template <typename T, unsigned D>
class Vector : public ScalarVector<T,D,Vector<T,D>> {
private:
  T _data[D];

  friend typename ScalarVector<T,D,Vector<T,D>>::type;

protected:
  const T* getComponentPointer(unsigned iDim) const {
    return &_data[iDim];
  }
  T* getComponentPointer(unsigned iDim) {
    return &_data[iDim];
  }

public:
  Vector():
    _data{} { }

  template <typename W, typename IMPL>
  Vector(const ScalarVector<W,D,IMPL>& rhs)
  {
    for (unsigned iDim=0; iDim < D; ++iDim) {
      _data[iDim] = rhs[iDim];
    }
  }

  Vector(const Vector& rhs)
  {
    for (unsigned iDim=0; iDim < D; ++iDim) {
      _data[iDim] = rhs[iDim];
    }
  }

  Vector(const T v[D])
  {
    std::memcpy(_data, v, D*sizeof(T));
  }

  Vector(const std::vector<T>& v):
    Vector(v.data())
  { }

  Vector(std::initializer_list<T> v)
  {
    OLB_PRECONDITION(v.size() == D);
    std::copy(v.begin(), v.end(), _data);
  }

  Vector(T scalar)
  {
    for (unsigned iDim=0; iDim < D; ++iDim) {
      _data[iDim] = scalar;
    }
  }

  Vector(T a, T b)
  {
    OLB_PRECONDITION(D == 2);
    _data[0] = a;
    _data[1] = b;
  }

  Vector(T a, T b, T c)
  {
    OLB_PRECONDITION(D == 3);
    _data[0] = a;
    _data[1] = b;
    _data[2] = c;
  }

  /// Construct with entries given by a lambda expression
  template <typename F, typename = decltype(std::declval<F&>()(std::size_t{0}))>
  Vector(F&& f)
  {
    for (unsigned iDim=0; iDim < D; ++iDim) {
      _data[iDim] = f(iDim);
    }
  }

  template <typename U, typename IMPL_>
  Vector& operator = (const GenericVector<U,D,IMPL_>& rhs)
  {
    for (unsigned iDim=0; iDim < D; ++iDim) {
      this->operator[](iDim) = rhs[iDim];
    }
    return *this;
  }

  Vector& operator = (const Vector<T,D>& rhs)
  {
    this->operator=<T,Vector<T,D>>(rhs);
    return *this;
  }

  const T* data() const
  {
    return _data;
  }

  T* data()
  {
    return _data;
  }

};


template <typename T, typename IMPL, typename IMPL_>
Vector<T,3> crossProduct3D(const ScalarVector<T,3,IMPL>& a, const ScalarVector<T,3,IMPL_>& b)
{
  return Vector<T,3>(
    a[1]*b[2] - a[2]*b[1],
    a[2]*b[0] - a[0]*b[2],
    a[0]*b[1] - a[1]*b[0]
  );
}

template <typename T, unsigned D, typename IMPL>
Vector<T,D> normalize(const ScalarVector<T,D,IMPL>& a, T scale = T{1})
{
  const T invScale = scale / norm(a);
  return Vector<T,D>([invScale,&a](unsigned iDim) -> T {
    return a[iDim] * invScale;
  });
}

template <typename T, unsigned D, typename IMPL>
Vector<T,D> abs(const ScalarVector<T,D,IMPL>& a)
{
  return Vector<T,D>([a](unsigned iDim) -> T {
    return abs(a[iDim]);
  });
}


template <typename T, unsigned D, typename U, typename IMPL>
utilities::meta::enable_if_arithmetic_t<U, Vector<T,D>>
operator+ (U a, const ScalarVector<T,D,IMPL>& b)
{
  return Vector<T,D>(b) += a;
}

template <typename T, unsigned D, typename U, typename IMPL>
utilities::meta::enable_if_arithmetic_t<U, Vector<T,D>>
operator+ (const ScalarVector<T,D,IMPL>& a, U b)
{
  return Vector<T,D>(a) += b;
}

template <typename T, unsigned D, typename IMPL, typename IMPL_>
Vector<T,D> operator+ (const ScalarVector<T,D,IMPL>& a, const ScalarVector<T,D,IMPL_>& b)
{
  return Vector<T,D>(a) += b;
}

template <typename T, typename W, unsigned D, typename IMPL, typename IMPL_>
Vector<decltype(T{}+W{}),D> operator+ (const ScalarVector<T,D,IMPL>& a, const ScalarVector<W,D,IMPL_>& b)
{
  Vector<decltype(T{}+W{}),D> result;
  for (unsigned iDim=0; iDim < D; ++iDim) {
    result[iDim] = a[iDim] + b[iDim];
  }
  return result;
}

template <typename T, unsigned D, typename U, typename IMPL>
utilities::meta::enable_if_arithmetic_t<U, Vector<T,D>>
operator- (U a, const ScalarVector<T,D,IMPL>& b)
{
  return Vector<T,D>(a) - b;
}

template <typename T, unsigned D, typename U, typename IMPL>
utilities::meta::enable_if_arithmetic_t<U, Vector<T,D>>
operator- (const ScalarVector<T,D,IMPL>& a, U b)
{
  return Vector<T,D>(a) -= b;
}

template <typename T, unsigned D, typename IMPL, typename IMPL_>
Vector<T,D> operator- (const ScalarVector<T,D,IMPL>& a, const ScalarVector<T,D,IMPL_>& b)
{
  return Vector<T,D>(a) -= b;
}

template <typename T, typename W, unsigned D, typename IMPL, typename IMPL_>
Vector<decltype(T{}-W{}),D> operator- (const ScalarVector<T,D,IMPL>& a, const ScalarVector<W,D,IMPL_>& b)
{
  Vector<decltype(T{}-W{}),D> result;
  for (unsigned iDim=0; iDim < D; ++iDim) {
    result[iDim] = a[iDim] - b[iDim];
  }
  return result;
}

template <typename T, unsigned D, typename U, typename IMPL>
utilities::meta::enable_if_arithmetic_t<U, Vector<decltype(T{}*U{}),D>>
operator* (U a, const ScalarVector<T,D,IMPL>& b)
{
  Vector<decltype(T{}*U{}),D> result(b);
  return result *= a;
}

template <typename T, unsigned D, typename U, typename IMPL>
utilities::meta::enable_if_arithmetic_t<U, Vector<T,D>>
operator* (const ScalarVector<T,D,IMPL>& a, U b)
{
  Vector<decltype(T{}*U{}),D> result(a);
  return result *= b;
}

/// Inner product
template <typename T, unsigned D, typename IMPL, typename IMPL_>
T operator* (const ScalarVector<T,D,IMPL>& a, const ScalarVector<T,D,IMPL_>& b)
{
  T scalarProduct{};
  for (unsigned iDim=0; iDim < D; ++iDim) {
    scalarProduct += a[iDim] * b[iDim];
  }
  return scalarProduct;
}

template <typename T, unsigned D, typename U, typename IMPL>
utilities::meta::enable_if_arithmetic_t<U, Vector<T,D>>
operator/ (const ScalarVector<T,D,IMPL>& a, U b)
{
  return Vector<T,D>(a) /= b;
}

template <typename T, unsigned D, typename IMPL, typename IMPL_>
Vector<T,D> maxv(const ScalarVector<T,D,IMPL>& v, const ScalarVector<T,D,IMPL_>& w)
{
  return Vector<T,D>([&v,&w](unsigned iDim) -> T {
    return std::max(v[iDim], w[iDim]);
  });
}


} // end namespace olb

#endif
