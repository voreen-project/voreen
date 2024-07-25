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
#include <type_traits>
#include <array>

#include "scalarVector.h"
#include "utilities/omath.h"
#include "utilities/oalgorithm.h"
#include "olbDebug.h"

namespace olb {


/// Plain old scalar vector
template <typename T, unsigned D>
class Vector : public ScalarVector<T,D,Vector<T,D>> {
private:
  std::array<T,D> _data;

  friend typename ScalarVector<T,D,Vector<T,D>>::type;

protected:
  const T* getComponentPointer(unsigned iDim) const any_platform
  {
    return &_data[iDim];
  }
  T* getComponentPointer(unsigned iDim) any_platform
  {
    return &_data[iDim];
  }

public:
  using value_t = T;

  Vector() any_platform:
    _data{}
  { }

  Vector(const Vector& rhs) any_platform:
    _data(rhs._data)
  { }

  Vector(Vector&& rhs) any_platform:
    _data(rhs._data)
  { }

  template <typename W, typename IMPL>
  Vector(const ScalarVector<W,D,IMPL>& rhs) any_platform
  {
    for (unsigned iDim=0; iDim < D; ++iDim) {
      _data[iDim] = rhs[iDim];
    }
  }

  template <typename... I,
            typename std::enable_if_t< (sizeof...(I) > 1
                                    && (meta::is_arithmetic<I>::value && ...)
                                    && sizeof...(I) == D), int> = 0>
  Vector(I... indexes) any_platform:
    _data{T(indexes)...}
  { }

  // Should be declared explicit
  template <typename U>
  Vector(const U* v) any_platform
  {
    for (unsigned iDim=0; iDim < D; ++iDim) {
      _data[iDim] = static_cast<T>(v[iDim]);
    }
  }

  // Should be declared explicit
  Vector(const std::vector<T>& v)
  {
    OLB_PRECONDITION(v.size() == D);
    for (unsigned iDim=0; iDim < D; ++iDim) {
      _data[iDim] = v[iDim];
    }
  }

  #ifndef  __CUDA_ARCH__
  Vector(std::initializer_list<T> v)
  {
    OLB_PRECONDITION(v.size() == D);
    std::copy(v.begin(), v.end(), _data.begin());
  }
  #endif

  Vector(T scalar) any_platform
  {
    for (unsigned iDim=0; iDim < D; ++iDim) {
      _data[iDim] = scalar;
    }
  }

  /// Construct with entries given by a lambda expression
  template <typename F, typename = decltype(std::declval<F&>()(std::size_t{0}))>
  Vector(F&& f) any_platform
  {
    for (unsigned iDim=0; iDim < D; ++iDim) {
      _data[iDim] = f(iDim);
    }
  }

  template <typename U, typename IMPL_>
  Vector& operator = (const GenericVector<U,D,IMPL_>& rhs) any_platform
  {
    for (unsigned iDim=0; iDim < D; ++iDim) {
      this->operator[](iDim) = rhs[iDim];
    }
    return *this;
  }

  Vector& operator = (const Vector& rhs) any_platform
  {
    for (unsigned iDim=0; iDim < D; ++iDim) {
      _data[iDim] = rhs[iDim];
    }
    return *this;
  }

  Vector& operator = (const Vector&& rhs) any_platform
  {
    for (unsigned iDim=0; iDim < D; ++iDim) {
      _data[iDim] = rhs[iDim];
    }
    return *this;
  }

  const T* data() const any_platform
  {
    return _data.data();
  }

  T* data() any_platform
  {
    return _data.data();
  }

  auto begin()
  {
    return _data.begin();
  }

  auto end()
  {
    return _data.end();
  }

  int getDim() const any_platform
  {
    return D;
  }

  Vector<T,D+1> withPrefix(T prefix) const any_platform
  {
    Vector<T,D+1> tmp;
    tmp[0] = prefix;
    for (unsigned iDim=0; iDim < D; ++iDim) {
      tmp[iDim+1] = _data[iDim];
    }
    return tmp;
  }

  std::size_t getNblock() const { return 1; }
  std::size_t getSerializableSize() const
  {
    return D * sizeof(T);
  };
  bool* getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode)
  {
    std::size_t currentBlock = 0;
    bool* dataPtr = nullptr;
    //registerVar(iBlock, sizeBlock, currentBlock, dataPtr, *_data.data(), D);
    if (iBlock == currentBlock) {
      sizeBlock = sizeof(T) * D;
      dataPtr = (bool*) _data.data();
    }
    currentBlock++;
    return dataPtr;
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
  using namespace util;
  return Vector<T,D>([&a](unsigned iDim) -> T {
    return abs(a[iDim]);
  });
}


template <typename T, unsigned D, typename U, typename IMPL>
meta::enable_if_arithmetic_t<U, Vector<T,D>>
operator+ (U a, const ScalarVector<T,D,IMPL>& b) any_platform
{
  return Vector<T,D>(b) += a;
}

template <typename T, unsigned D, typename U, typename IMPL>
meta::enable_if_arithmetic_t<U, Vector<T,D>>
operator+ (const ScalarVector<T,D,IMPL>& a, U b) any_platform
{
  return Vector<T,D>(a) += b;
}

template <typename T, unsigned D, typename IMPL, typename IMPL_>
Vector<T,D> operator+ (const ScalarVector<T,D,IMPL>& a, const ScalarVector<T,D,IMPL_>& b) any_platform
{
  return Vector<T,D>(a) += b;
}

template <typename T, typename W, unsigned D, typename IMPL, typename IMPL_>
Vector<decltype(T{}+W{}),D> operator+ (const ScalarVector<T,D,IMPL>& a, const ScalarVector<W,D,IMPL_>& b) any_platform
{
  Vector<decltype(T{}+W{}),D> result;
  for (unsigned iDim=0; iDim < D; ++iDim) {
    result[iDim] = a[iDim] + b[iDim];
  }
  return result;
}

template <typename T, unsigned D, typename U, typename IMPL>
meta::enable_if_arithmetic_t<U, Vector<T,D>>
operator- (U a, const ScalarVector<T,D,IMPL>& b) any_platform
{
  return Vector<T,D>(a) - b;
}

template <typename T, unsigned D, typename U, typename IMPL>
meta::enable_if_arithmetic_t<U, Vector<T,D>>
operator- (const ScalarVector<T,D,IMPL>& a, U b) any_platform
{
  return Vector<T,D>(a) -= b;
}

template <typename T, unsigned D, typename IMPL, typename IMPL_>
Vector<T,D> operator- (const ScalarVector<T,D,IMPL>& a, const ScalarVector<T,D,IMPL_>& b) any_platform
{
  return Vector<T,D>(a) -= b;
}

template <typename T, typename W, unsigned D, typename IMPL, typename IMPL_>
Vector<decltype(T{}-W{}),D> operator- (const ScalarVector<T,D,IMPL>& a, const ScalarVector<W,D,IMPL_>& b) any_platform
{
  Vector<decltype(T{}-W{}),D> result;
  for (unsigned iDim=0; iDim < D; ++iDim) {
    result[iDim] = a[iDim] - b[iDim];
  }
  return result;
}

template <typename T, unsigned D, typename U, typename IMPL>
meta::enable_if_arithmetic_t<U, Vector<decltype(T{}*U{}),D>>
operator* (U a, const ScalarVector<T,D,IMPL>& b) any_platform
{
  Vector<decltype(T{}*U{}),D> result(b);
  return result *= a;
}

template <typename T, unsigned D, typename U, typename IMPL>
meta::enable_if_arithmetic_t<U, Vector<T,D>>
operator* (const ScalarVector<T,D,IMPL>& a, U b) any_platform
{
  Vector<decltype(T{}*U{}),D> result(a);
  return result *= b;
}

/// Inner product
template <typename T, typename U, unsigned D, typename IMPL, typename IMPL_>
auto operator* (const ScalarVector<T,D,IMPL>& a, const ScalarVector<U,D,IMPL_>& b) any_platform
{
  decltype(T{}*U{}) scalarProduct{};
  for (unsigned iDim=0; iDim < D; ++iDim) {
    scalarProduct += a[iDim] * b[iDim];
  }
  return scalarProduct;
}

template <typename T, unsigned D, typename U, typename IMPL>
meta::enable_if_arithmetic_t<U, Vector<T,D>>
operator/ (const ScalarVector<T,D,IMPL>& a, U b) any_platform
{
  return Vector<T,D>(a) /= b;
}

template<typename T, unsigned D, typename U, typename IMPL>
meta::enable_if_arithmetic_t<U, bool>
operator< (U lhs, const ScalarVector<T,D,IMPL>& rhs) any_platform
{
  return Vector<U,D>(lhs) < rhs;
}

template<typename T, unsigned D, typename U, typename IMPL>
meta::enable_if_arithmetic_t<U, bool>
operator> (const ScalarVector<T,D,IMPL>& lhs, U rhs) any_platform
{
  return lhs > Vector<U,D>(rhs);
}

template<typename T, unsigned D, typename U, typename IMPL>
meta::enable_if_arithmetic_t<U, bool>
operator<= (U lhs, const ScalarVector<T,D,IMPL>& rhs) any_platform
{
  return Vector<U,D>(lhs) <= rhs;
}

template<typename T, unsigned D, typename U, typename IMPL>
meta::enable_if_arithmetic_t<U, bool>
operator>= (const ScalarVector<T,D,IMPL>& lhs, U rhs) any_platform
{
  return lhs >= Vector<U,D>(rhs);
}

template <typename T, unsigned D, typename IMPL, typename IMPL_>
Vector<T,D> maxv(const ScalarVector<T,D,IMPL>& v, const ScalarVector<T,D,IMPL_>& w)
{
  return Vector<T,D>([&v,&w](unsigned iDim) -> T {
    return util::max(v[iDim], w[iDim]);
  });
}


} // end namespace olb

#endif
