/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2020 Adrian Kummerlaender
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


#ifndef SCALAR_VECTOR_H
#define SCALAR_VECTOR_H

#include <vector>
#include <limits>
#include <iostream>
#include <type_traits>

#include "genericVector.h"
#include "meta.h"
#include "utilities/omath.h"
#include "utilities/oalgorithm.h"

namespace olb {


/// Vector of scalars
template<typename T, unsigned D, typename IMPL>
struct ScalarVector : public GenericVector<T,D,IMPL> {
  using type = GenericVector<T,D,IMPL>;

  static_assert(meta::is_arithmetic<T>::type::value, "T must be scalar type");

  ScalarVector() = default;
  ScalarVector(const ScalarVector&) = delete;
  ScalarVector(ScalarVector&&) = delete;
};

/// Squared euclidean vector norm
template<typename T, unsigned D, typename IMPL>
inline T norm_squared(const ScalarVector<T,D,IMPL>& a)
{
  T sqNorm{};
  for (unsigned iDim=0; iDim < D; ++iDim) {
    sqNorm += a[iDim] * a[iDim];
  }
  return sqNorm;
}

/// Euclidean vector norm
template<typename T, unsigned D, typename IMPL>
inline T norm(const ScalarVector<T,D,IMPL>& a)
{
  using namespace util;
  T sqNorm = norm_squared(a);
  return sqrt(sqNorm);
}

/// Returns true iff all components are within floating point error distance of 0
template<typename T, unsigned D, typename IMPL>
bool closeToZero(const ScalarVector<T,D,IMPL>& a)
{
  const T eps = std::numeric_limits<T>::epsilon();
  for (unsigned iDim=0; iDim < D; ++iDim) {
    if (util::fabs(a[iDim]) > eps) {
      return false;
    }
  }
  return true;
}

/// Copies data into a standard vector
template<typename T, unsigned D, typename IMPL>
std::vector<T> toStdVector(const ScalarVector<T,D,IMPL>& a)
{
  std::vector<T> v(D);
  for (unsigned iDim=0; iDim < D; ++iDim) {
    v[iDim] = a[iDim];
  }
  return v;
}

/// Returns true if all lhs components are smaller than rhs
template<typename T, unsigned D, typename U, typename IMPL, typename IMPL_>
inline bool operator< (const ScalarVector<T,D,IMPL>& lhs, const ScalarVector<U,D,IMPL_>& rhs)
{
  bool smaller = true;
  for (unsigned iDim=0; iDim < D; ++iDim) {
    smaller &= (lhs[iDim] < rhs[iDim]);
  }
  return smaller;
}

/// Returns true if all lhs components are greater than rhs
template<typename T, unsigned D, typename U, typename IMPL, typename IMPL_>
inline bool operator> (const ScalarVector<T,D,IMPL>& lhs, const ScalarVector<U,D,IMPL_>& rhs)
{
  return rhs < lhs;
}

/// Returns true if all lhs components are smaller or equal than rhs
template<typename T, unsigned D, typename U, typename IMPL, typename IMPL_>
inline bool operator<= (const ScalarVector<T,D,IMPL>& lhs, const ScalarVector<U,D,IMPL_>& rhs)
{
  bool smallerEq = true;
  for (unsigned iDim=0; iDim < D; ++iDim) {
    smallerEq &= (lhs[iDim] <= rhs[iDim]);
  }
  return smallerEq;
}

/// Returns true if all lhs components are smaller or equal than rhs
template<typename T, unsigned D, typename U, typename IMPL, typename IMPL_>
inline bool operator>= (const ScalarVector<T,D,IMPL>& lhs, const ScalarVector<U,D,IMPL_>& rhs)
{
  return rhs <= lhs;
}

/// Print vector entries to ostream in a human-readable fashion
template<typename T, unsigned D, typename IMPL>
inline std::ostream& operator << (std::ostream& os, const ScalarVector<T,D,IMPL>& o)
{
  if (D > 0) {
    os << "[";
    for (unsigned iDim=0; iDim < D-1; ++iDim) {
      os << o[iDim] << " ";
    }
    os << o[D-1]<<"]";
  }
  else {
    os << "[empty]";
  }
  return os;
}


}

#endif
