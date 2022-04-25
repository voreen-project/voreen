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


#ifndef GENERIC_VECTOR_H
#define GENERIC_VECTOR_H

#include <type_traits>

#include "meta.h"

namespace olb {


/// Generic vector of values supporting basic arithmetic
template<typename T, unsigned D, typename IMPL>
struct GenericVector {
  GenericVector() = default;
  GenericVector(const GenericVector&) = delete;
  GenericVector(GenericVector&&) = delete;

  GenericVector& operator = (GenericVector&& rhs) = delete;

  using value_type = T;
  static constexpr unsigned d = D;

  const T& operator [] (unsigned iDim) const any_platform
  {
    return *static_cast<const IMPL*>(this)->getComponentPointer(iDim);
  }

  T& operator [] (unsigned iDim) any_platform
  {
    return *static_cast<IMPL*>(this)->getComponentPointer(iDim);
  }

  template<typename IMPL_>
  IMPL& operator = (const GenericVector<T,D,IMPL_>& rhs) any_platform
  {
    for (unsigned iDim=0; iDim < D; ++iDim) {
      operator[](iDim) = rhs[iDim];
    }
    return *static_cast<IMPL*>(this);
  }

  template<typename U, typename IMPL_>
  IMPL& operator += (const GenericVector<U,D,IMPL_>& rhs) any_platform
  {
    for (unsigned iDim=0; iDim < D; ++iDim) {
      operator[](iDim) += rhs[iDim];
    }
    return *static_cast<IMPL*>(this);
  }

  template<typename IMPL_>
  IMPL& operator -= (const GenericVector<T,D,IMPL_>& rhs) any_platform
  {
    for (unsigned iDim=0; iDim < D; ++iDim) {
      operator[](iDim) -= rhs[iDim];
    }
    return *static_cast<IMPL*>(this);
  }

  template<typename U, typename IMPL_>
  IMPL& operator *= (const GenericVector<U,D,IMPL_>& rhs) any_platform
  {
    for (unsigned iDim=0; iDim < D; ++iDim) {
      operator[](iDim) *= rhs[iDim];
    }
    return *static_cast<IMPL*>(this);
  }

  template<typename U>
  meta::enable_if_arithmetic_t<U, IMPL&> operator += (const U& rhs) any_platform
  {
    for (unsigned iDim=0; iDim < D; ++iDim) {
      operator[](iDim) += rhs;
    }
    return *static_cast<IMPL*>(this);
  }

  template<typename U>
  meta::enable_if_arithmetic_t<U, IMPL&> operator -= (const U& rhs) any_platform
  {
    for (unsigned iDim=0; iDim < D; ++iDim) {
      operator[](iDim) -= rhs;
    }
    return *static_cast<IMPL*>(this);
  }

  template<typename U>
  meta::enable_if_arithmetic_t<U, IMPL&> operator *= (const U& rhs) any_platform
  {
    for (unsigned iDim=0; iDim < D; ++iDim) {
      operator[](iDim) *= rhs;
    }
    return *static_cast<IMPL*>(this);
  }

  template<typename U>
  meta::enable_if_arithmetic_t<U, IMPL&> operator /= (const U& rhs) any_platform
  {
    for (unsigned iDim=0; iDim < D; ++iDim) {
      operator[](iDim) /= rhs;
    }
    return *static_cast<IMPL*>(this);
  }


  template<typename IMPL_>
  bool operator == (const GenericVector<T,D,IMPL_>& rhs) const any_platform
  {
    bool isEqual = true;
    for (unsigned iDim=0; iDim < D; ++iDim) {
      isEqual &= operator[](iDim) == rhs[iDim];
    }
    return isEqual;
  }

  template<typename IMPL_>
  bool operator != (const GenericVector<T,D,IMPL_>& rhs) const any_platform
  {
    for (unsigned iDim=0; iDim < D; ++iDim) {
      if (operator[](iDim) != rhs[iDim]) {
        return true;
      }
    }
    return false;
  }

};


}

#endif
