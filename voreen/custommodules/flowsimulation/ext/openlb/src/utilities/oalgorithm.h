/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2020 Jan E. Marquardt
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

#ifndef OLB_OALGORITHM_H
#define OLB_OALGORITHM_H

#include <algorithm>

#include "core/meta.h"

namespace olb {

namespace util {

// Max
template< typename T > constexpr T max( T a, meta::id_t<T> b )
{
  return std::max(a, b);
}
template< typename T, class Compare > constexpr T max( T a, meta::id_t<T> b, Compare comp )
{
  return std::max(a, b, comp);
}
template< typename T > constexpr T max( std::initializer_list<T> ilist )
{
  return std::max(ilist);
}
template< typename T, class Compare > constexpr T max( std::initializer_list<T> ilist, Compare comp )
{
  return std::max(ilist, comp);
}

template <>
float max<float>(float x, float y) any_platform
{
  return std::fmax(x, y);
}

template <>
double max<double> (double x, double y) any_platform
{
  return std::fmax(x, y);
}

// Min
template< typename T > constexpr T min( T a, meta::id_t<T> b )
{
  return std::min(a, b);
}
template< typename T, class Compare > constexpr T min( T a, meta::id_t<T> b, Compare comp )
{
  return std::min(a, b, comp);
}
template< typename T > constexpr T min( std::initializer_list<T> ilist )
{
  return std::min(ilist);
}
template< typename T, class Compare > constexpr T min( std::initializer_list<T> ilist, Compare comp )
{
  return std::min(ilist, comp);
}

template <>
float min<float>(float x, float y) any_platform
{
  return std::fmin(x, y);
}

template <>
double min<double>(double x, double y) any_platform
{
  return std::fmin(x, y);
}


} // namespace util

} // namespace olb

#endif
