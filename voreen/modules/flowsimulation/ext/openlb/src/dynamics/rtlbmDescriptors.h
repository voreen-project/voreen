/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2017-2019 Albert Mink, Adrian Kummerlaender
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
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA 02110-1301, USA.
*/

/** \file
 *  -- header file
 */
#ifndef RTLBM_DESCRIPTORS_H
#define RTLBM_DESCRIPTORS_H

#include "latticeDescriptors.h"
#include "utilities/omath.h"

namespace olb {

namespace descriptors {

namespace tag {

struct RTLBM : public CATEGORY, public DESCRIPTOR_TAG { };

}


namespace rtlbm_data {

using utilities::Fraction;

template <unsigned D, unsigned Q>
platform_constant Fraction t[Q] = {};

template <unsigned D, unsigned Q>
platform_constant double norm_c[Q] = {};

template <>
platform_constant_definition Fraction t<3,7>[7] = {
  0,
  {1, 6}, {1, 6}, {1, 6},
  {1, 6}, {1, 6}, {1, 6}
};

template <>
platform_constant_definition double norm_c<3,7>[7] = {
  0.0,
  1.0, 1.0, 1.0,
  1.0, 1.0, 1.0
};

template <>
platform_constant_definition Fraction t<3,15>[15] = {
  0,
  {1, 15}, {1, 15}, {1, 15},
  {3, 40}, {3, 40}, {3, 40}, {3, 40},
  {1, 15}, {1, 15}, {1, 15},
  {3, 40}, {3, 40}, {3, 40}, {3, 40}
};

template <>
platform_constant_definition double norm_c<3,15>[15] = {
  0.0,
  1.0, 1.0, 1.0,
  1.73205080757, 1.73205080757, 1.73205080757, 1.73205080757,
  1.0, 1.0, 1.0,
  1.73205080757, 1.73205080757, 1.73205080757, 1.73205080757
};

template <>
platform_constant_definition Fraction t<3,27>[27] = {
  0,

  {1, 21}, {1, 21}, {1, 21},
  {4, 105}, {4, 105}, {4, 105},
  {4, 105}, {4, 105}, {4, 105},
  {9, 280}, {9, 280}, {9, 280}, {9, 280},

  {1, 21}, {1, 21}, {1, 21},
  {4, 105}, {4, 105}, {4, 105},
  {4, 105}, {4, 105}, {4, 105},
  {9, 280}, {9, 280}, {9, 280}, {9, 280}
};

template <>
platform_constant_definition double norm_c<3,27>[27] = {
  0.0,
  1.0, 1.0, 1.0,
  1.41421356237, 1.41421356237, 1.41421356237,
  1.41421356237, 1.41421356237, 1.41421356237,
  1.73205080757, 1.73205080757, 1.73205080757, 1.73205080757,
  1.0, 1.0, 1.0,
  1.41421356237, 1.41421356237, 1.41421356237,
  1.41421356237, 1.41421356237, 1.41421356237,
  1.73205080757, 1.73205080757, 1.73205080757, 1.73205080757
};

}

template <typename T, unsigned D, unsigned Q>
constexpr T t(unsigned iPop, tag::RTLBM)
{
  return rtlbm_data::t<D,Q>[iPop].template as<T>();
}

template <typename T, unsigned D, unsigned Q>
constexpr T norm_c(unsigned iPop, tag::RTLBM)
{
  return T{rtlbm_data::norm_c<D,Q>[iPop]};
}

template <typename T, typename DESCRIPTOR>
constexpr T norm_c(unsigned iPop)
{
  return norm_c<T, DESCRIPTOR::d, DESCRIPTOR::q>(iPop, typename DESCRIPTOR::category_tag());
}

}  // namespace descriptors

}  // namespace olb

#endif
