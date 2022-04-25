/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2019 Adrian Kummerlaender
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

#ifndef DESCRIPTOR_FUNCTION_H
#define DESCRIPTOR_FUNCTION_H

#include <type_traits>
#include <stdexcept>

#include "descriptorTag.h"

#include "utilities/fraction.h"
#include "core/meta.h"

namespace olb {

template <typename, unsigned> class Vector;

namespace descriptors {

/// \defgroup descriptor_interface Descriptor functions
/// \ingroup descriptor
//@{

/// \defgroup descriptor_interface_details Descriptor data
/// \ingroup descriptor_interface
//@{

namespace data {

using utilities::Fraction;

template <unsigned D, unsigned Q>
platform_constant int vicinity = {};

template <unsigned D, unsigned Q>
platform_constant int c[Q][D] = {};

template <unsigned D, unsigned Q>
platform_constant int opposite[Q] = {};

template <unsigned D, unsigned Q>
platform_constant Fraction t[Q] = {};

template <unsigned D, unsigned Q>
platform_constant Fraction cs2 = {};

template <unsigned D, unsigned Q>
platform_constant Fraction lambda_e = {};

template <unsigned D, unsigned Q>
platform_constant Fraction lambda_h = {};

}

template <unsigned D, unsigned Q>
constexpr int vicinity() any_platform
{
  return data::vicinity<D,Q>;
}

template <unsigned D, unsigned Q>
constexpr int c(unsigned iPop, unsigned iDim) any_platform
{
  return data::c<D,Q>[iPop][iDim];
}

template <unsigned D, unsigned Q>
constexpr Vector<int,D> c(unsigned iPop) any_platform
{
  return Vector<int,D>(data::c<D,Q>[iPop]);
}

template <unsigned D, unsigned Q>
constexpr int opposite(unsigned iPop) any_platform
{
  return data::opposite<D,Q>[iPop];
}

template <typename T, unsigned D, unsigned Q>
constexpr T t(unsigned iPop, tag::DEFAULT) any_platform
{
  return data::t<D,Q>[iPop].template as<T>();
}

template <typename T, unsigned D, unsigned Q>
constexpr T invCs2() any_platform
{
  return data::cs2<D,Q>.template inverseAs<T>();
}

template <typename T, unsigned D, unsigned Q>
constexpr T lambda_e() any_platform
{
  return data::lambda_e<D,Q>.template inverseAs<T>();
}

template <typename T, unsigned D, unsigned Q>
constexpr T lambda_h() any_platform
{
  return data::lambda_h<D,Q>.template inverseAs<T>();
}

//@}

template <typename DESCRIPTOR>
constexpr int d() any_platform
{
  return DESCRIPTOR::d;
}

struct dimension {
  template <typename DESCRIPTOR>
  constexpr auto operator()(meta::id<DESCRIPTOR> = meta::id<DESCRIPTOR>{}) const {
    return d<DESCRIPTOR>;
  }
};

template <typename DESCRIPTOR>
constexpr int q() any_platform
{
  return DESCRIPTOR::q;
}

template <typename DESCRIPTOR>
constexpr int vicinity() any_platform
{
  return vicinity<DESCRIPTOR::d, DESCRIPTOR::q>();
}

template <typename DESCRIPTOR>
constexpr int c(unsigned iPop, unsigned iDim) any_platform
{
  return c<DESCRIPTOR::d, DESCRIPTOR::q>(iPop, iDim);
}

template <typename DESCRIPTOR>
constexpr Vector<int,DESCRIPTOR::d> c(unsigned iPop) any_platform
{
  return c<DESCRIPTOR::d, DESCRIPTOR::q>(iPop);
}

template <typename DESCRIPTOR>
constexpr int opposite(unsigned iPop) any_platform
{
  return opposite<DESCRIPTOR::d, DESCRIPTOR::q>(iPop);
}

template <typename T, typename DESCRIPTOR>
constexpr T t(unsigned iPop) any_platform
{
  return t<T, DESCRIPTOR::d, DESCRIPTOR::q>(iPop, typename DESCRIPTOR::category_tag());
}

template <typename T, typename DESCRIPTOR>
constexpr T invCs2() any_platform
{
  return invCs2<T, DESCRIPTOR::d, DESCRIPTOR::q>();
}

template <typename T, typename DESCRIPTOR>
constexpr T lambda_e() any_platform
{
  return lambda_e<T, DESCRIPTOR::d, DESCRIPTOR::q>();
}

template <typename T, typename DESCRIPTOR>
constexpr T lambda_h() any_platform
{
  return lambda_h<T, DESCRIPTOR::d, DESCRIPTOR::q>();
}

/// Return array of population indices of DESCRIPTOR matching PREDICATE
template <typename DESCRIPTOR, typename PREDICATE>
constexpr auto filter_population_indices(PREDICATE predicate) any_platform {
  return meta::filter_index_sequence(predicate,
                                     std::make_index_sequence<DESCRIPTOR::q>());
}

//@}

}

}

#endif
