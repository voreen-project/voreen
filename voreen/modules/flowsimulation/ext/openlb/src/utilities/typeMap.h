/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Adrian Kummerlaender
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

#ifndef UTILITIES_TYPE_MAP_H
#define UTILITIES_TYPE_MAP_H

#include "core/meta.h"

namespace olb {

namespace meta {

template <typename...>
struct map_keys;

template <typename KEY, typename VALUE, typename... TAIL>
struct map_keys<KEY,VALUE,TAIL...> {
  static_assert(sizeof...(TAIL) % 2 == 0, "TAIL must be valid map size");
  using type = typename map_keys<TAIL...>::type::template push<KEY>;
};

template <>
struct map_keys<> {
  using type = list<>;
};

template <typename...>
struct map_values;

template <typename KEY, typename VALUE, typename... TAIL>
struct map_values<KEY,VALUE,TAIL...> {
  static_assert(sizeof...(TAIL) % 2 == 0, "TAIL must be valid map size");
  using type = typename map_values<TAIL...>::type::template push<VALUE>;
};

template <>
struct map_values<> {
  using type = list<>;
};

/// Map of types
/**
 * To be tidied up and ported onto newer type list in improvement/communicator
 **/
template <typename... KVs>
struct map {
  using keys_t   = typename map_keys<KVs...>::type;
  using values_t = typename map_values<KVs...>::type;

  template <typename KEY>
  using value = typename values_t::template get<(keys_t::template index<KEY>())>;

  static constexpr unsigned size = keys_t::size;
};

}

}

#endif
