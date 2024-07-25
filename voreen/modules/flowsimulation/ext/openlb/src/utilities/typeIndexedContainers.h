/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021-22 Adrian Kummerlaender, Julius Jessberger
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

#ifndef TYPE_INDEXED_CONTAINERS_H
#define TYPE_INDEXED_CONTAINERS_H

#include <vector>
#include <optional>

#include "core/meta.h"

namespace olb {

namespace utilities {

template<typename KEYS, typename VALUE>
class FixedTypeIndexedMap {
private:
  std::array<VALUE,KEYS::size> _index;

public:
  template<typename TYPE>
  static constexpr bool provides() {
    return KEYS::template contains<TYPE>();
  }

  template<typename TYPE>
  const VALUE& get() const {
    return _index[KEYS::template index<TYPE>()];
  }
  template<typename TYPE>
  VALUE& get() {
    return _index[KEYS::template index<TYPE>()];
  }

  template<typename TYPE>
  void set(VALUE value) {
    _index[KEYS::template index<TYPE>()] = value;
  }

};

/// (Time) efficient mapping between TYPEs and VALUEs
/**
 * CONTEXT should be set to the containing type / some one-of identifier
 * to restrict the index space to the types that are actually used in the
 * specific context.
 **/
template<typename VALUE, typename CONTEXT=void>
class TypeIndexedMap {
private:
  std::vector<std::optional<VALUE>> _index;

  static std::size_t next_index();
  template <typename TYPE>
  static std::size_t get_index();

public:
  TypeIndexedMap() {
    _index.resize(1);
  }

  template<typename TYPE>
  bool provides() const;

  template<typename TYPE>
  std::size_t index() const;

  template<typename TYPE>
  const VALUE& get() const;
  template<typename TYPE>
  VALUE& get();

  template<typename TYPE>
  void set(VALUE value);

};

template <typename VALUE, typename CONTEXT>
std::size_t TypeIndexedMap<VALUE,CONTEXT>::next_index()
{
  std::size_t curr;
  #ifdef PARALLEL_MODE_OMP
  #pragma omp critical
  #endif
  {
    static std::size_t index = 0;
    curr = index++;
  }
  return curr;
}

template <typename VALUE, typename CONTEXT>
template <typename TYPE>
std::size_t TypeIndexedMap<VALUE,CONTEXT>::get_index()
{
  static std::size_t index = next_index();
  return index;
}

template <typename VALUE, typename CONTEXT>
template <typename TYPE>
bool TypeIndexedMap<VALUE,CONTEXT>::provides() const
{
  const std::size_t index = get_index<TYPE>();
  if (index >= _index.size()) {
    return false;
  } else {
    return _index[index].has_value();
  }
}

template <typename VALUE, typename CONTEXT>
template <typename TYPE>
std::size_t TypeIndexedMap<VALUE,CONTEXT>::index() const
{
  return get_index<TYPE>();
}

template <typename VALUE, typename CONTEXT>
template <typename TYPE>
const VALUE& TypeIndexedMap<VALUE,CONTEXT>::get() const
{
  const std::size_t index = get_index<TYPE>();
  if (auto& value = _index[index]) {
    return value.value();
  } else {
    throw std::out_of_range("VALUE for TYPE doesn't exist as it was not set.");
  }
}

template <typename VALUE, typename CONTEXT>
template <typename TYPE>
VALUE& TypeIndexedMap<VALUE,CONTEXT>::get()
{
  const std::size_t index = get_index<TYPE>();
  if (auto& value = _index[index]) {
    return value.value();
  } else {
    throw std::out_of_range("VALUE for TYPE doesn't exist as it was not set.");
  }
}

template <typename VALUE, typename CONTEXT>
template <typename TYPE>
void TypeIndexedMap<VALUE,CONTEXT>::set(VALUE value)
{
  const std::size_t index = get_index<TYPE>();
  if (index >= _index.size()) {
    _index.resize(2*index);
  }
  _index[index] = value;
}


/// Mapping between KEYs and instances of type VALUEs
// Template argument MAP is something like the meta::map
template <typename MAP>
struct TypeIndexedTuple {
  template <typename... TYPES>
  using SharedPointerTuple = std::tuple<std::shared_ptr<TYPES>...>;

  /// An std::tuple with shared_ptrs to the values of MAP as types.
  using tuple_t = typename MAP::values_t::template decompose_into<SharedPointerTuple>;
  tuple_t tuple;

  /// Access by KEYS of MAP
  template <typename KEY>
  constexpr auto& get() {
    return std::get<(MAP::keys_t::template index<KEY>())>(tuple);
  }
  /// Access by KEYS of MAP
  template <typename KEY>
  constexpr auto const& get() const {
    return std::get<(MAP::keys_t::template index<KEY>())>(tuple);
  }
};

}

}

#endif
