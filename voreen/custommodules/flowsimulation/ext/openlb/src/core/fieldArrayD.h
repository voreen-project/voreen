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

#ifndef FIELD_ARRAY_D_H
#define FIELD_ARRAY_D_H

#include "serializer.h"
#include "columnVector.h"
#include "utilities/meta.h"
#include "dynamics/descriptorBase.h"
#include "core/vector.h"

#include <memory>
#include <tuple>
#include <typeindex>

namespace olb {

/// Vector storing a single field instance 
template<typename T, typename DESCRIPTOR, typename FIELD>
using FieldD = Vector<
  typename FIELD::template value_type<T>,
  DESCRIPTOR::template size<FIELD>()
>;

/// Vector-compatible pointer to a single field instance
template<typename T, typename DESCRIPTOR, typename FIELD>
using FieldPtr = VectorPtr<
  typename FIELD::template value_type<T>,
  DESCRIPTOR::template size<FIELD>()
>;

/// Read-only vector-compatible pointer to a single field instance
template<typename T, typename DESCRIPTOR, typename FIELD>
using ConstFieldPtr = ConstVectorPtr<
  typename FIELD::template value_type<T>,
  DESCRIPTOR::template size<FIELD>()
>;

/// SoA storage for instances of a single FIELD
template<typename T, typename DESCRIPTOR, typename FIELD>
struct FieldArrayD : public ColumnVector<typename FIELD::template value_type<T>,
                                         DESCRIPTOR::template size<FIELD>()> {
  using value_type = typename FIELD::template value_type<T>;

  template<typename TAG>
  FieldArrayD(std::size_t count, TAG tag):
    ColumnVector<value_type,DESCRIPTOR::template size<FIELD>()>(count, tag)
  { }

  FieldArrayD(std::size_t count):
    ColumnVector<value_type,DESCRIPTOR::template size<FIELD>()>(count)
  { }

  ConstFieldPtr<T,DESCRIPTOR,FIELD> getFieldPointer(std::size_t i) const {
    return ConstVectorPtr<value_type,DESCRIPTOR::template size<FIELD>()>(*this, i);
  }
  FieldPtr<T,DESCRIPTOR,FIELD> getFieldPointer(std::size_t i) {
    return VectorPtr<value_type,DESCRIPTOR::template size<FIELD>()>(*this, i);
  }

};

/// Helper for referring to arbitrary FieldArrayD instances
/**
 * Enables runtime-allocated field storage in BlockDynamicFieldsD*D
 **/
template<typename T, typename DESCRIPTOR>
class AnyFieldArrayD {
private:
  const std::type_index _field_type_index;
  std::unique_ptr<ColumnVectorBase> _field_array;

public:
  template<typename FIELD_ID>
  AnyFieldArrayD(FIELD_ID id, std::size_t count):
    _field_type_index(typeid(typename FIELD_ID::type)),
    _field_array(new FieldArrayD<T,DESCRIPTOR,typename FIELD_ID::type>(count))
  { }

  template<typename FIELD>
  const FieldArrayD<T,DESCRIPTOR,FIELD>& as() const {
    OLB_ASSERT(std::type_index(typeid(FIELD)) == _field_type_index, "Invalid cast");
    return static_cast<const FieldArrayD<T,DESCRIPTOR,FIELD>&>(*_field_array);
  }

  template<typename FIELD>
  FieldArrayD<T,DESCRIPTOR,FIELD>& as() {
    OLB_ASSERT(std::type_index(typeid(FIELD)) == _field_type_index, "Invalid cast");
    return static_cast<FieldArrayD<T,DESCRIPTOR,FIELD>&>(*_field_array);
  }

  const Serializable& asSerializable() const {
    return static_cast<const Serializable&>(*_field_array);
  }
  Serializable& asSerializable() {
    return static_cast<Serializable&>(*_field_array);
  }

};


/// Static storage for multiple FIELDS
template<typename T, typename DESCRIPTOR, typename... FIELDS>
class MultiFieldArrayD : public Serializable {
private:
  std::size_t _count;
  std::tuple<FieldArrayD<T,DESCRIPTOR,FIELDS>...> _data;

public:
  MultiFieldArrayD(std::size_t count):
    _count(count),
    // Trickery to construct each member of _data with `count`.
    // Uses the comma operator in conjunction with type dropping.
    _data((utilities::meta::void_t<FIELDS>(), count)...)
  { }

  /// Change number of rows
  /**
   * Drops the last (newCount-_count) rows when shrinking
   **/
  void resize(std::size_t newCount) {
    utilities::meta::swallow((get<FIELDS>().resize(newCount), 0)...);
    _count = newCount;
  }

  /// Swap contents of rows i and j
  void swap(std::size_t i, std::size_t j) {
    utilities::meta::swallow((get<FIELDS>().swap(i,j), 0)...);
  }

  template <typename FIELD>
  std::enable_if_t<utilities::meta::list_contains_item<FIELD,FIELDS...>::type::value, const FieldArrayD<T,DESCRIPTOR,FIELD>&>
  get() const {
    return std::get<descriptors::getIndexInFieldList<FIELD,FIELDS...>()>(_data);
  }

  template <typename FIELD>
  std::enable_if_t<!utilities::meta::list_contains_item<FIELD,FIELDS...>::type::value, const FieldArrayD<T,DESCRIPTOR,FIELD>&>
  get() const {
    throw std::invalid_argument("This MultiFieldArrayD does not provide FIELD.");
  }

  template <typename FIELD>
  std::enable_if_t<utilities::meta::list_contains_item<FIELD,FIELDS...>::type::value, FieldArrayD<T,DESCRIPTOR,FIELD>&>
  get() {
    return std::get<descriptors::getIndexInFieldList<FIELD,FIELDS...>()>(_data);
  }

  template <typename FIELD>
  std::enable_if_t<!utilities::meta::list_contains_item<FIELD,FIELDS...>::type::value, FieldArrayD<T,DESCRIPTOR,FIELD>&>
  get() {
    throw std::invalid_argument("This MultiFieldArrayD does not provide FIELD.");
  }

  /// Apply generic lambda expression to each FIELD of a cell
  template <typename F>
  void forFieldsAt(std::size_t idx, F f) {
    utilities::meta::swallow(
      (f(get<FIELDS>().getFieldPointer(idx), utilities::meta::id<FIELDS>{}), 0)...);
  }

  /// Number of data blocks for the serializable interface
  std::size_t getNblock() const override;
  /// Binary size for the serializer
  std::size_t getSerializableSize() const override;
  /// Return a pointer to the memory of the current block and its size for the serializable interface
  bool* getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode) override;

};

}

#endif
