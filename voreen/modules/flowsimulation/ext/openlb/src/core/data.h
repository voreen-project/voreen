/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2022 Adrian Kummerlaender
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

#ifndef CORE_DATA_H
#define CORE_DATA_H

#include "platform/platform.h"

#include "utilities/typeIndexedContainers.h"
#include "dynamics/context.h"

#include <optional>

namespace olb {

template <typename T, typename DESCRIPTOR, Platform PLATFORM> class ConcreteBlockMask;

/// Describe FieldArray of a FIELD in Data
/**
 * i.e. Data::getField<Array<FIELD>> will expose FieldArrayD<FIELD>
 **/
template <typename FIELD>
struct Array {
  template <typename T, typename DESCRIPTOR, Platform PLATFORM>
  using type = FieldArrayD<T,DESCRIPTOR,PLATFORM,FIELD>;
};

/// Describe paramaters of DYNAMICS in Data
/**
 * i.e. Data::getField<DynamicsParameters<DYNAMICS>> will expose DynamicsParametersD<DYNAMICS>
 **/
template <typename DYNAMICS>
struct DynamicsParameters {
  template <typename T, typename DESCRIPTOR, Platform PLATFORM>
  using type = DynamicsParametersD<T,DESCRIPTOR,PLATFORM,DYNAMICS>;
};

/// Describe mask of DYNAMICS in Data
/**
 * i.e. Data::getField<DynamicsMask<DYNAMICS>> will expose ConcreteBlockMask
 **/
template <typename DYNAMICS>
struct DynamicsMask {
  template <typename T, typename DESCRIPTOR, Platform PLATFORM>
  using type = ConcreteBlockMask<T,DESCRIPTOR,PLATFORM>;
};

template <typename CONTEXT>
struct TrivialParameters {
  template <typename T, typename DESCRIPTOR, Platform PLATFORM>
  using type = TrivialParametersD<T,DESCRIPTOR,PLATFORM,CONTEXT>;
};

/// Helper for referring to arbitrary data instances
/**
 * Enables runtime-allocated field storage in Data
 **/
template<typename T, typename DESCRIPTOR, Platform PLATFORM>
class AnyFieldType {
private:
  std::unique_ptr<Serializable> _data;
  std::function<void(ProcessingContext)> _setProcessingContext;

public:
  template<typename FIELD_TYPE, typename... ARGS>
  AnyFieldType(meta::id<FIELD_TYPE>, ARGS&&... args):
    _data(new typename FIELD_TYPE::template type<T,DESCRIPTOR,PLATFORM>(
      std::forward<decltype(args)>(args)...)),
    _setProcessingContext([&](ProcessingContext context) {
      static_cast<typename FIELD_TYPE::template type<T,DESCRIPTOR,PLATFORM>*>(_data.get())->setProcessingContext(context);
    })
  { }

  template<typename FIELD_TYPE>
  auto* as() {
    return static_cast<typename FIELD_TYPE::template type<T,DESCRIPTOR,PLATFORM>*>(_data.get());
  }

  Serializable* asSerializable() {
    return _data.get();
  }

  template<typename TYPE>
  std::optional<TYPE*> tryAs() {
    if (TYPE* ptr = dynamic_cast<TYPE*>(_data.get())) {
      return std::optional<TYPE*>(ptr);
    } else {
      return std::nullopt;
    }
  }

  void setProcessingContext(ProcessingContext context) {
    _setProcessingContext(context);
  }

};

/// Efficient indexing of dynamically allocated data fields
/**
 * May be overloaded to manage platform specific access structures
 **/
template<typename T, typename DESCRIPTOR, Platform PLATFORM>
class FieldTypeRegistry {
private:
  utilities::TypeIndexedMap<AnyFieldType<T,DESCRIPTOR,PLATFORM>*,FieldTypeRegistry> _index;

public:
  /// Returns true iff FIELD_TYPE is registered
  /**
   * FIELD_TYPE is registered if track<FIELD_TYPE> was called previously
   **/
  template <typename FIELD_TYPE>
  bool provides() const {
    return _index.template provides<FIELD_TYPE>();
  }

  /// Return pointer to FIELD_TYPE data
  /**
   * FIELD_TYPE must be registered
   **/
  template <typename FIELD_TYPE>
  AnyFieldType<T,DESCRIPTOR,PLATFORM>* get() {
    return _index.template get<FIELD_TYPE>();
  }

  /// Track newly allocated FIELD_TYPE
  template <typename FIELD_TYPE>
  void track(AnyFieldType<T,DESCRIPTOR,PLATFORM>* fieldType) {
    _index.template set<FIELD_TYPE>(fieldType);
  }

};


/// Storage of any FIELD_TYPE data on PLATFORM
/**
 * Manages access and serialization of dynamically allocated FIELD_TYPE data.
 * This is the core storage class used for ConcreteBlockLattice on any platform
 * supported by OpenLB.
 **/
template<typename T, typename DESCRIPTOR, Platform PLATFORM>
class Data final : public Serializable {
private:
  /// Map for dynamic-dispatch of field accesses
  std::map<std::type_index, std::unique_ptr<AnyFieldType<T,DESCRIPTOR,PLATFORM>>> _map;
  /// Efficient static-dispatch access to fields
  FieldTypeRegistry<T,DESCRIPTOR,PLATFORM> _registry;
  /// Map for ordering serialization of fields
  std::map<std::string, Serializable*> _serializationNames;

public:
  /// Returns true iff FIELD_TYPE is allocated and can be accessed
  template <typename FIELD_TYPE>
  bool provides() const {
    return _registry.template provides<FIELD_TYPE>();
  }

  /// Allocate and return FIELD_TYPE data
  /**
   * Must only be called once for each FIELD_TYPE
   **/
  template <typename FIELD_TYPE, typename... ARGS>
  auto& allocate(ARGS&&... args) {
    _map[typeid(FIELD_TYPE)] = std::make_unique<AnyFieldType<T,DESCRIPTOR,PLATFORM>>(
      meta::id<FIELD_TYPE>(), std::forward<decltype(args)>(args)...);

    AnyFieldType<T,DESCRIPTOR,PLATFORM>* newField = _map[typeid(FIELD_TYPE)].get();
    _registry.template track<FIELD_TYPE>(newField);
    return *newField->template as<FIELD_TYPE>();
  }

  /// Sets FIELD_TYPE serialization state to active
  /**
   * By default Data doesn't (de)serialize any fields but
   * ConcreteBlockLattice sets all descriptor-declared
   * field arrays to active.
   **/
  template <typename FIELD_TYPE>
  void setSerialization(bool active) {
    if (active) {
      _serializationNames[meta::name<FIELD_TYPE>()] = _registry.template get<FIELD_TYPE>()->asSerializable();
    } else {
      _serializationNames.erase(meta::name<FIELD_TYPE>());
    }
  }

  /// Return reference to data of FIELD_TYPE
  template <typename FIELD_TYPE>
  auto& get() {
    return *_registry.template get<FIELD_TYPE>()->template as<FIELD_TYPE>();
  }

  /// Return reference to data of Array<FIELD>, i.e. FieldArrayD<FIELD>
  /**
   * Most common field type, used for any lattice structured data in ConcreteBlockLattice
   **/
  template <typename FIELD>
  FieldArrayD<T,DESCRIPTOR,PLATFORM,FIELD>& getFieldArray() {
    return get<Array<FIELD>>();
  }

  void setProcessingContext(ProcessingContext context) {
    for (auto& [_, field] : _map) {
      field->setProcessingContext(context);
    }
  }

  /// Call f for each managed AnyFieldType
  template <typename F>
  void forEach(F f) {
    for (auto& [_, field] : _map) {
      f(*field);
    }
  }

  /// Call f for each managed AnyFieldType of TYPE
  /**
   * Used to e.g. iterate over all field arrays
   **/
  template <typename TYPE, typename F>
  void forEachCastable(F f) {
    forEach([&f](AnyFieldType<T,DESCRIPTOR,PLATFORM>& any) {
      if (auto casted = any.template tryAs<TYPE>()) {
        f(*casted);
      }
    });
  }

  /// Expose FieldTypeRegistry for device-support
  auto& getRegistry() {
    return _registry;
  }

  /// Number of data blocks for the serializable interface
  std::size_t getNblock() const override;
  /// Binary size for the serializer
  std::size_t getSerializableSize() const override;
  /// Return a pointer to the memory of the current block and its size for the serializable interface
  bool* getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode) override;

};

template<typename T, typename DESCRIPTOR, Platform PLATFORM>
std::size_t Data<T,DESCRIPTOR,PLATFORM>::getNblock() const
{
  std::size_t nBlock = 0;
  for (const auto& [_, serializable] : _serializationNames) {
    nBlock += serializable->getNblock();
  }
  return nBlock;
}

template<typename T, typename DESCRIPTOR, Platform PLATFORM>
std::size_t Data<T,DESCRIPTOR,PLATFORM>::getSerializableSize() const
{
  std::size_t size = 0;
  for (const auto& [_, serializable] : _serializationNames) {
    size += serializable->getSerializableSize();
  }
  return size;
}

template<typename T, typename DESCRIPTOR, Platform PLATFORM>
bool* Data<T,DESCRIPTOR,PLATFORM>::getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode)
{
  std::size_t currentBlock = 0;
  bool* dataPtr = nullptr;

  for (const auto& [_, serializable] : _serializationNames) {
    registerSerializableOfConstSize(iBlock, sizeBlock, currentBlock, dataPtr, *serializable, loadingMode);
  }

  return dataPtr;
}

}

#endif
