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

#ifndef GPU_CUDA_REGISTRY_H
#define GPU_CUDA_REGISTRY_H

#include "utilities/typeIndexedContainers.h"

#include <map>

namespace olb {

namespace gpu {

namespace cuda {

/// Mapping of TYPE in CONTEXT to runtime-fixed index
/**
 * Used for dynamic access to field arrays
 **/
template <typename CONTEXT, typename TYPE>
__constant__ std::size_t field_type_index;

/// Type-erased pointer to FieldArrayD device data
/**
 * Used for dynamic communication buffer (de)serialization
 **/
struct AnyDeviceFieldArrayD {
  void**   data;
  unsigned column_count;
  unsigned element_size;

  std::uint8_t* operator[](unsigned iD) __device__ {
    return reinterpret_cast<std::uint8_t*>(data[iD]);
  }
};

}

}

/// Maintain on-device structure for dynamic field access
/**
 * Used to enable access to fields that are dynamically allocated
 * instead of being declared by the descriptor.
 **/
template<typename T, typename DESCRIPTOR>
class FieldTypeRegistry<T,DESCRIPTOR,Platform::GPU_CUDA> {
private:
  /// Host-side field type index
  utilities::TypeIndexedMap<AnyFieldType<T,DESCRIPTOR,Platform::GPU_CUDA>*,FieldTypeRegistry> _index;
  /// Host-side version of gpu::cuda::AnyDeviceFieldArrayD
  struct FieldArrayPointer {
    gpu::cuda::device::unique_ptr<void*> data;
    const unsigned column_count;
    const unsigned element_size;
  };

  std::map<std::type_index,FieldArrayPointer> _fieldArrayPointers;

  std::vector<void**>                   _indexOnHost;
  gpu::cuda::device::unique_ptr<void**> _indexOnDevice;

  bool _modified;

public:
  FieldTypeRegistry():
    _indexOnHost(2, nullptr),
    _modified{true} { }

  template <typename FIELD_TYPE>
  bool provides() const {
    return _index.template provides<FIELD_TYPE>();
  }

  template <typename FIELD_TYPE>
  AnyFieldType<T,DESCRIPTOR,Platform::GPU_CUDA>* get() {
    return _index.template get<FIELD_TYPE>();
  }

  template <typename FIELD_TYPE>
  void track(AnyFieldType<T,DESCRIPTOR,Platform::GPU_CUDA>* fieldType) {
    _index.template set<FIELD_TYPE>(fieldType);

    // Copy host-side field index to device-side field index
    std::size_t index = _index.template index<FIELD_TYPE>();
    cudaMemcpyToSymbol(gpu::cuda::field_type_index<void,FIELD_TYPE>, &index, sizeof(std::size_t));
    gpu::cuda::device::check();
    if (index >= _indexOnHost.size()) {
      _indexOnHost.resize(2*index);
    }

    // Update device-side pointers to Array<FIELD> field types
    using ConcreteFieldType = typename FIELD_TYPE::template type<T,DESCRIPTOR,Platform::GPU_CUDA>;
    if constexpr (std::is_base_of_v<ColumnVectorBase, ConcreteFieldType>) {
      using field_t = typename ConcreteFieldType::field_t;
      auto& fieldArray = *fieldType->template as<FIELD_TYPE>();

      std::array<void*,DESCRIPTOR::template size<field_t>()> componentPointers;
      for (unsigned iD=0; iD < componentPointers.size(); ++iD) {
        componentPointers[iD] = fieldArray[iD].deviceData();
      }

      auto componentPointersOnDevice = gpu::cuda::device::malloc<void*>(DESCRIPTOR::template size<field_t>());
      gpu::cuda::device::copyToDevice(componentPointers.data(),
                                      componentPointersOnDevice,
                                      componentPointers.size()*sizeof(void*));

      _indexOnHost[index] = componentPointersOnDevice;
      _fieldArrayPointers.emplace(typeid(field_t), FieldArrayPointer{
        gpu::cuda::device::unique_ptr<void*>(componentPointersOnDevice),
        componentPointers.size(),
        sizeof(typename field_t::template value_type<T>)
      });
    }

    _modified = true;
  }

  void*** deviceData() {
    if (_modified) {
      _indexOnDevice = gpu::cuda::device::malloc<void**>(_indexOnHost.size());
      gpu::cuda::device::copyToDevice(_indexOnHost.data(),
                                      _indexOnDevice.get(),
                                      _indexOnHost.size()*sizeof(void**));
      _modified = false;
    }
    return _indexOnDevice.get();
  }

  template <typename FIELD>
  void refreshDeviceFieldArray(FieldArrayD<T,DESCRIPTOR,Platform::GPU_CUDA,FIELD>& fieldArray) {
    std::array<void*,DESCRIPTOR::template size<FIELD>()> componentPointers;
    for (unsigned iD=0; iD < componentPointers.size(); ++iD) {
      componentPointers[iD] = fieldArray[iD].deviceData();
    }
    auto& ptr = _fieldArrayPointers.at(typeid(FIELD));
    gpu::cuda::device::copyToDevice(componentPointers.data(),
                                    ptr.data.get(),
                                    componentPointers.size()*sizeof(void*));
  }

  gpu::cuda::AnyDeviceFieldArrayD deviceFieldArray(std::type_index field) {
    auto& ptr = _fieldArrayPointers.at(field);
    return gpu::cuda::AnyDeviceFieldArrayD{
      ptr.data.get(),
      ptr.column_count,
      ptr.element_size
    };
  }

  std::vector<gpu::cuda::AnyDeviceFieldArrayD> deviceFieldArrays(
    const std::vector<std::type_index>& fields) {
    std::vector<gpu::cuda::AnyDeviceFieldArrayD> deviceFields;
    deviceFields.reserve(fields.size());
    for (std::type_index field : fields) {
      deviceFields.emplace_back(deviceFieldArray(field));
    }
    return deviceFields;
  }

};

}

#endif
