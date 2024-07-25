/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2022 Adrian Kummerlaender
 *
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

#ifndef GPU_CUDA_COLUMN_H
#define GPU_CUDA_COLUMN_H

#include <memory>
#include <array>
#include <stdexcept>

#include "core/platform/platform.h"
#include "core/serializer.h"
#include "communication/communicatable.h"

#include "device.h"

#include <thrust/gather.h>
#include <thrust/scatter.h>

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

namespace olb {

namespace gpu {

namespace cuda {

/// Plain column for CUDA GPU targets
template<typename T>
class Column final : public AbstractColumn<T>
                   , public Serializable {
private:
  std::size_t _count;

  thrust::host_vector<T>   _hostData;
  thrust::device_vector<T> _deviceData;

public:
  using value_t = T;

  Column(std::size_t count):
    _count(count),
    _hostData(count),
    _deviceData(_hostData)
  { }

  Column(Column<T>&& rhs):
    _count(rhs._count),
    _hostData(std::move(rhs._hostData)),
    _deviceData(std::move(rhs._deviceData))
  { }

  virtual ~Column() = default;

  void resize(std::size_t newCount)
  {
    _hostData.resize(newCount);
    _deviceData.resize(newCount);
    _count = newCount;
  }

  const T& operator[](std::size_t i) const override
  {
    return _hostData[i];
  }

  T& operator[](std::size_t i) override
  {
    return _hostData[i];
  }

  std::size_t size() const
  {
    return _count;
  }

  const T* data() const
  {
    return _hostData.data();
  }

  T* data()
  {
    return _hostData.data();
  }

  const T* deviceData() const
  {
    return _deviceData.data().get();
  }

  T* deviceData()
  {
    return _deviceData.data().get();
  }

  void setProcessingContext(ProcessingContext);

  /// Number of data blocks for the serializable interface
  std::size_t getNblock() const override;
  /// Binary size for the serializer
  std::size_t getSerializableSize() const override;
  /// Return a pointer to the memory of the current block and its size for the serializable interface
  bool* getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode) override;

};


/// Virtual memory based cyclic column for usage in ColumnVector
/**
 * Column type used for propagatable population storage using
 * the virtual memory PS pattern.
 **/
template<typename T>
class CyclicColumn final : public AbstractCyclicColumn<T>
                         , public Serializable {
private:
  const std::ptrdiff_t _count;
  const std::size_t    _size;

  std::unique_ptr<T[]> _hostData;

  CUmemGenericAllocationHandle _handle;
  CUmemAllocationProp _prop{};
  CUmemAccessDesc _access{};
  CUdeviceptr _ptr;

  T* _deviceBase;
  T* _devicePopulation;

  std::ptrdiff_t _shift;

public:
  using value_t = T;

  CyclicColumn(std::size_t count):
    _count(device::getPageAlignedCount<T>(count)),
    _size(_count * sizeof(T)),
    _hostData(new T[_count] { }),
    _shift(0)
  {
    int device = device::get();

    _prop.type = CU_MEM_ALLOCATION_TYPE_PINNED;
    _prop.location.type = CU_MEM_LOCATION_TYPE_DEVICE;
    _prop.location.id = device;
    device::check(cuMemAddressReserve(&_ptr, 2 * _size, 0, 0, 0));

    // per-population handle until cuMemMap accepts non-zero offset
    device::check(cuMemCreate(&_handle, _size, &_prop, 0));
    device::check(cuMemMap(_ptr,         _size, 0, _handle, 0));
    device::check(cuMemMap(_ptr + _size, _size, 0, _handle, 0));

    _access.location.type = CU_MEM_LOCATION_TYPE_DEVICE;
    _access.location.id = device;
    _access.flags = CU_MEM_ACCESS_FLAGS_PROT_READWRITE;
    device::check(cuMemSetAccess(_ptr, 2 * _size, &_access, 1));

    _deviceBase = reinterpret_cast<T*>(_ptr);
    _devicePopulation = _deviceBase;
  }

  ~CyclicColumn() {
    device::check(cuMemUnmap(_ptr, 2 * _size));
    device::check(cuMemRelease(_handle));
    device::check(cuMemAddressFree(_ptr, 2 * _size));
  }

  const T& operator[](std::size_t i) const override
  {
    return _hostData[i];
  }

  T& operator[](std::size_t i) override
  {
    return _hostData[i];
  }

  const T* deviceData() const
  {
    return _devicePopulation;
  }

  T* deviceData()
  {
    return _devicePopulation;
  }

  std::size_t size() const
  {
    return _count;
  }

  void refresh()
  {
    _devicePopulation = _deviceBase + _shift;
  }

  void rotate(std::ptrdiff_t offset)
  {
    _shift -= offset;
    if (_shift >= _count) {
      _shift -= _count;
    }
    else if (_shift < 0) {
      _shift += _count;
    }
    refresh();
  }

  void setProcessingContext(ProcessingContext);

  /// Number of data blocks for the serializable interface
  std::size_t getNblock() const override;
  /// Binary size for the serializer
  std::size_t getSerializableSize() const override;
  /// Return a pointer to the memory of the current block and its size for the serializable interface
  bool* getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode) override;
  void postLoad() override { refresh(); }

};

}

}


/// Declare gpu::cuda::Column as the AbstractColumn implementation for GPU CUDA targets
template <typename T>
struct ImplementationOf<AbstractColumn<T>,Platform::GPU_CUDA> {
  using type = gpu::cuda::Column<T>;
};

/// Declare gpu::cuda::CyclicColumn as the AbstractCyclicColumn implementation for GPU CUDA targets
template <typename T>
struct ImplementationOf<AbstractCyclicColumn<T>,Platform::GPU_CUDA> {
  using type = gpu::cuda::CyclicColumn<T>;
};


/// Communicatable implementation for a single gpu::cuda::Column
/**
 * This is only implemented to match the interface of CPU communication.
 * The actual communication buffers for communicating block lattice data
 * are handled using custom kernel functions in the the Platform::GPU_CUDA
 * specialization of ConcreteBlockCommunicator.
 **/
template <typename T>
class ConcreteCommunicatable<gpu::cuda::Column<T>> final : public Communicatable {
private:
  gpu::cuda::Column<T>& _column;

public:
  ConcreteCommunicatable(gpu::cuda::Column<T>& column):
    _column{column} { }

  /// Get serialized size for data at locations `indices`
  std::size_t size(ConstSpan<CellID> indices) const
  {
    return indices.size() * sizeof(T);
  }

  /// Serialize data at locations `indices` to `buffer`
  std::size_t serialize(ConstSpan<CellID> indices,
                        std::uint8_t* buffer) const
  {
    thrust::gather(thrust::device,
                   thrust::device_pointer_cast(indices.begin()),
                   thrust::device_pointer_cast(indices.end()),
                   thrust::device_pointer_cast(_column.deviceData()),
                   thrust::device_pointer_cast(reinterpret_cast<T*>(buffer)));
    return indices.size() * sizeof(T);
  }

  /// Deserialize data at locations `indices` to `buffer`
  std::size_t deserialize(ConstSpan<CellID> indices,
                          const std::uint8_t* buffer)
  {
    thrust::scatter(thrust::device,
                    thrust::device_pointer_cast(reinterpret_cast<const T*>(buffer)),
                    thrust::device_pointer_cast(reinterpret_cast<const T*>(buffer) + indices.size()),
                    thrust::device_pointer_cast(indices.begin()),
                    thrust::device_pointer_cast(_column.deviceData()));
    return indices.size() * sizeof(T);
  }

};

/// Communicatable implementation for a single gpu::cuda::CyclicColumn
template <typename T>
class ConcreteCommunicatable<gpu::cuda::CyclicColumn<T>> final : public Communicatable {
private:
  gpu::cuda::CyclicColumn<T>& _column;

public:
  ConcreteCommunicatable(gpu::cuda::CyclicColumn<T>& column):
    _column{column} { }

  /// Get serialized size for data at locations `indices`
  std::size_t size(ConstSpan<CellID> indices) const
  {
    return indices.size() * sizeof(T);
  }

  /// Serialize data at locations `indices` to `buffer`
  std::size_t serialize(ConstSpan<CellID> indices,
                        std::uint8_t* buffer) const
  {
    thrust::gather(thrust::device,
                   thrust::device_pointer_cast(indices.begin()),
                   thrust::device_pointer_cast(indices.end()),
                   thrust::device_pointer_cast(_column.deviceData()),
                   thrust::device_pointer_cast(reinterpret_cast<T*>(buffer)));
    return indices.size() * sizeof(T);
  }

  /// Deserialize data at locations `indices` to `buffer`
  std::size_t deserialize(ConstSpan<CellID> indices,
                          const std::uint8_t* buffer)
  {
    thrust::scatter(thrust::device,
                    thrust::device_pointer_cast(reinterpret_cast<const T*>(buffer)),
                    thrust::device_pointer_cast(reinterpret_cast<const T*>(buffer) + indices.size()),
                    thrust::device_pointer_cast(indices.begin()),
                    thrust::device_pointer_cast(_column.deviceData()));
    return indices.size() * sizeof(T);
  }

};

}

#endif

#include "column.hh"
