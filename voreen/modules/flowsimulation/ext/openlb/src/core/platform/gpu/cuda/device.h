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

#ifndef GPU_CUDA_DEVICE_H
#define GPU_CUDA_DEVICE_H

#include <stdexcept>

#include <cuda.h>

namespace olb {

namespace gpu {

namespace cuda {

/// Basic wrappers of common CUDA functions
namespace device {

void check() {
  cudaError_t error = cudaGetLastError();
  if (error != cudaSuccess) {
    throw std::runtime_error(cudaGetErrorString(error));
  }
}

/// Check CUDA driver errors
void check(CUresult result) {
  if (result != CUDA_SUCCESS) {
    const char* description{};
    cuGetErrorString(result, &description);
    throw std::runtime_error(std::string(description));
  }
}

void synchronize() {
  cudaDeviceSynchronize();
  check();
}

int get() {
  int device;
  cudaGetDevice(&device);
  check();
  return device;
}

void copyToHost(void* src, void* dst, std::size_t count) {
  cudaMemcpy(dst, src, count, cudaMemcpyDeviceToHost);
  check();
}

void copyToDevice(void* src, void* dst, std::size_t count) {
  cudaMemcpy(dst, src, count, cudaMemcpyHostToDevice);
  check();
}

template <typename T>
T* malloc(std::size_t size) {
  T* ptr{};
  cudaMalloc(&ptr, size*sizeof(T));
  check();
  return ptr;
}

/// Returns GPU memory page size
std::size_t getDevicePageSize() {
  std::size_t granularity = 0;
  CUmemAllocationProp prop = {};
  prop.type = CU_MEM_ALLOCATION_TYPE_PINNED;
  prop.location.type = CU_MEM_LOCATION_TYPE_DEVICE;
  prop.location.id = device::get();
  check(cuMemGetAllocationGranularity(&granularity, &prop, CU_MEM_ALLOC_GRANULARITY_MINIMUM));
  return granularity;
}

/// Returns `count` rounded up to be a multiple of `getDevicePageSize`
template <typename T>
std::size_t getPageAlignedCount(std::size_t count)
{
  const std::size_t page_size = getDevicePageSize();
  const std::size_t size = ((count * sizeof(T) - 1) / page_size + 1) * page_size;
  const std::size_t volume = size / sizeof(T);

  if (size % page_size != 0) {
    throw std::invalid_argument("Buffer size must be multiple of PAGE_SIZE");
  }

  return volume;
};

template <typename T>
class unique_ptr {
private:
  T* _ptr;

public:
  unique_ptr():
    _ptr{nullptr} { }

  unique_ptr(T* ptr):
    _ptr{ptr} { }

  unique_ptr(unique_ptr&& rhs):
    _ptr{rhs.release()} { }

  ~unique_ptr() {
    if (_ptr != nullptr) {
      cudaFree(_ptr);
      check();
    }
  }

  unique_ptr& operator=(unique_ptr&& rhs) {
    if (_ptr != nullptr) {
      cudaFree(_ptr);
      check();
    }
    _ptr = rhs.release();
    return *this;
  }

  void reset(T* ptr) {
    operator=(unique_ptr<T>(ptr));
  }

  unique_ptr& operator=(T* ptr) {
    reset(ptr);
    return *this;
  }

  T* release() {
    T* ptr = _ptr;
    _ptr = nullptr;
    return ptr;
  }

  const T* get() const {
    return _ptr;
  }
  T* get() {
    return _ptr;
  }

};

class Stream {
private:
  cudaStream_t _stream;

public:
  Stream(unsigned int flags) {
    cudaStreamCreateWithFlags(&_stream, flags);
    check();
  }

  ~Stream() {
    cudaStreamDestroy(_stream);
    check();
  }

  cudaStream_t& get() {
    return _stream;
  }

  void synchronize() {
    cudaStreamSynchronize(_stream);
    check();
  }

};

}

}

}

}

#endif
