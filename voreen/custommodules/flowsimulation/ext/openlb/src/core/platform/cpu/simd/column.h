/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Adrian Kummerlaender
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

#ifndef CPU_SIMD_COLUMN_H
#define CPU_SIMD_COLUMN_H

#include <memory>
#include <array>
#include <cstring>
#include <stdexcept>

#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <asm/unistd_64.h>

#include "core/platform/platform.h"
#include "core/serializer.h"
#include "communication/communicatable.h"

#include "pack.h"

namespace olb {

namespace cpu {

namespace simd {

/// Plain column for SIMD CPU targets
template<typename T>
class Column final : public AbstractColumn<T>
                   , public Serializable {
private:
  std::size_t _count;
  std::unique_ptr<T[]> _data;

public:
  using value_t = T;

  Column(std::size_t count):
    _count(count),
    _data(new T[count] {})
  { }

  Column(Column<T>&& rhs):
    _count(rhs._count),
    _data(rhs._data.release())
  { }

  virtual ~Column() = default;

  void resize(std::size_t count)
  {
    std::unique_ptr<T[]> data = std::unique_ptr<T[]>(new T[count] { });
    std::copy(_data.get(), _data.get() + std::min(_count, count), data.get());
    _data.swap(data);
    _count = count;
  }

  const T& operator[](std::size_t i) const override
  {
    return _data[i];
  }

  T& operator[](std::size_t i) override
  {
    return _data[i];
  }

  std::size_t size() const
  {
    return _count;
  }

  const T* data() const
  {
    return _data.get();
  }

  T* data()
  {
    return _data.get();
  }

  void setProcessingContext(ProcessingContext) { };

  /// Number of data blocks for the serializable interface
  std::size_t getNblock() const override;
  /// Binary size for the serializer
  std::size_t getSerializableSize() const override;
  /// Return a pointer to the memory of the current block and its size for the serializable interface
  bool* getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode) override;

};


const int PROT_RW = PROT_READ | PROT_WRITE;

template <typename T>
std::size_t getPageAlignedCount(std::size_t count)
{
  const std::size_t page_size = sysconf(_SC_PAGESIZE);
  const std::size_t size = ((count * sizeof(T) - 1) / page_size + 1) * page_size;
  const std::size_t volume = size / sizeof(T);

  if (size % page_size != 0) {
    throw std::invalid_argument("Buffer size must be multiple of PAGE_SIZE");
  }

  return volume;
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

  std::uint8_t* _buffer;

  T* _base;
  T* _f;

  std::ptrdiff_t _shift;

public:
  using value_t = T;

  CyclicColumn(std::size_t count):
    _count(getPageAlignedCount<T>(count)),
    _size(_count * sizeof(T)),
    _shift(0)
  {
  #ifdef __NR_memfd_create
    // Open anonymous file for physical lattice memory
    // Manual call of "memfd_create("openlb", MFD_CLOEXEC)" in case GLIB is old
    const int shm_file = syscall(__NR_memfd_create, "openlb", MFD_CLOEXEC);
  #else
    std::string shm_path = "/openlb_block_XXXXXX";
    const int shm_name = mkstemp(const_cast<char*>(shm_path.data()));
    if (shm_name != -1) {
      throw std::runtime_error("Could not generate unique shared memory object name");
    }
    // Open shared memory object as physical lattice memory
    const int shm_file = shm_open(shm_path.c_str(), O_CREAT | O_RDWR | O_EXCL | O_CLOEXEC, S_IRUSR | S_IWUSR);
    shm_unlink(shm_path.c_str());
  #endif
    if (shm_file == -1) {
      throw std::runtime_error("Failed to create shared memory object");
    }

    // Resize to fit lattice populations
    if (ftruncate(shm_file, _size) == -1) {
      throw std::runtime_error("Failed to resize shared memory object");
    }

    // Allocate virtual address space for q times two consecutive lattices
    _buffer = static_cast<std::uint8_t*>(
      mmap(NULL, 2 * _size, PROT_NONE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0));

    // Map single physical lattice into virtual address space
    mmap(_buffer,         _size, PROT_RW, MAP_SHARED | MAP_FIXED, shm_file, 0);
    mmap(_buffer + _size, _size, PROT_RW, MAP_SHARED | MAP_FIXED, shm_file, 0);

    // Store base pointer for reference
    _base = reinterpret_cast<T*>(_buffer);
    // Initialize shiftable f pointer to be used for lattice access
    _f = _base;
  }

  ~CyclicColumn() {
    munmap(_buffer, 2 * _size);
  }

  const T& operator[](std::size_t i) const override
  {
    return _f[i];
  }

  T& operator[](std::size_t i) override
  {
    return _f[i];
  }

  std::size_t size() const
  {
    return _count;
  }

  void refresh()
  {
    _f = _base + _shift;
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

  void setProcessingContext(ProcessingContext) { };

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

/// Declare cpu::sisd::Column as the AbstractColumn implementation for CPU SISD targets
template <typename T>
struct ImplementationOf<AbstractColumn<T>,Platform::CPU_SIMD> {
  using type = cpu::simd::Column<T>;
};

/// Declare cpu::sisd::CyclicColumn as the AbstractCyclicColumn implementation for CPU SISD targets
template <typename T>
struct ImplementationOf<AbstractCyclicColumn<T>,Platform::CPU_SIMD> {
  using type = cpu::simd::CyclicColumn<T>;
};


template <typename T>
class ConcreteCommunicatable<cpu::simd::CyclicColumn<T>> final : public Communicatable {
private:
  cpu::simd::CyclicColumn<T>& _column;

public:
  ConcreteCommunicatable(cpu::simd::CyclicColumn<T>& column):
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
    if (meta::is_aligned<T>(buffer)) {
      T* target = reinterpret_cast<T*>(buffer);
      const typename cpu::simd::Pack<T>::index_t* offsets = indices.data();

      for (std::size_t i=0;
           i < (indices.size() / cpu::simd::Pack<T>::size) * cpu::simd::Pack<T>::size;
           i       += cpu::simd::Pack<T>::size,
           target  += cpu::simd::Pack<T>::size,
           offsets += cpu::simd::Pack<T>::size) {
        cpu::simd::store(target, {&_column[0], offsets});
      }
      for (std::size_t i=(indices.size() / cpu::simd::Pack<T>::size) * cpu::simd::Pack<T>::size;
           i < indices.size();
           ++i) {
        *(target++) = _column[indices[i]];
      }
    } else {
      std::uint8_t* target = buffer;
      for (CellID index : indices) {
        std::memcpy(target, reinterpret_cast<const void*>(&_column[index]), sizeof(T));
        target += sizeof(T);
      }
    }
    return indices.size() * sizeof(T);
  }

  /// Deserialize data at locations `indices` to `buffer`
  std::size_t deserialize(ConstSpan<CellID> indices,
                          const std::uint8_t* buffer)
  {
    if (meta::is_aligned<T>(buffer)) {
      const T* source = reinterpret_cast<const T*>(buffer);
      const typename cpu::simd::Pack<T>::index_t* offsets = indices.data();

      for (std::size_t i=0;
           i < (indices.size() / cpu::simd::Pack<T>::size) * cpu::simd::Pack<T>::size;
           i       += cpu::simd::Pack<T>::size,
           source  += cpu::simd::Pack<T>::size,
           offsets += cpu::simd::Pack<T>::size) {
        cpu::simd::store(&_column[0], cpu::simd::Pack<T>{source}, offsets);
      }
      for (std::size_t i=(indices.size() / cpu::simd::Pack<T>::size) * cpu::simd::Pack<T>::size;
           i < indices.size();
           ++i) {
        _column[indices[i]] = *(source++);
      }
    } else {
      const std::uint8_t* source = buffer;
      for (CellID index : indices) {
        std::memcpy(reinterpret_cast<void*>(&_column[index]), source, sizeof(T));
        source += sizeof(T);
      }
    }
    return indices.size() * sizeof(T);
  }

};


}

#endif

#include "column.hh"
