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

#ifndef GPU_CUDA_COLUMN_HH
#define GPU_CUDA_COLUMN_HH

#include "column.h"

#include <thrust/gather.h>

#include <cstring>

namespace olb {

namespace gpu {

namespace cuda {

template<typename T>
void Column<T>::setProcessingContext(ProcessingContext context)
{
  switch (context) {
  case ProcessingContext::Evaluation:
    device::copyToHost(_deviceData.data().get(), _hostData.data(), size()*sizeof(T));
    return;
  case ProcessingContext::Simulation:
    device::copyToDevice(_hostData.data(), _deviceData.data().get(), size()*sizeof(T));
    return;
  }
}

template<typename T>
std::size_t Column<T>::getNblock() const
{
  return 2;
}

template<typename T>
std::size_t Column<T>::getSerializableSize() const
{
  return _count * sizeof(T) + sizeof(std::size_t);
}

template<typename T>
bool* Column<T>::getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode)
{
  std::size_t currentBlock = 0;
  bool* dataPtr = nullptr;

  registerVar(iBlock, sizeBlock, currentBlock, dataPtr, _count);
  if (loadingMode && iBlock == 1) {
    resize(_count);
  }
  registerVar(iBlock, sizeBlock, currentBlock, dataPtr, *data(), _count);

  return dataPtr;
}


template<typename T>
void CyclicColumn<T>::setProcessingContext(ProcessingContext context)
{
  switch (context) {
  case ProcessingContext::Evaluation:
    device::copyToHost(_devicePopulation, _hostData.get(), _size);
    return;
  case ProcessingContext::Simulation:
    device::copyToDevice(_hostData.get(), _devicePopulation, _size);
    return;
  }
}

template<typename T>
std::size_t CyclicColumn<T>::getNblock() const
{
  return 3;
}

template<typename T>
std::size_t CyclicColumn<T>::getSerializableSize() const
{
  return _count * sizeof(T) + sizeof(std::ptrdiff_t) + sizeof(std::size_t);
}

template<typename T>
bool* CyclicColumn<T>::getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode)
{
  std::size_t currentBlock = 0;
  bool* dataPtr = nullptr;

  registerVar(iBlock, sizeBlock, currentBlock, dataPtr, _shift);
  registerVar(iBlock, sizeBlock, currentBlock, dataPtr, _count);
  registerVar(iBlock, sizeBlock, currentBlock, dataPtr, *_hostData.get(), _count);

  return dataPtr;
}

}

}

}

#endif
