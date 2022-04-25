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

#ifndef PLATFORM_GPU_CUDA_COMMUNICATOR_HH
#define PLATFORM_GPU_CUDA_COMMUNICATOR_HH

#include "communication/mpiRequest.h"
#include "communication/communicatable.h"
#include "communication/superCommunicationTagCoordinator.h"

#include "registry.h"

#include <thrust/device_vector.h>

namespace olb {

namespace gpu {

namespace cuda {

namespace kernel {

/// CUDA kernel for gathering FIELD data of lattice at indices into buffer
template <typename CONTEXT, typename FIELD>
void gather_field(CONTEXT lattice,
                  const CellID* indices, std::size_t nIndices,
                  typename FIELD::template value_type<typename CONTEXT::value_t>* buffer) __global__ {
  const CellID iIndex = blockIdx.x * blockDim.x + threadIdx.x;
  if (!(iIndex < nIndices)) {
    return;
  }
  auto* field = lattice.template getField<FIELD>();
  for (unsigned iD=0; iD < CONTEXT::descriptor_t::template size<FIELD>(); ++iD) {
    buffer[iD*nIndices+iIndex] = field[iD][indices[iIndex]];
  }
}

/// CUDA kernel for copying FIELD data of sourceLattice at sourceIndices into targetLattice at targetIndices
template <typename SOURCE, typename TARGET, typename FIELD>
void copy_field(SOURCE sourceLattice, TARGET targetLattice,
                const CellID* sourceIndices,
                const CellID* targetIndices,
                std::size_t nIndices) __global__ {
  const CellID iIndex = blockIdx.x * blockDim.x + threadIdx.x;
  if (!(iIndex < nIndices)) {
    return;
  }
  auto* source = sourceLattice.template getField<FIELD>();
  auto* target = targetLattice.template getField<FIELD>();
  for (unsigned iD=0; iD < SOURCE::descriptor_t::template size<FIELD>(); ++iD) {
    target[iD][targetIndices[iIndex]] = source[iD][sourceIndices[iIndex]];
  }
}

/// CUDA kernel for gathering fields at indices into buffer
__global__ void gather_any_fields(AnyDeviceFieldArrayD* fields, std::size_t nFields,
                                  const CellID* indices, std::size_t nIndices,
                                  std::uint8_t* buffer) {
  const CellID iIndex = blockIdx.x * blockDim.x + threadIdx.x;
  if (!(iIndex < nIndices)) {
    return;
  }
  for (unsigned iField=0; iField < nFields; ++iField) {
    auto& field = fields[iField];
    for (unsigned iD=0; iD < field.column_count; ++iD) {
      memcpy(buffer + (iD*nIndices + iIndex)*field.element_size,
             field[iD] + indices[iIndex]*field.element_size,
             field.element_size);
    }
    buffer += nIndices*field.column_count*field.element_size;
  }
}

/// CUDA kernel for copying sourceFields at sourceIndices to targetFields at targetIndices
/**
 * source and target fields may be of different block lattices but must represent the
 * same field types in the same sequence
 **/
__global__ void copy_any_fields(AnyDeviceFieldArrayD* sourceFields,
                                AnyDeviceFieldArrayD* targetFields,
                                std::size_t nFields,
                                const CellID* sourceIndices,
                                const CellID* targetIndices,
                                std::size_t nIndices) {
  const CellID iIndex = blockIdx.x * blockDim.x + threadIdx.x;
  if (!(iIndex < nIndices)) {
    return;
  }
  for (unsigned iField=0; iField < nFields; ++iField) {
    auto& sourceField = sourceFields[iField];
    auto& targetField = targetFields[iField];
    for (unsigned iD=0; iD < sourceField.column_count; ++iD) {
      memcpy(targetField[iD] + targetIndices[iIndex]*sourceField.element_size,
             sourceField[iD] + sourceIndices[iIndex]*sourceField.element_size,
             sourceField.element_size);
    }
  }
}

/// CUDA kernel for scattering FIELD data in buffer to indices in lattice
template <typename CONTEXT, typename FIELD>
void scatter_field(CONTEXT lattice,
                   const CellID* indices, std::size_t nIndices,
                   typename FIELD::template value_type<typename CONTEXT::value_t>* buffer) __global__ {
  const CellID iIndex = blockIdx.x * blockDim.x + threadIdx.x;
  if (!(iIndex < nIndices)) {
    return;
  }
  auto* field = lattice.template getField<FIELD>();
  for (unsigned iD=0; iD < CONTEXT::descriptor_t::template size<FIELD>(); ++iD) {
    field[iD][indices[iIndex]] = buffer[iD*nIndices+iIndex];
  }
}

/// CUDA kernel for scattering fields in buffer to indices in lattice
__global__ void scatter_any_fields(AnyDeviceFieldArrayD* fields, std::size_t nFields,
                                   const CellID* indices, std::size_t nIndices,
                                   std::uint8_t* buffer) {
  const CellID iIndex = blockIdx.x * blockDim.x + threadIdx.x;
  if (!(iIndex < nIndices)) {
    return;
  }
  for (unsigned iField=0; iField < nFields; ++iField) {
    auto& field = fields[iField];
    for (unsigned iD=0; iD < field.column_count; ++iD) {
      memcpy(field[iD] + indices[iIndex]*field.element_size,
             buffer + (iD*nIndices + iIndex)*field.element_size,
             field.element_size);
    }
    buffer += nIndices*field.column_count*field.element_size;
  }
}

}

/// Blocking gather of FIELD at given indices into buffer
template <typename FIELD, typename CONTEXT>
void gather_field(CONTEXT& lattice, const thrust::device_vector<CellID>& indices, std::uint8_t* buffer) {
  const auto block_size = 32;
  const auto block_count = (indices.size() + block_size - 1) / block_size;
  kernel::gather_field<CONTEXT,FIELD><<<block_count,block_size>>>(
    lattice,
    indices.data().get(), indices.size(),
    reinterpret_cast<typename FIELD::template value_type<typename CONTEXT::value_t>*>(buffer));
  device::check();
}

/// Non-blocking gather of FIELD at given indices into buffer
template <typename FIELD, typename CONTEXT>
void async_gather_field(cudaStream_t stream,
                        CONTEXT& lattice,
                        const thrust::device_vector<CellID>& indices,
                        std::uint8_t* buffer) {
  const auto block_size = 32;
  const auto block_count = (indices.size() + block_size - 1) / block_size;
  kernel::gather_field<CONTEXT,FIELD><<<block_count,block_size,0,stream>>>(
    lattice,
    indices.data().get(), indices.size(),
    reinterpret_cast<typename FIELD::template value_type<typename CONTEXT::value_t>*>(buffer));
  device::check();
}

/// Non-blocking copy of FIELD at given indices from sourceLattice to targetLattice
template <typename FIELD, typename SOURCE, typename TARGET>
void async_copy_field(cudaStream_t stream,
                      SOURCE& sourceLattice,
                      TARGET& targetLattice,
                      const thrust::device_vector<CellID>& sourceIndices,
                      const thrust::device_vector<CellID>& targetIndices) {
  const auto block_size = 32;
  const auto block_count = (sourceIndices.size() + block_size - 1) / block_size;
  kernel::copy_field<SOURCE,TARGET,FIELD><<<block_count,block_size,0,stream>>>(
    sourceLattice,
    targetLattice,
    sourceIndices.data().get(),
    targetIndices.data().get(),
    sourceIndices.size());
  device::check();
}


/// Blocking gather of fields at given indices into buffer
void gather_any_fields(thrust::device_vector<AnyDeviceFieldArrayD>& fields,
                       const thrust::device_vector<CellID>& indices,
                       std::uint8_t* buffer) {
  const auto block_size = 32;
  const auto block_count = (indices.size() + block_size - 1) / block_size;
  kernel::gather_any_fields<<<block_count,block_size>>>(
    fields.data().get(), fields.size(),
    indices.data().get(), indices.size(),
    buffer);
  device::check();
}

/// Non-blocking gather of fields at given indices into buffer
void async_gather_any_fields(cudaStream_t stream,
                             thrust::device_vector<AnyDeviceFieldArrayD>& fields,
                             const thrust::device_vector<CellID>& indices,
                             std::uint8_t* buffer) {
  const auto block_size = 32;
  const auto block_count = (indices.size() + block_size - 1) / block_size;
  kernel::gather_any_fields<<<block_count,block_size,0,stream>>>(
    fields.data().get(), fields.size(),
    indices.data().get(), indices.size(),
    buffer);
  device::check();
}

/// Non-blocking copy of fields at given indices from sourceIndices to targetIndices
void async_copy_any_fields(cudaStream_t stream,
                           thrust::device_vector<AnyDeviceFieldArrayD>& sourceFields,
                           thrust::device_vector<AnyDeviceFieldArrayD>& targetFields,
                           const thrust::device_vector<CellID>& sourceIndices,
                           const thrust::device_vector<CellID>& targetIndices) {
  const auto block_size = 32;
  const auto block_count = (sourceIndices.size() + block_size - 1) / block_size;
  kernel::copy_any_fields<<<block_count,block_size,0,stream>>>(
    sourceFields.data().get(),  targetFields.data().get(),  sourceFields.size(),
    sourceIndices.data().get(), targetIndices.data().get(), sourceIndices.size());
  device::check();
}

/// Blocking scatter of FIELD data in buffer to given indices
template <typename FIELD, typename CONTEXT>
void scatter_field(CONTEXT& lattice, const thrust::device_vector<CellID>& indices, std::uint8_t* buffer) {
  const auto block_size = 32;
  const auto block_count = (indices.size() + block_size - 1) / block_size;
  kernel::scatter_field<CONTEXT,FIELD><<<block_count,block_size>>>(
    lattice,
    indices.data().get(), indices.size(),
    reinterpret_cast<typename FIELD::template value_type<typename CONTEXT::value_t>*>(buffer));
  device::check();
}

/// Non-blocking scatter of FIELD data in buffer to given indices
template <typename FIELD, typename CONTEXT>
void async_scatter_field(cudaStream_t stream,
                         CONTEXT& lattice,
                         const thrust::device_vector<CellID>& indices,
                         std::uint8_t* buffer) {
  const auto block_size = 32;
  const auto block_count = (indices.size() + block_size - 1) / block_size;
  kernel::scatter_field<CONTEXT,FIELD><<<block_count,block_size,0,stream>>>(
    lattice,
    indices.data().get(), indices.size(),
    reinterpret_cast<typename FIELD::template value_type<typename CONTEXT::value_t>*>(buffer));
  device::check();
}

/// Blocking scatter of fields data in buffer to given indices
void scatter_any_fields(thrust::device_vector<AnyDeviceFieldArrayD>& fields,
                        const thrust::device_vector<CellID>& indices,
                        std::uint8_t* buffer) {
  const auto block_size = 32;
  const auto block_count = (indices.size() + block_size - 1) / block_size;
  kernel::scatter_any_fields<<<block_count,block_size>>>(
    fields.data().get(), fields.size(),
    indices.data().get(), indices.size(),
    buffer);
  device::check();
}

/// Non-blocking scatter of fields data in buffer to given indices
void async_scatter_any_fields(cudaStream_t stream,
                              thrust::device_vector<AnyDeviceFieldArrayD>& fields,
                              const thrust::device_vector<CellID>& indices,
                              std::uint8_t* buffer) {
  const auto block_size = 32;
  const auto block_count = (indices.size() + block_size - 1) / block_size;
  kernel::scatter_any_fields<<<block_count,block_size,0,stream>>>(
    fields.data().get(), fields.size(),
    indices.data().get(), indices.size(),
    buffer);
  device::check();
}

}

}

/// Wrapper for a local plain-copy block communication request
/**
 * Currently only supports GPU-internal communication.
 * This will need to be improved for node-level heterogeneous simulations.
 **/
template <typename T, typename DESCRIPTOR>
class ConcreteBlockCommunicator<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>>::CopyTask {
private:
  thrust::device_vector<gpu::cuda::AnyDeviceFieldArrayD> _sourceFields;
  thrust::device_vector<gpu::cuda::AnyDeviceFieldArrayD> _targetFields;

  const bool _onlyPopulationField;

  const thrust::device_vector<CellID> _targetCells;
  const thrust::device_vector<CellID> _sourceCells;

  ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>& _target;
  ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>& _source;

  std::unique_ptr<gpu::cuda::device::Stream> _stream;

public:
  CopyTask(
    const std::vector<std::type_index>& fields,
    const std::vector<CellID>& targetCells, ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>& target,
    const std::vector<CellID>& sourceCells, ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>& source):
    _sourceFields(source.getDataRegistry().deviceFieldArrays(fields)),
    _targetFields(target.getDataRegistry().deviceFieldArrays(fields)),
    _onlyPopulationField(fields.size() == 1 && fields[0] == typeid(descriptors::POPULATION)),
    _targetCells(targetCells),
    _sourceCells(sourceCells),
    _target(target),
    _source(source),
    _stream(std::make_unique<gpu::cuda::device::Stream>(cudaStreamNonBlocking))
  {
    OLB_ASSERT(_sourceCells.size() == _targetCells.size(),
               "Source cell count must match target cell count");
  }

  ~CopyTask()
  {
    wait();
  }

  void copy()
  {
    if (_onlyPopulationField) {
      gpu::cuda::DeviceContext<T,DESCRIPTOR> sourceLattice(_source);
      gpu::cuda::DeviceContext<T,DESCRIPTOR> targetLattice(_target);
      gpu::cuda::async_copy_field<descriptors::POPULATION>(_stream->get(), sourceLattice, targetLattice, _sourceCells, _targetCells);
    } else {
      gpu::cuda::async_copy_any_fields(_stream->get(), _sourceFields, _targetFields, _sourceCells, _targetCells);
    }
  }

  void wait()
  {
    _stream->synchronize();
  }

};

#ifdef PARALLEL_MODE_MPI

/// Wrapper for a non-blocking block propagation send request
template <typename T, typename DESCRIPTOR>
class ConcreteBlockCommunicator<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>>::SendTask {
private:
  thrust::device_vector<gpu::cuda::AnyDeviceFieldArrayD> _fields;
  const bool _onlyPopulationField;

  const thrust::device_vector<CellID> _cells;

  ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>& _source;

  gpu::cuda::device::unique_ptr<std::uint8_t> _buffer;
  std::unique_ptr<MpiSendRequest> _request;

  std::unique_ptr<gpu::cuda::device::Stream> _stream;

public:
  SendTask(MPI_Comm comm, int tag, int rank,
           const std::vector<std::type_index>& fields,
           const std::vector<CellID>& cells,
           ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>& block):
    _fields(block.getDataRegistry().deviceFieldArrays(fields)),
    _onlyPopulationField(fields.size() == 1 && fields[0] == typeid(descriptors::POPULATION)),
    _cells(cells),
    _source(block),
    _stream(std::make_unique<gpu::cuda::device::Stream>(cudaStreamNonBlocking))
  {
    std::size_t size = 0;
    for (auto& field : fields) {
      size += _source.getCommunicatable(field).size(cells);
    }
    _buffer = gpu::cuda::device::malloc<std::uint8_t>(size);
    _request = std::make_unique<MpiSendRequest>(
      _buffer.get(), size, rank, tag, comm);
  }

  ~SendTask()
  {
    _stream->synchronize();
    wait();
  }

  void prepare()
  {
    if (_onlyPopulationField) {
      gpu::cuda::DeviceContext<T,DESCRIPTOR> lattice(_source);
      gpu::cuda::async_gather_field<descriptors::POPULATION>(_stream->get(), lattice, _cells, _buffer.get());
    } else {
      gpu::cuda::async_gather_any_fields(_stream->get(), _fields, _cells, _buffer.get());
    }
  }

  void send()
  {
    _stream->synchronize();
    _request->start();
  }

  void wait()
  {
    _request->wait();
  }
};

/// Wrapper for a non-blocking block propagation receive request
template <typename T, typename DESCRIPTOR>
class ConcreteBlockCommunicator<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>>::RecvTask {
private:
  const int _tag;
  const int _rank;

  thrust::device_vector<gpu::cuda::AnyDeviceFieldArrayD> _fields;
  const bool _onlyPopulationField;

  const thrust::device_vector<CellID> _cells;

  ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>& _target;

  gpu::cuda::device::unique_ptr<std::uint8_t> _buffer;
  std::unique_ptr<MpiRecvRequest> _request;

  std::unique_ptr<gpu::cuda::device::Stream> _stream;

public:
  /// Manual replacement for std::reference_wrapper<RecvTask>
  /**
   * Used to track pending receive requests in std::set.
   *
   * This is a workaround for problematic external definition of
   * dependently-typed comparision operators for nested classes.
   * Reconsider as soon as depending on C++17 is allowed.
   **/
  class ref {
  private:
    RecvTask& _task;
  public:
    ref(std::unique_ptr<RecvTask>& task): _task(*task) { };

    RecvTask* operator->() const
    {
      return &_task;
    }

    bool operator <(const ref& rhs) const
    {
      return _task < rhs._task;
    }
  };

  RecvTask(MPI_Comm comm, int tag, int rank,
           const std::vector<std::type_index>& fields,
           const std::vector<CellID>& cells,
           ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>& block):
    _tag(tag),
    _rank(rank),
    _fields(block.getDataRegistry().deviceFieldArrays(fields)),
    _onlyPopulationField(fields.size() == 1 && fields[0] == typeid(descriptors::POPULATION)),
    _cells(cells),
    _target(block),
    _stream(std::make_unique<gpu::cuda::device::Stream>(cudaStreamNonBlocking))
  {
    std::size_t size = 0;
    for (auto& field : fields) {
      size += _target.getCommunicatable(field).size(cells);
    }
    _buffer = gpu::cuda::device::malloc<std::uint8_t>(size);
    _request = std::make_unique<MpiRecvRequest>(
      _buffer.get(), size, _rank, _tag, comm);
  }

  ~RecvTask()
  {
    wait();
  }

  bool operator<(const RecvTask& rhs) const
  {
    return  _rank  < rhs._rank
        || (_rank == rhs._rank && _tag < rhs._tag);
  }

  void receive()
  {
    _request->start();
  };

  bool isDone()
  {
    return _request->isDone();
  }

  void unpack()
  {
    if (_onlyPopulationField) {
      gpu::cuda::DeviceContext<T,DESCRIPTOR> lattice(_target);
      gpu::cuda::async_scatter_field<descriptors::POPULATION>(_stream->get(), lattice, _cells, _buffer.get());
    } else {
      gpu::cuda::async_scatter_any_fields(_stream->get(), _fields, _cells, _buffer.get());
    }
  }

  void wait()
  {
    _stream->synchronize();
  }

};

#endif // PARALLEL_MODE_MPI

template <typename T, typename DESCRIPTOR>
ConcreteBlockCommunicator<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>>::ConcreteBlockCommunicator(
  SuperLattice<T,DESCRIPTOR>& super,
  LoadBalancer<T>& loadBalancer,
#ifdef PARALLEL_MODE_MPI
  SuperCommunicationTagCoordinator<T>& tagCoordinator,
  MPI_Comm comm,
#endif
  int iC,
  const BlockCommunicationNeighborhood<T,DESCRIPTOR::d>& neighborhood):
  _iC(iC)
#ifdef PARALLEL_MODE_MPI
, _mpiCommunicator(comm)
#endif
{
#ifdef PARALLEL_MODE_MPI
  neighborhood.forRemoteNeighbors([&](int remoteC) {
    if (!neighborhood.getCellsOutboundTo(remoteC).empty()) {
      _sendTasks.emplace_back(std::make_unique<SendTask>(
        _mpiCommunicator, tagCoordinator.get(remoteC, loadBalancer.glob(_iC)),
        loadBalancer.rank(remoteC),
        neighborhood.getFieldsCommonWith(remoteC),
        neighborhood.getCellsOutboundTo(remoteC),
        super.template getBlock<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>>(_iC)));
    }
    if (!neighborhood.getCellsInboundFrom(remoteC).empty()) {
      _recvTasks.emplace_back(std::make_unique<RecvTask>(
        _mpiCommunicator, tagCoordinator.get(remoteC, loadBalancer.glob(_iC)),
        loadBalancer.rank(remoteC),
        neighborhood.getFieldsCommonWith(remoteC),
        neighborhood.getCellsInboundFrom(remoteC),
        super.template getBlock<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>>(_iC)));
    }
  });
#endif // PARALLEL_MODE_MPI

  neighborhood.forLocalNeighbors([&](int localC) {
    if (!neighborhood.getCellsInboundFrom(localC).empty()) {
      _copyTasks.emplace_back(std::make_unique<CopyTask>(
        neighborhood.getFieldsCommonWith(localC),
        neighborhood.getCellsInboundFrom(localC),   super.template getBlock<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>>(_iC),
        neighborhood.getCellsRequestedFrom(localC), super.template getBlock<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>>(loadBalancer.loc(localC))));
    }
  });
}

template <typename T, typename DESCRIPTOR>
void ConcreteBlockCommunicator<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>>::copy()
{
  for (auto& task : _copyTasks) {
    task->copy();
  }
}

#ifdef PARALLEL_MODE_MPI

template <typename T, typename DESCRIPTOR>
void ConcreteBlockCommunicator<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>>::receive()
{
  for (auto& task : _recvTasks) {
    task->receive();
  }
}

template <typename T, typename DESCRIPTOR>
void ConcreteBlockCommunicator<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>>::send()
{
  for (auto& task : _sendTasks) {
    task->prepare();
  }
  for (auto& task : _sendTasks) {
    task->send();
  }
}

template <typename T, typename DESCRIPTOR>
void ConcreteBlockCommunicator<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>>::unpack()
{
  std::set<typename RecvTask::ref> pending(_recvTasks.begin(), _recvTasks.end());
  while (!pending.empty()) {
    auto task_iterator = pending.begin();
    while (task_iterator != pending.end()) {
      auto& task = *task_iterator;
      if (task->isDone()) {
        task->unpack();
        task_iterator = pending.erase(task_iterator);
      }
      else {
        ++task_iterator;
      }
    }
  }
}

template <typename T, typename DESCRIPTOR>
void ConcreteBlockCommunicator<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>>::wait()
{
  for (auto& task : _recvTasks) {
    task->wait();
  }
  for (auto& task : _sendTasks) {
    task->wait();
  }
  for (auto& task : _copyTasks) {
    task->wait();
  }
}

#else // PARALLEL_MODE_MPI

template <typename T, typename DESCRIPTOR>
void ConcreteBlockCommunicator<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>>::wait()
{
  for (auto& task : _copyTasks) {
    task->wait();
  }
}

#endif

}

#endif
