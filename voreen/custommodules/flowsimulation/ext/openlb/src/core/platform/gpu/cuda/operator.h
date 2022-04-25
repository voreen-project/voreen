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

#ifndef GPU_CUDA_OPERATOR_H
#define GPU_CUDA_OPERATOR_H

#include "core/operator.h"

#include "mask.h"
#include "context.h"

namespace olb {

/// Implementations of GPU specifics
namespace gpu {

/// Implementations of Nvidia CUDA specifics
namespace cuda {

/// Masked application of DYNAMICS::apply for use in kernel::call_operators
template <typename DYNAMICS>
class MaskedCollision {
private:
  /// Pointer to on-device parameters structure to be passed to DYNAMICS::apply
  typename DYNAMICS::ParametersD* _parameters;
  /// Pointer to on-device mask array
  bool* _mask;

  template <typename T, typename DESCRIPTOR>
  CellStatistic<T> apply(DeviceContext<T,DESCRIPTOR>& lattice, CellID iCell) __device__ {
    DataOnlyCell<T,DESCRIPTOR> cell(lattice, iCell);
    return DYNAMICS().apply(cell, *_parameters);
  }

public:
  /// Constructor (commonly called on the host side)
  /**
   * See e.g. getFusedCollisionO
   **/
  MaskedCollision(typename DYNAMICS::ParametersD* parameters, bool* mask) any_platform:
    _parameters{parameters},
    _mask{mask}
  { }

  /// Chainable call operator for use in kernel::call_operators
  /**
   * Returns true iff MaskedCollision applies to iCell, enabling easy chaining
   * by a fold expression to yield a fused collision kernel.
   **/
  template <typename T, typename DESCRIPTOR>
  bool operator()(DeviceContext<T,DESCRIPTOR>& lattice, CellID iCell) __device__ {
    if (_mask[iCell]) {
      apply(lattice, iCell);
      return true;
    }
    return false;
  }

  /// Chainable call operator with statistics storage
  template <typename T, typename DESCRIPTOR>
  bool operator()(DeviceContext<T,DESCRIPTOR>& lattice, CellID iCell, CellStatistic<T>& statistic) __device__ {
    if (_mask[iCell]) {
      statistic = apply(lattice, iCell);
      return true;
    }
    return false;
  }

};

/// List-based application of DYNAMICS::apply for use in kernel::call_list_operators
template <typename DYNAMICS>
class ListedCollision {
private:
  typename DYNAMICS::ParametersD _parameters;

public:
  ListedCollision(typename DYNAMICS::ParametersD& parameters) __host__:
    _parameters{parameters}
  { }

  template <typename T, typename DESCRIPTOR>
  bool operator()(DeviceContext<T,DESCRIPTOR>& lattice, CellID iCell) __device__ {
    DataOnlyCell<T,DESCRIPTOR> cell(lattice, iCell);
    DYNAMICS().apply(cell, _parameters);
    return true;
  }

  template <typename T, typename DESCRIPTOR>
  bool operator()(DeviceContext<T,DESCRIPTOR>& lattice, CellID iCell, CellStatistic<T>& statistic) __device__ {
    DataOnlyCell<T,DESCRIPTOR> cell(lattice, iCell);
    statistic = DYNAMICS().apply(cell, _parameters);
    return true;
  }

};


/// Masked application of OPERATOR::apply
template <typename OPERATOR>
class MaskedPostProcessor {
private:
  /// Pointer to on-device mask array
  bool* _mask;

public:
  MaskedPostProcessor(bool* mask) any_platform:
    _mask{mask}
  { }

  template <typename T, typename DESCRIPTOR>
  bool operator()(DeviceBlockLattice<T,DESCRIPTOR>& lattice, CellID iCell) __device__ {
    if (_mask[iCell]) {
      Cell<T,DESCRIPTOR> cell(lattice, iCell);
      OPERATOR().apply(cell);
      return true;
    }
    return false;
  }

};

/// List-based application of OPERATOR::apply
/**
 * Most common approach to calling post processors on device data
 **/
template <typename OPERATOR>
struct ListedPostProcessor {
  template <typename T, typename DESCRIPTOR>
  bool operator()(DeviceBlockLattice<T,DESCRIPTOR>& lattice, CellID iCell) __device__ {
    Cell<T,DESCRIPTOR> cell(lattice, iCell);
    OPERATOR().apply(cell);
    return true;
  }

};

/// List-based application of OPERATOR::apply with parameters
template <typename OPERATOR>
class ListedPostProcessorWithParameters {
private:
  typename OPERATOR::ParametersD _parameters;

public:
  ListedPostProcessorWithParameters(typename OPERATOR::ParametersD& parameters) __host__:
    _parameters{parameters}
  { }

  template <typename T, typename DESCRIPTOR>
  bool operator()(DeviceBlockLattice<T,DESCRIPTOR>& lattice, CellID iCell) __device__ {
    Cell<T,DESCRIPTOR> cell(lattice, iCell);
    OPERATOR().apply(cell, _parameters);
    return true;
  }

};

/// CUDA kernels to execute collisions and post processors
namespace kernel {

/// CUDA kernel for applying purely local collision steps
template <typename CONTEXT, typename... OPERATORS>
void call_operators(CONTEXT lattice, bool* subdomain, OPERATORS... ops) __global__ {
  const CellID iCell = blockIdx.x * blockDim.x + threadIdx.x;
  if (!(iCell < lattice.getNcells()) || !subdomain[iCell]) {
    return;
  }
  (ops(lattice, iCell) || ... );
}

/// CUDA kernel for applying purely local collision steps while tracking statistics
/**
 * Statistics data is reduced by StatisticsPostProcessor
 **/
template <typename CONTEXT, typename... OPERATORS>
void call_operators_with_statistics(CONTEXT lattice, bool* subdomain, OPERATORS... ops) __global__ {
  const CellID iCell = blockIdx.x * blockDim.x + threadIdx.x;
  if (!(iCell < lattice.getNcells()) || !subdomain[iCell]) {
    return;
  }
  typename CONTEXT::value_t** statistic = lattice.template getField<descriptors::STATISTIC>();
  int* statisticGenerated = lattice.template getField<descriptors::STATISTIC_GENERATED>()[0];
  CellStatistic<typename CONTEXT::value_t> cellStatistic{-1, -1};
  (ops(lattice, iCell, cellStatistic) || ... );
  if (cellStatistic) {
    statisticGenerated[iCell] = 1;
    statistic[0][iCell] = cellStatistic.rho;
    statistic[1][iCell] = cellStatistic.uSqr;
  } else {
    statisticGenerated[iCell] = 0;
    statistic[0][iCell] = 0;
    statistic[1][iCell] = 0;
  }
}

/// CUDA kernel for applying generic OPERATORS with OperatorScope::PerCell or ListedCollision
template <typename CONTEXT, typename... OPERATORS>
void call_list_operators(CONTEXT lattice,
                         const CellID* indices, std::size_t nIndices,
                         OPERATORS... ops) __global__ {
  const std::size_t iIndex = blockIdx.x * blockDim.x + threadIdx.x;
  if (!(iIndex < nIndices)) {
    return;
  }
  (ops(lattice, indices[iIndex]) || ... );
}

/// CUDA kernel for applying ListedCollision
/**
 * Statistics data is reduced by StatisticsPostProcessor
 **/
template <typename CONTEXT, typename... OPERATORS>
void call_list_operators_with_statistics(CONTEXT lattice,
                                         const CellID* indices, std::size_t nIndices,
                                         OPERATORS... ops) __global__ {
  const std::size_t iIndex = blockIdx.x * blockDim.x + threadIdx.x;
  if (!(iIndex < nIndices)) {
    return;
  }
  typename CONTEXT::value_t** statistic = lattice.template getField<descriptors::STATISTIC>();
  int* statisticGenerated = lattice.template getField<descriptors::STATISTIC_GENERATED>()[0];
  CellStatistic<typename CONTEXT::value_t> cellStatistic{-1, -1};
  (ops(lattice, indices[iIndex], cellStatistic) || ... );
  if (cellStatistic) {
    statisticGenerated[indices[iIndex]] = 1;
    statistic[0][indices[iIndex]] = cellStatistic.rho;
    statistic[1][indices[iIndex]] = cellStatistic.uSqr;
  } else {
    statisticGenerated[indices[iIndex]] = 0;
    statistic[0][indices[iIndex]] = 0;
    statistic[1][indices[iIndex]] = 0;
  }
}

}

/// Apply masked collision operators to lattice
/**
 * ARGS are instances of MaskedCollision or DynamicDispatchCollision
 **/
template <typename CONTEXT, typename... ARGS>
void call_operators(CONTEXT& lattice, bool* subdomain, ARGS&&... args) {
  const auto block_size = 32;
  const auto block_count = (lattice.getNcells() + block_size - 1) / block_size;
  kernel::call_operators<CONTEXT,ARGS...><<<block_count,block_size>>>(
    lattice, subdomain, std::forward<decltype(args)>(args)...);
  device::check();
}

/// Apply masked collision operators to lattice (async)
template <typename CONTEXT, typename... ARGS>
void async_call_operators(cudaStream_t stream, CONTEXT& lattice, bool* subdomain, ARGS&&... args) {
  const auto block_size = 32;
  const auto block_count = (lattice.getNcells() + block_size - 1) / block_size;
  kernel::call_operators<CONTEXT,ARGS...><<<block_count,block_size,0,stream>>>(
    lattice, subdomain, std::forward<decltype(args)>(args)...);
  device::check();
}

/// Apply masked collision operators to lattice while tracking statistics
/**
 * ARGS are instances of MaskedCollision or DynamicDispatchCollision
 **/
template <typename CONTEXT, typename... ARGS>
void call_operators_with_statistics(CONTEXT& lattice, bool* subdomain, ARGS&&... args) {
  const auto block_size = 32;
  const auto block_count = (lattice.getNcells() + block_size - 1) / block_size;
  kernel::call_operators_with_statistics<CONTEXT,ARGS...><<<block_count,block_size>>>(
    lattice, subdomain, std::forward<decltype(args)>(args)...);
  device::check();
}

/// Apply masked collision operators to lattice while tracking statistics (async)
template <typename CONTEXT, typename... ARGS>
void async_call_operators_with_statistics(cudaStream_t stream, CONTEXT& lattice, bool* subdomain, ARGS&&... args) {
  const auto block_size = 32;
  const auto block_count = (lattice.getNcells() + block_size - 1) / block_size;
  kernel::call_operators_with_statistics<CONTEXT,ARGS...><<<block_count,block_size,0,stream>>>(
    lattice, subdomain, std::forward<decltype(args)>(args)...);
  device::check();
}

/// Helper for constructing fused collision operators
/**
 * This is a convenient way for improving performance by injecting
 * application knowledge. E.g. if the lattice contains primarily BGK
 * and BounceBack dynamics this can be declared using:
 *
 * \code{.cpp}
 * superLattice.forBlocksOnPlatform<Platform::GPU_CUDA>([](auto& block) {
 *   block.setCollisionO(
 *     gpu::cuda::getFusedCollisionO<T,DESCRIPTOR,
 *                                   BGKdynamics<T,DESCRIPTOR>,
 *                                   BounceBack<T,DESCRIPTOR>,
 *                                   BounceBackVelocity<T,DESCRIPTOR>>());
 * });
 * \endcode
 **/
template <typename T, typename DESCRIPTOR, typename... DYNAMICS>
std::function<void(ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>&)>
getFusedCollisionO() {
  return [](ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>& block) {
    bool* subdomain = block.template getData<CollisionSubdomainMask>().deviceData();
    DeviceContext<T,DESCRIPTOR> lattice(block);
    if (block.statisticsEnabled()) {
      call_operators_with_statistics(
        lattice,
        subdomain,
        MaskedCollision<DYNAMICS>{
          block.template getData<DynamicsParameters<DYNAMICS>>().deviceData(),
          block.template getData<DynamicsMask<DYNAMICS>>().deviceData()
        }...,
        DynamicDispatchCollision{});
    } else {
      call_operators(
        lattice,
        subdomain,
        MaskedCollision<DYNAMICS>{
          block.template getData<DynamicsParameters<DYNAMICS>>().deviceData(),
          block.template getData<DynamicsMask<DYNAMICS>>().deviceData()
        }...,
        DynamicDispatchCollision{});
    }
  };
}

/// Apply operators to listed cell indices
/**
 * Used to call post processors in ConcreteBlockO with OperatorScope::PerCell
 **/
template <typename CONTEXT, typename... ARGS>
void call_list_operators(CONTEXT& lattice,
                         const thrust::device_vector<CellID>& cells,
                         ARGS&&... args) {
  const auto block_size = 32;
  const auto block_count = (cells.size() + block_size - 1) / block_size;
  kernel::call_list_operators<CONTEXT,ARGS...><<<block_count, block_size>>>(
    lattice,
    cells.data().get(), cells.size(),
    std::forward<decltype(args)>(args)...);
  device::check();
}

/// Apply operators to listed cell indices (async version)
template <typename CONTEXT, typename... ARGS>
void async_call_list_operators(cudaStream_t stream,
                               CONTEXT& lattice,
                               const thrust::device_vector<CellID>& cells,
                               ARGS&&... args) {
  const auto block_size = 32;
  const auto block_count = (cells.size() + block_size - 1) / block_size;
  kernel::call_list_operators<<<block_count,block_size,0,stream>>>(
    lattice,
    cells.data().get(), cells.size(),
    std::forward<decltype(args)>(args)...);
  device::check();
}

/// Apply ListedCollision with statistics (async version)
template <typename CONTEXT, typename... ARGS>
void async_call_list_operators_with_statistics(cudaStream_t stream,
                                               CONTEXT& lattice,
                                               const thrust::device_vector<CellID>& cells,
                                               ARGS&&... args) {
  const auto block_size = 32;
  const auto block_count = (cells.size() + block_size - 1) / block_size;
  kernel::call_list_operators_with_statistics<<<block_count,block_size,0,stream>>>(
    lattice,
    cells.data().get(), cells.size(),
    std::forward<decltype(args)>(args)...);
  device::check();
}

}

}


/// Application of the collision step on a concrete CUDA block
template <typename T, typename DESCRIPTOR, typename DYNAMICS>
class ConcreteBlockCollisionO<T,DESCRIPTOR,Platform::GPU_CUDA,DYNAMICS> final
  : public BlockCollisionO<T,DESCRIPTOR,Platform::GPU_CUDA> {
private:
  std::unique_ptr<DYNAMICS> _dynamics;

  DynamicsParametersD<T,DESCRIPTOR,Platform::GPU_CUDA,DYNAMICS>* _parameters;
  ConcreteBlockMask<T,DESCRIPTOR,Platform::GPU_CUDA>*            _mask;

  gpu::cuda::Dynamics<T,DESCRIPTOR>** _dynamicsOfCells;

  gpu::cuda::device::unique_ptr<gpu::cuda::Dynamics<T,DESCRIPTOR>> _deviceDynamics;

  thrust::host_vector<CellID>   _hostCells;
  thrust::device_vector<CellID> _deviceCells;
  bool _modified;

  gpu::cuda::device::Stream _stream;

  /// Apply DYNAMICS using its mask and fall back to dynamic dispatch for others
  void applyDominant(ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>& block,
                     ConcreteBlockMask<T,DESCRIPTOR,Platform::GPU_CUDA>&    subdomain)
  {
    gpu::cuda::DeviceContext<T,DESCRIPTOR> lattice(block);
    if (block.statisticsEnabled()) {
      gpu::cuda::call_operators_with_statistics(
        lattice,
        subdomain.deviceData(),
        gpu::cuda::MaskedCollision<DYNAMICS>{_parameters->deviceData(), _mask->deviceData()},
        gpu::cuda::DynamicDispatchCollision{});
    } else {
      gpu::cuda::call_operators(
        lattice,
        subdomain.deviceData(),
        gpu::cuda::MaskedCollision<DYNAMICS>{_parameters->deviceData(), _mask->deviceData()},
        gpu::cuda::DynamicDispatchCollision{});
    }
  }

  /// Apply only DYNAMICS, do not apply others
  void applyIndividual(ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>& block,
                       ConcreteBlockMask<T,DESCRIPTOR,Platform::GPU_CUDA>&    subdomain)
  {
    gpu::cuda::DeviceContext<T,DESCRIPTOR> lattice(block);
    // Primitive heuristic for preferring mask-based to list-based dispatch
    if (_mask->weight() > 0.5*subdomain.weight()) {
      if (block.statisticsEnabled()) {
        gpu::cuda::async_call_operators_with_statistics(
          _stream.get(),
          lattice,
          subdomain.deviceData(),
          gpu::cuda::MaskedCollision<DYNAMICS>{_parameters->deviceData(), _mask->deviceData()});
      } else {
        gpu::cuda::async_call_operators(
          _stream.get(),
          lattice,
          subdomain.deviceData(),
          gpu::cuda::MaskedCollision<DYNAMICS>{_parameters->deviceData(), _mask->deviceData()});
      }

    // Use list of cell indices
    } else {
      // Update cell list from mask
      if (_modified) {
        _hostCells.clear();
        for (CellID iCell=0; iCell  < block.getNcells(); ++iCell) {
          if (_mask->operator[](iCell)) {
            _hostCells.push_back(iCell);
          }
        }
        _deviceCells = _hostCells;
        _modified = false;
      }

      if (block.statisticsEnabled()) {
        gpu::cuda::async_call_list_operators_with_statistics(
          _stream.get(),
          lattice,
          _deviceCells,
          gpu::cuda::ListedCollision<DYNAMICS>{_parameters->parameters});
      } else {
        gpu::cuda::async_call_list_operators(
          _stream.get(),
          lattice,
          _deviceCells,
          gpu::cuda::ListedCollision<DYNAMICS>{_parameters->parameters});
      }
    }
  }

public:
  ConcreteBlockCollisionO():
    _dynamics(new DYNAMICS()),
    _parameters(nullptr),
    _mask(nullptr),
    _modified(true),
    _stream(cudaStreamDefault)
  { }

  std::type_index id() const override
  {
    return typeid(DYNAMICS);
  }

  std::size_t weight() const override
  {
    return _mask->weight();
  }

  void set(CellID iCell, bool state, bool overlap) override
  {
    /// Only unmask cells that actually do something
    if constexpr (!std::is_same_v<DYNAMICS,NoDynamics<T,DESCRIPTOR>>) {
      if (!overlap) {
        _mask->set(iCell, state);
      }
    }
    if (state) {
      _dynamicsOfCells[iCell] = _deviceDynamics.get();
    }
    _modified = true;
  }

  Dynamics<T,DESCRIPTOR>* getDynamics() override
  {
    return _dynamics.get();
  }

  void setup(ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>& block) override
  {
    // Fetch pointers to DYNAMICS-specific parameter and mask data
    _parameters = &block.template getData<DynamicsParameters<DYNAMICS>>();
    _mask = &block.template getData<DynamicsMask<DYNAMICS>>();

    {
      // Construct on-device dynamics proxy for dynamic dispatch
      _deviceDynamics = gpu::cuda::device::malloc<gpu::cuda::ConcreteDynamics<T,DESCRIPTOR,DYNAMICS>>(1);
      gpu::cuda::kernel::construct_dynamics<T,DESCRIPTOR,DYNAMICS><<<1,1>>>(
        _deviceDynamics.get(),
        _parameters->deviceData());
      gpu::cuda::device::check();

      // Fetch pointer to on-device dynamic-dispatch field
      _dynamicsOfCells = block.template getField<gpu::cuda::DYNAMICS<T,DESCRIPTOR>>()[0].data();
    }
  }

  /// Apply collision to subdomain of block using strategy
  /**
   * The subdomain argument is currently assumed to be the core mask of BlockDynamicsMap
   **/
  void apply(ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>& block,
             ConcreteBlockMask<T,DESCRIPTOR,Platform::GPU_CUDA>&    subdomain,
             CollisionDispatchStrategy                              strategy) override
  {
    switch (strategy) {
    case CollisionDispatchStrategy::Dominant:
      return applyDominant(block, subdomain);
    case CollisionDispatchStrategy::Individual:
      return applyIndividual(block, subdomain);
    default:
      throw std::runtime_error("Invalid collision dispatch strategy");
    }
  }

};


/// Application of a cell-wise OPERATOR on a concrete CUDA block
template <typename T, typename DESCRIPTOR, typename OPERATOR>
class ConcreteBlockO<T,DESCRIPTOR,Platform::GPU_CUDA,OPERATOR,OperatorScope::PerCell> final
  : public BlockO<T,DESCRIPTOR,Platform::GPU_CUDA> {
private:
  thrust::host_vector<CellID>   _hostCells;
  thrust::device_vector<CellID> _deviceCells;
  bool _modified;

  gpu::cuda::device::Stream _stream;

public:
  ConcreteBlockO():
    _modified{false},
    _stream{cudaStreamDefault} { }

  std::type_index id() const override
  {
    return typeid(OPERATOR);
  }

  void set(CellID iCell, bool state) override
  {
    if (state) {
      _hostCells.push_back(iCell);
      _modified = true;
    }
  }

  void setup(ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>& block) override
  {
    _modified = false;
  }

  void apply(ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>& block) override
  {
    if (_hostCells.size() > 0) {
      if (_modified) {
        _deviceCells = _hostCells;
        _modified = false;
      }
      gpu::cuda::DeviceBlockLattice<T,DESCRIPTOR> lattice(block);
      gpu::cuda::async_call_list_operators(_stream.get(),
                                           lattice,
                                           _deviceCells,
                                           gpu::cuda::ListedPostProcessor<OPERATOR>{});
    }
  }

};


/// Application of a parametrized cell-wise OPERATOR on a concrete CUDA block
template <typename T, typename DESCRIPTOR, typename OPERATOR>
class ConcreteBlockO<T,DESCRIPTOR,Platform::GPU_CUDA,OPERATOR,OperatorScope::PerCellWithParameters> final
  : public BlockO<T,DESCRIPTOR,Platform::GPU_CUDA> {
private:
  thrust::host_vector<CellID>   _hostCells;
  thrust::device_vector<CellID> _deviceCells;
  bool _modified;

  gpu::cuda::device::Stream _stream;

  TrivialParametersD<T,DESCRIPTOR,Platform::GPU_CUDA,OPERATOR>* _parameters;

public:
  ConcreteBlockO():
    _modified{false},
    _stream{cudaStreamDefault} { }

  std::type_index id() const override
  {
    return typeid(OPERATOR);
  }

  void set(CellID iCell, bool state) override
  {
    if (state) {
      _hostCells.push_back(iCell);
      _modified = true;
    }
  }

  void setup(ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>& block) override
  {
    _modified = false;
    _parameters = &block.template getData<TrivialParameters<OPERATOR>>();
  }

  void apply(ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>& block) override
  {
    if (_hostCells.size() > 0) {
      if (_modified) {
        _deviceCells = _hostCells;
        _modified = false;
      }
      gpu::cuda::DeviceBlockLattice<T,DESCRIPTOR> lattice(block);
      gpu::cuda::async_call_list_operators(_stream.get(),
                                           lattice,
                                           _deviceCells,
                                           gpu::cuda::ListedPostProcessorWithParameters<OPERATOR>{_parameters->parameters});
    }
  }

};


/// Application of a block-wise OPERATOR on a concrete CUDA block
/**
 * e.g. StatisticsPostProcessor
 **/
template <typename T, typename DESCRIPTOR, typename OPERATOR>
class ConcreteBlockO<T,DESCRIPTOR,Platform::GPU_CUDA,OPERATOR,OperatorScope::PerBlock> final
  : public BlockO<T,DESCRIPTOR,Platform::GPU_CUDA> {
public:
  std::type_index id() const override
  {
    return typeid(OPERATOR);
  }

  void set(CellID iCell, bool state) override
  {
    throw std::logic_error("BlockO::set not supported for OperatorScope::PerBlock");
  }

  void setup(ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>& block) override
  {
    OPERATOR().setup(block);
  }

  void apply(ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>& block) override
  {
    OPERATOR().apply(block);
  }

};


}

#include "communicator.h"
#include "statistics.h"

#endif
