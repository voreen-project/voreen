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

#ifndef GPU_CUDA_DYNAMICS_H
#define GPU_CUDA_DYNAMICS_H

#include "dynamics/context.h"

namespace olb {

namespace gpu {

namespace cuda {

template <typename T, typename DESCRIPTOR> class DeviceContext;
template <typename T, typename DESCRIPTOR> class DataOnlyCell;

/// Virtual interface for device-side dynamically-dispatched dynamics access
template <typename T, typename DESCRIPTOR>
struct Dynamics {
  virtual CellStatistic<T> collide(DeviceContext<T,DESCRIPTOR> lattice, CellID iCell) __device__ = 0;

  virtual T    computeRho (DataOnlyCell<T,DESCRIPTOR>& cell              ) __device__ = 0;
  virtual void computeU   (DataOnlyCell<T,DESCRIPTOR>& cell,         T* u) __device__ = 0;
  virtual void computeJ   (DataOnlyCell<T,DESCRIPTOR>& cell,         T* j) __device__ = 0;
  virtual void computeRhoU(DataOnlyCell<T,DESCRIPTOR>& cell, T& rho, T* u) __device__ = 0;

  virtual void computeStress    (DataOnlyCell<T,DESCRIPTOR>& cell, T& rho, T* u, T* pi) __device__ = 0;
  virtual void computeAllMomenta(DataOnlyCell<T,DESCRIPTOR>& cell, T& rho, T* u, T* pi) __device__ = 0;

  virtual T computeEquilibrium(int iPop, T rho, T* u) __device__ = 0;

  virtual T getOmegaOrFallback(T fallback) __device__ = 0;

  void iniEquilibrium(DataOnlyCell<T,DESCRIPTOR>& cell, T rho, T* u) __device__ {
    for (unsigned iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      cell[iPop] = computeEquilibrium(iPop, rho, u);
    }
  };

};

/// On-device field mirroring BlockDynamicsMap
template <typename T, typename DESCRIPTOR>
struct DYNAMICS : public descriptors::TYPED_FIELD_BASE<Dynamics<T,DESCRIPTOR>*,1> { };

/// Implementation of gpu::cuda::Dynamics for concrete DYNAMICS
template <typename T, typename DESCRIPTOR, typename DYNAMICS>
class ConcreteDynamics final : public Dynamics<T,DESCRIPTOR> {
private:
  typename DYNAMICS::ParametersD* _parameters;

public:
  ConcreteDynamics(typename DYNAMICS::ParametersD* parameters) __device__:
    _parameters{parameters} {
  }

  CellStatistic<T> collide(DeviceContext<T,DESCRIPTOR> lattice, CellID iCell) override __device__ {
    DataOnlyCell<T,DESCRIPTOR> cell(lattice, iCell);
    return DYNAMICS().apply(cell, *_parameters);
  }

  T computeRho(DataOnlyCell<T,DESCRIPTOR>& cell) override __device__ {
    return DYNAMICS::MomentaF().computeRho(cell);
  }
  void computeU(DataOnlyCell<T,DESCRIPTOR>& cell, T* u) override __device__ {
    DYNAMICS::MomentaF().computeU(cell, u);
  }
  void computeJ(DataOnlyCell<T,DESCRIPTOR>& cell, T* j) override __device__ {
    DYNAMICS::MomentaF().computeJ(cell, j);
  }
  void computeRhoU(DataOnlyCell<T,DESCRIPTOR>& cell, T& rho, T* u) override __device__ {
    DYNAMICS::MomentaF().computeRhoU(cell, rho, u);
  }
  void computeStress(DataOnlyCell<T,DESCRIPTOR>& cell, T& rho, T* u, T* pi) override __device__ {
    DYNAMICS::MomentaF().computeStress(cell, rho, u, pi);
  }
  void computeAllMomenta(DataOnlyCell<T,DESCRIPTOR>& cell, T& rho, T* u, T* pi) override __device__ {
    DYNAMICS::MomentaF().computeAllMomenta(cell, rho, u, pi);
  }

  T getOmegaOrFallback(T fallback) override __device__ {
    if constexpr (DYNAMICS::ParametersD::fields_t::template contains<descriptors::OMEGA>()) {
      return _parameters->template get<descriptors::OMEGA>();
    } else {
      return fallback;
    }
  }

  T computeEquilibrium(int iPop, T rho, T* u) override __device__ {
    return DYNAMICS().computeEquilibrium(iPop, rho, u);
  }

};

/// Last node in a MaskedDynamics chain in kernel::call_operators
struct DynamicDispatchCollision {
  template <typename T, typename DESCRIPTOR>
  bool operator()(DeviceContext<T,DESCRIPTOR>& lattice, CellID iCell) __device__ {
    if (auto* collisionO = lattice.template getField<DYNAMICS<T,DESCRIPTOR>>()[0][iCell]) {
      collisionO->collide(lattice, iCell);
      return true;
    }
    return false;
  }

  template <typename T, typename DESCRIPTOR>
  bool operator()(DeviceContext<T,DESCRIPTOR>& lattice, CellID iCell, CellStatistic<T>& statistic) __device__ {
    if (auto* collisionO = lattice.template getField<DYNAMICS<T,DESCRIPTOR>>()[0][iCell]) {
      statistic = collisionO->collide(lattice, iCell);
      return true;
    }
    return false;
  }

};

namespace kernel {

/// CUDA kernel for constructing on-device ConcreteDynamics
template <typename T, typename DESCRIPTOR, typename DYNAMICS, typename PARAMETERS=typename DYNAMICS::ParametersD>
void construct_dynamics(void* target, PARAMETERS* parameters) __global__ {
  new (target) ConcreteDynamics<T,DESCRIPTOR,DYNAMICS>(parameters);
}

}

}

}


/// Representation of DynamicsParameters<DYNAMICS> for CUDA block lattice
template <typename T, typename DESCRIPTOR, typename DYNAMICS>
class DynamicsParametersD<T,DESCRIPTOR,Platform::GPU_CUDA,DYNAMICS> final
  : public AbstractDynamicsParameters<T,DESCRIPTOR>
  , public Serializable {
private:
  gpu::cuda::device::unique_ptr<typename DYNAMICS::ParametersD> _deviceParameters;

public:
  typename DYNAMICS::ParametersD parameters;

  DynamicsParametersD(std::size_t): // TODO: Implement more generic non-cellwise field allocation in Data
    _deviceParameters{gpu::cuda::device::malloc<typename DYNAMICS::ParametersD>(1)},
    parameters{}
  {
    gpu::cuda::device::copyToDevice(&parameters,
                                    _deviceParameters.get(),
                                    sizeof(typename DYNAMICS::ParametersD));
  }

  AbstractParameters<T,DESCRIPTOR>& asAbstract() override {
    return parameters;
  }

  void setProcessingContext(ProcessingContext context) override {
    if (context == ProcessingContext::Simulation) {
      gpu::cuda::device::copyToDevice(&parameters,
                                      _deviceParameters.get(),
                                      sizeof(typename DYNAMICS::ParametersD));
    }
  }

  typename DYNAMICS::ParametersD* deviceData() {
    return _deviceParameters.get();
  }

  /// Number of data blocks for the serializable interface
  std::size_t getNblock() const override;
  /// Binary size for the serializer
  std::size_t getSerializableSize() const override;
  /// Return a pointer to the memory of the current block and its size for the serializable interface
  bool* getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode) override;

};

template<typename T, typename DESCRIPTOR, typename DYNAMICS>
std::size_t DynamicsParametersD<T,DESCRIPTOR,Platform::GPU_CUDA,DYNAMICS>::getNblock() const
{
  return decltype(parameters)::fields_t::size;
}

template<typename T, typename DESCRIPTOR, typename DYNAMICS>
std::size_t DynamicsParametersD<T,DESCRIPTOR,Platform::GPU_CUDA,DYNAMICS>::getSerializableSize() const
{
  std::size_t size = 0;
  decltype(parameters)::fields_t::for_each([&size](auto field) {
    using field_t = typename decltype(field)::type;
    size += FieldD<T,DESCRIPTOR,field_t>{}.getSerializableSize();
  });
  return size;
}

template<typename T, typename DESCRIPTOR, typename DYNAMICS>
bool* DynamicsParametersD<T,DESCRIPTOR,Platform::GPU_CUDA,DYNAMICS>::getBlock(
  std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode)
{
  std::size_t currentBlock = 0;
  bool* dataPtr = nullptr;
  decltype(parameters)::fields_t::for_each([&](auto field) {
    using field_t = typename decltype(field)::type;
    if constexpr (DESCRIPTOR::template size<field_t>() == 1) {
      registerVar(iBlock, sizeBlock, currentBlock, dataPtr,
                  parameters.template get<field_t>(), loadingMode);
    } else {
      registerSerializableOfConstSize(iBlock, sizeBlock, currentBlock, dataPtr,
                                      parameters.template get<field_t>(), loadingMode);
    }
  });
  return dataPtr;
}

}

#endif
