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

#ifndef GPU_CUDA_CONTEXT_H
#define GPU_CUDA_CONTEXT_H

#include "registry.h"
#include "dynamics.h"

namespace olb {

template <typename T, typename DESCRIPTOR, Platform PLATFORM> class ConcreteBlockLattice;

namespace gpu {

namespace cuda {

/// Return pointers to on-device data for FIELD on lattice
template <typename T, typename DESCRIPTOR, typename FIELD>
static auto getDeviceFieldPointer(ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>& lattice) {
  auto& fieldArray = lattice.template getField<FIELD>();
  return Vector<typename FIELD::template value_type<T>*,DESCRIPTOR::template size<FIELD>()>([&](unsigned iD) {
    return fieldArray[iD].deviceData();
  });
}

/// Return tuple of pointers to on-device data for multiple FIELDS on lattice
template <typename T, typename DESCRIPTOR, typename... FIELDS>
static auto getDeviceFieldPointers(typename meta::list<FIELDS...>,
                                   ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>& lattice) {
  return std::make_tuple(getDeviceFieldPointer<T,DESCRIPTOR,FIELDS>(lattice)...);
}

/// Structure for passing pointers to on-device data into CUDA kernels
template <typename T, typename DESCRIPTOR>
class DeviceContext {
private:
  const std::size_t _nCells;

  /// Pointers to descriptor-declared field data
  decltype(getDeviceFieldPointers(typename DESCRIPTOR::fields_t(),
                                  std::declval<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>&>()))
  _staticFieldsD;
  /// Pointer to type-erased field data
  /**
   * Maintained by FieldTypeRegistry<T,DESCRIPTOR,Platform::GPU_CUDA>
   **/
  void*** _customFieldsD;

public:
  using value_t = T;
  using descriptor_t = DESCRIPTOR;

  DeviceContext(ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>& lattice) __host__:
    _nCells(lattice.getNcells()),
    _staticFieldsD(getDeviceFieldPointers(typename DESCRIPTOR::fields_t(), lattice)),
    _customFieldsD(lattice.getDataRegistry().deviceData())
  { }

  std::size_t getNcells() const any_platform {
    return _nCells;
  }

  template <typename FIELD>
  auto getField() __device__ {
    if constexpr (DESCRIPTOR::template provides<FIELD>()) {
      return std::get<(DESCRIPTOR::fields_t::template index<FIELD>())>(_staticFieldsD).data();
    } else {
      return reinterpret_cast<typename FIELD::template value_type<T>**>(
        _customFieldsD[field_type_index<void,Array<FIELD>>]);
    }
  }

};

/// Pointer to row of a D-dimensional field
template <typename T, typename DESCRIPTOR, typename FIELD>
class FieldPtr : public ScalarVector<typename FIELD::template value_type<T>,
                                     DESCRIPTOR::template size<FIELD>(),
                                     FieldPtr<T,DESCRIPTOR,FIELD>> {
private:
  DeviceContext<T,DESCRIPTOR>& _data;
  std::size_t _index;

  friend typename ScalarVector<typename FIELD::template value_type<T>,
                               DESCRIPTOR::template size<FIELD>(),
                               FieldPtr>::type;

protected:
  const typename FIELD::template value_type<T>* getComponentPointer(unsigned iDim) const __device__
  {
    return _data.template getField<FIELD>()[iDim] + _index;
  }
  typename FIELD::template value_type<T>* getComponentPointer(unsigned iDim) __device__
  {
    return _data.template getField<FIELD>()[iDim] + _index;
  }

public:
  FieldPtr(DeviceContext<T,DESCRIPTOR>& data, std::size_t index) __device__:
    _data(data),
    _index(index) { }

  FieldPtr(FieldPtr<T,DESCRIPTOR,FIELD>&& rhs) __device__:
    _data(rhs._data),
    _index(rhs._index) { }

  template <typename U, typename IMPL>
  FieldPtr<T,DESCRIPTOR,FIELD>& operator=(const GenericVector<U,DESCRIPTOR::template size<FIELD>(),IMPL>& rhs) __device__
  {
    for (unsigned iDim=0; iDim < DESCRIPTOR::template size<FIELD>(); ++iDim) {
      this->operator[](iDim) = rhs[iDim];
    }
    return *this;
  }

};

/// Device-side implementation of the data-only Cell concept for collision steps
template <typename T, typename DESCRIPTOR>
class DataOnlyCell {
protected:
  DeviceContext<T,DESCRIPTOR>& _data;
  CellID _iCell;

public:
  using value_t = T;
  using descriptor_t = DESCRIPTOR;

  DataOnlyCell(DeviceContext<T,DESCRIPTOR>& data, CellID iCell) __device__:
    _data(data),
    _iCell(iCell)
  { }

  template <typename FIELD>
  typename FIELD::template value_type<T> getFieldComponent(unsigned iD) __device__ {
    return _data.template getField<FIELD>()[iD][_iCell];
  }

  value_t& operator[](int iPop) __device__ {
    return _data.template getField<descriptors::POPULATION>()[iPop][_iCell];
  }

  template <typename FIELD>
  auto getField() const __device__ {
    auto fieldArray = _data.template getField<FIELD>();
    if constexpr (descriptor_t::template size<FIELD>() == 1) {
      return fieldArray[0][_iCell];
    } else {
      return FieldD<value_t,descriptor_t,FIELD>([&](unsigned iD) {
        return fieldArray[iD][_iCell];
      });
    }
  }

  template <typename FIELD>
  FieldPtr<T,DESCRIPTOR,FIELD> getFieldPointer() __device__ {
    return FieldPtr<T,DESCRIPTOR,FIELD>(_data, _iCell);
  }

  template <typename FIELD>
  void setField(FieldD<value_t,descriptor_t,FIELD>&& value) __device__ {
    auto fieldArray = _data.template getField<FIELD>();
    for (unsigned iD=0; iD < descriptor_t::template size<FIELD>(); ++iD) {
      fieldArray[iD][_iCell] = value[iD];
    }
  }

  template <typename FIELD>
  void setField(const FieldD<value_t,descriptor_t,FIELD>& value) __device__ {
    auto fieldArray = _data.template getField<FIELD>();
    for (unsigned iD=0; iD < descriptor_t::template size<FIELD>(); ++iD) {
      fieldArray[iD][_iCell] = value[iD];
    }
  }

};


template <typename T, typename DESCRIPTOR>
class Cell;

/// Device-side view of a block lattice
/**
 * Used for non-local operators such as post processors
 **/
template <typename T, typename DESCRIPTOR>
class DeviceBlockLattice final : public DeviceContext<T,DESCRIPTOR> {
private:
  LatticeR<DESCRIPTOR::d> _projection;
  int _padding;

public:
  DeviceBlockLattice(ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>& lattice) __host__:
    DeviceContext<T,DESCRIPTOR>(lattice),
    _padding{lattice.getPadding()}
  {
    auto size = lattice.getExtent() + 2*_padding;
    if constexpr (DESCRIPTOR::d == 3) {
      _projection = {size[1]*size[2], size[2], 1};
    } else {
      _projection = {size[1], 1};
    }
  }

  CellID getCellId(LatticeR<DESCRIPTOR::d> loc) const __device__ {
    return (loc+_padding) * _projection;
  }

  CellDistance getNeighborDistance(LatticeR<DESCRIPTOR::d> dir) const __device__ {
    return dir * _projection;
  }

  Cell<T,DESCRIPTOR> get(CellID iCell) __device__ {
    return Cell<T,DESCRIPTOR>(*this, iCell);
  }

  Cell<T,DESCRIPTOR> get(LatticeR<DESCRIPTOR::d> loc) __device__ {
    return get(getCellId(loc));
  }

};

/// Device-side implementation of the Cell concept for post processors
/**
 * Adds neighborhood and dynamically-dispatched momenta access to DataOnlyCell
 **/
template <typename T, typename DESCRIPTOR>
class Cell final : public DataOnlyCell<T,DESCRIPTOR> {
protected:
  DeviceBlockLattice<T,DESCRIPTOR>& getBlock() __device__ {
    return static_cast<DeviceBlockLattice<T,DESCRIPTOR>&>(this->_data);
  }

public:
  Cell(DeviceBlockLattice<T,DESCRIPTOR>& data, CellID iCell) __device__:
    DataOnlyCell<T,DESCRIPTOR>(data, iCell)
  { }

  Cell<T,DESCRIPTOR> neighbor(LatticeR<DESCRIPTOR::d> offset) __device__ {
    return {getBlock(), CellID(this->_iCell + getBlock().getNeighborDistance(offset))};
  }

  Dynamics<T,DESCRIPTOR>& getDynamics() __device__ {
    return *this->template getField<DYNAMICS<T,DESCRIPTOR>>();
  }

  T computeRho() __device__ {
    return getDynamics().computeRho(*this);
  }
  void computeU(T* u) __device__ {
    getDynamics().computeU(*this, u);
  }
  void computeJ(T* j) __device__ {
    getDynamics().computeJ(*this, j);
  }
  void computeRhoU(T& rho, T* u) __device__ {
    getDynamics().computeRhoU(*this, rho, u);
  }
  void computeStress(T* pi) __device__ {
    T rho, u[DESCRIPTOR::d] { };
    getDynamics().computeRhoU(*this, rho, u);
    getDynamics().computeStress(*this, rho, u, pi);
  }
  void computeAllMomenta(T& rho, T* u, T* pi) __device__ {
    getDynamics().computeAllMomenta(*this, rho, u, pi);
  }

  void iniEquilibrium(T rho, T* u) __device__ {
    getDynamics().iniEquilibrium(*this, rho, u);
  }

};

}

}

template <typename T, typename DESCRIPTOR, typename CONTEXT>
class TrivialParametersD<T,DESCRIPTOR,Platform::GPU_CUDA,CONTEXT> final : public Serializable {
private:
  gpu::cuda::device::unique_ptr<typename CONTEXT::ParametersD> _deviceParameters;

public:
  typename CONTEXT::ParametersD parameters;

  TrivialParametersD(std::size_t): // TODO: Implement more generic non-cellwise field allocation in Data
    _deviceParameters{gpu::cuda::device::malloc<typename CONTEXT::ParametersD>(1)},
    parameters{}
  { }

  /// Return abstract interface to host-side parameters
  typename CONTEXT::ParametersD& asAbstract() {
    return parameters;
  }

  void setProcessingContext(ProcessingContext context) {
    if (context == ProcessingContext::Simulation) {
      gpu::cuda::device::copyToDevice(&parameters,
                                      _deviceParameters.get(),
                                      sizeof(typename CONTEXT::ParametersD));
    }
  }

  /// Number of data blocks for the serializable interface
  std::size_t getNblock() const override;
  /// Binary size for the serializer
  std::size_t getSerializableSize() const override;
  /// Return a pointer to the memory of the current block and its size for the serializable interface
  bool* getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode) override;

};

template <typename T, typename DESCRIPTOR, typename CONTEXT>
std::size_t TrivialParametersD<T,DESCRIPTOR,Platform::GPU_CUDA,CONTEXT>::getNblock() const
{
  return 1;
}

template <typename T, typename DESCRIPTOR, typename CONTEXT>
std::size_t TrivialParametersD<T,DESCRIPTOR,Platform::GPU_CUDA,CONTEXT>::getSerializableSize() const
{
  return sizeof(typename CONTEXT::ParametersD);
}

template <typename T, typename DESCRIPTOR, typename CONTEXT>
bool* TrivialParametersD<T,DESCRIPTOR,Platform::GPU_CUDA,CONTEXT>::getBlock(
  std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode)
{
  std::size_t currentBlock = 0;
  bool* dataPtr = nullptr;
  registerVar(iBlock, sizeBlock, currentBlock, dataPtr,
              parameters, loadingMode);
  return dataPtr;
}

}

#endif
