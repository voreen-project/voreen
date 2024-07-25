/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Adrian Kummerlaender
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

#ifndef SIMD_OPERATOR_H
#define SIMD_OPERATOR_H

#include "mask.h"

#include "core/meta.h"
#include "core/operator.h"

#include "core/platform/cpu/cell.h"

namespace olb {

template <typename T>
struct CellStatistic<cpu::simd::Pack<T>> {
  cpu::simd::Pack<T> rho;
  cpu::simd::Pack<T> uSqr;

  operator bool() const {
    return true;
  }
};

namespace cpu {

/// Implementations of vector CPU specifics
namespace simd {

/// SIMD-specific pointer to a pack of rows of a D-dimensional field
template <typename T, unsigned D>
class FieldPtr : public ScalarVector<Pack<T>,D,FieldPtr<T,D>> {
private:
  ColumnVector<Column<T>,D>& _data;
  std::size_t _index;

  friend typename ScalarVector<Pack<T>,D,FieldPtr>::type;

  std::array<Pack<T>,D> _packs;

protected:
  const Pack<T>* getComponentPointer(unsigned iDim) const
  {
    return &_packs[iDim];
  }
  Pack<T>* getComponentPointer(unsigned iDim)
  {
    return &_packs[iDim];
  }

public:
  FieldPtr(ColumnVector<Column<T>,D>& columns, std::size_t index):
    _data(columns),
    _index(index) {
    meta::call_n_times<D>([&](unsigned iD) {
      _packs[iD] = &_data[iD][_index];
    });
  }

  FieldPtr(FieldPtr<T,D>&& rhs):
    _data(rhs._data),
    _index(rhs._index),
    _packs(rhs._packs) { }

  template <typename U, typename IMPL>
  FieldPtr<T,D>& operator=(const GenericVector<U,D,IMPL>& rhs)
  {
    for (unsigned iDim=0; iDim < D; ++iDim) {
      this->operator[](iDim) = rhs[iDim];
    }
    return *this;
  }

};

/// Implementation of the Cell concept for vectorized collision operators
/**
 * Any fields that are potentially changed by a specific collision operator
 * must be declared as RW_FIELDS to mask any writeback operations.
 *
 * Read-only fields require no special consideration.
 **/
template <typename T, typename DESCRIPTOR, typename V, typename... RW_FIELDS>
class Cell {
private:
  using rw_fields = meta::list<RW_FIELDS...>;

  ConcreteBlockLattice<T,DESCRIPTOR,Platform::CPU_SIMD>& _lattice;
  std::size_t _iCell;
  Mask<T>& _mask;

  /// Storage for packs of fields that are both read and written
  /**
   * Commonly this is only descriptors::POPULATION
   **/
  std::tuple<FieldD<V,DESCRIPTOR,RW_FIELDS>...> _fields;

public:
  using value_t = V;
  using descriptor_t = DESCRIPTOR;

  /// Load r/w fields into SIMD packs
  Cell(ConcreteBlockLattice<T,DESCRIPTOR,Platform::CPU_SIMD>& lattice, std::size_t iCell, Mask<T>& mask):
    _lattice(lattice),
    _iCell(iCell),
    _mask(mask)
  {
    rw_fields::for_each([&](auto field) {
      using FIELD = typename decltype(field)::type;
      auto& array = _lattice.template getField<FIELD>();
      auto& pack = std::get<(rw_fields::template index<FIELD>())>(_fields);
      //meta::call_n_times<(DESCRIPTOR::template size<FIELD>())>([&](unsigned iD) {
      for (unsigned iD=0; iD < DESCRIPTOR::template size<FIELD>(); ++iD) {
        pack[iD] = cpu::simd::Pack<T>(&array[iD][_iCell]);
      }
    });
  }

  /// Store modified r/w fields back into lattice taking into account the mask
  ~Cell()
  {
    rw_fields::for_each([&](auto field) {
      using FIELD = typename decltype(field)::type;
      auto& array = _lattice.template getField<FIELD>();
      auto& pack = std::get<(rw_fields::template index<FIELD>())>(_fields);
      //meta::call_n_times<(DESCRIPTOR::template size<FIELD>())>([&](unsigned iD) {
      for (unsigned iD=0; iD < DESCRIPTOR::template size<FIELD>(); ++iD) {
        cpu::simd::maskstore(&array[iD][_iCell], _mask, pack[iD]);
      }
    });
  }

  /// Return reference to iPop population pack
  V& operator[](unsigned iPop) {
    return std::get<(rw_fields::template index<descriptors::POPULATION>())>(_fields)[iPop];
  }

  /// Return pack-valued copy of FIELD
  template <typename FIELD>
  auto getField() const {
    if constexpr (rw_fields::template contains<FIELD>()) {
      if constexpr (DESCRIPTOR::template size<FIELD>() == 1) {
        return std::get<(rw_fields::template index<FIELD>())>(_fields)[0];
      } else {
        return std::get<(rw_fields::template index<FIELD>())>(_fields);
      }
    } else {
      auto& fieldArray = _lattice.template getField<FIELD>();
      if constexpr (DESCRIPTOR::template size<FIELD>() == 1) {
        return cpu::simd::Pack<T>(&fieldArray[0][_iCell]);
      } else {
        return FieldD<V,DESCRIPTOR,FIELD>([&](unsigned iD) {
          return &fieldArray[iD][_iCell];
        });
      }
    }
    __builtin_unreachable();
  }

  /// Set compoents of FIELD from pack-valued vector
  template <typename FIELD>
  void setField(FieldD<V,DESCRIPTOR,FIELD>&& value) {
    if constexpr (rw_fields::template contains<FIELD>()) {
      std::get<(rw_fields::template index<FIELD>())>(_fields) = value;
    } else {
      auto& array = _lattice.template getField<FIELD>();
      for (unsigned iD=0; iD < DESCRIPTOR::template size<FIELD>(); ++iD) {
        cpu::simd::maskstore(&array[iD][_iCell], _mask, value[iD]);
      }
    }
  }

  /// Return reference to pack-valued interim storage vector of r/w field
  template <typename FIELD>
  std::enable_if_t<rw_fields::template contains<FIELD>(), FieldD<V,DESCRIPTOR,FIELD>&>
  getFieldPointer() {
    return std::get<(rw_fields::template index<FIELD>())>(_fields);
  }

  /// Return pack-valued copy of non r/w field
  template <typename FIELD>
  std::enable_if_t<!rw_fields::template contains<FIELD>(), FieldD<V,DESCRIPTOR,FIELD>>
  getFieldPointer() {
    return getField<FIELD>();
  }

  /// Return reference to pack-valued interim storage component of r/w field
  template <typename FIELD>
  std::enable_if_t<rw_fields::template contains<FIELD>(),V&>
  getFieldComponent(unsigned iD) {
    return std::get<(rw_fields::template index<FIELD>())>(_fields)[iD];
  }

  /// Return pack-valued copy of non r/w field component
  template <typename FIELD>
  std::enable_if_t<!rw_fields::template contains<FIELD>(),V>
  getFieldComponent(unsigned iD) {
    return &_lattice.template getField<FIELD>()[iD][_iCell];
  }

};


/// Implementation of cpu::Dynamics for concrete DYNAMICS on SIMD blocks
template <typename T, typename DESCRIPTOR, typename DYNAMICS>
class ConcreteDynamics final : public cpu::Dynamics<T,DESCRIPTOR,Platform::CPU_SIMD> {
private:
  typename DYNAMICS::ParametersD* _parameters;

public:
  ConcreteDynamics(typename DYNAMICS::ParametersD* parameters):
    _parameters{parameters} {
  }

  CellStatistic<T> collide(cpu::Cell<T,DESCRIPTOR,Platform::CPU_SIMD>& cell) override {
    return DYNAMICS().apply(cell, *_parameters);
  }

  T computeRho(cpu::Cell<T,DESCRIPTOR,Platform::CPU_SIMD>& cell) override {
    return typename DYNAMICS::MomentaF().computeRho(cell);
  }
  void computeU(cpu::Cell<T,DESCRIPTOR,Platform::CPU_SIMD>& cell, T* u) override {
    typename DYNAMICS::MomentaF().computeU(cell, u);
  }
  void computeJ(cpu::Cell<T,DESCRIPTOR,Platform::CPU_SIMD>& cell, T* j) override {
    typename DYNAMICS::MomentaF().computeJ(cell, j);
  }
  void computeRhoU(cpu::Cell<T,DESCRIPTOR,Platform::CPU_SIMD>& cell, T& rho, T* u) override {
    typename DYNAMICS::MomentaF().computeRhoU(cell, rho, u);
  }
  void computeStress(cpu::Cell<T,DESCRIPTOR,Platform::CPU_SIMD>& cell, T& rho, T* u, T* pi) override {
    typename DYNAMICS::MomentaF().computeStress(cell, rho, u, pi);
  }
  void computeAllMomenta(cpu::Cell<T,DESCRIPTOR,Platform::CPU_SIMD>& cell, T& rho, T* u, T* pi) override {
    typename DYNAMICS::MomentaF().computeAllMomenta(cell, rho, u, pi);
  }

  T getOmegaOrFallback(T fallback) override {
    if constexpr (DYNAMICS::ParametersD::fields_t::template contains<descriptors::OMEGA>()) {
      return _parameters->template get<descriptors::OMEGA>();
    } else {
      return fallback;
    }
    __builtin_unreachable();
  }

  T computeEquilibrium(int iPop, T rho, T* u) override {
    return DYNAMICS().computeEquilibrium(iPop, rho, u);
  }

};

}

}


/// Application of the collision step on a concrete SIMD block
/**
 * ConcreteBlockCollisionO::apply allows for applying DYNAMICS
 * on a given cell index or an entire block. The latter option
 * accepts a ConcreteBlockMask instance describing the subset
 * of cells for which DYNAMICS is to be vectorized.
 **/
template <typename T, typename DESCRIPTOR, typename DYNAMICS>
class ConcreteBlockCollisionO<T,DESCRIPTOR,Platform::CPU_SIMD,DYNAMICS> final
  : public BlockCollisionO<T,DESCRIPTOR,Platform::CPU_SIMD> {
private:
  std::unique_ptr<DYNAMICS> _dynamics;
  std::unique_ptr<cpu::Dynamics<T,DESCRIPTOR,Platform::CPU_SIMD>> _concreteDynamics;

  DynamicsParametersD<T,DESCRIPTOR,Platform::CPU_SIMD,DYNAMICS>* _parameters;
  ConcreteBlockMask<T,DESCRIPTOR,Platform::CPU_SIMD>*            _mask;

  cpu::Dynamics<T,DESCRIPTOR,Platform::CPU_SIMD>** _dynamicsOfCells;

  /// Helper for dispatching collision on non-DYNAMICS cells
  void applyOther(ConcreteBlockLattice<T,DESCRIPTOR,Platform::CPU_SIMD>& block,
                  std::size_t                                            iCell)
  {
    cpu::Cell<T,DESCRIPTOR,Platform::CPU_SIMD> cell(block, iCell);
    if (auto cellStatistic = _dynamicsOfCells[iCell]->collide(cell)) {
      block.getStatistics().incrementStats(cellStatistic.rho, cellStatistic.uSqr);
    }
  }

  /// Apply collision on cell range [iCell,iCell+pack_size) of block
  void apply(ConcreteBlockLattice<T,DESCRIPTOR,Platform::CPU_SIMD>& block,
             ConcreteBlockMask<T,DESCRIPTOR,Platform::CPU_SIMD>&    subdomain,
             ConcreteBlockMask<T,DESCRIPTOR,Platform::CPU_SIMD>&    mask,
             typename DYNAMICS::ParametersD&                        parameters,
             std::size_t                                            iCell)
  {
    if constexpr (dynamics::is_vectorizable_v<DYNAMICS>) {
      if (cpu::simd::Mask<T> m = {mask.raw(), iCell}) {
        cpu::simd::Cell<T,DESCRIPTOR,cpu::simd::Pack<T>,descriptors::POPULATION> cell(block, iCell, m);
        auto cellStatistic = DYNAMICS().apply(cell, parameters);
        auto& statistics = block.getStatistics();
        for (unsigned i=0; i < cpu::simd::Pack<T>::size; ++i) {
          if (mask[iCell+i]) {
            if (cellStatistic.rho[i] != T{-1}) {
              statistics.incrementStats(cellStatistic.rho[i], cellStatistic.uSqr[i]);
            }
          } else if (subdomain[iCell+i]) {
            applyOther(block, iCell+i);
          }
        }
      } else {
        for (std::size_t i=iCell; i < iCell+cpu::simd::Pack<T>::size; ++i) {
          if (subdomain[i]) {
            applyOther(block, i);
          }
        }
      }
    }
  }

public:
  ConcreteBlockCollisionO():
    _dynamics(new DYNAMICS()),
    _parameters(nullptr),
    _mask(nullptr)
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
      _dynamicsOfCells[iCell] = _concreteDynamics.get();
    }
  }

  Dynamics<T,DESCRIPTOR>* getDynamics() override
  {
    return _dynamics.get();
  }

  void setup(ConcreteBlockLattice<T,DESCRIPTOR,Platform::CPU_SIMD>& block) override
  {
    _parameters = &block.template getData<DynamicsParameters<DYNAMICS>>();
    _mask = &block.template getData<DynamicsMask<DYNAMICS>>();
    if constexpr (dynamics::has_parametrized_momenta_v<DYNAMICS>) {
      _dynamics->setMomentaParameters(&_parameters->parameters);
    }

    _concreteDynamics.reset(new cpu::simd::ConcreteDynamics<T,DESCRIPTOR,DYNAMICS>(&_parameters->parameters));
    // Fetch pointer to concretized dynamic-dispatch field
    _dynamicsOfCells = block.template getField<cpu::DYNAMICS<T,DESCRIPTOR,Platform::CPU_SIMD>>()[0].data();
  }

  /// Apply collision on entire block
  void apply(ConcreteBlockLattice<T,DESCRIPTOR,Platform::CPU_SIMD>& block,
             ConcreteBlockMask<T,DESCRIPTOR,Platform::CPU_SIMD>&    subdomain,
             CollisionDispatchStrategy                              strategy) override
  {
    if (strategy != CollisionDispatchStrategy::Dominant) {
      throw std::runtime_error("Platform::CPU_SIMD currently only support CollisionDispatchStrategy::Dominant");
    }

    auto& mask = *_mask;
    if constexpr (dynamics::is_vectorizable_v<DYNAMICS>) {
      // Ensure that serialized mask storage is up-to-date
      mask.setProcessingContext(ProcessingContext::Simulation);
      #ifdef PARALLEL_MODE_OMP
      #pragma omp parallel for schedule(static)
      #endif
      for (std::size_t iCell=0; iCell < block.getNcells(); iCell += cpu::simd::Pack<T>::size) {
        apply(block, subdomain, mask, _parameters->parameters, iCell);
      }
    } else { // Fallback for non-vectorizable collision operators
      #ifdef PARALLEL_MODE_OMP
      #pragma omp parallel for schedule(static)
      #endif
      for (std::size_t iCell=0; iCell < block.getNcells(); ++iCell) {
        if (mask[iCell]) {
          cpu::Cell<T,DESCRIPTOR,Platform::CPU_SIMD> cell(block, iCell);
          if (auto cellStatistic = DYNAMICS().apply(cell, _parameters->parameters)) {
            block.getStatistics().incrementStats(cellStatistic.rho, cellStatistic.uSqr);
          }
        } else if (subdomain[iCell]) {
          applyOther(block, iCell);
        }
      }
    }
  }

};


/// Application of a cell-wise OPERATOR on a concrete vector CPU block
template <typename T, typename DESCRIPTOR, typename OPERATOR>
class ConcreteBlockO<T,DESCRIPTOR,Platform::CPU_SIMD,OPERATOR,OperatorScope::PerCell> final
  : public BlockO<T,DESCRIPTOR,Platform::CPU_SIMD> {
private:
  std::vector<CellID> _cells;

public:
  std::type_index id() const override
  {
    return typeid(OPERATOR);
  }

  void set(CellID iCell, bool state) override
  {
    if (state) {
      _cells.emplace_back(iCell);
    }
  }

  void setup(ConcreteBlockLattice<T,DESCRIPTOR,Platform::CPU_SIMD>& block) override
  { }

  void apply(ConcreteBlockLattice<T,DESCRIPTOR,Platform::CPU_SIMD>& block) override
  {
    if (_cells.size() > 0) {
      cpu::Cell<T,DESCRIPTOR,Platform::CPU_SIMD> cell(block, 0);
      for (CellID iCell : _cells) {
        cell.setCellId(iCell);
        OPERATOR().apply(cell);
      }
    }
  }

};


template <typename T, typename DESCRIPTOR, typename OPERATOR>
class ConcreteBlockO<T,DESCRIPTOR,Platform::CPU_SIMD,OPERATOR,OperatorScope::PerCellWithParameters> final
  : public BlockO<T,DESCRIPTOR,Platform::CPU_SIMD> {
private:
  std::vector<CellID> _cells;

  TrivialParametersD<T,DESCRIPTOR,Platform::CPU_SIMD,OPERATOR>* _parameters;

public:
  std::type_index id() const override
  {
    return typeid(OPERATOR);
  }

  void set(CellID iCell, bool state) override
  {
    if (state) {
      _cells.emplace_back(iCell);
    }
  }

  void setup(ConcreteBlockLattice<T,DESCRIPTOR,Platform::CPU_SIMD>& block) override
  {
    _parameters = &block.template getData<TrivialParameters<OPERATOR>>();
  }

  void apply(ConcreteBlockLattice<T,DESCRIPTOR,Platform::CPU_SIMD>& block) override
  {
    if (_cells.size() > 0) {
      cpu::Cell<T,DESCRIPTOR,Platform::CPU_SIMD> cell(block, 0);
      for (CellID iCell : _cells) {
        cell.setCellId(iCell);
        OPERATOR().apply(cell, _parameters->parameters);
      }
    }
  }

};


/// Application of a block-wise OPERATOR on a concrete vector CPU block
/**
 * e.g. StatisticsPostProcessor
 **/
template <typename T, typename DESCRIPTOR, typename OPERATOR>
class ConcreteBlockO<T,DESCRIPTOR,Platform::CPU_SIMD,OPERATOR,OperatorScope::PerBlock> final
  : public BlockO<T,DESCRIPTOR,Platform::CPU_SIMD> {
public:
  std::type_index id() const override
  {
    return typeid(OPERATOR);
  }

  void set(CellID iCell, bool state) override
  {
    throw std::logic_error("BlockO::set not supported for OperatorScope::PerBlock");
  }

  void setup(ConcreteBlockLattice<T,DESCRIPTOR,Platform::CPU_SIMD>& block) override
  {
    OPERATOR().setup(block);
  }

  void apply(ConcreteBlockLattice<T,DESCRIPTOR,Platform::CPU_SIMD>& block) override
  {
    OPERATOR().apply(block);
  }

};

}

#endif
