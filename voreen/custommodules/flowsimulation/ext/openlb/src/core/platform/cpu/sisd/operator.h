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

#ifndef SISD_OPERATOR_H
#define SISD_OPERATOR_H

#include "mask.h"
#include "core/operator.h"
#include "core/platform/cpu/cell.h"

namespace olb {

namespace cpu {

/// Implementations of scalar CPU specifics
namespace sisd {

/// Implementation of cpu::Dynamics for concrete DYNAMICS on SISD blocks
template <typename T, typename DESCRIPTOR, typename DYNAMICS>
class ConcreteDynamics final : public cpu::Dynamics<T,DESCRIPTOR,Platform::CPU_SISD> {
private:
  typename DYNAMICS::ParametersD* _parameters;

public:
  ConcreteDynamics(typename DYNAMICS::ParametersD* parameters):
    _parameters{parameters} {
  }

  CellStatistic<T> collide(cpu::Cell<T,DESCRIPTOR,Platform::CPU_SISD>& cell) override {
    return DYNAMICS().apply(cell, *_parameters);
  }

  T computeRho(cpu::Cell<T,DESCRIPTOR,Platform::CPU_SISD>& cell) override {
    return typename DYNAMICS::MomentaF().computeRho(cell);
  }
  void computeU(cpu::Cell<T,DESCRIPTOR,Platform::CPU_SISD>& cell, T* u) override {
    typename DYNAMICS::MomentaF().computeU(cell, u);
  }
  void computeJ(cpu::Cell<T,DESCRIPTOR,Platform::CPU_SISD>& cell, T* j) override {
    typename DYNAMICS::MomentaF().computeJ(cell, j);
  }
  void computeRhoU(cpu::Cell<T,DESCRIPTOR,Platform::CPU_SISD>& cell, T& rho, T* u) override {
    typename DYNAMICS::MomentaF().computeRhoU(cell, rho, u);
  }
  void computeStress(cpu::Cell<T,DESCRIPTOR,Platform::CPU_SISD>& cell, T& rho, T* u, T* pi) override {
    typename DYNAMICS::MomentaF().computeStress(cell, rho, u, pi);
  }
  void computeAllMomenta(cpu::Cell<T,DESCRIPTOR,Platform::CPU_SISD>& cell, T& rho, T* u, T* pi) override {
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

/// Application of the collision step on a concrete SISD block
/**
 * ConcreteBlockCollisionO::apply allows for applying DYNAMICS
 * on a given cell index or an entire block. The latter option
 * accepts a ConcreteBlockMask instance describing the subset
 * of cells for which DYNAMICS is to be called directly
 * instead of via roundtrip.
 **/
template <typename T, typename DESCRIPTOR, typename DYNAMICS>
class ConcreteBlockCollisionO<T,DESCRIPTOR,Platform::CPU_SISD,DYNAMICS> final
  : public BlockCollisionO<T,DESCRIPTOR,Platform::CPU_SISD> {
private:
  std::unique_ptr<DYNAMICS> _dynamics;
  std::unique_ptr<cpu::Dynamics<T,DESCRIPTOR,Platform::CPU_SISD>> _concreteDynamics;

  DynamicsParametersD<T,DESCRIPTOR,Platform::CPU_SISD,DYNAMICS>* _parameters;
  ConcreteBlockMask<T,DESCRIPTOR,Platform::CPU_SISD>*            _mask;

  cpu::Dynamics<T,DESCRIPTOR,Platform::CPU_SISD>** _dynamicsOfCells;

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

  void setup(ConcreteBlockLattice<T,DESCRIPTOR,Platform::CPU_SISD>& block) override
  {
    _parameters = &block.template getData<DynamicsParameters<DYNAMICS>>();
    _mask = &block.template getData<DynamicsMask<DYNAMICS>>();
    if constexpr (dynamics::has_parametrized_momenta_v<DYNAMICS>) {
      _dynamics->setMomentaParameters(&_parameters->parameters);
    }

    _concreteDynamics.reset(new cpu::sisd::ConcreteDynamics<T,DESCRIPTOR,DYNAMICS>(&_parameters->parameters));
    // Fetch pointer to concretized dynamic-dispatch field
    _dynamicsOfCells = block.template getField<cpu::DYNAMICS<T,DESCRIPTOR,Platform::CPU_SISD>>()[0].data();
  }

  /// Apply collision on subdomain of block
  /**
   * Loop excludes overlap areas of block as collisions are never applied there.
   * This assumes that `subdomain` is the core mask of BlockDynamicsMap.
   **/
  void apply(ConcreteBlockLattice<T,DESCRIPTOR,Platform::CPU_SISD>& block,
             ConcreteBlockMask<T,DESCRIPTOR,Platform::CPU_SISD>&    subdomain,
             CollisionDispatchStrategy                              strategy) override
  {
    if (strategy != CollisionDispatchStrategy::Dominant) {
      throw std::runtime_error("Platform::CPU_SISD currently only support CollisionDispatchStrategy::Dominant");
    }

    auto& parameters = _parameters->parameters;
    auto& mask = *_mask;
    if constexpr (DESCRIPTOR::d == 3) {
      #ifdef PARALLEL_MODE_OMP
      #pragma omp parallel for schedule(dynamic,1)
      #endif
      for (int iX=0; iX < block.getNx(); ++iX) {
        auto& statistics = block.getStatistics();
        for (int iY=0; iY < block.getNy(); ++iY) {
          std::size_t iCell = block.getCellId(iX,iY,0);
          for (int iZ=0; iZ < block.getNz(); ++iZ) {
            cpu::Cell<T,DESCRIPTOR,Platform::CPU_SISD> cell(block, iCell);
            if (auto cellStatistic = mask[iCell] ? DYNAMICS().apply(cell, parameters)
                                                 : _dynamicsOfCells[iCell]->collide(cell)) {
              statistics.incrementStats(cellStatistic.rho, cellStatistic.uSqr);
            }
            iCell += 1;
          }
        }
      }
    } else {
      #ifdef PARALLEL_MODE_OMP
      #pragma omp parallel for schedule(dynamic,1)
      #endif
      for (int iX=0; iX < block.getNx(); ++iX) {
        auto& statistics = block.getStatistics();
        std::size_t iCell = block.getCellId(iX,0);
        for (int iY=0; iY < block.getNy(); ++iY) {
          cpu::Cell<T,DESCRIPTOR,Platform::CPU_SISD> cell(block, iCell);
          if (auto cellStatistic = mask[iCell] ? DYNAMICS().apply(cell, parameters)
                                               : _dynamicsOfCells[iCell]->collide(cell)) {
            statistics.incrementStats(cellStatistic.rho, cellStatistic.uSqr);
          }
          iCell += 1;
        }
      }
    }
  }

};


/// Application of a cell-wise OPERATOR on a concrete scalar CPU block
template <typename T, typename DESCRIPTOR, typename OPERATOR>
class ConcreteBlockO<T,DESCRIPTOR,Platform::CPU_SISD,OPERATOR,OperatorScope::PerCell> final
  : public BlockO<T,DESCRIPTOR,Platform::CPU_SISD> {
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

  void setup(ConcreteBlockLattice<T,DESCRIPTOR,Platform::CPU_SISD>& block) override
  { }

  void apply(ConcreteBlockLattice<T,DESCRIPTOR,Platform::CPU_SISD>& block) override
  {
    if (_cells.size() > 0) {
      cpu::Cell<T,DESCRIPTOR,Platform::CPU_SISD> cell(block, 0);
      for (CellID iCell : _cells) {
        cell.setCellId(iCell);
        OPERATOR().apply(cell);
      }
    }
  }

};


template <typename T, typename DESCRIPTOR, typename OPERATOR>
class ConcreteBlockO<T,DESCRIPTOR,Platform::CPU_SISD,OPERATOR,OperatorScope::PerCellWithParameters> final
  : public BlockO<T,DESCRIPTOR,Platform::CPU_SISD> {
private:
  std::vector<CellID> _cells;

  TrivialParametersD<T,DESCRIPTOR,Platform::CPU_SISD,OPERATOR>* _parameters;

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

  void setup(ConcreteBlockLattice<T,DESCRIPTOR,Platform::CPU_SISD>& block) override
  {
    _parameters = &block.template getData<TrivialParameters<OPERATOR>>();
  }

  void apply(ConcreteBlockLattice<T,DESCRIPTOR,Platform::CPU_SISD>& block) override
  {
    if (_cells.size() > 0) {
      cpu::Cell<T,DESCRIPTOR,Platform::CPU_SISD> cell(block, 0);
      for (CellID iCell : _cells) {
        cell.setCellId(iCell);
        OPERATOR().apply(cell, _parameters->parameters);
      }
    }
  }

};


/// Application of a block-wise OPERATOR on a concrete scalar CPU block
/**
 * e.g. StatisticsPostProcessor
 **/
template <typename T, typename DESCRIPTOR, typename OPERATOR>
class ConcreteBlockO<T,DESCRIPTOR,Platform::CPU_SISD,OPERATOR,OperatorScope::PerBlock> final
  : public BlockO<T,DESCRIPTOR,Platform::CPU_SISD> {
public:
  std::type_index id() const override
  {
    return typeid(OPERATOR);
  }

  void set(CellID iCell, bool state) override
  {
    throw std::logic_error("BlockO::set not supported for OperatorScope::PerBlock");
  }

  void setup(ConcreteBlockLattice<T,DESCRIPTOR,Platform::CPU_SISD>& block) override
  {
    OPERATOR().setup(block);
  }

  void apply(ConcreteBlockLattice<T,DESCRIPTOR,Platform::CPU_SISD>& block) override
  {
    OPERATOR().apply(block);
  }

};

}

#endif
