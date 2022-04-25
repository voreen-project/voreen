/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2019 Adrian Kummerlaender
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

#ifndef BLOCK_DYNAMICS_MAP_H
#define BLOCK_DYNAMICS_MAP_H

#include <functional>
#include <typeindex>
#include <algorithm>

#include "dynamics/dynamics.h"
#include "platform/cpu/cell.h"

namespace olb {

template<typename T, typename DESCRIPTOR> struct AbstractCollisionO;
template<typename T, typename DESCRIPTOR, Platform PLATFORM> struct BlockCollisionO;
template<typename T, typename DESCRIPTOR, Platform PLATFORM> class ConcreteBlockLattice;
template<typename T, typename DESCRIPTOR, Platform PLATFORM> class ConcreteBlockMask;
template<typename T, typename DESCRIPTOR, Platform PLATFORM, typename DYNAMICS> class ConcreteBlockCollisionO;

/// Concrete CPU dynamics for legacy dynamics
/**
 * Required to enable interaction between generic operators / non-legacy
 * post processors and legacy dynamics. Note that this is only supported
 * on CPU targets.
 **/
template <typename T, typename DESCRIPTOR, Platform PLATFORM>
class LegacyConcreteDynamics final : public cpu::Dynamics<T,DESCRIPTOR,PLATFORM> {
private:
  BlockLattice<T,DESCRIPTOR>& _lattice;
  Dynamics<T,DESCRIPTOR>* _dynamics;

public:
  LegacyConcreteDynamics(ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>& lattice,
                         Dynamics<T,DESCRIPTOR>* dynamics):
    _lattice(lattice),
    _dynamics{dynamics} {
  }

  CellStatistic<T> collide(cpu::Cell<T,DESCRIPTOR,PLATFORM>& cell) override {
    auto abstractCell = _lattice.get(cell.getCellId());
    return _dynamics->collide(abstractCell);
  }

  T computeRho(cpu::Cell<T,DESCRIPTOR,PLATFORM>& cell) override {
    auto abstractCell = _lattice.get(cell.getCellId());
    return _dynamics->computeRho(abstractCell);
  }
  void computeU(cpu::Cell<T,DESCRIPTOR,PLATFORM>& cell, T* u) override {
    auto abstractCell = _lattice.get(cell.getCellId());
    _dynamics->computeU(abstractCell, u);
  }
  void computeJ(cpu::Cell<T,DESCRIPTOR,PLATFORM>& cell, T* j) override {
    auto abstractCell = _lattice.get(cell.getCellId());
    _dynamics->computeJ(abstractCell, j);
  }
  void computeRhoU(cpu::Cell<T,DESCRIPTOR,PLATFORM>& cell, T& rho, T* u) override {
    auto abstractCell = _lattice.get(cell.getCellId());
    _dynamics->computeRhoU(abstractCell, rho, u);
  }
  void computeStress(cpu::Cell<T,DESCRIPTOR,PLATFORM>& cell, T& rho, T* u, T* pi) override {
    auto abstractCell = _lattice.get(cell.getCellId());
    _dynamics->computeStress(abstractCell, rho, u, pi);
  }
  void computeAllMomenta(cpu::Cell<T,DESCRIPTOR,PLATFORM>& cell, T& rho, T* u, T* pi) override {
    auto abstractCell = _lattice.get(cell.getCellId());
    _dynamics->computeAllMomenta(abstractCell, rho, u, pi);
  }

  T getOmegaOrFallback(T fallback) override {
    throw std::runtime_error("getOmegaOrFallback not supported for legacy dynamics");
  }

  T computeEquilibrium(int iPop, T rho, T* u) override {
    return _dynamics->computeEquilibrium(iPop, rho, u);
  }

};

/// Concrete collision operator for legacy dynamics
/**
 * To be removed once both all OpenLB-included dynamics are adapted
 * and legacy support was included in at least one release.
 *
 * Legacy dynamics are only supported on CPU targets.
 **/
template <typename T, typename DESCRIPTOR, Platform PLATFORM>
class LegacyBlockCollisionO final : public BlockCollisionO<T,DESCRIPTOR,PLATFORM> {
private:
  ConcreteBlockMask<T,DESCRIPTOR,PLATFORM> _mask;
  Dynamics<T,DESCRIPTOR>** const _legacyDynamics;
  ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>* _block;

  /// Map of concrete dynamics to support legacy dynamics in the context of
  /// generic operators / non-legacy post processors
  std::map<Dynamics<T,DESCRIPTOR>*,
           std::unique_ptr<cpu::Dynamics<T,DESCRIPTOR,PLATFORM>>> _concreteDynamics;
  cpu::Dynamics<T,DESCRIPTOR,PLATFORM>** _dynamicsOfCells;

  /// Resolve legacy dynamics at iCell into matching LegacyConcreteDynamics
  cpu::Dynamics<T,DESCRIPTOR,PLATFORM>* resolve(CellID iCell)
  {
    auto iter = _concreteDynamics.find(_legacyDynamics[iCell]);
    if (iter == _concreteDynamics.end()) {
      iter = std::get<0>(_concreteDynamics.emplace(
        _legacyDynamics[iCell],
        new LegacyConcreteDynamics<T,DESCRIPTOR,PLATFORM>(*_block,
                                                          _legacyDynamics[iCell])));
    }
    return std::get<1>(*iter).get();
  }

public:
  LegacyBlockCollisionO(std::size_t count, Dynamics<T,DESCRIPTOR>** dynamics):
    _mask{count},
    _legacyDynamics{dynamics}
  { }

  std::type_index id() const override
  {
    return typeid(LegacyBlockCollisionO);
  }

  std::size_t weight() const override
  {
    return _mask.weight();
  }

  void set(CellID iCell, bool state, bool overlap) override
  {
    if (!overlap) {
      _mask.set(iCell, state);
    }
    if (state) {
      _dynamicsOfCells[iCell] = resolve(iCell);
    }
  }

  Dynamics<T,DESCRIPTOR>* getDynamics() override
  {
    return nullptr;
  }

  void setup(ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>& block) override
  {
    // Remember block in order to construct LegacyConcreteDynamics
    _block = &block;
    // Fetch pointer to concretized dynamic-dispatch field
    _dynamicsOfCells = block.template getField<cpu::DYNAMICS<T,DESCRIPTOR,PLATFORM>>()[0].data();
  }

  /// Apply collision on subdomain of block
  /**
   * Loop excludes overlap areas of block as collisions are never applied there.
   * This assumes that `subdomain` is the core mask of BlockDynamicsMap.
   **/
  void apply(ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>& block,
             ConcreteBlockMask<T,DESCRIPTOR,PLATFORM>&    subdomain,
             CollisionDispatchStrategy                    strategy) override
  {
    if (strategy != CollisionDispatchStrategy::Dominant) {
      throw std::runtime_error("LegacyBlockCollisionO only supports CollisionDispatchStrategy::Dominant");
    }

    if constexpr (DESCRIPTOR::d == 3) {
      auto cell = block.get(0,0,0);
      #ifdef PARALLEL_MODE_OMP
      #pragma omp parallel for schedule(dynamic,1) firstprivate(cell)
      #endif
      for (int iX=0; iX < block.getNx(); ++iX) {
        for (int iY=0; iY < block.getNy(); ++iY) {
          for (int iZ=0; iZ < block.getNz(); ++iZ) {
            CellID iCell = block.getCellId(iX,iY,iZ);
            if (_mask[iCell]) {
              cell.setCellId(iCell);
              if (auto cellStatistic = _legacyDynamics[iCell]->collide(cell)) {
                block.getStatistics().incrementStats(cellStatistic.rho, cellStatistic.uSqr);
              }
            } else {
              cpu::Cell<T,DESCRIPTOR,PLATFORM> cell(block, iCell);
              if (auto cellStatistic = _dynamicsOfCells[iCell]->collide(cell)) {
                block.getStatistics().incrementStats(cellStatistic.rho, cellStatistic.uSqr);
              }
            }
          }
        }
      }
    } else {
      auto cell = block.get(0,0);
      #ifdef PARALLEL_MODE_OMP
      #pragma omp parallel for firstprivate(cell)
      #endif
      for (int iX=0; iX < block.getNx(); ++iX) {
        for (int iY=0; iY < block.getNy(); ++iY) {
          CellID iCell = block.getCellId(iX,iY);
          if (_mask[iCell]) {
            cell.setCellId(iCell);
            if (auto cellStatistic = _legacyDynamics[iCell]->collide(cell)) {
              block.getStatistics().incrementStats(cellStatistic.rho, cellStatistic.uSqr);
            }
          } else {
            cpu::Cell<T,DESCRIPTOR,PLATFORM> cell(block, iCell);
            if (auto cellStatistic = _dynamicsOfCells[iCell]->collide(cell)) {
              block.getStatistics().incrementStats(cellStatistic.rho, cellStatistic.uSqr);
            }
          }
        }
      }
    }
  }
};

/// Factory for instances of a specific Dynamics type
/**
 * Factory callables for DYNAMICS and ConcreteBlockCollisionO<DYNAMICS> are constructed
 * at DynamicsPromise construction time. Recipients accepting such _promised dynamics_
 * are not obligated to actually _realize_ this promise.
 *
 * This structure is needed to bridge the gap between high level virtual interfaces and
 * efficient platform-specific implementations with full type knowledge.
 * i.e. DynamicsPromise can carry full dynamics type information through the virtual
 * layer surrounding concrete platform-specialized block lattices (needed as virtual
 * template methods are not supported by C++).
 **/
template <typename T, typename DESCRIPTOR>
class DynamicsPromise {
protected:
  std::type_index _id;
  std::function<Dynamics<T,DESCRIPTOR>*()>                   _dynamicsConstructor;
  std::function<AbstractCollisionO<T,DESCRIPTOR>*(Platform)> _operatorConstructor;

public:
  template <typename DYNAMICS>
  DynamicsPromise(meta::id<DYNAMICS> id = meta::id<DYNAMICS>{}):
    _id(typeid(DYNAMICS)),
    _dynamicsConstructor([]() -> Dynamics<T,DESCRIPTOR>* {
      return new DYNAMICS();
    }),
    _operatorConstructor([](Platform platform) -> AbstractCollisionO<T,DESCRIPTOR>* {
      static_assert(dynamics::is_generic_v<T,DESCRIPTOR,DYNAMICS>,
                    "Promised DYNAMICS must be generic");
      switch (platform) {
      #ifdef PLATFORM_CPU_SISD
      case Platform::CPU_SISD:
        return new ConcreteBlockCollisionO<T,DESCRIPTOR,Platform::CPU_SISD,DYNAMICS>();
      #endif
      #ifdef PLATFORM_CPU_SIMD
      case Platform::CPU_SIMD:
        return new ConcreteBlockCollisionO<T,DESCRIPTOR,Platform::CPU_SIMD,DYNAMICS>();
      #endif
      #ifdef PLATFORM_GPU_CUDA
      case Platform::GPU_CUDA:
        return new ConcreteBlockCollisionO<T,DESCRIPTOR,Platform::GPU_CUDA,DYNAMICS>();
      #endif
      default:
        throw std::invalid_argument("Invalid PLATFORM");
      }
    })
  { }

  /// Returns type index of the promised DYNAMICS
  std::type_index id() const {
    return _id;
  }

  /// Returns new instance of the promised DYNAMICS
  Dynamics<T,DESCRIPTOR>* realize() {
    return _dynamicsConstructor();
  };

  /// Returns new instance of the collision operator for promised DYNAMICS
  template <Platform PLATFORM>
  BlockCollisionO<T,DESCRIPTOR,PLATFORM>* realize() {
    return static_cast<BlockCollisionO<T,DESCRIPTOR,PLATFORM>*>(
      _operatorConstructor(PLATFORM));
  };

};

template <typename DYNAMICS>
DynamicsPromise(meta::id<DYNAMICS>) -> DynamicsPromise<typename DYNAMICS::value_t,
                                                       typename DYNAMICS::descriptor_t>;


/// Mask describing the subdomain on which to apply the collision step
/**
 * By default this is the entire non-overlap area of the block lattice
 **/
struct CollisionSubdomainMask {
  template <typename T, typename DESCRIPTOR, Platform PLATFORM>
  using type = ConcreteBlockMask<T,DESCRIPTOR,PLATFORM>;
};

/// Map between cell indices and concrete dynamics
/**
 * Central class for managing and applying dynamics of / to a block lattice.
 *
 * Maintains both a map of cell indices to dynamics and
 * of dynamics to a mask of all assigned cell indices.
 *
 * Automatically determines the dominant dynamics and applies
 * the collision step to all cells via the dominant dynamics'
 * ConcreteBlockCollisionO.
 **/
template<typename T, typename DESCRIPTOR, Platform PLATFORM>
class BlockDynamicsMap {
private:
  ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>& _lattice;

  std::map<std::type_index,
           std::unique_ptr<BlockCollisionO<T,DESCRIPTOR,PLATFORM>>> _map;

  std::unique_ptr<Dynamics<T,DESCRIPTOR>*[]>                 _dynamicsOfCells;
  std::unique_ptr<BlockCollisionO<T,DESCRIPTOR,PLATFORM>*[]> _operatorOfCells;

  /// Subdomain on which to apply collisions
  ConcreteBlockMask<T,DESCRIPTOR,PLATFORM>& _coreMask;
  /// Pointer to collision operator with highest cell fraction
  BlockCollisionO<T,DESCRIPTOR,PLATFORM>* _dominantCollisionO;

  /// Returns (Dynamics, CollisionO, Mask) tuple for promise
  BlockCollisionO<T,DESCRIPTOR,PLATFORM>& resolve(DynamicsPromise<T,DESCRIPTOR>&& promise)
  {
    auto iter = _map.find(promise.id());
    if (iter == _map.end()) {
      iter = _map.emplace(std::piecewise_construct,
                          std::forward_as_tuple(promise.id()),
                          std::forward_as_tuple(promise.template realize<PLATFORM>())).first;
      iter->second->setup(_lattice);
    }
    return *(iter->second);
  }

  /// Assigns dynamics to cell index and updates block masks
  void set(std::size_t iCell, BlockCollisionO<T,DESCRIPTOR,PLATFORM>& collisionO)
  {
    if (_operatorOfCells[iCell] != nullptr) {
      _map.at(_operatorOfCells[iCell]->id())->set(iCell, false, !_coreMask[iCell]);
    }
    collisionO.set(iCell, true, !_coreMask[iCell]);
    _operatorOfCells[iCell] = &collisionO;
    if (Dynamics<T,DESCRIPTOR>* dynamics = collisionO.getDynamics()) {
      _dynamicsOfCells[iCell] = dynamics;
    }
  }

public:
  /// Constructor for a BlockDynamicsMap
  BlockDynamicsMap(ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>& lattice):
    _lattice(lattice),
    _dynamicsOfCells(new Dynamics<T,DESCRIPTOR>*                [_lattice.getNcells()] { nullptr }),
    _operatorOfCells(new BlockCollisionO<T,DESCRIPTOR,PLATFORM>*[_lattice.getNcells()] { nullptr }),
    _coreMask(lattice.template getData<CollisionSubdomainMask>()),
    _dominantCollisionO(nullptr)
  {
    _lattice.forCoreSpatialLocations([&](LatticeR<DESCRIPTOR::d> lattice) {
      _coreMask.set(_lattice.getCellId(lattice), true);
    });

    if constexpr (isPlatformCPU(PLATFORM)) {
      /// Setup support for legacy dynamics
      using LegacyO = LegacyBlockCollisionO<T,DESCRIPTOR,PLATFORM>;
      _map[typeid(LegacyO)] = std::make_unique<LegacyO>(_lattice.getNcells(),
                                                        _dynamicsOfCells.get());
      _map[typeid(LegacyO)]->setup(_lattice);
    }
  }

  /// Returns a pointer to dynamics matching the promise
  /**
   * Promise is only realized if dynamics were not previously constructed
   **/
  Dynamics<T,DESCRIPTOR>* get(DynamicsPromise<T,DESCRIPTOR>&& promise)
  {
    return resolve(std::forward<decltype(promise)>(promise)).getDynamics();
  }

  const Dynamics<T,DESCRIPTOR>* get(std::size_t iCell) const
  {
    return _dynamicsOfCells[iCell];
  }

  Dynamics<T,DESCRIPTOR>* get(std::size_t iCell)
  {
    return _dynamicsOfCells[iCell];
  }

  /// Assigns promised dynamics to cell index iCell
  void set(std::size_t iCell, DynamicsPromise<T,DESCRIPTOR>&& promise)
  {
    set(iCell, resolve(std::forward<decltype(promise)>(promise)));
    _dominantCollisionO = nullptr;
  }

  /// Assigns dynamics pointer to cell index iCell
  /**
   * Dynamics at the pointed address must either be legacy or previously
   * constructed from a promise.
   **/
  void set(std::size_t iCell, Dynamics<T,DESCRIPTOR>* dynamics)
  {
    auto iter = _map.find(dynamics->id());
    if (iter != _map.end()) { // Non-legacy dynamics
      set(iCell, *iter->second);
    } else { // Legacy dynamics
      iter = _map.find(typeid(LegacyBlockCollisionO<T,DESCRIPTOR,PLATFORM>));
      if (iter != _map.end()) {
        // Order is important here as LegacyBlockCollisionO uses the dynamics
        // pointer to construct matching LegacyConcreteDynamics internally.
        _dynamicsOfCells[iCell] = dynamics;
        set(iCell, *iter->second);
      } else {
        throw std::runtime_error("Legacy dynamics not supported on this platform");
      }
    }
    _dominantCollisionO = nullptr;
  }

  /// Executes local collision step for entire non-overlap area of lattice
  /**
   * When using the default CollisionDispatchStrategy::Dominant, the most frequently
   * assigned dynamics resp. collision operator is applied and responsible for also
   * colliding cells excluded by its masks using their respective dynamics.
   *
   * Correspondingly, the alternative CollisionDispatchStrategy::Individual applies
   * all dynamics separately using a list-based approach.
   **/
  void collide(CollisionDispatchStrategy strategy)
  {
    switch (strategy) {
    case CollisionDispatchStrategy::Dominant:
      if (!_dominantCollisionO) {
        _dominantCollisionO = std::max_element(_map.begin(),
                                               _map.end(),
                                               [](const auto& lhs, const auto& rhs) -> bool {
                                                 return lhs.second->weight() < rhs.second->weight();
                                               })->second.get();
      }
      _dominantCollisionO->apply(_lattice, _coreMask, strategy);
      break;

    case CollisionDispatchStrategy::Individual:
      for (auto& [id, collisionO] : _map) {
        if (collisionO->weight() > 0) {
          collisionO->apply(_lattice, _coreMask, strategy);
        }
      }
      #ifdef PLATFORM_GPU_CUDA
      gpu::cuda::device::synchronize();
      #endif
      break;

    default:
      throw std::runtime_error("Invalid collision dispatch strategy");
      break;
    }
  }

  /// Returns a human-readable string listing all managed dynamics and their assigned fraction of cells
  std::string describe()
  {
    std::stringstream out;
    for (auto& [_, collisionO] : _map) {
      out << collisionO->getDynamics()->getName() << ", "
          << static_cast<double>(collisionO.weight()) / (_coreMask.weight())
          << std::endl;
    }
    return out.str();
  }

};


}

#endif
