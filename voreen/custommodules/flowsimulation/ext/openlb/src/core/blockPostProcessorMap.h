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

#ifndef BLOCK_POST_PROCESSOR_MAP_H
#define BLOCK_POST_PROCESSOR_MAP_H

#include <functional>
#include <typeindex>

#include "operator.h"
#include "postProcessing.h"

namespace olb {


/// Factory for instances of a specific POST_PROCESSOR type
/**
 * Factory callable for ConcreteBlockO<POST_PROCESSOR> is constructed at
 * PostProcessorPromise construction time. Recipients accepting such
 * _promised post processors_ are not obligated to actually _realize_
 * this promise.
 *
 * Analogously to DynamicsPromise, this structure is needed to bridge the gap between
 * high level virtual interfaces and efficient platform-specific implementations with
 * full type knowledge.
 **/
template <typename T, typename DESCRIPTOR>
class PostProcessorPromise {
protected:
  std::type_index _id;
  int _priority;
  OperatorScope _scope;
  std::function<AbstractBlockO<T,DESCRIPTOR>*(Platform)> _constructor;

public:
  template <typename POST_PROCESSOR>
  PostProcessorPromise(meta::id<POST_PROCESSOR> id = meta::id<POST_PROCESSOR>{}):
    _id(typeid(POST_PROCESSOR)),
    _priority(POST_PROCESSOR().getPriority()),
    _scope(POST_PROCESSOR::scope),
    _constructor([](Platform platform) -> AbstractBlockO<T,DESCRIPTOR>* {
      switch (platform) {
      #ifdef PLATFORM_CPU_SISD
      case Platform::CPU_SISD:
        return new ConcreteBlockO<T,DESCRIPTOR,Platform::CPU_SISD,POST_PROCESSOR,POST_PROCESSOR::scope>();
      #endif
      #ifdef PLATFORM_CPU_SIMD
      case Platform::CPU_SIMD:
        return new ConcreteBlockO<T,DESCRIPTOR,Platform::CPU_SIMD,POST_PROCESSOR,POST_PROCESSOR::scope>();
      #endif
      #ifdef PLATFORM_GPU_CUDA
      case Platform::GPU_CUDA:
        return new ConcreteBlockO<T,DESCRIPTOR,Platform::GPU_CUDA,POST_PROCESSOR,POST_PROCESSOR::scope>();
      #endif
      default:
        throw std::invalid_argument("Invalid PLATFORM");
      }
    })
  { }

  /// Returns type index of the promised POST_PROCESSOR
  std::type_index id() const {
    return _id;
  }

  int priority() const {
    return _priority;
  }

  OperatorScope scope() const {
    return _scope;
  }

  template <Platform PLATFORM>
  BlockO<T,DESCRIPTOR,PLATFORM>* realize() {
    return static_cast<BlockO<T,DESCRIPTOR,PLATFORM>*>(
      _constructor(PLATFORM));
  };

};

template <typename PP>
PostProcessorPromise(meta::id<PP>) -> PostProcessorPromise<typename PP::value_t,
                                                           typename PP::descriptor_t>;

/// Block operator for supporting legacy post processor in the new operator-centric framework
template <typename T, typename DESCRIPTOR, Platform PLATFORM>
class LegacyBlockPostProcessorO final : public BlockO<T,DESCRIPTOR,PLATFORM> {
private:
  /// List of legacy post processors
  std::vector<std::unique_ptr<PostProcessor<T,DESCRIPTOR>>> _postProcessors;

public:
  std::type_index id() const override
  {
    return typeid(LegacyBlockPostProcessorO);
  }

  void set(CellID iCell, bool state) override
  {
    throw std::logic_error("Invalid legacy post processor setter");
  }

  void setup(ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>& block) override { }

  void apply(ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>& block) override
  {
    #ifdef PLATFORM_GPU_CUDA
    if constexpr (PLATFORM == Platform::GPU_CUDA) {
      if (!_postProcessors.empty()) {
        throw std::runtime_error("Legacy post processors not supported on GPU_CUDA");
      }
    }
    #else // CPU_* platform
    #ifdef PARALLEL_MODE_OMP
    if (static_cast<std::ptrdiff_t>(_postProcessors.size()) >= omp.get_size()) {
      #pragma omp parallel for schedule(dynamic)
      for (auto& postProcessor : _postProcessors) {
        if constexpr (DESCRIPTOR::d == 3) {
          postProcessor->processSubDomain(block, 0, block.getNx()-1, 0, block.getNy()-1, 0, block.getNz()-1);
        } else {
          postProcessor->processSubDomain(block, 0, block.getNx()-1, 0, block.getNy()-1);
        }
      }
      // Disable OpenMP-parallel execution of single post processors (e.g. statistics, bouzidi (hacky))
      // To be done correctly when implementing BlockPostProcessorMap
    } else {
      #endif
      for (auto& postProcessor : _postProcessors) {
        if constexpr (DESCRIPTOR::d == 3) {
          postProcessor->processSubDomain(block, 0, block.getNx()-1, 0, block.getNy()-1, 0, block.getNz()-1);
        } else {
          postProcessor->processSubDomain(block, 0, block.getNx()-1, 0, block.getNy()-1);
        }
      }
      #ifdef PARALLEL_MODE_OMP
    }
    #endif
    #endif
  }

  void add(PostProcessor<T,DESCRIPTOR>* postProcessor)
  {
    _postProcessors.emplace_back(postProcessor);
  }

};

/// Map of post processors of a single priority and stage
/**
 * Maintained in ConcreteBlockLattice
 **/
template<typename T, typename DESCRIPTOR, Platform PLATFORM>
class BlockPostProcessorMap {
private:
  /// Lattice reference to be passed to BlockO::setup and used for apply
  ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>& _lattice;

  /// Type indexed map of BlockO instances
  /**
   * Only add new BlockO instances to this map via resolve
   **/
  std::map<std::type_index,
           std::unique_ptr<BlockO<T,DESCRIPTOR,PLATFORM>>> _map;

  /// BlockO for managing legacy post processors
  /**
   * Will be removed once all post processors are refactored into operators
   **/
  LegacyBlockPostProcessorO<T,DESCRIPTOR,PLATFORM> _legacyPostProcessors;

  /// Resolve promised post processor into concrete BlockO instance
  BlockO<T,DESCRIPTOR,PLATFORM>& resolve(PostProcessorPromise<T,DESCRIPTOR>&& promise)
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

public:
  BlockPostProcessorMap(ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>* lattice):
    _lattice{*lattice}
  { }

  void addLegacy(PostProcessor<T,DESCRIPTOR>* postProcessor)
  {
    _legacyPostProcessors.add(postProcessor);
  }

  void add(std::size_t iCell, PostProcessorPromise<T,DESCRIPTOR>&& promise)
  {
    resolve(std::forward<decltype(promise)>(promise)).set(iCell, true);
  }

  void add(PostProcessorPromise<T,DESCRIPTOR>&& promise)
  {
    resolve(std::forward<decltype(promise)>(promise));
  }

  /// Apply all managed post processors to lattice
  /**
   * All post processors within a single BlockPostProcessorMap should be expected
   * to be executed in parallel as far as feasible. E.g. (non-)legacy post processors
   * on CPU blocks are parallelized using OpenMP if enabled while per-cell operators
   * on GPU blocks are parallelized both on the level of kernels and by running
   * multiple kernels asynchronously in separate streams.
   *
   * In general, interdependent post processors should be separated into different
   * priorities (not doing so also works in many cases but leads to non-deterministic
   * data relationships).
   **/
  void apply()
  {
    _legacyPostProcessors.apply(_lattice);

    #ifdef PARALLEL_MODE_OMP
    if (_map.size() <= 1) {
      // prevent nested OMP for statistics post processor (hacky, to be reimplemented)
      for (auto& [_, postProcessor] : _map) {
        postProcessor->apply(_lattice);
      }
    } else {
      #pragma omp parallel
      {
        #pragma omp single
        {
          for (auto& [_, postProcessor] : _map) {
            auto* ptr = postProcessor.get();
            #pragma omp task
            {
              ptr->apply(_lattice);
            }
          }
        }
      }
    }
    #else
    for (auto& [_, postProcessor] : _map) {
      postProcessor->apply(_lattice);
    }
    #endif

    #ifdef PLATFORM_GPU_CUDA
    if constexpr (PLATFORM == Platform::GPU_CUDA) {
      gpu::cuda::device::synchronize();
    }
    #endif
  }

};


}

#endif
