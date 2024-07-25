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

#ifndef DYNAMICS_CONTEXT_D_H
#define DYNAMICS_CONTEXT_D_H

#include "core/fieldParametersD.h"

namespace olb {


/// Base of DynamicsParametersD
/**
 * Used for platform-agnostic access to dynamics parameters.
 **/
template <typename T, typename DESCRIPTOR>
struct AbstractDynamicsParameters {
  virtual AbstractParameters<T,DESCRIPTOR>& asAbstract() = 0;

  virtual void setProcessingContext(ProcessingContext context) = 0;
};

/// Concrete storage of DYNAMICS::ParametersD in olb::Data
template <typename T, typename DESCRIPTOR, Platform PLATFORM, typename DYNAMICS>
struct DynamicsParametersD final : public AbstractDynamicsParameters<T,DESCRIPTOR>
                                 , public Serializable {
  typename DYNAMICS::ParametersD parameters;

  DynamicsParametersD(std::size_t): // TODO: Implement more generic non-cellwise field allocation in Data
    parameters{}
  { }

  /// Return abstract interface to host-side parameters
  AbstractParameters<T,DESCRIPTOR>& asAbstract() override {
    return parameters;
  }

  void setProcessingContext(ProcessingContext context) override { }

  /// Number of data blocks for the serializable interface
  std::size_t getNblock() const override;
  /// Binary size for the serializer
  std::size_t getSerializableSize() const override;
  /// Return a pointer to the memory of the current block and its size for the serializable interface
  bool* getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode) override;

};

template <typename T, typename DESCRIPTOR, Platform PLATFORM, typename DYNAMICS>
std::size_t DynamicsParametersD<T,DESCRIPTOR,PLATFORM,DYNAMICS>::getNblock() const
{
  return decltype(parameters)::fields_t::size;
}

template <typename T, typename DESCRIPTOR, Platform PLATFORM, typename DYNAMICS>
std::size_t DynamicsParametersD<T,DESCRIPTOR,PLATFORM,DYNAMICS>::getSerializableSize() const
{
  std::size_t size = 0;
  decltype(parameters)::fields_t::for_each([&size](auto field) {
    using field_t = typename decltype(field)::type;
    size += FieldD<T,DESCRIPTOR,field_t>{}.getSerializableSize();
  });
  return size;
}

template <typename T, typename DESCRIPTOR, Platform PLATFORM, typename DYNAMICS>
bool* DynamicsParametersD<T,DESCRIPTOR,PLATFORM,DYNAMICS>::getBlock(
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

template <typename T, typename DESCRIPTOR, Platform PLATFORM, typename CONTEXT>
struct TrivialParametersD final : public Serializable {
  typename CONTEXT::ParametersD parameters;

  TrivialParametersD(std::size_t): // TODO: Implement more generic non-cellwise field allocation in Data
    parameters{}
  { }

  /// Return abstract interface to host-side parameters
  typename CONTEXT::ParametersD& asAbstract() {
    return parameters;
  }

  void setProcessingContext(ProcessingContext context) { }

  /// Number of data blocks for the serializable interface
  std::size_t getNblock() const override;
  /// Binary size for the serializer
  std::size_t getSerializableSize() const override;
  /// Return a pointer to the memory of the current block and its size for the serializable interface
  bool* getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode) override;

};

template <typename T, typename DESCRIPTOR, Platform PLATFORM, typename CONTEXT>
std::size_t TrivialParametersD<T,DESCRIPTOR,PLATFORM,CONTEXT>::getNblock() const
{
  return 1;
}

template <typename T, typename DESCRIPTOR, Platform PLATFORM, typename CONTEXT>
std::size_t TrivialParametersD<T,DESCRIPTOR,PLATFORM,CONTEXT>::getSerializableSize() const
{
  return sizeof(typename CONTEXT::ParametersD);
}

template <typename T, typename DESCRIPTOR, Platform PLATFORM, typename CONTEXT>
bool* TrivialParametersD<T,DESCRIPTOR,PLATFORM,CONTEXT>::getBlock(
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
