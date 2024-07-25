/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2020 Adrian Kummerlaender
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

#ifndef SUPER_COMMUNICATOR_HH
#define SUPER_COMMUNICATOR_HH

#include "superCommunicator.h"
#include "superCommunicationTagCoordinator.hh"
#include "blockCommunicationNeighborhood.hh"

#include <algorithm>
#include <stdexcept>

namespace olb {

// *INDENT-OFF*

template <typename T, typename SUPER>
SuperCommunicator<T,SUPER>::SuperCommunicator(
  SUPER& super):
  _super(super)
#ifdef PARALLEL_MODE_MPI
, _tagCoordinator(super.getLoadBalancer())
#endif
{
#ifdef PARALLEL_MODE_MPI
  if (MPI_Comm_dup(MPI_COMM_WORLD, &_neighborhoodComm) != MPI_SUCCESS) {
    throw std::runtime_error("Unable to duplicate MPI communicator");
  }
  if (MPI_Comm_dup(MPI_COMM_WORLD, &_communicatorComm) != MPI_SUCCESS) {
    throw std::runtime_error("Unable to duplicate MPI communicator");
  }
#endif

  auto& cuboidGeometry = _super.getCuboidGeometry();
  auto& load = _super.getLoadBalancer();

  for (int iC = 0; iC < load.size(); ++iC) {
    _blockNeighborhoods.emplace_back(  cuboidGeometry
                                     , load, load.glob(iC)
                                     , _super.getOverlap()
#ifdef PARALLEL_MODE_MPI
                                     , _neighborhoodComm
#endif
                                     );
  }
}

template <typename T, typename SUPER>
SuperCommunicator<T,SUPER>::~SuperCommunicator()
{
#ifdef PARALLEL_MODE_MPI
  MPI_Comm_free(&_neighborhoodComm);
  MPI_Comm_free(&_communicatorComm);
#endif
}

template <typename T, typename SUPER>
void SuperCommunicator<T,SUPER>::exchangeRequests()
{
  auto& load = _super.getLoadBalancer();

  for (int iC = 0; iC < load.size(); ++iC) {
    _blockNeighborhoods[iC].maintain();
  }

  for (int iC = 0; iC < load.size(); ++iC) {
    for (std::type_index field : _fieldsRequested) {
      _blockNeighborhoods[iC].setFieldAvailability(
        field, _super.getBlock(iC).hasCommunicatable(field));
    }
  }

  for (int iC = 0; iC < load.size(); ++iC) {
    _blockNeighborhoods[iC].forLocalNeighbors([&](int localC) {
      _blockNeighborhoods[iC].setFieldsAvailability(localC, _super.getBlock(load.loc(localC)));
    });
  }

#ifdef PARALLEL_MODE_MPI
  _tagCoordinator.template coordinate<SUPER::d>(_blockNeighborhoods);

  for (int iC = 0; iC < load.size(); ++iC) {
    _blockNeighborhoods[iC].send(_tagCoordinator);
  }
  for (int iC = 0; iC < load.size(); ++iC) {
    _blockNeighborhoods[iC].receive(_tagCoordinator);
  }
  for (int iC = 0; iC < load.size(); ++iC) {
    _blockNeighborhoods[iC].wait();
  }
#endif

  _blockCommunicators.clear();
  _blockCommunicators.resize(load.size());
  for (int iC = 0; iC < load.size(); ++iC) {
    auto* block = &_super.getBlock(iC);
    _blockCommunicators[iC] = callUsingConcretePlatform<typename SUPER::block_t>(
      block->getPlatform(),
      block,
      [&](auto* concreteBlock) -> std::unique_ptr<BlockCommunicator> {
        return std::make_unique<ConcreteBlockCommunicator<std::remove_reference_t<decltype(*concreteBlock)>>>(
          _super,
          load,
#ifdef PARALLEL_MODE_MPI
          _tagCoordinator,
          _communicatorComm,
#endif
          iC,
          _blockNeighborhoods[iC]);
      });
  }

#ifdef PARALLEL_MODE_MPI
  for (int iC = 0; iC < load.size(); ++iC) {
    _blockNeighborhoods[iC].forRemoteNeighbors([&](int remoteC) {
      _remoteCuboidNeighborhood.emplace(remoteC);
    });
  }
#endif

  _ready = true;
}

template <typename T, typename SUPER>
void SuperCommunicator<T,SUPER>::requestCell(LatticeR<SUPER::d+1> latticeR)
{
  _blockNeighborhoods[latticeR[0]].requestCell(latticeR.data()+1);
  _ready = false;
  _enabled = true;
}

template <typename T, typename SUPER>
void SuperCommunicator<T,SUPER>::requestOverlap(int width)
{
  auto& load = _super.getLoadBalancer();
  for (int iC = 0; iC < load.size(); ++iC) {
    _blockNeighborhoods[iC].requestOverlap(width);
  }
  _ready = false;
  _enabled = true;
}

template <typename T, typename SUPER>
void SuperCommunicator<T,SUPER>::communicate()
{
  if (!_enabled) {
    return;
  }
  if (!_ready) {
    throw std::logic_error("Requests must be re-exchanged after any changes");
  }

#ifdef PLATFORM_GPU_CUDA
  gpu::cuda::device::synchronize();
#endif

  auto& load = _super.getLoadBalancer();
#ifdef PARALLEL_MODE_MPI
  for (int iC = 0; iC < load.size(); ++iC) {
    _blockCommunicators[iC]->receive();
  }
  for (int iC = 0; iC < load.size(); ++iC) {
    _blockCommunicators[iC]->send();
  }
  for (int iC = 0; iC < load.size(); ++iC) {
    _blockCommunicators[iC]->unpack();
  }
  for (int iC = 0; iC < load.size(); ++iC) {
    _blockCommunicators[iC]->copy();
  }
  for (int iC = 0; iC < load.size(); ++iC) {
    _blockCommunicators[iC]->wait();
  }
#else // PARALLEL_MODE_MPI
  for (int iC = 0; iC < load.size(); ++iC) {
    _blockCommunicators[iC]->copy();
  }
  for (int iC = 0; iC < load.size(); ++iC) {
    _blockCommunicators[iC]->wait();
  }
#endif
}

template <typename T, typename SUPER>
const std::set<int>& SuperCommunicator<T,SUPER>::getRemoteCuboids() const
{
  if (_ready) {
    return _remoteCuboidNeighborhood;
  } else {
    throw std::logic_error("Requests must be re-exchanged after any changes");
  }
}

// *INDENT-ON*

}

#endif
