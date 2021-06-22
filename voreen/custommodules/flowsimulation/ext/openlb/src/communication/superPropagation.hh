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

#ifndef SUPER_PROPAGATION_HH
#define SUPER_PROPAGATION_HH

#include "superPropagation.h"
#include "core/olbDebug.h"
#include "core/superLattice2D.h"
#include "core/superLattice3D.h"

#include <algorithm>

namespace olb {


#ifdef PARALLEL_MODE_MPI

template <typename T>
class SuperCommunicationTagCoordinator<T>::ChannelId {
private:
  const int _iC;
  const int _jC;

public:
  ChannelId(int iC, int jC):
    _iC(std::min(iC,jC)),
    _jC(std::max(iC,jC)) { }

  bool operator==(const ChannelId& rhs) const {
    return _iC == rhs._iC
        && _jC == rhs._jC;
  }

  bool operator<(const ChannelId& rhs) const {
    return  _iC  < rhs._iC
        || (_iC == rhs._iC && _jC < rhs._jC);
  }

};

template <typename T>
SuperCommunicationTagCoordinator<T>::SuperCommunicationTagCoordinator(LoadBalancer<T>& loadBalancer):
  _loadBalancer(loadBalancer) { }

template <typename T>
template <typename DESCRIPTOR>
void SuperCommunicationTagCoordinator<T>::init(
  std::vector<BlockStaticPopulationPropagationNeighborhood<T,DESCRIPTOR>>& neighborhood)
{
  OLB_PRECONDITION(!_initialized);

  for (int iC = 0; iC < _loadBalancer.size(); ++iC) {
    for (const auto partner : neighborhood[iC].getOutCells()) {
      const int               jC    = partner.first;
      const std::vector<int>& cells = partner.second;
      if (!cells.empty() && !_loadBalancer.isLocal(jC)) {
        _tags[_loadBalancer.rank(jC)][{_loadBalancer.glob(iC),jC}] = -1;
      }
    }
  }

  for (auto rank=_tags.begin(); rank != _tags.end(); ++rank) {
    for (auto tag=(*rank).second.begin(); tag != (*rank).second.end(); ++tag) {
      (*tag).second = std::distance((*rank).second.begin(), tag);
    }
  }

  _initialized = true;
}

template <typename T>
int SuperCommunicationTagCoordinator<T>::get(int iC, int jC) {
  return _tags[_loadBalancer.rank(iC)][{iC,jC}];
}

#endif // PARALLEL_MODE_MPI

template <typename T, typename DESCRIPTOR>
SuperPropagationCommunicator<T,DESCRIPTOR>::SuperPropagationCommunicator(
  SuperLattice<T,DESCRIPTOR>& superLattice):
  _superLattice(superLattice)
#ifdef PARALLEL_MODE_MPI
  , _tagCoordinator(superLattice.getLoadBalancer())
#endif
{
#ifdef PARALLEL_MODE_MPI
  MPI_Comm_dup(MPI_COMM_WORLD, &_mpi_propagation_comm);
#endif
}

template <typename T, typename DESCRIPTOR>
void SuperPropagationCommunicator<T,DESCRIPTOR>::init() {
  OLB_PRECONDITION(!_initialized);

  auto& load = _superLattice.getLoadBalancer();

  for (int iC = 0; iC < load.size(); ++iC) {
    _blockNeighborhoods.emplace_back(_superLattice, load.glob(iC), 1);
  }

#ifdef PARALLEL_MODE_MPI
  for (int iC = 0; iC < load.size(); ++iC) {
    _blockNeighborhoods[iC].send();
  }
  for (int iC = 0; iC < load.size(); ++iC) {
    _blockNeighborhoods[iC].receive();
  }
  for (int iC = 0; iC < load.size(); ++iC) {
    _blockNeighborhoods[iC].wait();
  }

  _tagCoordinator.template init<DESCRIPTOR>(_blockNeighborhoods);
#endif

  for (int iC = 0; iC < load.size(); ++iC) {
    _blockCommunicators.emplace_back(
      _superLattice,
#ifdef PARALLEL_MODE_MPI
      _tagCoordinator,
      _mpi_propagation_comm,
#endif
      load.glob(iC),
      _blockNeighborhoods[iC].getOutCells(),
      _blockNeighborhoods[iC].getInCells(),
      _blockNeighborhoods[iC].getRemoteOutCells());
  }

  _initialized = true;
}

template <typename T, typename DESCRIPTOR>
void SuperPropagationCommunicator<T,DESCRIPTOR>::propagate() {
  OLB_PRECONDITION(_initialized);

  auto& load = _superLattice.getLoadBalancer();
#ifdef PARALLEL_MODE_MPI
  for (int iC = 0; iC < load.size(); ++iC) {
    _blockCommunicators[iC].receive();
  }
  for (int iC = 0; iC < load.size(); ++iC) {
    _blockCommunicators[iC].send();
  }
  for (int iC = 0; iC < load.size(); ++iC) {
    _blockCommunicators[iC].unpack();
  }
  for (int iC = 0; iC < load.size(); ++iC) {
    _blockCommunicators[iC].copy();
  }
  for (int iC = 0; iC < load.size(); ++iC) {
    _blockCommunicators[iC].wait();
  }
#else // PARALLEL_MODE_MPI
  for (int iC = 0; iC < load.size(); ++iC) {
    _blockCommunicators[iC].copy();
  }
#endif
}


}

#endif
