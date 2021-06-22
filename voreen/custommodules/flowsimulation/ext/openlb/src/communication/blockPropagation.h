/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2020 Adrian Kummerlaender, Mathias J. Krause
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

#ifndef BLOCK_PROPAGATION_H
#define BLOCK_PROPAGATION_H

#include "communication/mpiManager.h"
#include "communication/loadBalancer.h"
#include "utilities/aliases.h"

#include <map>
#include <memory>


/// Disable to switch back to old block propagation
/**
 * The interim implementation improves performance but is still work in progress.
 **/
//#define NEW_INTERIM_BLOCK_PROPAGATION

namespace olb {

template<typename T> class SuperCommunicationTagCoordinator;


template<typename T, typename DESCRIPTOR>
class BlockStaticPopulationPropagationNeighborhood {
private:
  LoadBalancer<T>& _loadBalancer;

  const int _iC;

  const int _overlap;
  const int _extendedOverlap;

  std::map<int, std::vector<int>> _localInCells;
  std::map<int, std::vector<int>> _localOutCells;

  std::map<int, std::vector<int>> _remoteOutCells;

#ifdef PARALLEL_MODE_MPI
  std::map<int, singleton::MpiNonBlockingHelper> _helpers;
#endif

public:
  BlockStaticPopulationPropagationNeighborhood(
    SuperLattice2D<T,DESCRIPTOR>& superLattice,
    int iC,
    int overlap
  );

  BlockStaticPopulationPropagationNeighborhood(
    SuperLattice3D<T,DESCRIPTOR>& superLattice,
    int iC,
    int overlap
  );

#ifdef PARALLEL_MODE_MPI
  void send();
  void receive();
  void wait();
#endif

  const std::map<int, std::vector<int>>& getOutCells() const;
  const std::map<int, std::vector<int>>& getInCells() const;
  const std::map<int, std::vector<int>>& getRemoteOutCells() const;

};


template<typename T, typename DESCRIPTOR>
class BlockStaticPopulationCommunicator {
private:
  const int _iC;
#ifdef PARALLEL_MODE_MPI
  MPI_Comm& _comm;
#endif

  class CopyTask;
#ifdef PARALLEL_MODE_MPI
  class SendTask;
  class RecvTask;
#endif

  std::vector<CopyTask> _copyTasks;
#ifdef PARALLEL_MODE_MPI
  std::vector<SendTask> _sendTasks;
  std::vector<RecvTask> _recvTasks;
#endif

public:
  BlockStaticPopulationCommunicator(
    SuperLattice<T,DESCRIPTOR>& superLattice,
#ifdef PARALLEL_MODE_MPI
    SuperCommunicationTagCoordinator<T>& tagCoordinator,
    MPI_Comm& comm,
#endif
    int iC,
    const std::map<int, std::vector<int>>& outCells,
    const std::map<int, std::vector<int>>& inCells,
    const std::map<int, std::vector<int>>& remoteOutCells
  );

  BlockStaticPopulationCommunicator(const BlockStaticPopulationCommunicator& rhs) = delete;
  BlockStaticPopulationCommunicator(BlockStaticPopulationCommunicator&& rhs) = default;

  void copy();

#ifdef PARALLEL_MODE_MPI
  void receive();
  void send();
  void unpack();
  void wait();
#endif

};


}

#endif
