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

#ifndef BLOCK_PROPAGATION_HH
#define BLOCK_PROPAGATION_HH

#include <map>
#include <set>

#include "blockPropagation.h"

#include "core/superLattice2D.h"
#include "core/superLattice3D.h"

namespace olb {


template <typename T, typename DESCRIPTOR>
BlockStaticPopulationPropagationNeighborhood<T,DESCRIPTOR>::BlockStaticPopulationPropagationNeighborhood(
  SuperLattice2D<T,DESCRIPTOR>& superLattice,
  int iC,
  int overlap
):
  _loadBalancer(superLattice.getLoadBalancer()),
  _iC(iC),
  _overlap(overlap),
  _extendedOverlap(superLattice.getOverlap())
{
  auto& cuboidGeometry = superLattice.getCuboidGeometry();

  const int nX  = cuboidGeometry.get(_iC).getNx();
  const int nY  = cuboidGeometry.get(_iC).getNy();

  BlockStructure2D extendedBlock(nX, nY, _extendedOverlap);

  for (int iX=0; iX < nX+2*_overlap; ++iX) {
    for (int iY=0; iY < nY+2*_overlap; ++iY) {
      if (iX < _overlap || iX > nX + _overlap - 1 ||
          iY < _overlap || iY > nY + _overlap - 1) {
        int latticeR[3];
        latticeR[0] = _iC;
        latticeR[1] = iX - _overlap;
        latticeR[2] = iY - _overlap;

        T physR[2];
        cuboidGeometry.getPhysR(physR, latticeR);

        int remoteLatticeR[3];
        if (cuboidGeometry.getLatticeR(remoteLatticeR, physR)) {
          _localInCells[remoteLatticeR[0]].emplace_back(extendedBlock.getCellId(
            latticeR[1] + _extendedOverlap,
            latticeR[2] + _extendedOverlap
          ));

          Cuboid2D<T>& remoteCuboid = cuboidGeometry.get(remoteLatticeR[0]);
          const int remoteNx = remoteCuboid.getNx();
          const int remoteNy = remoteCuboid.getNy();
          BlockStructure2D remoteExtendedBlock(remoteNx, remoteNy, _extendedOverlap);

          _remoteOutCells[remoteLatticeR[0]].emplace_back(remoteExtendedBlock.getCellId(
            remoteLatticeR[1] + _extendedOverlap,
            remoteLatticeR[2] + _extendedOverlap
          ));

#ifdef PARALLEL_MODE_MPI
          _helpers[remoteLatticeR[0]].allocate(1);
#endif
        }
      }
    }
  }
}

template <typename T, typename DESCRIPTOR>
BlockStaticPopulationPropagationNeighborhood<T,DESCRIPTOR>::BlockStaticPopulationPropagationNeighborhood(
  SuperLattice3D<T,DESCRIPTOR>& superLattice,
  int iC,
  int overlap
):
  _loadBalancer(superLattice.getLoadBalancer()),
  _iC(iC),
  _overlap(overlap),
  _extendedOverlap(superLattice.getOverlap())
{
  auto& cuboidGeometry = superLattice.getCuboidGeometry();  

  const int nX  = cuboidGeometry.get(_iC).getNx();
  const int nY  = cuboidGeometry.get(_iC).getNy();
  const int nZ  = cuboidGeometry.get(_iC).getNz();

  BlockStructure3D extendedBlock(nX, nY, nZ, _extendedOverlap);

  for (int iX=0; iX < nX+2*_overlap; ++iX) {
    for (int iY=0; iY < nY+2*_overlap; ++iY) {
      for (int iZ=0; iZ < nZ+2*_overlap; ++iZ) {
        if (iX < _overlap || iX > nX + _overlap - 1 ||
            iY < _overlap || iY > nY + _overlap - 1 ||
            iZ < _overlap || iZ > nZ + _overlap - 1) {
          int latticeR[4];
          latticeR[0] = _iC;
          latticeR[1] = iX - _overlap;
          latticeR[2] = iY - _overlap;
          latticeR[3] = iZ - _overlap;

          T physR[3];
          cuboidGeometry.getPhysR(physR, latticeR);

          int remoteLatticeR[4];
          if (cuboidGeometry.getLatticeR(remoteLatticeR, physR)) {
            _localInCells[remoteLatticeR[0]].emplace_back(extendedBlock.getCellId(
              latticeR[1] + _extendedOverlap,
              latticeR[2] + _extendedOverlap,
              latticeR[3] + _extendedOverlap
            ));

            Cuboid3D<T>& remoteCuboid = cuboidGeometry.get(remoteLatticeR[0]);
            const int remoteNx = remoteCuboid.getNx();
            const int remoteNy = remoteCuboid.getNy();
            const int remoteNz = remoteCuboid.getNz();
            BlockStructure3D remoteExtendedBlock(remoteNx, remoteNy, remoteNz, _extendedOverlap);

            _remoteOutCells[remoteLatticeR[0]].emplace_back(remoteExtendedBlock.getCellId(
              remoteLatticeR[1] + _extendedOverlap,
              remoteLatticeR[2] + _extendedOverlap,
              remoteLatticeR[3] + _extendedOverlap
            ));

#ifdef PARALLEL_MODE_MPI
            _helpers[remoteLatticeR[0]].allocate(1);
#endif
          }
        }
      }
    }
  }
}

#ifdef PARALLEL_MODE_MPI

template <typename T, typename DESCRIPTOR>
void BlockStaticPopulationPropagationNeighborhood<T,DESCRIPTOR>::send()
{
  for (auto& partner : _remoteOutCells) {
    const int         iC    = partner.first;
    std::vector<int>& cells = partner.second;
    singleton::mpi().iSend(
      cells.data(),
      cells.size(),
      _loadBalancer.rank(iC),
      _helpers[iC].get_mpiRequest(0),
      _iC);
  }
}

template <typename T, typename DESCRIPTOR>
void BlockStaticPopulationPropagationNeighborhood<T,DESCRIPTOR>::receive()
{
  for (const auto& partner : _localInCells) {
    const int iC = partner.first;
    _localOutCells[iC].resize(
      singleton::mpi().probeReceiveSize(_loadBalancer.rank(iC), MPI_INT, iC));

    singleton::mpi().receive(
      _localOutCells[iC].data(),
      _localOutCells[iC].size(),
      _loadBalancer.rank(iC),
      iC);
  }
}

template <typename T, typename DESCRIPTOR>
void BlockStaticPopulationPropagationNeighborhood<T,DESCRIPTOR>::wait()
{
  for (const auto& partner : _remoteOutCells) {
    const int iC = partner.first;
    _helpers[iC].wait(0);
  }
}

#endif // PARALLEL_MODE_MPI

template <typename T, typename DESCRIPTOR>
const std::map<int, std::vector<int>>&
BlockStaticPopulationPropagationNeighborhood<T,DESCRIPTOR>::getOutCells() const
{
  return _localOutCells;
}

template <typename T, typename DESCRIPTOR>
const std::map<int, std::vector<int>>&
BlockStaticPopulationPropagationNeighborhood<T,DESCRIPTOR>::getInCells() const
{
  return _localInCells;
}

template <typename T, typename DESCRIPTOR>
const std::map<int, std::vector<int>>&
BlockStaticPopulationPropagationNeighborhood<T,DESCRIPTOR>::getRemoteOutCells() const
{
  return _remoteOutCells;
}


/// Wrapper for a local plain-copy block propagation request
template<typename T, typename DESCRIPTOR>
class BlockStaticPopulationCommunicator<T,DESCRIPTOR>::CopyTask {
private:
  const std::vector<int>& _targetCells;
  const std::vector<int>& _sourceCells;

  BlockStaticPopulationD<T,DESCRIPTOR>& _targetPopulation;
  const BlockStaticPopulationD<T,DESCRIPTOR>& _sourcePopulation;

public:
  CopyTask(
    const std::vector<int>& targetCells,       BlockStaticPopulationD<T,DESCRIPTOR>& targetPopulation,
    const std::vector<int>& sourceCells, const BlockStaticPopulationD<T,DESCRIPTOR>& sourcePopulation):
    _targetCells(targetCells),
    _sourceCells(sourceCells),
    _targetPopulation(targetPopulation),
    _sourcePopulation(sourcePopulation)
  {
    OLB_ASSERT(_sourceCells.size() == _targetCells.size(), "Source cell count must match target cell count");
  }

  void copy()
  {
    for (unsigned iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      for (unsigned i=0; i < _sourceCells.size(); ++i) {
        T* dst = _targetPopulation.getPopulationPointer(iPop, _targetCells[i]);
        const T* src = _sourcePopulation.getPopulationPointer(iPop, _sourceCells[i]);
        *dst = *src;
      }
    }
  };
};

#ifdef PARALLEL_MODE_MPI

/// Wrapper for a non-blocking block propagation send request
template<typename T, typename DESCRIPTOR>
class BlockStaticPopulationCommunicator<T,DESCRIPTOR>::SendTask {
private:
  MPI_Comm& _comm;

  const int _tag;
  const int _rank;
  const std::vector<int>& _cells;

  BlockStaticPopulationD<T,DESCRIPTOR>& _population;

  std::unique_ptr<T[]> _buffer;
  singleton::MpiNonBlockingHelper _helper;

public:
  SendTask(MPI_Comm& comm, int tag, int rank, const std::vector<int>& cells, BlockStaticPopulationD<T,DESCRIPTOR>& population):
    _comm(comm),
    _tag(tag),
    _rank(rank),
    _cells(cells),
    _population(population),
    _buffer(new T[_cells.size() * DESCRIPTOR::q])
  {
    _helper.allocate(1);

    singleton::mpi().sendInit(
      _buffer.get(),
      _cells.size() * DESCRIPTOR::q,
      _rank,
      _helper.get_mpiRequest(0),
      _tag,
      _comm);
  }

  void send()
  {
    T* bufferPop = _buffer.get();
    for (unsigned iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      for (int iCell : _cells) {
        *(bufferPop++) = *_population.getPopulationPointer(iPop, iCell);
      }
    }

    _helper.start(0);
  };

  void wait()
  {
    _helper.wait(0);
  }
};

/// Wrapper for a non-blocking block propagation receive request
template<typename T, typename DESCRIPTOR>
class BlockStaticPopulationCommunicator<T,DESCRIPTOR>::RecvTask {
private:
  MPI_Comm& _comm;

  const int _tag;
  const int _rank;
  const std::vector<int>& _cells;

  BlockStaticPopulationD<T,DESCRIPTOR>& _population;

  std::unique_ptr<T[]> _buffer;
  singleton::MpiNonBlockingHelper _helper;

public:
  /// Manual replacement for std::reference_wrapper<RecvTask>
  /**
   * Used to track pending receive requests in std::set.
   *
   * This is a workaround for problematic external definition of
   * dependently-typed comparision operators for nested classes.
   * Reconsider as soon as depending on C++17 is allowed.
   **/
  class ref {
  private:
    RecvTask& _task;
  public:
    ref(RecvTask& task): _task(task) { };

    RecvTask* operator->() const
    {
      return &_task;
    }

    bool operator <(const ref& rhs) const
    {
      return _task < rhs._task;
    }
  };

  RecvTask(MPI_Comm& comm, int tag, int rank, const std::vector<int>& cells, BlockStaticPopulationD<T,DESCRIPTOR>& population):
    _comm(comm),
    _tag(tag),
    _rank(rank),
    _cells(cells),
    _population(population),
    _buffer(new T[_cells.size() * DESCRIPTOR::q])
  {
    _helper.allocate(1);

    singleton::mpi().recvInit(
      _buffer.get(),
      _cells.size() * DESCRIPTOR::q,
      rank,
      _helper.get_mpiRequest(0),
      _tag,
      _comm);
  }

  bool operator<(const RecvTask& rhs) const
  {
    return  _rank  < rhs._rank
        || (_rank == rhs._rank && _tag < rhs._tag);
  }

  void receive()
  {
    _helper.start(0);
  };

  bool isDone()
  {
    return _helper.isDone(0);
  }

  void unpack()
  {
    T* bufferPop = _buffer.get();
    for (unsigned iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      for (int iCell : _cells) {
        *_population.getPopulationPointer(iPop, iCell) = *(bufferPop++);
      }
    }
  }
};

#endif // PARALLEL_MODE_MPI

template <typename T, typename DESCRIPTOR>
BlockStaticPopulationCommunicator<T,DESCRIPTOR>::BlockStaticPopulationCommunicator(
  SuperLattice<T,DESCRIPTOR>& superLattice,
#ifdef PARALLEL_MODE_MPI
  SuperCommunicationTagCoordinator<T>& tagCoordinator,
  MPI_Comm& comm,
#endif
  int iC,
  const std::map<int, std::vector<int>>& outCells,
  const std::map<int, std::vector<int>>& inCells,
  const std::map<int, std::vector<int>>& remoteOutCells):
  _iC(iC)
#ifdef PARALLEL_MODE_MPI
  , _comm(comm)
#endif
{
  auto& loadBalancer = superLattice.getLoadBalancer();
  auto& population   = superLattice.getExtendedBlockLattice(loadBalancer.loc(_iC)).getStaticPopulationD();

#ifdef PARALLEL_MODE_MPI
  for (const auto& partner : outCells) {
    const int               iC    = partner.first;
    const std::vector<int>& cells = partner.second;
    if (!cells.empty() && !loadBalancer.isLocal(iC)) {
      _sendTasks.emplace_back(
        _comm, tagCoordinator.get(iC, _iC),
        loadBalancer.rank(iC), cells, population);
    }
  }
#endif // PARALLEL_MODE_MPI

  for (const auto& partner : inCells) {
    const int               iC    = partner.first;
    const std::vector<int>& cells = partner.second;
    if (!cells.empty()) {
      if (loadBalancer.isLocal(iC)) {
        auto& remotePopulation = superLattice.getExtendedBlockLattice(loadBalancer.loc(iC)).getStaticPopulationD();
        _copyTasks.emplace_back(cells, population, remoteOutCells.at(iC), remotePopulation);
      }
      else {
#ifdef PARALLEL_MODE_MPI
        _recvTasks.emplace_back(
          _comm, tagCoordinator.get(iC, _iC),
          loadBalancer.rank(iC), cells, population);
#endif // PARALLEL_MODE_MPI
      }
    }
  }
}

template <typename T, typename DESCRIPTOR>
void BlockStaticPopulationCommunicator<T,DESCRIPTOR>::copy()
{
  for (auto& task : _copyTasks) {
    task.copy();
  }
}

#ifdef PARALLEL_MODE_MPI

template <typename T, typename DESCRIPTOR>
void BlockStaticPopulationCommunicator<T,DESCRIPTOR>::receive()
{
  for (auto& task : _recvTasks) {
    task.receive();
  }
}

template <typename T, typename DESCRIPTOR>
void BlockStaticPopulationCommunicator<T,DESCRIPTOR>::send()
{
  for (auto& task : _sendTasks) {
    task.send();
  }
}

template <typename T, typename DESCRIPTOR>
void BlockStaticPopulationCommunicator<T,DESCRIPTOR>::unpack()
{
  std::set<typename RecvTask::ref> pending(_recvTasks.begin(), _recvTasks.end());
  while (!pending.empty()) {
    auto task_iterator = pending.begin();
    while (task_iterator != pending.end()) {
      auto& task = *task_iterator;
      if (task->isDone()) {
        task->unpack();
        task_iterator = pending.erase(task_iterator);
      }
      else {
        ++task_iterator;
      }
    }
  }
}

template <typename T, typename DESCRIPTOR>
void BlockStaticPopulationCommunicator<T,DESCRIPTOR>::wait()
{
  for (auto& task : _sendTasks) {
    task.wait();
  }
}

#endif // PARALLEL_MODE_MPI


}

#endif
