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

#ifndef BLOCK_COMMUNICATOR_HH
#define BLOCK_COMMUNICATOR_HH

#include "blockCommunicator.h"

#include "mpiRequest.h"
#include "communicatable.h"
#include "superCommunicationTagCoordinator.h"

namespace olb {

/// Wrapper for a local plain-copy block communication request
template <typename BLOCK>
class ConcreteBlockCommunicator<BLOCK>::CopyTask {
private:
  const std::vector<CellID>& _targetCells;
  const std::vector<CellID>& _sourceCells;

  MultiConcreteCommunicatable<BLOCK> _target;
  MultiConcreteCommunicatable<BLOCK> _source;

  std::unique_ptr<std::uint8_t[]> _buffer;

public:
  CopyTask(
    const std::vector<std::type_index>& fields,
    const std::vector<CellID>& targetCells, BLOCK& target,
    const std::vector<CellID>& sourceCells, BLOCK& source):
    _targetCells(targetCells),
    _sourceCells(sourceCells),
    _target(target, fields),
    _source(source, fields),
    _buffer(new std::uint8_t[_source.size(_sourceCells)] { })
  {
    OLB_ASSERT(_sourceCells.size() == _targetCells.size(),
               "Source cell count must match target cell count");
  }

  void copy()
  {
    _source.serialize(_sourceCells, _buffer.get());
    _target.deserialize(_targetCells, _buffer.get());
  };
};

#ifdef PARALLEL_MODE_MPI

/// Wrapper for a non-blocking block propagation send request
template <typename BLOCK>
class ConcreteBlockCommunicator<BLOCK>::SendTask {
private:
  const std::vector<CellID>& _cells;

  MultiConcreteCommunicatable<BLOCK> _source;

  std::unique_ptr<std::uint8_t[]> _buffer;
  MpiSendRequest _request;

public:
  SendTask(MPI_Comm comm, int tag, int rank,
           const std::vector<std::type_index>& fields,
           const std::vector<CellID>& cells,
           BLOCK& block):
    _cells(cells),
    _source(block, fields),
    _buffer(new std::uint8_t[_source.size(_cells)] { }),
    _request(_buffer.get(), _source.size(_cells),
             rank, tag, comm)
  { }

  void send()
  {
    _source.serialize(_cells, _buffer.get());
    _request.start();
  }

  void wait()
  {
    _request.wait();
  }
};

/// Wrapper for a non-blocking block propagation receive request
template <typename BLOCK>
class ConcreteBlockCommunicator<BLOCK>::RecvTask {
private:
  const int _tag;
  const int _rank;
  const std::vector<CellID>& _cells;

  MultiConcreteCommunicatable<BLOCK> _target;

  std::unique_ptr<std::uint8_t[]> _buffer;
  MpiRecvRequest _request;

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

  RecvTask(MPI_Comm comm, int tag, int rank,
           const std::vector<std::type_index>& fields,
           const std::vector<CellID>& cells,
           BLOCK& block):
    _tag(tag),
    _rank(rank),
    _cells(cells),
    _target(block, fields),
    _buffer(new std::uint8_t[_target.size(_cells)] { }),
    _request(_buffer.get(), _target.size(_cells),
             _rank, _tag, comm)
  { }

  bool operator<(const RecvTask& rhs) const
  {
    return  _rank  < rhs._rank
        || (_rank == rhs._rank && _tag < rhs._tag);
  }

  void receive()
  {
    _request.start();
  };

  bool isDone()
  {
    return _request.isDone();
  }

  void unpack()
  {
    _target.deserialize(_cells, _buffer.get());
  }
};

#endif // PARALLEL_MODE_MPI

template <typename BLOCK>
template <typename T, typename SUPER>
ConcreteBlockCommunicator<BLOCK>::ConcreteBlockCommunicator(
  SUPER& super,
  LoadBalancer<T>& loadBalancer,
#ifdef PARALLEL_MODE_MPI
  SuperCommunicationTagCoordinator<T>& tagCoordinator,
  MPI_Comm comm,
#endif
  int iC,
  const BlockCommunicationNeighborhood<T,SUPER::d>& neighborhood):
  _iC(iC)
#ifdef PARALLEL_MODE_MPI
, _mpiCommunicator(comm)
#endif
{
#ifdef PARALLEL_MODE_MPI
  neighborhood.forRemoteNeighbors([&](int remoteC) {
    if (!neighborhood.getCellsOutboundTo(remoteC).empty()) {
      _sendTasks.emplace_back(_mpiCommunicator, tagCoordinator.get(remoteC, loadBalancer.glob(_iC)),
                              loadBalancer.rank(remoteC),
                              neighborhood.getFieldsCommonWith(remoteC),
                              neighborhood.getCellsOutboundTo(remoteC),
                              super.template getBlock<BLOCK>(_iC));
    }
    if (!neighborhood.getCellsInboundFrom(remoteC).empty()) {
      _recvTasks.emplace_back(_mpiCommunicator, tagCoordinator.get(remoteC, loadBalancer.glob(_iC)),
                              loadBalancer.rank(remoteC),
                              neighborhood.getFieldsCommonWith(remoteC),
                              neighborhood.getCellsInboundFrom(remoteC),
                              super.template getBlock<BLOCK>(_iC));
    }
  });
#endif // PARALLEL_MODE_MPI

  neighborhood.forLocalNeighbors([&](int localC) {
    if (!neighborhood.getCellsInboundFrom(localC).empty()) {
      _copyTasks.emplace_back(neighborhood.getFieldsCommonWith(localC),
                              neighborhood.getCellsInboundFrom(localC),   super.template getBlock<BLOCK>(_iC),
                              neighborhood.getCellsRequestedFrom(localC), super.template getBlock<BLOCK>(loadBalancer.loc(localC)));
    }
  });
}

template <typename BLOCK>
void ConcreteBlockCommunicator<BLOCK>::copy()
{
  for (auto& task : _copyTasks) {
    task.copy();
  }
}

#ifdef PARALLEL_MODE_MPI

template <typename BLOCK>
void ConcreteBlockCommunicator<BLOCK>::receive()
{
  for (auto& task : _recvTasks) {
    task.receive();
  }
}

template <typename BLOCK>
void ConcreteBlockCommunicator<BLOCK>::send()
{
  for (auto& task : _sendTasks) {
    task.send();
  }
}

template <typename BLOCK>
void ConcreteBlockCommunicator<BLOCK>::unpack()
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

template <typename BLOCK>
void ConcreteBlockCommunicator<BLOCK>::wait()
{
  for (auto& task : _sendTasks) {
    task.wait();
  }
}

#else // PARALLEL_MODE_MPI

template <typename BLOCK>
void ConcreteBlockCommunicator<BLOCK>::wait()
{ }

#endif

}

#endif
