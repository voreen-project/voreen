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

#ifndef BLOCK_COMMUNICATION_NEIGHBORHOOD_HH
#define BLOCK_COMMUNICATION_NEIGHBORHOOD_HH

#include <stdexcept>
#include <algorithm>

#include "blockCommunicationNeighborhood.h"

#include "utilities/aliases.h"

namespace olb {

// *INDENT-OFF*

template <typename T, unsigned D>
BlockCommunicationNeighborhood<T,D>::BlockCommunicationNeighborhood(
    CuboidGeometry<T,D>& cuboidGeometry
  , LoadBalancer<T>& loadBalancer
  , int iC
  , int padding
#ifdef PARALLEL_MODE_MPI
  , MPI_Comm comm
#endif
):
    _cuboidGeometry(cuboidGeometry)
  , _loadBalancer(loadBalancer)
  , _iC(iC)
  , _padding(padding)
#ifdef PARALLEL_MODE_MPI
  , _neighborhoodComm(comm)
#endif
{
  // Ensure that any neighboring cuboids are interacted with during request negotiation
  requestOverlap(1);
  for (auto& [_, cells] : _cellsInboundFrom) {
    cells.clear();
  }
  for (auto& [_, cells] : _cellsRequestedFrom) {
    cells.clear();
  }
}

template <typename T, unsigned D>
void BlockCommunicationNeighborhood<T,D>::requestCell(LatticeR<D> latticeR)
{
  auto& cuboid = _cuboidGeometry.get(_iC);
  BlockStructureD<D> paddedBlock(cuboid.getExtent(), _padding);
  if (paddedBlock.isPadding(latticeR)) {
    T physR[D];
    // Read physR using global cuboid geometry instead of local cuboid
    // to resolve periodic boundaries
    _cuboidGeometry.getPhysR(physR, latticeR.withPrefix(_iC));
    int remoteLatticeR[D+1];
    if (_cuboidGeometry.getLatticeR(remoteLatticeR, physR)) {
      _cellsInboundFrom[remoteLatticeR[0]].emplace_back(paddedBlock.getCellId(latticeR));

      Cuboid<T,D>& remoteCuboid = _cuboidGeometry.get(remoteLatticeR[0]);
      BlockStructureD<D> remotePaddedBlock(remoteCuboid.getExtent(), _padding);

      _cellsRequestedFrom[remoteLatticeR[0]].emplace_back(
        remotePaddedBlock.getCellId(remoteLatticeR+1));
    }
  } else {
    throw std::logic_error("Requested cell is outside of available padding");
  }
}

template <typename T, unsigned D>
void BlockCommunicationNeighborhood<T,D>::requestOverlap(int width)
{
  if (width < 1) {
    throw std::logic_error("Overlap requests must have width >= 1");
  }
  if (width > _padding) {
    throw std::logic_error("Requested overlap exceeds available padding");
  }

  auto& cuboid = _cuboidGeometry.get(_iC);
  BlockStructureD<D> overlapBlock(cuboid.getExtent(), width);

  overlapBlock.forSpatialLocations([&](LatticeR<D> latticeR) {
    if (overlapBlock.isPadding(latticeR)) {
      requestCell(latticeR);
    }
  });
}

template <typename T, unsigned D>
void BlockCommunicationNeighborhood<T,D>::setFieldAvailability(std::type_index field, bool available)
{
  auto iter = std::find(_fieldsRequested.begin(), _fieldsRequested.end(), field);
  if (iter != _fieldsRequested.end()) {
    _fieldsAvailable[iter - _fieldsRequested.begin()] = available;
  }
}

template <typename T, unsigned D>
template <typename BLOCK>
void BlockCommunicationNeighborhood<T,D>::setFieldsAvailability(int iC, BLOCK& block)
{
  auto& fieldsCommonWith = _fieldsCommonWith[iC];
  fieldsCommonWith.clear();
  for (unsigned iField=0; iField < _fieldsRequested.size(); ++iField) {
    auto field = _fieldsRequested[iField];
    if (block.hasCommunicatable(field) && _fieldsAvailable[iField]) {
      fieldsCommonWith.emplace_back(field);
    }
  }
}

template <typename T, unsigned D>
void BlockCommunicationNeighborhood<T,D>::maintain()
{
  for (auto& [iC, _] : _cellsInboundFrom) {
    auto& iCells = _cellsInboundFrom[iC];
    auto& oCells = _cellsRequestedFrom[iC];
    // Compute sorted indices w.r.t. cell IDs in iCells
    std::vector<std::ptrdiff_t> p(iCells.size());
    std::iota(p.begin(), p.end(), 0);
    std::sort(p.begin(), p.end(), [&iCells](auto i, auto j) { return iCells[i] < iCells[j]; });
    // Reduce indices to only include unique values of iCells
    auto pu = std::unique(p.begin(), p.end(), [&iCells](auto i, auto j) { return iCells[i] == iCells[j]; });
    p.erase(pu, p.end());
    // Copy cell IDs at computed indices to new in and out lists
    {
      std::vector<CellID> buffer(p.size());
      std::transform(p.begin(), p.end(), buffer.begin(), [&iCells](auto i) { return iCells[i]; });
      iCells = buffer;
    }
    {
      std::vector<CellID> buffer(p.size());
      std::transform(p.begin(), p.end(), buffer.begin(), [&oCells](auto i) { return oCells[i]; });
      oCells = buffer;
    }
  }
}

#ifdef PARALLEL_MODE_MPI

template <typename T, unsigned D>
void BlockCommunicationNeighborhood<T,D>::send(SuperCommunicationTagCoordinator<T>& coordinator)
{
  for (auto& [iC, cells] : _cellsRequestedFrom) {
    if (!_loadBalancer.isLocal(iC)) {
      _fieldRequests[iC] = std::make_unique<MpiSendRequest>(
        _fieldsAvailable.get(), _fieldsRequested.size(),
        _loadBalancer.rank(iC), coordinator.get(iC, _iC, 0), _neighborhoodComm);
      _fieldRequests[iC]->start();

      _cellsRequests[iC] = std::make_unique<MpiSendRequest>(
        cells.data(), cells.size(),
        _loadBalancer.rank(iC), coordinator.get(iC, _iC, 1), _neighborhoodComm);
      _cellsRequests[iC]->start();
    }
  }
}

template <typename T, unsigned D>
void BlockCommunicationNeighborhood<T,D>::receive(SuperCommunicationTagCoordinator<T>& coordinator)
{
  forRemoteNeighbors([&](int iC) {
    std::unique_ptr<bool[]> fieldsAvailable(new bool[_fieldsRequested.size()] { });
    auto& fieldsCommonWith = _fieldsCommonWith[iC];
    fieldsCommonWith.clear();

    singleton::mpi().receive(
      fieldsAvailable.get(), _fieldsRequested.size(),
      _loadBalancer.rank(iC),
      coordinator.get(iC, _iC, 0),
      _neighborhoodComm);

    for (unsigned iField=0; iField < _fieldsRequested.size(); ++iField) {
      if (fieldsAvailable[iField] && _fieldsAvailable[iField]) {
        fieldsCommonWith.emplace_back(_fieldsRequested[iField]);
      }
    }

    std::size_t newOutboundCount = singleton::mpi().probeReceiveSize(_loadBalancer.rank(iC), MPI_UNSIGNED,
                                                                     coordinator.get(iC, _iC, 1),
                                                                     _neighborhoodComm);
    _cellsOutboundTo[iC].resize(newOutboundCount);

    singleton::mpi().receive(
      _cellsOutboundTo[iC].data(),
      _cellsOutboundTo[iC].size(),
      _loadBalancer.rank(iC),
      coordinator.get(iC, _iC, 1),
      _neighborhoodComm);
  });
}

template <typename T, unsigned D>
void BlockCommunicationNeighborhood<T,D>::wait()
{
  forRemoteNeighbors([&](int iC) {
    _fieldRequests[iC]->wait();
    _cellsRequests[iC]->wait();
  });
}

#endif // PARALLEL_MODE_MPI

template <typename T, unsigned D>
const std::vector<std::type_index>&
BlockCommunicationNeighborhood<T,D>::getFieldsCommonWith(int iC) const
{
  return _fieldsCommonWith.at(iC);
}

template <typename T, unsigned D>
const std::vector<CellID>&
BlockCommunicationNeighborhood<T,D>::getCellsOutboundTo(int iC) const
{
  return _cellsOutboundTo.at(iC);
}

template <typename T, unsigned D>
const std::vector<CellID>&
BlockCommunicationNeighborhood<T,D>::getCellsInboundFrom(int iC) const
{
  return _cellsInboundFrom.at(iC);
}

template <typename T, unsigned D>
const std::vector<CellID>&
BlockCommunicationNeighborhood<T,D>::getCellsRequestedFrom(int iC) const
{
  return _cellsRequestedFrom.at(iC);
}

// *INDENT-ON*

}

#endif
