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

#ifndef BLOCK_COMMUNICATION_NEIGHBORHOOD_H
#define BLOCK_COMMUNICATION_NEIGHBORHOOD_H

#include "core/blockStructure.h"
#include "communication/mpiRequest.h"
#include "communication/loadBalancer.h"
#include "communication/communicatable.h"
#include "utilities/aliases.h"

#include <map>
#include <memory>

namespace olb {

// *INDENT-OFF*

template<typename T> class SuperCommunicationTagCoordinator;

/// Configurable overlap communication neighborhood of a block
/**
 * Managed by SuperCommunicator
 **/
template<typename T, unsigned D>
class BlockCommunicationNeighborhood {
private:
  CuboidGeometry<T,D>& _cuboidGeometry;
  LoadBalancer<T>& _loadBalancer;

  const int _iC;
  const int _padding;

  std::vector<std::type_index> _fieldsRequested;
  std::unique_ptr<bool[]>      _fieldsAvailable;
  std::map<int, std::vector<std::type_index>> _fieldsCommonWith;

  std::map<int, std::vector<CellID>> _cellsInboundFrom;
  std::map<int, std::vector<CellID>> _cellsOutboundTo;
  std::map<int, std::vector<CellID>> _cellsRequestedFrom;

#ifdef PARALLEL_MODE_MPI
  MPI_Comm _neighborhoodComm;

  std::map<int, std::unique_ptr<MpiSendRequest>> _fieldRequests;
  std::map<int, std::unique_ptr<MpiSendRequest>> _cellsRequests;
#endif

public:
  BlockCommunicationNeighborhood(  CuboidGeometry<T,D>& cuboidGeometry
                                 , LoadBalancer<T>& loadBalancer
                                 , int iC
                                 , int padding
#ifdef PARALLEL_MODE_MPI
                                 , MPI_Comm comm
#endif
                                 );

  /// Request field and provides local availability
  template <typename FIELD>
  void requestField() {
    if (std::find(_fieldsRequested.begin(), _fieldsRequested.end(), typeid(FIELD)) == _fieldsRequested.end()) {
      _fieldsRequested.emplace_back(typeid(FIELD));
      _fieldsAvailable.reset(static_cast<bool*>(
        std::realloc(static_cast<void*>(_fieldsAvailable.release()),
                     _fieldsRequested.size() * sizeof(bool))));
    }
  }
  /// Request individual cell for communication
  void requestCell(LatticeR<D> latticeR);
  /// Request all cells in overlap of size width for communication
  void requestOverlap(int width);

  /// Update local availability of previously requested field
  void setFieldAvailability(std::type_index field, bool available);
  /// Update outbound availabilities for locally available neighbor block
  /**
   * Used for local non-MPI exchange of field availabilities
   **/
  template <typename BLOCK>
  void setFieldsAvailability(int iC, BLOCK& block);

  void maintain();

#ifdef PARALLEL_MODE_MPI
  void send(SuperCommunicationTagCoordinator<T>&);
  void receive(SuperCommunicationTagCoordinator<T>&);
  void wait();

  /// Calls f(iC) for every non-local cuboid ID iC
  template <typename F>
  void forRemoteNeighbors(F f) const {
    for (const auto& [iC, _] : _cellsRequestedFrom) {
      if (!_loadBalancer.isLocal(iC)) {
        f(iC);
      }
    }
  }
#endif

  /// Calls f(iC) for every local cuboid ID iC
  template <typename F>
  void forLocalNeighbors(F f) const {
    for (const auto& [iC, _] : _cellsRequestedFrom) {
      if (_loadBalancer.isLocal(iC)) {
        f(iC);
      }
    }
  }

  const std::vector<std::type_index>& getFieldsCommonWith(int iC) const;

  const std::vector<CellID>& getCellsOutboundTo(int iC) const;
  const std::vector<CellID>& getCellsInboundFrom(int iC) const;
  const std::vector<CellID>& getCellsRequestedFrom(int iC) const;

};

// *INDENT-ON*

}

#endif
