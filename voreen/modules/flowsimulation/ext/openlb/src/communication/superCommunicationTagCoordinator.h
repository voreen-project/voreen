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

#ifndef SUPER_COMMUNICATION_TAG_COORDINATOR_H
#define SUPER_COMMUNICATION_TAG_COORDINATOR_H

#include <map>

#include "communication/loadBalancer.h"
#include "blockCommunicationNeighborhood.h"

namespace olb {

#ifdef PARALLEL_MODE_MPI

/// Communication-free negotation of unique tags for inter-cuboid communication
/**
 * The cuboid IDs themselves can not be used safely as the MPI standard
 * only supports tag sizes up to 32768.
 *
 * This class provides unique tags for the pair-wise communication of cuboids
 * between their owning ranks, supporting pair-wise neighborhood sizes up to
 * 32768 (depending on number of used groups) and "arbitrarily" large global
 * cuboid counts.
 **/
template <typename T>
class SuperCommunicationTagCoordinator {
private:
  LoadBalancer<T>& _loadBalancer;

  class ChannelId;

  std::map<int,std::map<ChannelId,int>> _tags;

public:
  SuperCommunicationTagCoordinator(LoadBalancer<T>& loadBalancer);

  /// Generate unique tags for given block neighborhood system
  template <unsigned D>
  void coordinate(std::vector<BlockCommunicationNeighborhood<T,D>>& neighborhood);

  /// Returns unique tag for communication between cuboids iC and jC
  /**
   * iGroup can be used as a namespace substitute. It simply offsets the tags
   * group-times by the pair-wise neighborhood size, correspondingly reducing
   * said neighborhood's maximum supported size.
   **/
  int get(int iC, int jC, int iGroup=0);

};

#endif // PARALLEL_MODE_MPI

}

#endif
