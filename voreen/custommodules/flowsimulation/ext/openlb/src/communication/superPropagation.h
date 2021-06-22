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

#ifndef SUPER_PROPAGATION_H
#define SUPER_PROPAGATION_H

#include <map>
#include <vector>

#include "blockPropagation.h"
#include "utilities/aliases.h"
#include "communication/loadBalancer.h"

namespace olb {


#ifdef PARALLEL_MODE_MPI

template <typename T>
class SuperCommunicationTagCoordinator {
private:
  LoadBalancer<T>& _loadBalancer;

  class ChannelId;

  std::map<int,std::map<ChannelId,int>> _tags;
  bool _initialized = false;

public:
  SuperCommunicationTagCoordinator(LoadBalancer<T>& loadBalancer);

  template <typename DESCRIPTOR>
  void init(std::vector<BlockStaticPopulationPropagationNeighborhood<T,DESCRIPTOR>>& neighborhood);

  int get(int iC, int jC);

};

#endif // PARALLEL_MODE_MPI


template <typename T, typename DESCRIPTOR>
class SuperPropagationCommunicator {
private:
  SuperLattice<T,DESCRIPTOR>& _superLattice;

#ifdef PARALLEL_MODE_MPI
  SuperCommunicationTagCoordinator<T> _tagCoordinator;
  MPI_Comm _mpi_propagation_comm;
#endif

  std::vector<BlockStaticPopulationPropagationNeighborhood<T,DESCRIPTOR>> _blockNeighborhoods;
  std::vector<BlockStaticPopulationCommunicator<T,DESCRIPTOR>>            _blockCommunicators;

  bool _initialized = false;

public:
  SuperPropagationCommunicator(SuperLattice<T,DESCRIPTOR>& superLattice);

  void init();
  void propagate();

};


}

#endif
