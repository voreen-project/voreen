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

#ifndef PLATFORM_GPU_CUDA_COMMUNICATOR_H
#define PLATFORM_GPU_CUDA_COMMUNICATOR_H

#include <map>
#include <set>
#include <typeindex>

#include "communication/blockCommunicator.h"

namespace olb {

template <typename T, typename DESCRIPTOR>
class SuperLattice;

template <typename T, typename DESCRIPTOR>
class ConcreteBlockCommunicator<ConcreteBlockLattice<T,DESCRIPTOR,Platform::GPU_CUDA>>
  final : public BlockCommunicator {
private:
  const int _iC;
#ifdef PARALLEL_MODE_MPI
  MPI_Comm _mpiCommunicator;
#endif

  class CopyTask;
#ifdef PARALLEL_MODE_MPI
  class SendTask;
  class RecvTask;
#endif

  std::vector<std::unique_ptr<CopyTask>> _copyTasks;
#ifdef PARALLEL_MODE_MPI
  std::vector<std::unique_ptr<SendTask>> _sendTasks;
  std::vector<std::unique_ptr<RecvTask>> _recvTasks;
#endif

public:
  ConcreteBlockCommunicator(SuperLattice<T,DESCRIPTOR>& super,
                            LoadBalancer<T>& loadBalancer,
#ifdef PARALLEL_MODE_MPI
                            SuperCommunicationTagCoordinator<T>& tagCoordinator,
                            MPI_Comm comm,
#endif
                            int iC,
                            const BlockCommunicationNeighborhood<T,DESCRIPTOR::d>& neighborhood);

  void copy() override;

#ifdef PARALLEL_MODE_MPI
  void receive() override;
  void send() override;
  void unpack() override;
#endif

  void wait() override;

};


}

#endif

#include "communicator.hh"
