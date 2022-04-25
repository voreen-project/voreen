/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Adrian Kummerlaender
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

#ifndef MPI_REQUEST_H
#define MPI_REQUEST_H
#ifdef PARALLEL_MODE_MPI

#include "mpi.h"
#include "mpiManager.h"

namespace olb {

/// Basic wrapper around a single MPI_Request
class MpiRequest {
protected:
  MPI_Request _request;
  MPI_Status  _status;

public:
  MpiRequest():
    _request{}, _status{} { };

  void start() {
    MPI_Start(&_request);
  }

  void wait() {
    MPI_Wait(&_request, &_status);
  }

  bool isDone() {
    int done;
    MPI_Test(&_request, &done, MPI_STATUS_IGNORE);
    return done;
  }
};

/// Non-blocking MPI send request
class MpiSendRequest : public MpiRequest {
public:
  template <typename T>
  MpiSendRequest(T* buffer, std::size_t size,
                 int rank, int tag, MPI_Comm communicator)
  {
    singleton::mpi().sendInit(
      buffer, size,
      rank,
      &this->_request,
      tag,
      communicator);
  }

  ~MpiSendRequest() {
    //MPI_Request_free(&this->_request);
  }
};

/// Non-blocking MPI receive request
class MpiRecvRequest : public MpiRequest {
public:
  template <typename T>
  MpiRecvRequest(T* buffer, std::size_t size,
                 int rank, int tag, MPI_Comm communicator)
  {
    singleton::mpi().recvInit(
      buffer, size,
      rank,
      &this->_request,
      tag,
      communicator);
  }

  ~MpiRecvRequest() {
    //MPI_Request_free(&this->_request);
  }
};

}

#endif
#endif
