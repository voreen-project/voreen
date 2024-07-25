/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007 The OpenLB project
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

/** \file
 * Wrapper functions that simplify the use of MPI
 */

#ifndef MPI_MANAGER_H
#define MPI_MANAGER_H

#ifdef PARALLEL_MODE_MPI
#include "mpi.h"
#include <vector>
#include <memory>
#include "core/blockData.h"
#endif
#include <string>
#include "io/ostreamManager.h"
#include "utilities/aDiff.h"

namespace olb {

namespace singleton {

#ifdef PARALLEL_MODE_MPI

/// Helper class for non blocking MPI communication

class MpiNonBlockingHelper {
private:
  /// Size of the vector _mpiRequest/_mpiStatus
  unsigned _size;
  /// vector of MPI_Request
  std::unique_ptr<MPI_Request[]> _mpiRequest;
  /// vector of MPI_Status
  std::unique_ptr<MPI_Status[]> _mpiStatus;
public:
  MpiNonBlockingHelper();
  ~MpiNonBlockingHelper() = default;

  MpiNonBlockingHelper(MpiNonBlockingHelper&& rhs) = default;
  MpiNonBlockingHelper(const MpiNonBlockingHelper&) = delete;
  MpiNonBlockingHelper& operator=(const MpiNonBlockingHelper&) = delete;

  /// Allocates memory
  void allocate(unsigned i);
  /// Reset
  void free();

  /// Returns the size of the vector _mpiRequest/_mpiStatus
  unsigned get_size() const;

  /// Get the specified request object
  MPI_Request* get_mpiRequest(int i=0) const;
  /// Get the specified status object
  MPI_Status* get_mpiStatus(int i=0) const;

  void start(int i);
  void wait(int i);
  bool isDone(int i);

  /// Swap method
  void swap(MpiNonBlockingHelper& rhs);
};

/// Wrapper functions that simplify the use of MPI

class MpiManager {
public:
  /// Initializes the mpi manager
  void init(int *argc, char ***argv, bool verbose=true);
  /// Returns the number of processes
  int getSize() const;
  /// Returns the process ID
  int getRank() const;
  /// Returns process ID of main processor
  int bossId() const;
  /// Tells whether current processor is main processor
  bool isMainProcessor() const;
  /// Returns universal MPI-time in seconds
  double getTime() const;

  /// Synchronizes the processes
  void barrier(MPI_Comm comm = MPI_COMM_WORLD);

  /// Sends data at *buf, blocking
  template <typename T>
  void send(T *buf, int count, int dest, int tag = 0, MPI_Comm comm = MPI_COMM_WORLD);
  template <typename T,unsigned DIM>
  void send(util::ADf<T,DIM> *buf, int count, int dest, int tag = 0, MPI_Comm comm = MPI_COMM_WORLD);

  /// Initialize persistent non-blocking send
  template <typename T>
  void sendInit(T *buf, int count, int dest, MPI_Request* request, int tag = 0, MPI_Comm comm = MPI_COMM_WORLD);
  template <typename T,unsigned DIM>
  void sendInit(util::ADf<T,DIM> *buf, int count, int dest, MPI_Request* request, int tag = 0, MPI_Comm comm = MPI_COMM_WORLD);

  /// Sends data at *buf, non blocking
  template <typename T>
  void iSend(T *buf, int count, int dest, MPI_Request* request, int tag = 0, MPI_Comm comm = MPI_COMM_WORLD);
  template <typename T,unsigned DIM>
  void iSend(util::ADf<T,DIM> *buf, int count, int dest, MPI_Request* request, int tag = 0, MPI_Comm comm = MPI_COMM_WORLD);

  /// Sends data at *buf, non blocking and buffered
  template <typename T>
  void ibSend(T *buf, int count, int dest, MPI_Request* request, int tag = 0, MPI_Comm comm = MPI_COMM_WORLD);
  template <typename T,unsigned DIM>
  void ibSend(util::ADf<T,DIM> *buf, int count, int dest, MPI_Request* request, int tag = 0, MPI_Comm comm = MPI_COMM_WORLD);

  /// Probe size of incoming message
  std::size_t probeReceiveSize(int source, MPI_Datatype type, int tag = 0, MPI_Comm comm = MPI_COMM_WORLD);

  /// Receives data at *buf, blocking
  template <typename T>
  void receive(T *buf, int count, int source, int tag = 0, MPI_Comm comm = MPI_COMM_WORLD);
  template <typename T,unsigned DIM>
  void receive(util::ADf<T,DIM> *buf, int count, int source, int tag = 0, MPI_Comm comm = MPI_COMM_WORLD);

  /// Initialize persistent non-blocking receive
  template <typename T>
  void recvInit(T *buf, int count, int dest, MPI_Request* request, int tag = 0, MPI_Comm comm = MPI_COMM_WORLD);
  template <typename T,unsigned DIM>
  void recvInit(util::ADf<T,DIM> *buf, int count, int dest, MPI_Request* request, int tag = 0, MPI_Comm comm = MPI_COMM_WORLD);

  /// Receives data at *buf, non blocking
  template <typename T>
  void iRecv(T *buf, int count, int source, MPI_Request* request, int tag = 0, MPI_Comm comm = MPI_COMM_WORLD);
  template <typename T,unsigned DIM>
  void iRecv(util::ADf<T,DIM> *buf, int count, int source, MPI_Request* request, int tag = 0, MPI_Comm comm = MPI_COMM_WORLD);

  /// Send and receive data between two partners
  template <typename T>
  void sendRecv(T *sendBuf, T *recvBuf, int count, int dest, int source, int tag = 0,
                MPI_Comm comm = MPI_COMM_WORLD);
  template <typename T,unsigned DIM>
  void sendRecv(util::ADf<T,DIM> *sendBuf, util::ADf<T,DIM> *recvBuf, int count,
                int dest, int source, int tag = 0, MPI_Comm comm = MPI_COMM_WORLD);

  /// Sends data to master processor
  template <typename T>
  void sendToMaster(T* sendBuf, int sendCount, bool iAmRoot, MPI_Comm comm = MPI_COMM_WORLD);

  /// Scatter data from one processor over multiple processors
  template <typename T>
  void scatterV(T *sendBuf, T *recvBuf, int* sendCounts, int root = 0, MPI_Comm comm = MPI_COMM_WORLD);

  /// Gather data from multiple processors to one processor
  template <typename T>
  void gatherV(T* sendBuf, T* recvBuf, int *recvCounts, int root = 0, MPI_Comm comm = MPI_COMM_WORLD);


  /// Broadcast data from one processor to multiple processors
  template <typename T>
  void bCast(T* sendBuf, int sendCount, int root = 0, MPI_Comm comm = MPI_COMM_WORLD);
  template <typename T,unsigned DIM>
  void bCast(util::ADf<T,DIM>* sendBuf, int sendCount, int root = 0, MPI_Comm comm = MPI_COMM_WORLD);
  template <typename T,unsigned DIM>
  void bCast(BlockData<2,util::ADf<T,DIM>,util::ADf<T,DIM>>& sendData, int root = 0, MPI_Comm comm = MPI_COMM_WORLD);

  /// Broadcast data when root is unknown to other processors
  template <typename T>
  void bCastThroughMaster(T* sendBuf, int sendCount, bool iAmRoot, MPI_Comm comm = MPI_COMM_WORLD);
  template <typename T,unsigned DIM>
  void bCastThroughMaster(util::ADf<T,DIM>* sendBuf, int sendCount, bool iAmRoot, MPI_Comm comm = MPI_COMM_WORLD);

  /// Special case for broadcasting strings. Memory handling is automatic.
  void bCast(std::string& message, int root = 0);
  /// Special case for broadcasting BlockData2D
  void bCast(BlockData<2,double,double>& sendData, int root = 0, MPI_Comm comm = MPI_COMM_WORLD);
  /// Special case for broadcasting BlockData2D
  void bCast(BlockData<2,float,float>& sendData, int root = 0, MPI_Comm comm = MPI_COMM_WORLD);

  /// Reduction operation toward one processor
  template <typename T>
  void reduce(T& sendVal, T& recvVal, MPI_Op op, int root = 0, MPI_Comm = MPI_COMM_WORLD);
  template <typename T,unsigned DIM>
  void reduce(util::ADf<T,DIM>& sendVal, util::ADf<T,DIM>& recvVal,
              MPI_Op op, int root = 0, MPI_Comm = MPI_COMM_WORLD);
  template <typename T,unsigned DIM>
  void reduce(BlockData<2,util::ADf<T,DIM>,util::ADf<T,DIM>>& sendVal,
              BlockData<2,util::ADf<T,DIM>,util::ADf<T,DIM>>& recvVal,
              MPI_Op op, int root = 0, MPI_Comm comm = MPI_COMM_WORLD);

  /// Element-per-element reduction of a vector of data
  template <typename T>
  void reduceVect(std::vector<T>& sendVal, std::vector<T>& recvVal,
                  MPI_Op op, int root = 0, MPI_Comm comm = MPI_COMM_WORLD);

  /// Reduction operation, followed by a broadcast
  template <typename T>
  void reduceAndBcast(T& reductVal, MPI_Op op, int root = 0, MPI_Comm comm = MPI_COMM_WORLD);
  template <typename T,unsigned DIM>
  void reduceAndBcast(util::ADf<T,DIM>& reductVal, MPI_Op op, int root = 0, MPI_Comm comm = MPI_COMM_WORLD);

  /// Complete a non-blocking MPI operation
  void wait(MPI_Request* request, MPI_Status* status);

  /// Complete a series of non-blocking MPI operations
  void waitAll(MpiNonBlockingHelper& mpiNbHelper);

private:
  /// Implementation code for Scatter
  template <typename T>
  void scatterv_impl(T *sendBuf, int* sendCounts, int* displs,
                     T* recvBuf, int recvCount, int root, MPI_Comm comm);

  /// Implementation code for Gather
  template <typename T>
  void gatherv_impl(T* sendBuf, int sendCount, T* recvBuf, int* recvCounts, int* displs,
                    int root, MPI_Comm comm);
private:
  MpiManager();
  ~MpiManager();
private:
  int numTasks, taskId;
  bool ok;
  mutable OstreamManager clout;

  friend MpiManager& mpi();
};

#else

class MpiManager {
public:
  /// Initializes the mpi manager
  void init(int *argc, char ***argv, bool verbose=false) { }
  /// Returns the number of processes
  int getSize() const
  {
    return 1;
  }
  /// Returns the process ID
  int getRank() const
  {
    return 0;
  }
  /// Returns process ID of main processor
  int bossId() const
  {
    return 0;
  }
  /// Tells whether current processor is main processor
  bool isMainProcessor() const
  {
    return true;
  }

  /// Synchronizes the processes
  void barrier() const {};

  friend MpiManager& mpi();
};

#endif

inline MpiManager& mpi()
{
  static MpiManager instance;
  return instance;
}

#ifdef PARALLEL_MODE_MPI

MpiManager::MpiManager() : ok(false), clout(std::cout,"MpiManager")
{ }

MpiManager::~MpiManager()
{
  if (ok) {
    MPI_Finalize();
    ok = false;
  }
}

void MpiManager::init(int *argc, char ***argv, bool verbose)
{

  int ok1 = MPI_Init(argc, argv);
  int ok2 = MPI_Comm_rank(MPI_COMM_WORLD, &taskId);
  int ok3 = MPI_Comm_size(MPI_COMM_WORLD, &numTasks);
  int ok4 = MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL);
  ok = (ok1 == MPI_SUCCESS && ok2 == MPI_SUCCESS && ok3 == MPI_SUCCESS && ok4 == MPI_SUCCESS);
  if (verbose) {
    clout << "Sucessfully initialized, numThreads=" << getSize() << std::endl;
  }
}

int MpiManager::getSize() const
{
  return numTasks;
}

int MpiManager::getRank() const
{
  return taskId;
}

int MpiManager::bossId() const
{
  return 0;
}

bool MpiManager::isMainProcessor() const
{
  return bossId() == getRank();
}

double MpiManager::getTime() const
{
  if (!ok) {
    return 0.;
  }
  return MPI_Wtime();
}

void MpiManager::barrier(MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Barrier(comm);
}

template <>
void MpiManager::send<bool>(bool *buf, int count, int dest, int tag, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Send(static_cast<void*>(buf), count, MPI_BYTE, dest, tag, comm);
}

template <>
void MpiManager::send<char>(char *buf, int count, int dest, int tag, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Send(static_cast<void*>(buf), count, MPI_CHAR, dest, tag, comm);
}

template <>
void MpiManager::send<std::uint8_t>(std::uint8_t *buf, int count, int dest, int tag, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Send(static_cast<void*>(buf), count, MPI_BYTE, dest, tag, comm);
}

template <>
void MpiManager::send<int>(int *buf, int count, int dest, int tag, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Send(static_cast<void*>(buf), count, MPI_INT, dest, tag, comm);
}

template <>
void MpiManager::send<float>(float *buf, int count, int dest, int tag, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Send(static_cast<void*>(buf), count, MPI_FLOAT, dest, tag, comm);
}

template <>
void MpiManager::send<double>(double *buf, int count, int dest, int tag, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Send(static_cast<void*>(buf), count, MPI_DOUBLE, dest, tag, comm);
}

template <>
void MpiManager::sendInit<double>(double *buf, int count, int dest, MPI_Request* request, int tag, MPI_Comm comm)
{
  if (ok) {
    MPI_Send_init(buf, count, MPI_DOUBLE, dest, tag, comm, request);
  }
}

template <>
void MpiManager::sendInit<std::size_t>(std::size_t *buf, int count, int dest, MPI_Request* request, int tag, MPI_Comm comm)
{
  if (ok) {
    MPI_Send_init(buf, count, MPI_UNSIGNED_LONG, dest, tag, comm, request);
  }
}

template <>
void MpiManager::sendInit<std::uint32_t>(std::uint32_t *buf, int count, int dest, MPI_Request* request, int tag, MPI_Comm comm)
{
  if (ok) {
    MPI_Send_init(buf, count, MPI_UNSIGNED, dest, tag, comm, request);
  }
}

template <>
void MpiManager::sendInit<std::uint8_t>(std::uint8_t *buf, int count, int dest, MPI_Request* request, int tag, MPI_Comm comm)
{
  if (ok) {
    MPI_Send_init(buf, count, MPI_BYTE, dest, tag, comm, request);
  }
}

template <>
void MpiManager::sendInit<int>(int *buf, int count, int dest, MPI_Request* request, int tag, MPI_Comm comm)
{
  if (ok) {
    MPI_Send_init(buf, count, MPI_INT, dest, tag, comm, request);
  }
}

template <>
void MpiManager::sendInit<bool>(bool *buf, int count, int dest, MPI_Request* request, int tag, MPI_Comm comm)
{
  if (ok) {
    MPI_Send_init(static_cast<void*>(buf), count, MPI_BYTE, dest, tag, comm, request);
  }
}

template <>
void MpiManager::iSend<bool>
(bool *buf, int count, int dest, MPI_Request* request, int tag, MPI_Comm comm)
{
  if (ok) {
    MPI_Isend(static_cast<void*>(buf), count, MPI_BYTE, dest, tag, comm, request);
  }
}

template <>
void MpiManager::iSend<char>
(char *buf, int count, int dest, MPI_Request* request, int tag, MPI_Comm comm)
{
  if (ok) {
    MPI_Isend(static_cast<void*>(buf), count, MPI_CHAR, dest, tag, comm, request);
  }
}

template <>
void MpiManager::iSend<std::uint8_t>
(std::uint8_t *buf, int count, int dest, MPI_Request* request, int tag, MPI_Comm comm)
{
  if (ok) {
    MPI_Isend(static_cast<void*>(buf), count, MPI_BYTE, dest, tag, comm, request);
  }
}

template <>
void MpiManager::iSend<int>
(int *buf, int count, int dest, MPI_Request* request, int tag, MPI_Comm comm)
{
  if (ok) {
    MPI_Isend(static_cast<void*>(buf), count, MPI_INT, dest, tag, comm, request);
  }
}

template <>
void MpiManager::iSend<std::size_t>
(std::size_t *buf, int count, int dest, MPI_Request* request, int tag, MPI_Comm comm)
{
  if (ok) {
    MPI_Isend(static_cast<void*>(buf), count, MPI_UNSIGNED_LONG, dest, tag, comm, request);
  }
}

template <>
void MpiManager::iSend<std::uint32_t>
(std::uint32_t *buf, int count, int dest, MPI_Request* request, int tag, MPI_Comm comm)
{
  if (ok) {
    MPI_Isend(static_cast<void*>(buf), count, MPI_UNSIGNED, dest, tag, comm, request);
  }
}

template <>
void MpiManager::iSend<float>
(float *buf, int count, int dest, MPI_Request* request, int tag, MPI_Comm comm)
{
  if (ok) {
    MPI_Isend(static_cast<void*>(buf), count, MPI_FLOAT, dest, tag, comm, request);
  }
}

template <>
void MpiManager::iSend<double>
(double *buf, int count, int dest, MPI_Request* request, int tag, MPI_Comm comm)
{
  if (ok) {
    MPI_Isend(static_cast<void*>(buf), count, MPI_DOUBLE, dest, tag, comm, request);
  }
}

template <>
void MpiManager::iSend<long double>
(long double *buf, int count, int dest, MPI_Request* request, int tag, MPI_Comm comm)
{
  if (ok) {
    MPI_Isend(static_cast<void*>(buf), count, MPI_LONG_DOUBLE, dest, tag, comm, request);
  }
}


template <>
void MpiManager::ibSend<bool>
(bool *buf, int count, int dest, MPI_Request* request, int tag, MPI_Comm comm)
{
  if (ok) {
    MPI_Ibsend(static_cast<void*>(buf), count, MPI_BYTE, dest, tag, comm, request);
  }
}

template <>
void MpiManager::ibSend<char>
(char *buf, int count, int dest, MPI_Request* request, int tag, MPI_Comm comm)
{
  if (ok) {
    MPI_Ibsend(static_cast<void*>(buf), count, MPI_CHAR, dest, tag, comm, request);
  }
}

template <>
void MpiManager::ibSend<int>
(int *buf, int count, int dest, MPI_Request* request, int tag, MPI_Comm comm)
{
  if (ok) {
    MPI_Ibsend(static_cast<void*>(buf), count, MPI_INT, dest, tag, comm, request);
  }
}

template <>
void MpiManager::ibSend<float>
(float *buf, int count, int dest, MPI_Request* request, int tag, MPI_Comm comm)
{
  if (ok) {
    MPI_Ibsend(static_cast<void*>(buf), count, MPI_FLOAT, dest, tag, comm, request);
  }
}

template <>
void MpiManager::ibSend<double>
(double *buf, int count, int dest, MPI_Request* request, int tag, MPI_Comm comm)
{
  if (ok) {
    MPI_Ibsend(static_cast<void*>(buf), count, MPI_DOUBLE, dest, tag, comm, request);
  }
}

std::size_t MpiManager::probeReceiveSize(int source, MPI_Datatype type, int tag, MPI_Comm comm)
{
  MPI_Status status;
  MPI_Probe(source, tag, comm, &status);
  int requestSize;
  MPI_Get_count(&status, type, &requestSize);
  if (requestSize == MPI_UNDEFINED) {
    throw std::runtime_error("MPI_UNDEFINED");
  }
  return requestSize;
}

template <>
void MpiManager::receive<bool>(bool *buf, int count, int source, int tag, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Status status;
  MPI_Recv(static_cast<void*>(buf), count, MPI_BYTE, source, tag, comm, &status);
}


template <>
void MpiManager::receive<char>(char *buf, int count, int source, int tag, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Status status;
  MPI_Recv(static_cast<void*>(buf), count, MPI_CHAR, source, tag, comm, &status);
}

template <>
void MpiManager::receive<std::uint8_t>(std::uint8_t *buf, int count, int source, int tag, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Status status;
  MPI_Recv(static_cast<std::uint8_t*>(buf), count, MPI_BYTE, source, tag, comm, &status);
}

template <>
void MpiManager::receive<int>(int *buf, int count, int source, int tag, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Status status;
  MPI_Recv(static_cast<void*>(buf), count, MPI_INT, source, tag, comm, &status);
}

template <>
void MpiManager::receive<std::size_t>(std::size_t *buf, int count, int source, int tag, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Status status;
  MPI_Recv(static_cast<void*>(buf), count, MPI_UNSIGNED_LONG, source, tag, comm, &status);
}

template <>
void MpiManager::receive<std::uint32_t>(std::uint32_t *buf, int count, int source, int tag, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Status status;
  MPI_Recv(static_cast<void*>(buf), count, MPI_UNSIGNED, source, tag, comm, &status);
}

template <>
void MpiManager::receive<float>(float *buf, int count, int source, int tag, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Status status;
  MPI_Recv(static_cast<void*>(buf), count, MPI_FLOAT, source, tag, comm, &status);
}

template <>
void MpiManager::receive<double>(double *buf, int count, int source, int tag, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Status status;
  MPI_Recv(static_cast<void*>(buf), count, MPI_DOUBLE, source, tag, comm, &status);
}

template <>
void MpiManager::receive<long double>(long double *buf, int count, int source, int tag, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Status status;
  MPI_Recv(static_cast<void*>(buf), count, MPI_LONG_DOUBLE, source, tag, comm, &status);
}

template <>
void MpiManager::sendToMaster<bool>(bool* sendBuf, int sendCount, bool iAmRoot, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  if (iAmRoot && !isMainProcessor()) {
    send(sendBuf, sendCount, 0);
  }
  if (isMainProcessor() && !iAmRoot) {
    receive(sendBuf, sendCount, MPI_ANY_SOURCE);
  }
}

template <>
void MpiManager::sendToMaster<char>(char* sendBuf, int sendCount, bool iAmRoot, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  if (iAmRoot && !isMainProcessor()) {
    send(sendBuf, sendCount, 0);
  }
  if (isMainProcessor() && !iAmRoot) {
    receive(sendBuf, sendCount, MPI_ANY_SOURCE);
  }
}

template <>
void MpiManager::sendToMaster<int>(int* sendBuf, int sendCount, bool iAmRoot, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  if (iAmRoot && !isMainProcessor()) {
    send(sendBuf, sendCount, 0);
  }
  if (isMainProcessor() && !iAmRoot) {
    receive(sendBuf, sendCount, MPI_ANY_SOURCE);
  }
}

template <>
void MpiManager::sendToMaster<float>(float* sendBuf, int sendCount, bool iAmRoot, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  if (iAmRoot && !isMainProcessor()) {
    send(sendBuf, sendCount, 0);
  }
  if (isMainProcessor() && !iAmRoot) {
    receive(sendBuf, sendCount, MPI_ANY_SOURCE);
  }
}

template <>
void MpiManager::sendToMaster<double>(double* sendBuf, int sendCount, bool iAmRoot, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  if (iAmRoot && !isMainProcessor()) {
    send(sendBuf, sendCount, 0);
  }
  if (isMainProcessor() && !iAmRoot) {
    receive(sendBuf, sendCount, MPI_ANY_SOURCE);
  }
}

template <>
void MpiManager::recvInit<double>(double *buf, int count, int dest, MPI_Request* request, int tag, MPI_Comm comm)
{
  if (ok) {
    MPI_Recv_init(buf, count, MPI_DOUBLE, dest, tag, comm, request);
  }
}

template <>
void MpiManager::recvInit<int>(int *buf, int count, int dest, MPI_Request* request, int tag, MPI_Comm comm)
{
  if (ok) {
    MPI_Recv_init(buf, count, MPI_INT, dest, tag, comm, request);
  }
}

template <>
void MpiManager::recvInit<std::uint8_t>(std::uint8_t *buf, int count, int dest, MPI_Request* request, int tag, MPI_Comm comm)
{
  if (ok) {
    MPI_Recv_init(buf, count, MPI_BYTE, dest, tag, comm, request);
  }
}

template <>
void MpiManager::iRecv<bool>(bool *buf, int count, int source, MPI_Request* request, int tag, MPI_Comm comm)
{
  if (ok) {
    MPI_Irecv(static_cast<void*>(buf), count, MPI_BYTE, source, tag, comm, request);
  }
}

template <>
void MpiManager::iRecv<char>(char *buf, int count, int source, MPI_Request* request, int tag, MPI_Comm comm)
{
  if (ok) {
    MPI_Irecv(static_cast<void*>(buf), count, MPI_CHAR, source, tag, comm, request);
  }
}

template <>
void MpiManager::iRecv<int>(int *buf, int count, int source, MPI_Request* request, int tag, MPI_Comm comm)
{
  if (ok) {
    MPI_Irecv(static_cast<void*>(buf), count, MPI_INT, source, tag, comm, request);
  }
}

template <>
void MpiManager::iRecv<float>(float *buf, int count, int source, MPI_Request* request, int tag, MPI_Comm comm)
{
  if (ok) {
    MPI_Irecv(static_cast<void*>(buf), count, MPI_FLOAT, source, tag, comm, request);
  }
}

template <>
void MpiManager::iRecv<double>(double *buf, int count, int source, MPI_Request* request, int tag, MPI_Comm comm)
{
  if (ok) {
    MPI_Irecv(static_cast<void*>(buf), count, MPI_DOUBLE, source, tag, comm, request);
  }
}

template <>
void MpiManager::sendRecv<bool>
(bool *sendBuf, bool *recvBuf, int count, int dest, int source, int tag, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Status status;
  MPI_Sendrecv(static_cast<void*>(sendBuf),
               count,
               MPI_BYTE, dest, tag,
               static_cast<void*>(recvBuf),
               count,
               MPI_BYTE, source, tag, comm, &status);
}

template <>
void MpiManager::sendRecv<char>
(char *sendBuf, char *recvBuf, int count, int dest, int source, int tag, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Status status;
  MPI_Sendrecv(static_cast<void*>(sendBuf),
               count,
               MPI_CHAR, dest, tag,
               static_cast<void*>(recvBuf),
               count,
               MPI_CHAR, source, tag, comm, &status);
}

template <>
void MpiManager::sendRecv<int>
(int *sendBuf, int *recvBuf, int count, int dest, int source, int tag, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Status status;
  MPI_Sendrecv(static_cast<void*>(sendBuf),
               count,
               MPI_INT, dest, tag,
               static_cast<void*>(recvBuf),
               count,
               MPI_INT, source, tag, comm, &status);
}

template <>
void MpiManager::sendRecv<float>
(float *sendBuf, float *recvBuf, int count, int dest, int source, int tag, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Status status;
  MPI_Sendrecv(static_cast<void*>(sendBuf),
               count,
               MPI_FLOAT, dest, tag,
               static_cast<void*>(recvBuf),
               count,
               MPI_FLOAT, source, tag, comm, &status);
}

template <>
void MpiManager::sendRecv<long>
(long *sendBuf, long *recvBuf, int count, int dest, int source, int tag, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Status status;
  MPI_Sendrecv(static_cast<void*>(sendBuf),
               count,
               MPI_LONG, dest, tag,
               static_cast<void*>(recvBuf),
               count,
               MPI_LONG, source, tag, comm, &status);
}

template <>
void MpiManager::sendRecv<double>
(double *sendBuf, double *recvBuf, int count, int dest, int source, int tag, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Status status;
  MPI_Sendrecv(static_cast<void*>(sendBuf),
               count,
               MPI_DOUBLE, dest, tag,
               static_cast<void*>(recvBuf),
               count,
               MPI_DOUBLE, source, tag, comm, &status);
}

template <>
void MpiManager::sendRecv<long double>
(long double *sendBuf, long double *recvBuf, int count, int dest, int source, int tag, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Status status;
  MPI_Sendrecv(static_cast<void*>(sendBuf),
               count,
               MPI_LONG_DOUBLE, dest, tag,
               static_cast<void*>(recvBuf),
               count,
               MPI_LONG_DOUBLE, source, tag, comm, &status);
}

template <>
void MpiManager::scatterv_impl<bool>(bool* sendBuf, int* sendCounts, int* displs,
                                     bool* recvBuf, int recvCount, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Scatterv(static_cast<void*>(sendBuf),
               sendCounts, displs, MPI_BYTE,
               static_cast<void*>(recvBuf),
               recvCount, MPI_BYTE, root, comm);
}

template <>
void MpiManager::scatterv_impl<char>(char* sendBuf, int* sendCounts, int* displs,
                                     char* recvBuf, int recvCount, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Scatterv(static_cast<void*>(sendBuf),
               sendCounts, displs, MPI_CHAR,
               static_cast<void*>(recvBuf),
               recvCount, MPI_CHAR, root, comm);
}

template <>
void MpiManager::scatterv_impl<int>(int *sendBuf, int* sendCounts, int* displs,
                                    int* recvBuf, int recvCount, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Scatterv(static_cast<void*>(sendBuf),
               sendCounts, displs, MPI_INT,
               static_cast<void*>(recvBuf),
               recvCount, MPI_INT, root, comm);
}

template <>
void MpiManager::scatterv_impl<float>(float *sendBuf, int* sendCounts, int* displs,
                                      float* recvBuf, int recvCount, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Scatterv(static_cast<void*>(sendBuf),
               sendCounts, displs, MPI_FLOAT,
               static_cast<void*>(recvBuf),
               recvCount, MPI_FLOAT, root, comm);
}

template <>
void MpiManager::scatterv_impl<double>(double *sendBuf, int* sendCounts, int* displs,
                                       double* recvBuf, int recvCount, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Scatterv(static_cast<void*>(sendBuf),
               sendCounts, displs, MPI_DOUBLE,
               static_cast<void*>(recvBuf),
               recvCount, MPI_DOUBLE, root, comm);
}

template <>
void MpiManager::gatherv_impl<bool>(bool* sendBuf, int sendCount,
                                    bool* recvBuf, int* recvCounts, int* displs,
                                    int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Gatherv(static_cast<void*>(sendBuf), sendCount, MPI_BYTE,
              static_cast<void*>(recvBuf), recvCounts, displs, MPI_BYTE,
              root, comm);
}

template <>
void MpiManager::gatherv_impl<char>(char* sendBuf, int sendCount,
                                    char* recvBuf, int* recvCounts, int* displs,
                                    int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Gatherv(static_cast<void*>(sendBuf), sendCount, MPI_CHAR,
              static_cast<void*>(recvBuf), recvCounts, displs, MPI_CHAR,
              root, comm);
}

template <>
void MpiManager::gatherv_impl<int>(int* sendBuf, int sendCount,
                                   int* recvBuf, int* recvCounts, int* displs,
                                   int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Gatherv(static_cast<void*>(sendBuf), sendCount, MPI_INT,
              static_cast<void*>(recvBuf), recvCounts, displs, MPI_INT,
              root, comm);
}

template <>
void MpiManager::gatherv_impl<float>(float* sendBuf, int sendCount,
                                     float* recvBuf, int* recvCounts, int* displs,
                                     int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Gatherv(static_cast<void*>(sendBuf), sendCount, MPI_FLOAT,
              static_cast<void*>(recvBuf), recvCounts, displs, MPI_FLOAT,
              root, comm);
}

template <>
void MpiManager::gatherv_impl<double>(double* sendBuf, int sendCount,
                                      double* recvBuf, int* recvCounts, int* displs,
                                      int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Gatherv(static_cast<void*>(sendBuf), sendCount, MPI_DOUBLE,
              static_cast<void*>(recvBuf), recvCounts, displs, MPI_DOUBLE,
              root, comm);
}

template <>
void MpiManager::bCast<bool>(bool* sendBuf, int sendCount, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Bcast(static_cast<void*>(sendBuf),
            sendCount, MPI_BYTE, root, comm);
}

template <>
void MpiManager::bCast<char>(char* sendBuf, int sendCount, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Bcast(static_cast<void*>(sendBuf),
            sendCount, MPI_CHAR, root, comm);
}

template <>
void MpiManager::bCast<unsigned char>(unsigned char* sendBuf, int sendCount, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Bcast(static_cast<void*>(sendBuf),
            sendCount, MPI_UNSIGNED_CHAR, root, comm);
}

template <>
void MpiManager::bCast<int>(int* sendBuf, int sendCount, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Bcast(static_cast<void*>(sendBuf),
            sendCount, MPI_INT, root, comm);
}

template <>
void MpiManager::bCast<unsigned long>(unsigned long* sendBuf, int sendCount, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Bcast(static_cast<void*>(sendBuf),
            sendCount, MPI_UNSIGNED_LONG, root, comm);
}

template <>
void MpiManager::bCast<float>(float* sendBuf, int sendCount, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Bcast(static_cast<void*>(sendBuf),
            sendCount, MPI_FLOAT, root, comm);
}

template <>
void MpiManager::bCast<double>(double* sendBuf, int sendCount, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Bcast(static_cast<void*>(sendBuf),
            sendCount, MPI_DOUBLE, root, comm);
}

template <>
void MpiManager::bCast<std::string>(std::string* sendBuf, int sendCount, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  int length = (int) sendBuf->size();
  MPI_Bcast(static_cast<void*>(&length), 1, MPI_INT, root, comm);
  char* buffer = new char[length+1];
  if (getRank()==root) {
    std::copy(sendBuf->c_str(), sendBuf->c_str()+length+1, buffer);
  }
  MPI_Bcast(static_cast<void*>(buffer), length+1, MPI_CHAR, root, comm);
  if (getRank()!=root) {
    *sendBuf = buffer;
  }
  delete [] buffer;
}

void MpiManager::bCast(BlockData<2,double,double>& sendData, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  for (unsigned iD=0; iD < sendData.getSize(); ++iD) {
    MPI_Bcast(static_cast<void*>(sendData.getColumn(iD).data()),
              sendData.getNcells(), MPI_DOUBLE, root, comm);
  }
}

void MpiManager::bCast(BlockData<2,float,float>& sendData, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  for (unsigned iD=0; iD < sendData.getSize(); ++iD) {
    MPI_Bcast(static_cast<void*>(sendData.getColumn(iD).data()),
              sendData.getNcells(), MPI_FLOAT, root, comm);
  }
}

template <>
void MpiManager::bCastThroughMaster<bool>(bool* sendBuf, int sendCount, bool iAmRoot, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  if (iAmRoot && !isMainProcessor()) {
    send(sendBuf, sendCount, 0);
  }
  if (isMainProcessor() && !iAmRoot) {
    receive(sendBuf, sendCount, MPI_ANY_SOURCE);
  }
  bCast(sendBuf, sendCount, 0);
}

template <>
void MpiManager::bCastThroughMaster<char>(char* sendBuf, int sendCount, bool iAmRoot, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  if (iAmRoot && !isMainProcessor()) {
    send(sendBuf, sendCount, 0);
  }
  if (isMainProcessor() && !iAmRoot) {
    receive(sendBuf, sendCount, MPI_ANY_SOURCE);
  }
  bCast(sendBuf, sendCount, 0);
}

template <>
void MpiManager::bCastThroughMaster<int>(int* sendBuf, int sendCount, bool iAmRoot, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  if (iAmRoot && !isMainProcessor()) {
    send(sendBuf, sendCount, 0);
  }
  if (isMainProcessor() && !iAmRoot) {
    receive(sendBuf, sendCount, MPI_ANY_SOURCE);
  }
  bCast(sendBuf, sendCount, 0);
}

template <>
void MpiManager::bCastThroughMaster<float>(float* sendBuf, int sendCount, bool iAmRoot, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  if (iAmRoot && !isMainProcessor()) {
    send(sendBuf, sendCount, 0);
  }
  if (isMainProcessor() && !iAmRoot) {
    receive(sendBuf, sendCount, MPI_ANY_SOURCE);
  }
  bCast(sendBuf, sendCount, 0);
}

template <>
void MpiManager::bCastThroughMaster<double>(double* sendBuf, int sendCount, bool iAmRoot, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  if (iAmRoot && !isMainProcessor()) {
    send(sendBuf, sendCount, 0);
  }
  if (isMainProcessor() && !iAmRoot) {
    receive(sendBuf, sendCount, MPI_ANY_SOURCE);
  }
  bCast(sendBuf, sendCount, 0);
}

template <>
void MpiManager::reduce<bool>(bool& sendVal, bool& recvVal,  MPI_Op op, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Reduce(static_cast<void*>(&sendVal),
             static_cast<void*>(&recvVal), 1, MPI_BYTE, op, root, comm);
}

template <>
void MpiManager::reduce<char>(char& sendVal, char& recvVal,  MPI_Op op, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Reduce(static_cast<void*>(&sendVal),
             static_cast<void*>(&recvVal), 1, MPI_CHAR, op, root, comm);
}

template <>
void MpiManager::reduce<int>(int& sendVal, int& recvVal,  MPI_Op op, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Reduce(static_cast<void*>(&sendVal),
             static_cast<void*>(&recvVal), 1, MPI_INT, op, root, comm);
}

template <>
void MpiManager::reduce<float>(float& sendVal, float& recvVal,  MPI_Op op, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Reduce(static_cast<void*>(&sendVal),
             static_cast<void*>(&recvVal), 1, MPI_FLOAT, op, root, comm);
}

template <>
void MpiManager::reduce<double>(double& sendVal, double& recvVal,  MPI_Op op, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Reduce(static_cast<void*>(&sendVal),
             static_cast<void*>(&recvVal), 1, MPI_DOUBLE, op, root, comm);
}


template <>
void MpiManager::reduce<BlockData<2,double,int> >(BlockData<2,double,int>& sendVal, BlockData<2,double,int>& recvVal,  MPI_Op op, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  for (unsigned iD=0; iD < sendVal.getSize(); ++iD) {
    MPI_Reduce(static_cast<void*>(sendVal.getColumn(iD).data()),
               static_cast<void*>(recvVal.getColumn(iD).data()),
               sendVal.getNcells(), MPI_DOUBLE, op, root, comm);
  }
}

template <>
void MpiManager::reduce<BlockData<2,double,double> >(BlockData<2,double,double>& sendVal, BlockData<2,double,double>& recvVal,  MPI_Op op, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  for (unsigned iD=0; iD < sendVal.getSize(); ++iD) {
    MPI_Reduce(static_cast<void*>(sendVal.getColumn(iD).data()),
               static_cast<void*>(recvVal.getColumn(iD).data()),
               sendVal.getNcells(), MPI_DOUBLE, op, root, comm);
  }
}

template <>
void MpiManager::reduce<BlockData<2,float,float> >(BlockData<2,float,float>& sendVal, BlockData<2,float,float>& recvVal,  MPI_Op op, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  for (unsigned iD=0; iD < sendVal.getSize(); ++iD) {
    MPI_Reduce(static_cast<void*>(sendVal.getColumn(iD).data()),
               static_cast<void*>(recvVal.getColumn(iD).data()),
               sendVal.getNcells(), MPI_FLOAT, op, root, comm);
  }
}

/*template <>
void MpiManager::reduceVect<bool>(std::vector<bool>& sendVal, std::vector<bool>& recvVal,
                                  MPI_Op op, int root, MPI_Comm comm)
{
  if (!ok) return;
  MPI_Reduce(static_cast<void*>(&(sendVal[0])),
             static_cast<void*>(&(recvVal[0])),
             sendVal.size(), MPI_BYTE, op, root, comm);
}
*/
template <>
void MpiManager::reduceVect<char>(std::vector<char>& sendVal, std::vector<char>& recvVal,
                                  MPI_Op op, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Reduce(static_cast<void*>(&(sendVal[0])),
             static_cast<void*>(&(recvVal[0])),
             sendVal.size(), MPI_CHAR, op, root, comm);
}

template <>
void MpiManager::reduceVect<int>(std::vector<int>& sendVal, std::vector<int>& recvVal,
                                 MPI_Op op, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Reduce(static_cast<void*>(&(sendVal[0])),
             static_cast<void*>(&(recvVal[0])),
             sendVal.size(), MPI_INT, op, root, comm);
}

template <>
void MpiManager::reduceVect<float>(std::vector<float>& sendVal, std::vector<float>& recvVal,
                                   MPI_Op op, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Reduce(static_cast<void*>(&(sendVal[0])),
             static_cast<void*>(&(recvVal[0])),
             sendVal.size(), MPI_FLOAT, op, root, comm);
}

template <>
void MpiManager::reduceVect<double>(std::vector<double>& sendVal, std::vector<double>& recvVal,
                                    MPI_Op op, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  MPI_Reduce(static_cast<void*>(&(sendVal[0])),
             static_cast<void*>(&(recvVal[0])),
             sendVal.size(), MPI_DOUBLE, op, root, comm);
}

template <>
void MpiManager::reduceAndBcast<bool>(bool& reductVal, MPI_Op op, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  char recvVal;
  MPI_Reduce(static_cast<void*>(&reductVal), static_cast<void*>(&recvVal), 1, MPI_BYTE, op, root, comm);
  reductVal = recvVal;
  MPI_Bcast(&reductVal, 1, MPI_BYTE, root, comm);

}

template <>
void MpiManager::reduceAndBcast<char>(char& reductVal, MPI_Op op, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  char recvVal;
  MPI_Reduce(&reductVal, &recvVal, 1, MPI_CHAR, op, root, comm);
  reductVal = recvVal;
  MPI_Bcast(&reductVal, 1, MPI_CHAR, root, comm);

}

template <>
void MpiManager::reduceAndBcast<int>(int& reductVal, MPI_Op op, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  int recvVal;
  MPI_Reduce(&reductVal, &recvVal, 1, MPI_INT, op, root, comm);
  reductVal = recvVal;
  MPI_Bcast(&reductVal, 1, MPI_INT, root, comm);

}

template <>
void MpiManager::reduceAndBcast<float>(float& reductVal, MPI_Op op, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  float recvVal;
  MPI_Reduce(&reductVal, &recvVal, 1, MPI_FLOAT, op, root, comm);
  reductVal = recvVal;
  MPI_Bcast(&reductVal, 1, MPI_FLOAT, root, comm);

}

template <>
void MpiManager::reduceAndBcast<double>(double& reductVal, MPI_Op op, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  double recvVal;
  MPI_Reduce(&reductVal, &recvVal, 1, MPI_DOUBLE, op, root, comm);
  reductVal = recvVal;
  MPI_Bcast(&reductVal, 1, MPI_DOUBLE, root, comm);

}

template <>
void MpiManager::reduceAndBcast<long double>(long double& reductVal, MPI_Op op, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  long double recvVal;
  MPI_Reduce(&reductVal, &recvVal, 1, MPI_LONG_DOUBLE, op, root, comm);
  reductVal = recvVal;
  MPI_Bcast(&reductVal, 1, MPI_LONG_DOUBLE, root, comm);

}

template <>
void MpiManager::reduceAndBcast<long>(long& reductVal, MPI_Op op, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  long recvVal;
  MPI_Reduce(&reductVal, &recvVal, 1, MPI_LONG, op, root, comm);
  reductVal = recvVal;
  MPI_Bcast(&reductVal, 1, MPI_LONG, root, comm);

}

template <>
void MpiManager::reduceAndBcast<unsigned long>(unsigned long& reductVal, MPI_Op op, int root, MPI_Comm comm)
{
  if (!ok) {
    return;
  }
  unsigned long recvVal;
  MPI_Reduce(&reductVal, &recvVal, 1, MPI_UNSIGNED_LONG, op, root, comm);
  reductVal = recvVal;
  MPI_Bcast(&reductVal, 1, MPI_UNSIGNED_LONG, root, comm);

}

void MpiManager::wait(MPI_Request* request, MPI_Status* status)
{
  if (!ok) {
    return;
  }
  MPI_Wait(request, status);
}

void MpiManager::waitAll(MpiNonBlockingHelper& mpiNbHelper)
{
  if (!ok || mpiNbHelper.get_size() == 0) {
    return;
  }
  MPI_Waitall(mpiNbHelper.get_size(), mpiNbHelper.get_mpiRequest(), mpiNbHelper.get_mpiStatus());
}


MpiNonBlockingHelper::MpiNonBlockingHelper():
  _size(0)
{ }

void MpiNonBlockingHelper::swap(MpiNonBlockingHelper& rhs)
{
  std::swap(_size, rhs._size);
  std::swap(_mpiRequest, rhs._mpiRequest);
  std::swap(_mpiStatus, rhs._mpiStatus);
}

void MpiNonBlockingHelper::allocate(unsigned n)
{
  free();
  _mpiRequest.reset(new MPI_Request[n] { });
  _mpiStatus.reset(new MPI_Status[n] { });
  _size = n;
}

void MpiNonBlockingHelper::free()
{
  _size = 0;
}

unsigned MpiNonBlockingHelper::get_size() const
{
  return _size;
}

MPI_Request* MpiNonBlockingHelper::get_mpiRequest(int i) const
{
  OLB_PRECONDITION(size_t(i) < _size);
  return &_mpiRequest[i];
}

MPI_Status* MpiNonBlockingHelper::get_mpiStatus(int i) const
{
  OLB_PRECONDITION(size_t(i) < _size);
  return &_mpiStatus[i];
}

void MpiNonBlockingHelper::start(int i)
{
  MPI_Start(get_mpiRequest(i));
}

void MpiNonBlockingHelper::wait(int i)
{
  MPI_Wait(get_mpiRequest(i), get_mpiStatus(i));
}

bool MpiNonBlockingHelper::isDone(int i)
{
  int done;
  MPI_Test(get_mpiRequest(i), &done, MPI_STATUS_IGNORE);
  return done;
}

#endif  // PARALLEL_MODE_MPI

}  // namespace singleton

}  // namespace olb


#endif  // MPI_MANAGER_H
