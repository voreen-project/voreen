/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007, 2014 Mathias J. Krause
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

/** \file
 * The description of a 3D cuboid neighbourhood -- header file.
 */


#ifndef CUBOID_NEIGHBOURHOOD_3D_H
#define CUBOID_NEIGHBOURHOOD_3D_H

#include "communication/mpiManager.h"
#include "communication/superStructure3D.h"
#include "core/blockStructure3D.h"

#include <vector>


/// All OpenLB code is contained in this namespace.
namespace olb {

template<typename T> class SuperStructure3D;


/// Single 3D cuboid neighbourhoods are the basic component of a
/// 3D communicator
/** For each cuboid a cuboid neighbourhood is defined. It stores the
 * needed cell coordinates (Cell3D) of other cuboids. Futhermore this
 * class provides basic initialization and communication methods
 * for the class Communicator3D.
 *
 * WARNING: For unstructured grids there is an interpolation needed
 * for the method buffer_outData which is not yet implemented!
 *
 * This class is not intended to be derived from.
 */
template<typename T>
struct Cell3D {
  Cell3D() : latticeR {0,0,0,0} {};
  // local position latticeR: iC, iX, iY, iZ;
  int latticeR[4];
  std::size_t latticeCellId;
  // global position physR: x, y, z;
  T physR[3];

  bool operator==(Cell3D const& rhs) const
  {
      return latticeR == rhs.latticeR;
  };

  /// Copy constructor
  Cell3D(Cell3D const& rhs) = default;

  /// Copy assignment
  Cell3D& operator=(Cell3D const& rhs) = default;

};


template<typename T>
class CuboidNeighbourhood3D {
private:
  /// Reference to the super structure
  SuperStructure3D<T>& _superStructure;
  /// Cuboid ID
  const int _iCglob;
  /// Number of cubboids in the structure
  const int _nC;
  /// Delta of the cuboid
  const T _deltaC;
  /// Number of data to be transfered
  const std::size_t _nData;
  /// Size of underlying data type
  const std::size_t _nDataType;

  /// Cuboid representing overlapping _iCglob
  Cuboid3D<T> _iCuboid;
  /// Block representing _iCuboid for local cell ID calculation
  BlockStructure3D _iCblock;

  /// Internal needed cells
  std::vector<Cell3D<T> > _inCells;
  /// External  needed cells
  std::vector<Cell3D<T> > _outCells;
  /// Internal needed neighbour cuboid and number of cells
  std::vector<int> _inC;
  std::vector<int> _inN;
  /// External needed neighbour cuboid and number of cells
  std::vector<int> _outC;
  std::vector<int> _outN;
  /// Buffer for the internal needed data
  std::uint8_t** _inData;
  std::size_t* _inDataSize;
  /// Buffer for the external needed data
  std::uint8_t** _outData;
  std::size_t* _outDataSize;
  /// Buffer for the internal needed data
  T **_inDataCoordinates;
  /// Buffer for the external needed data
  T **_outDataCoordinates;
  /// Indecates if there was an initialization done
  bool _initInCNdone;
  int* _tempInCN;
  /// Indecates if there was an initialization done
  bool _initOutCNdone;
  /// Handels the MPI communication
#ifdef PARALLEL_MODE_MPI
  singleton::MpiNonBlockingHelper _mpiNbHelper;
#endif
public:
  /// Constructor
  //CuboidNeighbourhood2D() {};
  /// Constructor
  CuboidNeighbourhood3D(SuperStructure3D<T>& superStructure, int iC);
  /// Copy construction
  CuboidNeighbourhood3D(CuboidNeighbourhood3D<T> const& rhs);
  /// Copy assignment
  CuboidNeighbourhood3D operator=(CuboidNeighbourhood3D<T> rhs);
  /// Destructor
  ~CuboidNeighbourhood3D();

  /// Read only access to _inCells
  Cell3D<T> const& get_inCell(int i) const;
  /// Returns the number of cells in _inCells
  int get_inCellsSize() const;
  /// Read only access to _inC
  int const& get_inC(int i) const;
  /// Returns the number of cells in _inC
  int get_inCsize() const;
  /// Read and write access to **_inData
  std::uint8_t** get_inData();
  /// Read and write access to **_outData
  std::uint8_t** get_outData();

  std::size_t getLocalCellId(int iX, int iY, int iZ) const;
  /// Adds a cell to the vector _inCells
  void add_inCell(Cell3D<T> cell);
  /// Adds a cell to the vector _outCells
  void add_outCell(Cell3D<T> cell);
  /// Adds a cell to the vector _inCells
  ///  if the cell is not already there and
  ///  if there is another cuboid which can deliver the information
  void add_inCell(int iX, int iY, int iZ);
  /// Adds all cells with the distance overlap*_delta to
  /// the vector _inCells
  void add_inCells(int overlap);
  /// Initializes _inC and _inN
  void init_inCN();
  /// Initializes _outC and _outN
  void init_outCN();
  /// Initialization Helper
#ifdef PARALLEL_MODE_MPI
  void bufSend_inCells(singleton::MpiNonBlockingHelper& helper);
#else
  void bufSend_inCells() {}
#endif
  /// Initialization Helper
  void recWrite_outCells();
  /// Finishes a communication step
  void finish_comm();

  /// Buffers data to be send
  // WARNING: Here is interpolation needed if globX, globY
  // are not integers. This needs to be fixed if one will
  // use unstructured grids.
  void buffer_outData();
  /// Sends data (pure mpi method)
  void send_outData();
  /// Receives data (pure mpi method)
  void receive_inData();
  /// Writes all data to the corresponding lattice cells
  void write_inData();
  /// Resets neighbourhood after initialization (init_inCN
  /// and init_outCN)
  void reset();
};

}  // namespace olb

#endif
