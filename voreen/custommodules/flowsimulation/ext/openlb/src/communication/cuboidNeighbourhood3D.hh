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
 * The description of a 3D cuboid neighbourhood -- generic implementation.
 */


#ifndef CUBOID_NEIGHBOURHOOD_3D_HH
#define CUBOID_NEIGHBOURHOOD_3D_HH

#include "communication/mpiManager.h"
#include <vector>
#include <string>
#include "communication/cuboidNeighbourhood3D.h"
#include "geometry/cuboidGeometry3D.h"
#include "dynamics/dynamics.h"
#include "core/cell.h"
#include "communication/superStructure3D.h"


namespace olb {


//////////////// Class CuboidNeighbourhood3D //////////////////

template<typename T>
CuboidNeighbourhood3D<T>::CuboidNeighbourhood3D(SuperStructure3D<T>& superStructure, int iC):
  _superStructure(superStructure),
  _iCglob(iC),
  _nC(_superStructure.getCuboidGeometry().getNc()),
  _deltaC(superStructure.getCuboidGeometry().get(iC).getDeltaR()),
  _nData(_superStructure.getDataSize()),
  _nDataType(_superStructure.getDataTypeSize()),
  _iCuboid(_superStructure.getCuboidGeometry().get(_iCglob)),
  _iCblock(_iCuboid.getNx(), _iCuboid.getNy(), _iCuboid.getNz(), _superStructure.getOverlap())
{
  _initInCNdone  = false;
  _initOutCNdone = false;
}

template<typename T>
CuboidNeighbourhood3D<T>::CuboidNeighbourhood3D(CuboidNeighbourhood3D<T> const& rhs):
  _superStructure(rhs._superStructure),
  _iCglob(rhs._iCglob),
  _nC(rhs._nC),
  _deltaC(rhs._deltaC),
  _nData(rhs._nData),
  _nDataType(rhs._nDataType),
  _iCuboid(rhs._iCuboid),
  _iCblock(rhs._iCblock)
{
  _inCells       = rhs._inCells;
  _outCells      = rhs._outCells;
  _inC           = rhs._inC;
  _inN           = rhs._inN;
  _outC          = rhs._outC;
  _outN          = rhs._outN;

  _initInCNdone  = false;
  _initOutCNdone = false;
}

template<typename T>
CuboidNeighbourhood3D<T> CuboidNeighbourhood3D<T>::operator= (
  CuboidNeighbourhood3D<T> rhs )
{
  CuboidNeighbourhood3D<T> tmp(rhs);
  return tmp;
}

template<typename T>
CuboidNeighbourhood3D<T>::~CuboidNeighbourhood3D<T>()
{
  reset();
}


template<typename T>
Cell3D<T> const& CuboidNeighbourhood3D<T>::get_inCell(int i) const
{
  return _inCells[i];
}

template<typename T>
int CuboidNeighbourhood3D<T>::get_inCellsSize() const
{
  return _inCells.size();
}

template<typename T>
int const& CuboidNeighbourhood3D<T>::get_inC(int i) const
{
  return _inC[i];
}

template<typename T>
int CuboidNeighbourhood3D<T>::get_inCsize() const
{
  return _inC.size();
}

template<typename T>
std::uint8_t** CuboidNeighbourhood3D<T>::get_inData()
{
  return _inData;
}

template<typename T>
std::uint8_t** CuboidNeighbourhood3D<T>::get_outData()
{
  return _outData;
}


template<typename T>
std::size_t CuboidNeighbourhood3D<T>::getLocalCellId(int iX, int iY, int iZ) const
{
  const int overlap = _superStructure.getOverlap();
  return _iCblock.getCellId(iX+overlap, iY+overlap, iZ+overlap);
}

template<typename T>
void CuboidNeighbourhood3D<T>::add_inCell(Cell3D<T> cell)
{
  cell.latticeCellId = getLocalCellId(cell.latticeR[1], cell.latticeR[2], cell.latticeR[3]); 
  _inCells.emplace_back(cell);
}

template<typename T>
void CuboidNeighbourhood3D<T>::add_outCell(Cell3D<T> cell)
{
  cell.latticeCellId = getLocalCellId(cell.latticeR[1], cell.latticeR[2], cell.latticeR[3]); 
  _outCells.emplace_back(cell);
}

template<typename T>
void CuboidNeighbourhood3D<T>::add_inCell(int iX, int iY, int iZ)
{
  Cell3D<T> found;
  found.latticeR[0] = _iCglob;
  found.latticeR[1] = iX;
  found.latticeR[2] = iY;
  found.latticeR[3] = iZ;

  auto& cuboidGeometry = _superStructure.getCuboidGeometry();
  cuboidGeometry.getPhysR(found.physR, found.latticeR);

  std::vector<T> tmp(found.physR,found.physR + 3);
  if (cuboidGeometry.getC(tmp, found.latticeR[0]) ) {
    for (unsigned i=0; i < _inCells.size(); i++) {
      if (_inCells[i] == found) {
        return;
      }
    }
    add_inCell(found);
  }
}

template<typename T>
void CuboidNeighbourhood3D<T>::add_inCells(int overlap)
{
  auto& cuboidGeometry = _superStructure.getCuboidGeometry();
  const int nX  = cuboidGeometry.get(_iCglob).getNx();
  const int nY  = cuboidGeometry.get(_iCglob).getNy();
  const int nZ  = cuboidGeometry.get(_iCglob).getNz();

  for (int iX=0; iX < nX+2*overlap; ++iX) {
    for (int iY=0; iY < nY+2*overlap; ++iY) {
      for (int iZ=0; iZ < nZ+2*overlap; ++iZ) {
        if (iX < overlap || iX > nX + overlap - 1 ||
            iY < overlap || iY > nY + overlap - 1 ||
            iZ < overlap || iZ > nZ + overlap - 1) {
          Cell3D<T> found;
          found.latticeR[0] = _iCglob;
          found.latticeR[1] = iX - overlap;
          found.latticeR[2] = iY - overlap;
          found.latticeR[3] = iZ - overlap;

          _superStructure.getCuboidGeometry().getPhysR(found.physR, found.latticeR);
          std::vector<T> tmp(found.physR,found.physR + 3);
          if (_superStructure.getCuboidGeometry().getC(tmp, found.latticeR[0]) ) {
            add_inCell(found);
          }
        }
      }
    }
  }
}

template<typename T>
void CuboidNeighbourhood3D<T>::init_inCN()
{
  _inC.clear();
  _inN.clear();

  _inData = new std::uint8_t*[_nC];
  _inDataSize = new std::size_t[_nC];
  _inDataCoordinates = new T*[_nC];
  _tempInCN = new int[_nC];
  for (int i=0; i<_nC; i++) {
    _tempInCN[i]=0;
  }

  for (unsigned i=0; i<_inCells.size(); i++) {
    _tempInCN[_inCells[i].latticeR[0]]++;
  }
  for (int i=0; i<_nC; i++) {
    if (_tempInCN[i]!=0) {
      _inC.push_back(i);
      _inN.push_back(_tempInCN[i]);
      _inDataSize[i] = _tempInCN[i];
      _inData[i] = new std::uint8_t[_tempInCN[i]*_nData*_nDataType];
      _inDataCoordinates[i] = new T[_tempInCN[i]*3];
    }
    else {
      _inData[i] = NULL;
      _inDataCoordinates[i] = NULL;
    }
  }

#ifdef PARALLEL_MODE_MPI
  _mpiNbHelper.allocate(_nC);
  for (int i=0; i<_nC; i++) {
    int dRank = _superStructure.getLoadBalancer().rank(i);
    singleton::mpi().iSend(&_tempInCN[i] , 1, dRank, _mpiNbHelper.get_mpiRequest(i), _iCglob);
  }
#endif

  _initInCNdone = true;
}

template<typename T>
void CuboidNeighbourhood3D<T>::init_outCN()
{
  _outC.clear();
  _outN.clear();
  _outData = new std::uint8_t*[_nC];
  _outDataSize = new std::size_t[_nC];
  _outDataCoordinates = new T*[_nC];

  std::vector<int> temp(_nC,0);

  for (unsigned i=0; i<_outCells.size(); i++) {
    temp[_outCells[i].latticeR[0]]++;
  }

  for (int i=0; i < _nC; i++) {
#ifdef PARALLEL_MODE_MPI
    int sRank = _superStructure.getLoadBalancer().rank(i);
    singleton::mpi().receive(&temp[i], 1, sRank, i);
#endif
    if (temp[i] != 0) {
      _outC.push_back(i);
      _outN.push_back(temp[i]);
    }
    _outDataSize[i] = temp[i];
    _outData[i] = new std::uint8_t[temp[i]*_nData*_nDataType];
    _outDataCoordinates[i] = new T[temp[i]*3];
  }

  _initOutCNdone = true;
}

#ifdef PARALLEL_MODE_MPI
template<typename T>
void CuboidNeighbourhood3D<T>::bufSend_inCells(singleton::MpiNonBlockingHelper& helper)
{
  helper.free();

  std::vector<int> temp(_nC,0);
  for (unsigned i=0; i < _inCells.size(); i++) {
    int iC = _inCells[i].latticeR[0];
    if (singleton::mpi().getRank() != _superStructure.getLoadBalancer().rank(iC)) {
      _inDataCoordinates[iC][3*temp[iC]]   = _inCells[i].physR[0];
      _inDataCoordinates[iC][3*temp[iC]+1] = _inCells[i].physR[1];
      _inDataCoordinates[iC][3*temp[iC]+2] = _inCells[i].physR[2];
      temp[iC]++;
    }
  }

  helper.allocate(_inC.size());
  for (unsigned iC=0; iC<_inC.size(); iC++) {
    int dRank = _superStructure.getLoadBalancer().rank(_inC[iC]);
    singleton::mpi().iSend(_inDataCoordinates[_inC[iC]],
                           _inN[iC]*3, dRank, helper.get_mpiRequest(iC), _iCglob);
  }
}
#endif

template<typename T>
void CuboidNeighbourhood3D<T>::recWrite_outCells()
{
#ifdef PARALLEL_MODE_MPI
  for (unsigned iC=0; iC < _outC.size(); iC++) {
    int sRank = _superStructure.getLoadBalancer().rank(_outC[iC]);
    singleton::mpi().receive(_outDataCoordinates[_outC[iC]], _outN[iC]*3, sRank, _outC[iC]);
    if (singleton::mpi().getRank() != sRank) {
      Cell3D<T> found;
      for (int i=0; i < _outN[iC]; i++) {
        found.physR[0] = _outDataCoordinates[_outC[iC]][3*i];
        found.physR[1] = _outDataCoordinates[_outC[iC]][3*i+1];
        found.physR[2] = _outDataCoordinates[_outC[iC]][3*i+2];
        _superStructure.getCuboidGeometry().getLatticeR(found.latticeR, found.physR);
        found.latticeR[0] = _outC[iC];
        add_outCell(found);
      }
    }
  }
#endif
}

template<typename T>
void CuboidNeighbourhood3D<T>::finish_comm()
{
#ifdef PARALLEL_MODE_MPI
  singleton::mpi().waitAll(_mpiNbHelper);
#endif
}

template<typename T>
void CuboidNeighbourhood3D<T>::buffer_outData()
{
  const int iCloc = _superStructure.getLoadBalancer().loc(_iCglob);
  int iC = -1;
  for (unsigned iData=0; iData < _nData; ++iData) {
    std::vector<int> iCommCell(_nC,0);
    for (unsigned i=0; i < _outCells.size(); ++i, ++iCommCell[iC]) {
      iC = _outCells[i].latticeR[0];
      std::uint8_t* outPop = _outData[iC] + (iData*_outDataSize[iC] + iCommCell[iC])*_nDataType;
      std::uint8_t* locPop = _superStructure(iCloc, _outCells[i].latticeCellId, iData);
      std::memcpy(outPop, locPop, _nDataType);
    }
  }
}

template<typename T>
void CuboidNeighbourhood3D<T>::send_outData()
{
#ifdef PARALLEL_MODE_MPI
  for (unsigned iC=0; iC<_outC.size(); iC++) {
    int dRank = _superStructure.getLoadBalancer().rank(_outC[iC]);
   singleton::mpi().iSend(_outData[_outC[iC]],
                          _outN[iC]*_nData*_nDataType,
                          dRank, _mpiNbHelper.get_mpiRequest(iC), _iCglob);
  }
#endif
}

template<typename T>
void CuboidNeighbourhood3D<T>::receive_inData()
{
#ifdef PARALLEL_MODE_MPI
  for (unsigned iC=0; iC<_inC.size(); iC++) {
    int sRank = _superStructure.getLoadBalancer().rank(_inC[iC]);
    singleton::mpi().receive(_inData[_inC[iC]], _inN[iC]*_nData*_nDataType, sRank,_inC[iC]);
  }
#endif
}

template<typename T>
void CuboidNeighbourhood3D<T>::write_inData()
{
  const int iCloc = _superStructure.getLoadBalancer().loc(_iCglob);
  int iC = -1;
  for (unsigned iData=0; iData < _nData; ++iData) {
    std::vector<int> iCommCell(_nC,0);
    for (unsigned i=0; i < _inCells.size(); ++i, ++iCommCell[iC]) {
      iC = _inCells[i].latticeR[0];
      std::uint8_t* inPop  = _inData[iC] + (iData*_inDataSize[iC] + iCommCell[iC])*_nDataType;
      std::uint8_t* locPop = _superStructure(iCloc, _inCells[i].latticeCellId, iData);
      std::memcpy(locPop, inPop, _nDataType);
    }
  }
}

template<typename T>
void CuboidNeighbourhood3D<T>::reset()
{

  if (_initInCNdone) {
#ifdef PARALLEL_MODE_MPI
    for (int iC=0; iC<_nC; iC++) {
      delete[] _inData[iC];
      delete[] _inDataCoordinates[iC];
    }
#endif
    delete[] _inData;
    delete[] _inDataSize;
    delete[] _inDataCoordinates;
    delete[] _tempInCN;
    _initInCNdone = false;
  }
  if (_initOutCNdone) {
    for (int iC=0; iC<_nC; iC++) {
      delete[] _outData[iC];
      delete[] _outDataCoordinates[iC];
    }
    delete[] _outData;
    delete[] _outDataSize;
    delete[] _outDataCoordinates;
#ifdef PARALLEL_MODE_MPI
    _mpiNbHelper.free();
#endif
    _initOutCNdone = false;
  }
  _inCells.clear();
  _outCells.clear();
  _inC.clear();
  _outC.clear();
  _inN.clear();
  _outN.clear();
}

}  // namespace olb

#endif
