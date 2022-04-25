/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Lukas Baron, Mathias J. Krause, Albert Mink
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

#ifndef BLOCK_BASE_F_2D_HH
#define BLOCK_BASE_F_2D_HH

#include "blockBaseF2D.h"

namespace olb {


template <typename T>
BlockF2D<T>::BlockF2D(BlockStructureD<2>& blockStructure, int targetDim)
  : GenericF<T,int>(targetDim,2), _blockStructure(&blockStructure) { }

  //added from old
template <typename T>
BlockF2D<T>::BlockF2D(int targetDim)
  : GenericF<T,int>(targetDim,2), _blockStructure(nullptr) { }

template <typename T>
BlockStructureD<2>& BlockF2D<T>::getBlockStructure() //const
{
  return *_blockStructure;
}

  //added from old
template <typename T>
void BlockF2D<T>::setBlockStructure(BlockStructureD<2>* blockStructure)
{
  _blockStructure = blockStructure;
}

template <typename T,typename BaseType>
BlockDataF2D<T,BaseType>::BlockDataF2D(BlockData<2,T,BaseType>& blockData)
  : BlockF2D<T>(blockData, blockData.getSize()),
    _blockData(&blockData),
    _owning(false)
{ }

template <typename T,typename BaseType>
BlockDataF2D<T,BaseType>::BlockDataF2D(BlockF2D<BaseType>& f)
  : BlockF2D<T>(f.getBlockStructure(), f.getTargetDim()),
    _blockData(new BlockData<2,T,BaseType>(f)),
    _owning(true)
{ }

template <typename T,typename BaseType>
BlockDataF2D<T,BaseType>::BlockDataF2D(int nx, int ny, int size)
// hacky solution to both managing BlockData2D using std::unique_ptr and
// passing it down the line to the base class
  : BlockF2D<T>(*(new BlockData<2,T,BaseType>({{nx, ny}, 0}, size)), size),
    _blockData(static_cast<BlockData<2,T,BaseType>*>(&(this->getBlockStructure()))),
    _owning(true)
{ }

template <typename T,typename BaseType>
BlockDataF2D<T,BaseType>::~BlockDataF2D()
{
  if (_owning) {
    delete _blockData;
  }
}

template <typename T,typename BaseType>
BlockData<2,T,BaseType>& BlockDataF2D<T,BaseType>::getBlockData()
{
  return *_blockData;
}

template <typename T, typename BaseType>
bool BlockDataF2D<T,BaseType>::operator() (BaseType output[], const int input[])
{
  for (int iDim = 0; iDim < this->getTargetDim(); ++iDim) {
    output[iDim] = _blockData->get(input, iDim);
  }
  return true;
}


template <typename T>
BlockIdentity2D<T>::BlockIdentity2D(BlockF2D<T>& f)
  : BlockF2D<T>(f.getBlockStructure(),f.getTargetDim() ), _f(f)
{
  this->getName() = _f.getName();
  std::swap( _f._ptrCalcC, this->_ptrCalcC );
}

template <typename T>
bool BlockIdentity2D<T>::operator()(T output[], const int input[])
{
  return _f(output,input);
}


template <typename T>
BlockExtractComponentF2D<T>::BlockExtractComponentF2D(BlockF2D<T>& f, int extractDim)
  : BlockF2D<T>(f.getBlockStructure(),1 ), _f(f), _extractDim(extractDim)
{
  this->getName() = _f.getName();
}

template <typename T>
int BlockExtractComponentF2D<T>::getExtractDim()
{
  return _extractDim;
}

template <typename T>
bool BlockExtractComponentF2D<T>::operator()(T output[], const int input[])
{
  std::vector<T> outTmp(_f.getTargetDim(), T{});
  _f(outTmp.data(), input);
  output[0] = outTmp[_extractDim];
  return true;
}


template <typename T>
BlockExtractComponentIndicatorF2D<T>::BlockExtractComponentIndicatorF2D(
  BlockF2D<T>& f, int extractDim, BlockIndicatorF2D<T>& indicatorF)
  : BlockExtractComponentF2D<T>(f, extractDim),
    _indicatorF(indicatorF)
{
  this->getName() = f.getName();
}

template <typename T>
bool BlockExtractComponentIndicatorF2D<T>::operator()(T output[], const int input[])
{
  output[0] = T{};
  if (_indicatorF(input)) {
    return BlockExtractComponentF2D<T>::operator()(output, input);
  }
  return true;
}


template <typename T>
BlockExtractIndicatorF2D<T>::BlockExtractIndicatorF2D(
  BlockF2D<T>& f, BlockIndicatorF2D<T>& indicatorF)
  : BlockF2D<T>(f.getBlockStructure(), f.getTargetDim()),
    _f(f),
    _indicatorF(indicatorF)
{
  this->getName() = f.getName();
}

template <typename T>
bool BlockExtractIndicatorF2D<T>::operator()(T output[], const int input[])
{
  for (int i = 0; i < this->getTargetDim(); ++i) {
    output[i] = T{};
  }
  if (_indicatorF(input)) {
    _f(output, input);
  }
  return true;
}


template <typename T, typename DESCRIPTOR>
BlockLatticeF2D<T,DESCRIPTOR>::BlockLatticeF2D
(BlockLattice<T,DESCRIPTOR>& blockStructure, int targetDim)
  : BlockF2D<T>(blockStructure, targetDim), _blockLattice(blockStructure)
{ }
/*
template <typename T, typename DESCRIPTOR>
BlockLatticeF2D<T,DESCRIPTOR>::BlockLatticeF2D(BlockLatticeF2D<T,DESCRIPTOR> const& rhs)
  : BlockF2D<T>(rhs.getBlockStructure(), rhs.getTargetDim() ), _blockLattice(rhs.getBlock())
{ }

template <typename T, typename DESCRIPTOR>
BlockLatticeF2D<T,DESCRIPTOR>& BlockLatticeF2D<T,DESCRIPTOR>::operator=(BlockLatticeF2D<T,DESCRIPTOR> const& rhs)
{
  BlockLatticeF2D<T,DESCRIPTOR> tmp(rhs);
  return tmp;
}
*/
template <typename T, typename DESCRIPTOR>
BlockLattice<T,DESCRIPTOR>& BlockLatticeF2D<T, DESCRIPTOR>::getBlock()
{
  return _blockLattice;
}


template <typename T, typename DESCRIPTOR>
BlockLatticeIdentity2D<T,DESCRIPTOR>::BlockLatticeIdentity2D(
  BlockLatticeF2D<T,DESCRIPTOR>& f)
  : BlockLatticeF2D<T,DESCRIPTOR>(f.getBlock(),f.getTargetDim()),
    _f(f)
{
  this->getName() = _f.getName();
  std::swap( _f._ptrCalcC, this->_ptrCalcC );
}

template <typename T, typename DESCRIPTOR>
bool BlockLatticeIdentity2D<T,DESCRIPTOR>::operator()(T output[], const int input[])
{
  return _f(output,input);
}


template <typename T, typename DESCRIPTOR>
BlockLatticePhysF2D<T,DESCRIPTOR>::BlockLatticePhysF2D
(BlockLattice<T,DESCRIPTOR>& blockLattice, const UnitConverter<T,DESCRIPTOR>& converter, int targetDim)
  : BlockLatticeF2D<T,DESCRIPTOR>(blockLattice, targetDim), _converter(converter)
{ }

template <typename T, typename DESCRIPTOR, typename TDESCRIPTOR>
BlockLatticeThermalPhysF2D<T,DESCRIPTOR,TDESCRIPTOR>::BlockLatticeThermalPhysF2D
(BlockLattice<T,TDESCRIPTOR>& blockLattice, const ThermalUnitConverter<T,DESCRIPTOR,TDESCRIPTOR>& converter, int targetDim)
  : BlockLatticeF2D<T,TDESCRIPTOR>(blockLattice, targetDim), _converter(converter)
{ }




} // end namespace olb

#endif
