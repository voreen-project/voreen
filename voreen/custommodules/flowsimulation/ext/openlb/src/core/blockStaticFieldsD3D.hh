/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2019 Adrian Kummerlaender
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

#ifndef BLOCK_STATIC_FIELDS_D_3D_HH
#define BLOCK_STATIC_FIELDS_D_3D_HH

#include "blockStaticFieldsD3D.h"

namespace olb {


template<typename T, typename DESCRIPTOR>
BlockStaticFieldsD3D<T,DESCRIPTOR>::BlockStaticFieldsD3D(int nX, int nY, int nZ):
  BlockStructure3D(nX, nY, nZ),
  _data(this->getNcells())
{ }

template <typename T, typename DESCRIPTOR>
template <typename FIELD>
const FieldArrayD<T,DESCRIPTOR,FIELD>& BlockStaticFieldsD3D<T,DESCRIPTOR>::getFieldArray() const
{
  return _data.template get<FIELD>();
}

template <typename T, typename DESCRIPTOR>
template <typename FIELD>
FieldArrayD<T,DESCRIPTOR,FIELD>& BlockStaticFieldsD3D<T,DESCRIPTOR>::getFieldArray()
{
  return _data.template get<FIELD>();
}

template <typename T, typename DESCRIPTOR>
template <typename FIELD>
FieldD<T,DESCRIPTOR,FIELD> BlockStaticFieldsD3D<T,DESCRIPTOR>::getField(std::size_t iCell) const
{
  return FieldD<T,DESCRIPTOR,FIELD>([this,iCell](unsigned iDim) -> T {
    return getFieldArray<FIELD>()[iDim][iCell];
  });
}

template <typename T, typename DESCRIPTOR>
template <typename FIELD>
void BlockStaticFieldsD3D<T,DESCRIPTOR>::setField(std::size_t iCell, const FieldD<T,DESCRIPTOR,FIELD>& v)
{
  auto& field = getFieldArray<FIELD>();
  for (unsigned iDim=0; iDim < DESCRIPTOR::template size<FIELD>(); ++iDim) {
    field[iDim][iCell] = v[iDim];
  }
}

template <typename T, typename DESCRIPTOR>
template <typename FIELD>
const typename FIELD::template value_type<T>&
BlockStaticFieldsD3D<T,DESCRIPTOR>::getFieldComponent(std::size_t iCell, unsigned iDim) const
{
  return getFieldArray<FIELD>()[iDim][iCell];
}

template <typename T, typename DESCRIPTOR>
template <typename FIELD>
typename FIELD::template value_type<T>&
BlockStaticFieldsD3D<T,DESCRIPTOR>::getFieldComponent(std::size_t iCell, unsigned iDim)
{
  return getFieldArray<FIELD>()[iDim][iCell];
}

template<typename T, typename DESCRIPTOR>
std::size_t BlockStaticFieldsD3D<T,DESCRIPTOR>::getNblock() const
{
  return 3 + _data.getNblock();
}

template<typename T, typename DESCRIPTOR>
std::size_t BlockStaticFieldsD3D<T,DESCRIPTOR>::getSerializableSize() const
{
  return 3 * sizeof(int) + _data.getSerializableSize();
}

template<typename T, typename DESCRIPTOR>
bool* BlockStaticFieldsD3D<T,DESCRIPTOR>::getBlock(
  std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode)
{
  std::size_t currentBlock = 0;
  bool* dataPtr = nullptr;

  registerVar(iBlock, sizeBlock, currentBlock, dataPtr, this->_nx);
  registerVar(iBlock, sizeBlock, currentBlock, dataPtr, this->_ny);
  registerVar(iBlock, sizeBlock, currentBlock, dataPtr, this->_nz);
  registerSerializableOfConstSize(iBlock, sizeBlock, currentBlock, dataPtr, _data, loadingMode);

  return dataPtr;
}


}

#endif
