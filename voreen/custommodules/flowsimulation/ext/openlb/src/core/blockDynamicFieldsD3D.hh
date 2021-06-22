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

#ifndef BLOCK_DYNAMIC_FIELDS_D_3D_HH
#define BLOCK_DYNAMIC_FIELDS_D_3D_HH

#include "blockDynamicFieldsD3D.h"
#include "utilities/meta.h"

#include <stdexcept>
#include <tuple>

namespace olb {


template<typename T, typename DESCRIPTOR>
BlockDynamicFieldsD3D<T,DESCRIPTOR>::BlockDynamicFieldsD3D(int nX, int nY, int nZ):
  BlockStructure3D(nX, nY, nZ)
{ }

template <typename T, typename DESCRIPTOR>
template <typename FIELD>
const FieldArrayD<T,DESCRIPTOR,FIELD>& BlockDynamicFieldsD3D<T,DESCRIPTOR>::getFieldArray() const
{
  auto field = _fields.find(std::type_index(typeid(FIELD)));
  if (field == _fields.end()) {
    throw std::out_of_range("Dynamic FIELD doesn't exist as it was not set.");
  } else {
    return field->second.template as<FIELD>();
  }
}

template <typename T, typename DESCRIPTOR>
template <typename FIELD>
FieldArrayD<T,DESCRIPTOR,FIELD>& BlockDynamicFieldsD3D<T,DESCRIPTOR>::getFieldArray()
{
  auto field = _fields.find(std::type_index(typeid(FIELD)));
  if (field == _fields.end()) {
    auto result = _fields.emplace(std::piecewise_construct,
                                  std::forward_as_tuple(std::type_index(typeid(FIELD))),
                                  std::forward_as_tuple(utilities::meta::id<FIELD>(), this->getNcells()));
    return (*result.first).second.template as<FIELD>();
  } else {
    return field->second.template as<FIELD>();
  }
}

template <typename T, typename DESCRIPTOR>
template <typename FIELD>
FieldD<T,DESCRIPTOR,FIELD> BlockDynamicFieldsD3D<T,DESCRIPTOR>::getField(std::size_t iCell) const
{
  return FieldD<T,DESCRIPTOR,FIELD>([this,iCell](unsigned iDim) -> T {
    return getFieldArray<FIELD>()[iDim][iCell];
  });
}

template <typename T, typename DESCRIPTOR>
template <typename FIELD>
void BlockDynamicFieldsD3D<T,DESCRIPTOR>::setField(std::size_t iCell, const FieldD<T,DESCRIPTOR,FIELD>& v)
{
  auto& field = getFieldArray<FIELD>();
  for (unsigned iDim=0; iDim < DESCRIPTOR::template size<FIELD>(); ++iDim) {
    field[iDim][iCell] = v[iDim];
  }
}

template <typename T, typename DESCRIPTOR>
template <typename FIELD>
const typename FIELD::template value_type<T>&
BlockDynamicFieldsD3D<T,DESCRIPTOR>::getFieldComponent(std::size_t iCell, unsigned iDim) const
{
  return getFieldArray<FIELD>()[iDim][iCell];
}

template <typename T, typename DESCRIPTOR>
template <typename FIELD>
typename FIELD::template value_type<T>&
BlockDynamicFieldsD3D<T,DESCRIPTOR>::getFieldComponent(std::size_t iCell, unsigned iDim)
{
  return getFieldArray<FIELD>()[iDim][iCell];
}

template <typename T, typename DESCRIPTOR>
template <typename FIELD>
void BlockDynamicFieldsD3D<T,DESCRIPTOR>::setSerializationName(const std::string& name) {
  getFieldArray<FIELD>(); // ensure that FIELD exists
  _serialization_names[name] = &_fields.at(typeid(FIELD));
}

template<typename T, typename DESCRIPTOR>
std::size_t BlockDynamicFieldsD3D<T,DESCRIPTOR>::getNblock() const
{
  std::size_t nBlock = 3;
  for (const auto& pair : _serialization_names) {
    nBlock += std::get<1>(pair)->asSerializable().getNblock();
  }
  return nBlock;
}

template<typename T, typename DESCRIPTOR>
std::size_t BlockDynamicFieldsD3D<T,DESCRIPTOR>::getSerializableSize() const
{
  std::size_t size = 3 * sizeof(int);
  for (const auto& pair : _serialization_names) {
    size += std::get<1>(pair)->asSerializable().getSerializableSize();
  }
  return size;
}

template<typename T, typename DESCRIPTOR>
bool* BlockDynamicFieldsD3D<T,DESCRIPTOR>::getBlock(
  std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode)
{
  std::size_t currentBlock = 0;
  bool* dataPtr = nullptr;

  registerVar(iBlock, sizeBlock, currentBlock, dataPtr, this->_nx);
  registerVar(iBlock, sizeBlock, currentBlock, dataPtr, this->_ny);
  registerVar(iBlock, sizeBlock, currentBlock, dataPtr, this->_nz);
  for (auto& pair : _serialization_names) {
    auto& dynamicFieldSerializable = std::get<1>(pair)->asSerializable();
    registerSerializableOfConstSize(iBlock, sizeBlock, currentBlock, dataPtr, dynamicFieldSerializable, loadingMode);
  }

  return dataPtr;
}


}

#endif
