/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Adrian Kummerlaender
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

#ifndef BLOCK_LATTICE_GHOST_HH
#define BLOCK_LATTICE_GHOST_HH

#include "blockLatticeGhost.h"

namespace olb {


template<typename T, typename DESCRIPTOR>
BlockLatticeGhost<T,DESCRIPTOR>::BlockLatticeGhost(Vector<int,DESCRIPTOR::d> size, int padding)
  : BlockStructure<DESCRIPTOR>(size, padding),
    _unstructuredFieldsD()
{ }

template<typename T, typename DESCRIPTOR>
template<typename FIELD>
auto& BlockLatticeGhost<T,DESCRIPTOR>::getField(meta::id<FIELD>)
{
  static_assert(descriptors::is_unstructured_field<FIELD>::value,
                "Ghost block only provides unstructured fields");
  return *_unstructuredFieldsD.template get<FIELD>();
}

template<typename T, typename DESCRIPTOR>
template<typename FIELD>
const auto& BlockLatticeGhost<T,DESCRIPTOR>::getField(meta::id<FIELD>) const
{
  static_assert(descriptors::is_unstructured_field<FIELD>::value,
                "Ghost block only provides unstructured fields");
  return *_unstructuredFieldsD.template get<FIELD>();
}

template<typename T, typename DESCRIPTOR>
std::size_t BlockLatticeGhost<T,DESCRIPTOR>::getNblock() const
{
  return 1
         + _unstructuredFieldsD.getNblock();
}

template<typename T, typename DESCRIPTOR>
std::size_t BlockLatticeGhost<T,DESCRIPTOR>::getSerializableSize() const
{
  return sizeof(BlockStructure<DESCRIPTOR>)
         + _unstructuredFieldsD.getSerializableSize();
}

template<typename T, typename DESCRIPTOR>
bool* BlockLatticeGhost<T,DESCRIPTOR>::getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode)
{
  std::size_t currentBlock = 0;
  bool* dataPtr = nullptr;

  registerVar(iBlock, sizeBlock, currentBlock, dataPtr, static_cast<BlockStructure<DESCRIPTOR>&>(*this));
  registerSerializableOfConstSize (iBlock, sizeBlock, currentBlock, dataPtr, _unstructuredFieldsD, loadingMode);

  return dataPtr;
}

template<typename T, typename DESCRIPTOR>
std::size_t BlockLatticeGhost<T,DESCRIPTOR>::communicatableSize(const std::vector<std::type_index>& fields,
                                                                const std::vector<CellID>& cells) const
{
  return _unstructuredFieldsD.communicatableSize(fields, cells);
}

template<typename T, typename DESCRIPTOR>
std::size_t BlockLatticeGhost<T,DESCRIPTOR>::communicatableSerialize(const std::vector<std::type_index>& fields,
                                                                     const std::vector<CellID>& cells,
                                                                     std::uint8_t* buffer) const
{
  return _unstructuredFieldsD.communicatableSerialize(fields, cells, buffer);
}

template<typename T, typename DESCRIPTOR>
std::size_t BlockLatticeGhost<T,DESCRIPTOR>::communicatableDeserialize(const std::vector<std::type_index>& fields,
                                                                       const std::vector<CellID>& cells,
                                                                       const std::uint8_t* buffer)
{
  return _unstructuredFieldsD.communicatableDeserialize(fields, cells, buffer);
}


}

#endif
