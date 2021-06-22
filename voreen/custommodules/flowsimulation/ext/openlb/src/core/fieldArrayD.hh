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

#ifndef FIELD_ARRAY_D_HH
#define FIELD_ARRAY_D_HH

#include "fieldArrayD.h"


namespace olb {


template<typename T, typename DESCRIPTOR, typename... FIELDS>
std::size_t MultiFieldArrayD<T,DESCRIPTOR,FIELDS...>::getNblock() const
{
  return descriptors::getFieldListSize<DESCRIPTOR::d,DESCRIPTOR::q,FIELDS...>();
}

template<typename T, typename DESCRIPTOR, typename... FIELDS>
std::size_t MultiFieldArrayD<T,DESCRIPTOR,FIELDS...>::getSerializableSize() const
{
  return _count * descriptors::getFieldListSize<DESCRIPTOR::d,DESCRIPTOR::q,FIELDS...>() * sizeof(T);
}

template<typename T, typename DESCRIPTOR, typename... FIELDS>
bool* MultiFieldArrayD<T,DESCRIPTOR,FIELDS...>::getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode)
{
  std::size_t currentBlock = 0;
  bool* dataPtr = nullptr;

  utilities::meta::tuple_for_each(_data, [&](auto& field) {
    this->registerSerializableOfConstSize(iBlock, sizeBlock, currentBlock, dataPtr, field, loadingMode);
  });

  return dataPtr;
}


}

#endif
