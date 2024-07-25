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


template <typename T, typename DESCRIPTOR, Platform PLATFORM, typename... FIELDS>
template <typename FIELD>
const typename FIELD::template value_type<T>&
MultiFieldArrayD<T,DESCRIPTOR,PLATFORM,FIELDS...>::getFieldComponent(std::size_t iCell, unsigned iDim) const
{
  return get<FIELD>()[iDim][iCell];
}

template <typename T, typename DESCRIPTOR, Platform PLATFORM, typename... FIELDS>
template <typename FIELD>
typename FIELD::template value_type<T>&
MultiFieldArrayD<T,DESCRIPTOR,PLATFORM,FIELDS...>::getFieldComponent(std::size_t iCell, unsigned iDim)
{
  return get<FIELD>()[iDim][iCell];
}

template <typename T, typename DESCRIPTOR, Platform PLATFORM, typename... FIELDS>
template <typename F>
void MultiFieldArrayD<T,DESCRIPTOR,PLATFORM,FIELDS...>::forFields(F f) const
{
  fields_t::for_each([&](auto field) {
    f(get(field));
  });
}

template <typename T, typename DESCRIPTOR, Platform PLATFORM, typename... FIELDS>
template <typename F>
void MultiFieldArrayD<T,DESCRIPTOR,PLATFORM,FIELDS...>::forFields(F f)
{
  fields_t::for_each([&](auto field) {
    f(get(field));
  });
}

template <typename T, typename DESCRIPTOR, Platform PLATFORM, typename... FIELDS>
template <typename F>
void MultiFieldArrayD<T,DESCRIPTOR,PLATFORM,FIELDS...>::forFieldsAt(std::size_t idx, F f)
{
  fields_t::for_each([&](auto field) {
    f(get(field).getFieldPointer(idx), field);
  });
}

template <typename T, typename DESCRIPTOR, Platform PLATFORM, typename... FIELDS>
std::size_t MultiFieldArrayD<T,DESCRIPTOR,PLATFORM,FIELDS...>::getNblock() const
{
  return (get<FIELDS>().getNblock() + ... + 0);
}

template <typename T, typename DESCRIPTOR, Platform PLATFORM, typename... FIELDS>
std::size_t MultiFieldArrayD<T,DESCRIPTOR,PLATFORM,FIELDS...>::getSerializableSize() const
{
  return (get<FIELDS>().getSerializableSize() + ... + 0);
}

template <typename T, typename DESCRIPTOR, Platform PLATFORM, typename... FIELDS>
bool* MultiFieldArrayD<T,DESCRIPTOR,PLATFORM,FIELDS...>::getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode)
{
  std::size_t currentBlock = 0;
  bool* dataPtr = nullptr;

  meta::tuple_for_each(_static, [&](auto& field) {
    registerSerializableOfConstSize(iBlock, sizeBlock, currentBlock, dataPtr, field, loadingMode);
  });

  return dataPtr;
}


}

#endif
