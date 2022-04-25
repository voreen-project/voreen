/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2022 Nicolas Hafen, Mathias J. Krause
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

#ifndef PARTICLE_DESCRIPTOR_UTILITIES_H
#define PARTICLE_DESCRIPTOR_UTILITIES_H

namespace olb {

namespace descriptors {


/// Traversal of nested field contents for output and initialization
template <typename F, typename T, typename DESCRIPTOR, typename FIELDS>
struct access_field_content;

//Recursive retrieval of type names
template <typename F, typename T, typename DESCRIPTOR, typename ARRAYTYPE, typename HEAD, typename ...FIELDS>
void getFieldContent( F f, ARRAYTYPE& multiFieldArray, std::size_t iP )
{
  auto& field = multiFieldArray.template get<HEAD>();
  auto fieldContent = field.getRowPointer(iP);
  int fieldSize = fieldContent.getSize();
  std::stringstream fieldContentStream;   //Necessary until c++20
  fieldContentStream << fieldContent;     //Necessary until c++20
  if constexpr (!sizeof...(FIELDS)) {
    bool resetField = f( typeid(HEAD), fieldSize, fieldContentStream.str() );
    if (resetField) {
      for (int i=0; i<fieldSize; ++i) {
        auto singleEntry = fieldContent[i];
        bool isArithmetic = std::is_arithmetic<decltype(singleEntry)>::value;
        if (isArithmetic) {
          fieldContent[i] = 0;
        }
      }
    }
  }
  else {
    bool resetField = f( typeid(HEAD), fieldSize, fieldContentStream.str() );
    if (resetField) {
      for (int i=0; i<fieldSize; ++i) {
        auto singleEntry = fieldContent[i];
        bool isArithmetic = std::is_arithmetic<decltype(singleEntry)>::value;
        if (isArithmetic) {
          fieldContent[i] = 0;
        }
      }
    }
    getFieldContent<F,T,DESCRIPTOR,ARRAYTYPE,FIELDS...>( f, multiFieldArray, iP );
  }
}

//Recursive traversing and retrieval of nested type names (level two)
template <typename F, typename T, typename DESCRIPTOR, typename HEAD, typename ...FIELDS>
void getFieldContentL2( F f, DynamicFieldGroupsD<T, typename DESCRIPTOR::fieldList>& dynamicFieldGroups, std::size_t iP )
{
  if constexpr (!sizeof...(FIELDS)) {
    auto& fieldL2 = dynamicFieldGroups.template get<HEAD>();
    access_field_content<F,T,HEAD,typename HEAD::fieldList>::fields( f, fieldL2, iP );
  }
  else {
    auto& fieldL2 = dynamicFieldGroups.template get<HEAD>();
    access_field_content<F,T,HEAD,typename HEAD::fieldList>::fields( f, fieldL2, iP );
    getFieldContentL2<F,T, DESCRIPTOR, FIELDS...>( f, dynamicFieldGroups, iP );
  }
}

//Recursive traversing and retrieval of nested type names (level one)
template <typename F, typename T, typename DESCRIPTOR, typename... FIELDS>
struct access_field_content<F,T,DESCRIPTOR, meta::list<FIELDS...>> {
  using arrayType = MultiFieldArrayD<T,DESCRIPTOR,Platform::CPU_SISD,FIELDS...>; //ensure fixed type during recursion
  static constexpr void fields( F f, MultiFieldArrayD<T,DESCRIPTOR,Platform::CPU_SISD,FIELDS...>& multiFieldArray,
                                std::size_t iP )
  {
    return getFieldContent<F,T,DESCRIPTOR,arrayType,FIELDS...>( f, multiFieldArray, iP );
  }
  static constexpr void fieldsL2( F f, DynamicFieldGroupsD<T, typename DESCRIPTOR::fieldList>& dynamicFieldGroups,
                                  std::size_t iP )
  {
    return getFieldContentL2<F,T,DESCRIPTOR,FIELDS...>( f, dynamicFieldGroups, iP );
  }
};

} //namespace descriptors

} //namespace olb


#endif
