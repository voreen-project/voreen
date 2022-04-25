/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Nicolas Hafen, Mathias J. Krause
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

#ifndef CONTAINER_H
#define CONTAINER_H

namespace olb {

/**
 *  Container is a std::vector inspirated data wrapper that allows for simple
 *  content manipulation of its owned data. FieldArrayD, AnyFieldArrayD, MultiFieldArrayD 
 *  and DynamicFieldGroupsD do handle crucial data operations but do not provide
 *  simple access and manipulation functions as e.g. provided in std:: containers
 *  such as size and capacity handling. This Container class applicable to a provided
 *  FIELD_ARRAY_TYPE adds this functionality. 
 *  It is also intended to be a crucial part in the parallization refactoring of the
 *  particle system by replacing the currently used std::deque container.
 */

// *INDENT-OFF*


template<typename T, typename DESCRIPTOR, typename FIELD_ARRAY_TYPE>
class Container{
private:
  FIELD_ARRAY_TYPE _data;
  using size_type = std::size_t;
  /// Size
  size_type _size;
public:
  Container():
    _data(FIELD_ARRAY_TYPE(size_type(0))),_size(size_type(0)) { }

  Container( size_type count ):
    _data(FIELD_ARRAY_TYPE(count)),_size(count) { }
 

  //// Element access (according to std::vector)

  constexpr auto& at();                               //TODO: implement

  FIELD_ARRAY_TYPE& data(){ return _data; }

  //// Capacity (according to std::vector)

  constexpr bool empty();                             //TODO: implement

  constexpr size_type size(){ return _size;}

  constexpr size_type max_size();                     //TODO: implement

  constexpr void reserve( size_type new_capacity );   //TODO: implement

  constexpr size_type capacity(){ return _data.count(); }

  constexpr void shrink_to_fit(){ _data.resize(_size); }

  //// Modifiers (according to std::vector)
 
  constexpr void clear();                             //TODO: implement

  constexpr void push_back(){
    if ( this->capacity()==this->size() ){
      size_type new_capacity = 2*this->capacity();
      if (this->capacity()==0){ new_capacity=1; }
      this->resize( new_capacity );
    }
    _size++;
  }

  constexpr void pop_back(){ _size--; }

  void resize( size_type new_capacity ){ _data.resize(new_capacity); } 

  void swapElements(std::size_t i, std::size_t j)
  {
    _data.swap(i, j);
  }

  constexpr size_type getExtent(){ return _size; } //Analogy to lattice type (necessary for generic iteration)

};



// *INDENT-ON*


}

#endif
