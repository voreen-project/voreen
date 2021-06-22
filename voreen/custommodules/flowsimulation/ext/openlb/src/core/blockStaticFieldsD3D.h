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

#ifndef BLOCK_STATIC_FIELDS_D_3D_H
#define BLOCK_STATIC_FIELDS_D_3D_H

#include "blockStructure3D.h"
#include "blockData3D.h"
#include "serializer.h"
#include "fieldArrayD.h"

namespace olb {


/// Storage of data for fields declared by DESCRIPTOR
template<typename T, typename DESCRIPTOR>
class BlockStaticFieldsD3D : public BlockStructure3D, public Serializable {
private:
  /// Type list of all non-tag descriptor field
  template<typename... FIELDS>
  using data_field_list = utilities::meta::field_list_filter_t<descriptors::is_data_field, FIELDS...>;

  template<typename... FIELDS>
  using CurriedMultiFieldArrayD = MultiFieldArrayD<T,DESCRIPTOR,FIELDS...>;

  /// Declare _data to be a MultiFieldArrayD containing all non-tag descriptor fields
  typename DESCRIPTOR::template decomposeInto<data_field_list>
                     ::template decomposeInto<CurriedMultiFieldArrayD> _data;

public:
  BlockStaticFieldsD3D(int nX, int nY, int nZ);
  BlockStaticFieldsD3D(): BlockStaticFieldsD3D(1,1,1) { };

  /// Return read-only reference to the array containing FIELD
  template <typename FIELD>
  const FieldArrayD<T,DESCRIPTOR,FIELD>& getFieldArray() const;
  /// Return reference to the array containing FIELD
  template <typename FIELD>
  FieldArrayD<T,DESCRIPTOR,FIELD>& getFieldArray();

  /// Return copy of FIELD data for cell iCell
  template <typename FIELD>
  FieldD<T,DESCRIPTOR,FIELD> getField(std::size_t iCell) const;
  /// Set FIELD data at cell iCell
  template <typename FIELD>
  void setField(std::size_t iCell, const FieldD<T,DESCRIPTOR,FIELD>& v);

  template <typename FIELD>
  const typename FIELD::template value_type<T>& getFieldComponent(std::size_t iCell, unsigned iDim) const;
  template <typename FIELD>
  typename FIELD::template value_type<T>& getFieldComponent(std::size_t iCell, unsigned iDim);

  template <typename FIELD>
  ConstFieldPtr<T,DESCRIPTOR,FIELD> getFieldPointer(std::size_t iCell) const {
    return getFieldArray<FIELD>().getFieldPointer(iCell);
  }
  template <typename FIELD>
  FieldPtr<T,DESCRIPTOR,FIELD> getFieldPointer(std::size_t iCell) {
    return getFieldArray<FIELD>().getFieldPointer(iCell);
  }

  template <typename F>
  void forFieldsAt(std::size_t iCell, F f) {
    this->_data.forFieldsAt(iCell, f);
  }

  /// Number of data blocks for the serializable interface
  std::size_t getNblock() const override;
  /// Binary size for the serializer
  std::size_t getSerializableSize() const override;
  /// Return a pointer to the memory of the current block and its size for the serializable interface
  bool* getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode) override;

};


}

#endif
