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

#ifndef BLOCK_DYNAMIC_FIELDS_D_2D_H
#define BLOCK_DYNAMIC_FIELDS_D_2D_H

#include <unordered_map>
#include <map>
#include <string>

#include "blockStructure2D.h"
#include "serializer.h"
#include "fieldArrayD.h"

namespace olb {


/// Runtime-allocated field storage
template<typename T, typename DESCRIPTOR>
class BlockDynamicFieldsD2D : public BlockStructure2D, public Serializable {
private:
  std::unordered_map<std::type_index, AnyFieldArrayD<T,DESCRIPTOR>> _fields;
  std::map<std::string, AnyFieldArrayD<T,DESCRIPTOR>*> _serialization_names;

public:
  BlockDynamicFieldsD2D(int nX, int nY);
  BlockDynamicFieldsD2D(): BlockDynamicFieldsD2D(1,1) { };

  /// Return read-only reference to the array containing FIELD
  template <typename FIELD>
  const FieldArrayD<T,DESCRIPTOR,FIELD>& getFieldArray() const;
  /// Return reference to the array containing FIELD
  template <typename FIELD>
  FieldArrayD<T,DESCRIPTOR,FIELD>& getFieldArray();

  template <typename FIELD>
  FieldD<T,DESCRIPTOR,FIELD> getField(std::size_t iCell) const;
  template <typename FIELD>
  void setField(std::size_t iCell, const FieldD<T,DESCRIPTOR,FIELD>& v);

  template <typename FIELD>
  const typename FIELD::template value_type<T>&
  getFieldComponent(std::size_t iCell, unsigned iDim) const;
  template <typename FIELD>
  typename FIELD::template value_type<T>&
  getFieldComponent(std::size_t iCell, unsigned iDim);

  template <typename FIELD>
  ConstFieldPtr<T,DESCRIPTOR,FIELD> getFieldPointer(std::size_t iCell) const {
    return getFieldArray<FIELD>().getFieldPointer(iCell);
  }
  template <typename FIELD>
  FieldPtr<T,DESCRIPTOR,FIELD> getFieldPointer(std::size_t iCell) {
    return getFieldArray<FIELD>().getFieldPointer(iCell);
  }

  /// Set unique name for FIELD to be included in serialization
  /**
   * The name by itself doesn't matter as long as it is unique with respect to the
   * instance of BlockDynamicFieldsD.  Only used to establish order between fields
   * in the serialization stream.
   *
   * MUST be called prior to Serializable::load / Serializable::save
   **/
  template <typename FIELD>
  void setSerializationName(const std::string& name);

  std::size_t getNblock() const override;
  std::size_t getSerializableSize() const override;
  bool* getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode) override;

};


}

#endif
