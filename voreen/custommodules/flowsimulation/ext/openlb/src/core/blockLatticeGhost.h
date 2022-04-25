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

#ifndef BLOCK_LATTICE_GHOST_H
#define BLOCK_LATTICE_GHOST_H

#include "utilities/aliases.h"

#include "blockStructure.h"
#include "unstructuredFieldsD.h"
#include "serializer.h"

namespace olb {


template<typename T, typename DESCRIPTOR>
class BlockLatticeGhost final : public BlockStructure<DESCRIPTOR>
                              , public FieldsCommunicatable
                              , public Serializable {
private:
  /// Unstructured field data
  UnstructuredFieldsD<T,DESCRIPTOR> _unstructuredFieldsD;

public:
  BlockLatticeGhost(Vector<int,DESCRIPTOR::d> size, int padding=0);
  ~BlockLatticeGhost() = default;

  BlockLatticeGhost(const BlockLatticeGhost& rhs) = delete;
  BlockLatticeGhost(BlockLatticeGhost&&) = delete;
  BlockLatticeGhost& operator=(BlockLatticeGhost<T,DESCRIPTOR>&&) = delete;
  BlockLatticeGhost& operator=(const BlockLatticeGhost<T,DESCRIPTOR>&) = delete;

  template<typename FIELD>
  auto& getField(meta::id<FIELD> field = meta::id<FIELD>());
  template<typename FIELD>
  const auto& getField(meta::id<FIELD> field = meta::id<FIELD>()) const;

  std::size_t getNblock() const override;
  std::size_t getSerializableSize() const override;
  bool* getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode) override;

private:
  /// Get serialized size of `fields` at indices `cells`
  std::size_t communicatableSize(const std::vector<std::type_index>& fields,
                                 const std::vector<CellID>& cells) const override;
  /// Serialize data of `fields` at indices `cells` to `buffer`
  std::size_t communicatableSerialize(const std::vector<std::type_index>& fields,
                                      const std::vector<CellID>& cells,
                                      std::uint8_t* buffer) const override;
  /// Deserialize data of `fields` at indices `cells` to `buffer`
  std::size_t communicatableDeserialize(const std::vector<std::type_index>& fields,
                                const std::vector<CellID>& cells,
                                const std::uint8_t* buffer) override;

};


}

#endif
