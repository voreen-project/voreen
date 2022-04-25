/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2015 Mathias J. Krause
 *                2021 Adrian Kummerlaender
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

#ifndef BLOCK_DATA_H
#define BLOCK_DATA_H

#include <vector>

#include "blockStructure.h"
#include "serializer.h"
#include "utilities/aliases.h"
#include "dynamics/descriptorTag.h"
#include "communication/communicatable.h"

#include "platform/column.h"

namespace olb {

template<unsigned D, typename T, typename U>
class BlockData : public BlockStructureD<D>
                , public Serializable {
protected:
  const unsigned _size;
  std::vector<cpu::sisd::Column<U>> _data;
  ConcreteCommunicatable<std::vector<cpu::sisd::Column<U>>> _communicatable;

public:
  static constexpr Platform platform = Platform::CPU_SISD;

  struct DUMMY_FIELD : public descriptors::DESCRIPTOR_TAG { };

  BlockData(Cuboid<T,D>& cuboid, int overlap=0, int size=1);
  BlockData(BlockStructureD<D>&& block, int size=1);
  BlockData(BlockF<U,D>& blockF);
  BlockData(BlockData<D,T,U>&&) = default;
  virtual ~BlockData() = default;

  bool operator() (T output[], const int input[]);

  Column<U>& getColumn(unsigned iD);

  U& get(std::size_t iCell, int iD=0);
  U& get(LatticeR<D> latticeR, int iD=0);
  U get(LatticeR<D> latticeR, int iD=0) const;

  Platform getPlatform() const {
    return platform;
  }

  bool hasCommunicatable(std::type_index field) const {
    return field == typeid(DUMMY_FIELD);
  }
  auto& getCommunicatable(std::type_index field) {
    OLB_ASSERT(field == typeid(DUMMY_FIELD),
               "BlockData only offers DUMMY_FIELD for communication");
    return _communicatable;
  }

  unsigned getSize() const;

  std::size_t getNblock() const override;
  std::size_t getSerializableSize() const override;
  bool* getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode) override;


};

/// Curried BlockData template for use in callUsingConcretePlatform
template<unsigned D, typename T, typename U>
struct ConcretizableBlockData {

using base_t = BlockData<D,T,U>;

template <Platform PLATFORM>
using type = BlockData<D,T,U>;

};

}

#endif
