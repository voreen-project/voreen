/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2015 Mathias J. Krause, Benjamin Förster
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

/** \file
 * Dynamics for a generic 2D super data -- header file.
 */

#ifndef SUPER_DATA_H
#define SUPER_DATA_H

#include <vector>
#include <memory>

#include "blockData.h"
#include "utilities/aliases.h"
#include "communication/superStructure.h"
#include "communication/superCommunicator.h"
#include "dynamics/descriptorTag.h"


namespace olb {

template<unsigned D, typename T, typename U>
class SuperData : public SuperStructure<T,D> {
protected:
  /// Dimension of the data field
  const std::size_t _size;
  /// Vector of BlockData
  std::vector<std::unique_ptr<BlockData<D,T,U>>> _block;
  /// Inter-block communicator
  std::unique_ptr<SuperCommunicator<T,SuperData>> _communicator;

public:
  constexpr static unsigned d = D;

  using block_t = ConcretizableBlockData<D,T,U>;

  SuperData(CuboidGeometry<T,D>& cuboidGeometry,
            LoadBalancer<T>& loadBalancer,
            int overlap = 2,
            int size = 1);
  virtual ~SuperData() = default;
  SuperData(SuperF<D,T,U>& rhs);

  const BlockData<D,T,U>& getBlock(int iC) const;
  BlockData<D,T,U>& getBlock(int iC);

  template <typename BLOCK = BlockData<D,T,U>>
  BLOCK& getBlock(int iC);
  template <typename BLOCK = BlockData<D,T,U>>
  const BLOCK& getBlock(int iC) const;

  /// Communicate overlaps
  void communicate() override;

  /// Read only access to the dim of the data of the super structure
  int getDataSize() const;
  /// Read only access to the data type dim of the data of the super structure
  int getDataTypeSize() const;

};

}

#endif
