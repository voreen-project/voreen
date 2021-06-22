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

#ifndef BLOCK_STATIC_POPULATION_D_3D_H
#define BLOCK_STATIC_POPULATION_D_3D_H

#include "blockStructure3D.h"
#include "fieldArrayD.h"
#include "serializer.h"

#include "geometry/cuboidGeometry3D.h"
#include "communication/mpiManager.h"


namespace olb {


/// Storage of lattice populations with implicit propagation
/**
 * Projection between linear SoA storage and spatial embedding is
 * provided by BlockStructure3D.
 **/
template<typename T, typename DESCRIPTOR>
class BlockStaticPopulationD3D : public BlockStructure3D, public Serializable {
private:
  const std::size_t _size;

  FieldArrayD<T,DESCRIPTOR,descriptors::POPULATION> _data;

  std::ptrdiff_t _shift[DESCRIPTOR::q];
  std::size_t _remainder[DESCRIPTOR::q];

  T* _population[DESCRIPTOR::q][2];

public:
  BlockStaticPopulationD3D(int nX, int nY, int nZ);
  BlockStaticPopulationD3D() : BlockStaticPopulationD3D(1,1,1) { };

  T* getPopulationPointer(unsigned iPop, std::size_t iCell);
  const T* getPopulationPointer(unsigned iPop, std::size_t iCell) const;

  /// Update control structure to propagate populations
  void shift();
  /// Re-initialize control structure
  /**
   * Updates _population to reflect structure encoded by _shift
   **/
  void refreshControlStructure();

  /// Number of data blocks for the serializable interface
  std::size_t getNblock() const override;
  /// Binary size for the serializer
  std::size_t getSerializableSize() const override;
  /// Return a pointer to the memory of the current block and its size for the serializable interface
  bool* getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode) override;

};


}

#endif
