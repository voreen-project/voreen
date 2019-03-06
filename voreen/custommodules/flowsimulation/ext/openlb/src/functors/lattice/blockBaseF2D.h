/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2013 Albert Mink, Lukas Baron, Mathias J. Krause
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


#ifndef BLOCK_BASE_F_2D_H
#define BLOCK_BASE_F_2D_H

#include "functors/genericF.h"
#include "core/blockStructure2D.h"
#include "core/blockData2D.h"
#include "core/blockLatticeStructure2D.h"
#include "core/unitConverter.h"


/** Note: Throughout the whole source code directory genericFunctions, the
 *  template parameters for i/o dimensions are:
 *           F: S^m -> T^n  (S=source, T=target)
 */

namespace olb {

/// represents all functors that operate on a cuboid in general, mother class of BlockLatticeF, ..
template <typename T>
class BlockF2D : public GenericF<T,int> {
protected:
  BlockF2D(BlockStructure2D& blockStructure, int targetDim);
  BlockF2D(int targetDim);
  BlockStructure2D* _blockStructure;
public:
  // virtual due to BlockReduction2D1D
  virtual BlockStructure2D& getBlockStructure();
  void setBlockStructure(BlockStructure2D* blockStructure);

  // not used anymore? blockData2D computes min/max as well
  /// computes min/maxValue of blockStructure, [iDim]
  //  std::vector<T> getMinValue();
  //  std::vector<T> getMaxValue();

  BlockF2D<T>& operator-(BlockF2D<T>& rhs);
  BlockF2D<T>& operator+(BlockF2D<T>& rhs);
  BlockF2D<T>& operator*(BlockF2D<T>& rhs);
  BlockF2D<T>& operator/(BlockF2D<T>& rhs);
};

/// BlockDataF2D can store data of any BlockFunctor2D
template <typename T,typename BaseType>
class BlockDataF2D : public BlockF2D<T> {
protected:
  BlockDataF2D(BlockData2D<T,BaseType>& blockData);
  /// to store functor data, constuctor creates _blockData with functor data
  BlockDataF2D(BlockF2D<BaseType>& f);
  /// used for BlockReduction2D1D to build BlockData2D by the constructor
  BlockDataF2D(int nx, int ny, int size=1);
  BlockData2D<T,BaseType>& _blockData;
public:
  /// destructor is called if object was not created by passing a blockData
  ~BlockDataF2D() override;
  /// returns _blockData
  BlockData2D<T,BaseType>& getBlockData();
  /// access to _blockData via its get()
  bool operator() (T output[], const int input[]) override;
private:
  /// flag whether _blockData was allocated with new
  bool _isConstructed;
};

/// identity functor
template <typename T>
class BlockIdentity2D final : public BlockF2D<T> {
protected:
  BlockF2D<T>& _f;
public:
  BlockIdentity2D(BlockF2D<T>& f);
  // access operator should not delete f, since f still has the identity as child
  bool operator() (T output[], const int input[]) override;
};

/// represents all functors that operate on a Lattice in general, e.g. getVelocity(), getForce(), getPressure()
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticeF2D : public BlockF2D<T> {
protected:
  BlockLatticeF2D(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice, int targetDim);
  BlockLatticeStructure2D<T,DESCRIPTOR>& _blockLattice;
public:
  BlockLatticeStructure2D<T,DESCRIPTOR>& getBlockLattice();
};

/// represents all functors that operate on a Lattice with output in Phys, e.g. physVelocity(), physForce(), physPressure()
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticePhysF2D : public BlockLatticeF2D<T,DESCRIPTOR> {
protected:
  BlockLatticePhysF2D(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice,
                      const UnitConverter<T,DESCRIPTOR>& converter, int targetDim);
  const UnitConverter<T,DESCRIPTOR>& _converter;
};


/// represents all thermal functors that operate on a Lattice with output in Phys, e.g. physTemperature(), physHeatFlux()
template <typename T, template <typename U> class DESCRIPTOR, template <typename V> class ThermalDESCRIPTOR>
class BlockLatticeThermalPhysF2D : public BlockLatticeF2D<T,ThermalDESCRIPTOR> {
protected:
  BlockLatticeThermalPhysF2D(BlockLatticeStructure2D<T,ThermalDESCRIPTOR>& blockLattice,
                             const ThermalUnitConverter<T,DESCRIPTOR,ThermalDESCRIPTOR>& converter, int targetDim);
  const ThermalUnitConverter<T,DESCRIPTOR,ThermalDESCRIPTOR>& _converter;
};


} // end namespace olb

#endif
