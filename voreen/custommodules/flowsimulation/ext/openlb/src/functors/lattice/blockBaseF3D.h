/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Lukas Baron, Mathias J. Krause, Albert Mink
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

#ifndef BLOCK_BASE_F_3D_H
#define BLOCK_BASE_F_3D_H

#include "functors/genericF.h"
#include "core/blockData3D.h"
#include "core/blockStructure3D.h"
#include "core/blockLatticeStructure3D.h"
#include "core/unitConverter.h"

/** Note: Throughout the whole source code directory genericFunctions, the
 *  template parameters for i/o dimensions are:
 *           F: S^m -> T^n  (S=source, T=target)
 */

namespace olb {


/// represents all functors that operate on a cuboid in general, mother class of BlockLatticeF, ..
template <typename T>
class BlockF3D : public GenericF<T,int> {
protected:
  BlockF3D(BlockStructure3D& blockStructure, int targetDim);
  BlockStructure3D& _blockStructure;
public:
  /// virtual destructor for defined behaviour
  ~BlockF3D() override {};
  virtual BlockStructure3D& getBlockStructure() const;

  std::vector<T> getMinValue();
  std::vector<T> getMaxValue();

  BlockF3D<T>& operator-(BlockF3D<T>& rhs);
  BlockF3D<T>& operator+(BlockF3D<T>& rhs);
  BlockF3D<T>& operator*(BlockF3D<T>& rhs);
  BlockF3D<T>& operator/(BlockF3D<T>& rhs);
};

/// BlockDataF3D can store data of any BlockFunctor3D
template <typename T,typename BaseType>
class BlockDataF3D : public BlockF3D<BaseType> {
protected:
  BlockDataF3D(int nx, int ny, int nz, int size=1);
  BlockData3D<T,BaseType>& _blockData;
public:
  /// Constructor
  BlockDataF3D(BlockData3D<T,BaseType>& blockData);
  /// to store functor data, constuctor creates _blockData with functor data
  BlockDataF3D(BlockF3D<BaseType>& f);
  /// destructor is called if object was not created by passing a blockData
  ~BlockDataF3D() override;
  /// returns _blockData
  BlockData3D<T,BaseType>& getBlockData();
  /// access to _blockData via its get()
  bool operator() (BaseType output[], const int input[]) override;
private:
  /// flag whether _blockData was allocated with new
  bool _isConstructed;
};

/// identity functor
template <typename T>
class BlockIdentity3D final : public BlockF3D<T> {
protected:
  BlockF3D<T>& _f;
public:
  BlockIdentity3D(BlockF3D<T>& f);
  /// Assignment Operator
  //BlockIdentity3D<T> operator=(BlockF3D<T>& rhs);
  // access operator should not delete f, since f still has the identity as child
  bool operator() (T output[], const int input[]) override;
};

/// functor to extract one component
template <typename T>
class BlockExtractComponentF3D : public BlockF3D<T> {
protected:
  BlockF3D<T>& _f;
  int _extractDim;
public:
  BlockExtractComponentF3D(BlockF3D<T>& f, int extractDim);
  int getExtractDim();
  bool operator() (T output[], const int input[]);
};


/// functor to extract one component
template <typename T>
class BlockDotProductF3D : public BlockF3D<T> {
protected:
  BlockF3D<T>& _f;
  T _vector[];
public:
  BlockDotProductF3D(BlockF3D<T>& f, T vector[]);
  bool operator() (T output[], const int input[]);
};


/// represents all functors that operate on a Lattice in general, e.g. getVelocity(), getForce(), getPressure()
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticeF3D : public BlockF3D<T> {
protected:
  BlockLatticeF3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, int targetDim);
  BlockLatticeStructure3D<T,DESCRIPTOR>& _blockLattice;
public:
  /// Copy Constructor
  //BlockLatticeF3D(BlockLatticeF3D<T,DESCRIPTOR> const& rhs);
  /// Assignment Operator
  //BlockLatticeF3D<T,DESCRIPTOR>& operator=(BlockLatticeF3D<T,DESCRIPTOR> const& rhs);

  BlockLatticeStructure3D<T,DESCRIPTOR>& getBlockLattice();
};

/// represents all functors that operate on a Lattice with output in Phys, e.g. physVelocity(), physForce(), physPressure()
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticePhysF3D : public BlockLatticeF3D<T,DESCRIPTOR> {
protected:
  const UnitConverter<T,DESCRIPTOR>& _converter;
  BlockLatticePhysF3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
                      const UnitConverter<T,DESCRIPTOR>& converter, int targetDim);
public:
  /// Copy Constructor
  //BlockLatticePhysF3D(BlockLatticePhysF3D<T,DESCRIPTOR> const& rhs);
  /// Assignment Operator
  //BlockLatticePhysF3D<T,DESCRIPTOR>& operator=(BlockLatticePhysF3D<T,DESCRIPTOR> const& rhs);
};

/// represents all thermal functors that operate on a Lattice with output in Phys, e.g. physTemperature(), physHeatFlux()
template <typename T, template <typename U> class DESCRIPTOR, template <typename V> class ThermalDESCRIPTOR>
class BlockLatticeThermalPhysF3D : public BlockLatticeF3D<T,ThermalDESCRIPTOR> {
protected:
  BlockLatticeThermalPhysF3D(BlockLatticeStructure3D<T,ThermalDESCRIPTOR>& blockLattice,
                             const ThermalUnitConverter<T,DESCRIPTOR,ThermalDESCRIPTOR>& converter, int targetDim);
  const ThermalUnitConverter<T,DESCRIPTOR,ThermalDESCRIPTOR>& _converter;
};


} // end namespace olb

#endif
