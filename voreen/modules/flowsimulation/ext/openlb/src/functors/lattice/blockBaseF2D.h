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

#ifndef BLOCK_BASE_F_2D_H
#define BLOCK_BASE_F_2D_H

#include "functors/genericF.h"
#include "core/blockData.h"
#include "core/blockStructure.h"
#include "core/unitConverter.h"

#include <fstream>
#include <memory>

/* Note: Throughout the whole source code directory genericFunctions, the
 *  template parameters for i/o dimensions are:
 *           F: S^m -> T^n  (S=source, T=target)
 */

namespace olb {

template <typename T> class BlockIndicatorF2D;


/// represents all functors that operate on a cuboid in general, mother class of BlockLatticeF, ..
template <typename T>
class BlockF2D : public GenericF<T,int> {
protected:
  BlockF2D(BlockStructureD<2>& blockStructure, int targetDim);
  BlockF2D(int targetDim); //added from old
  BlockStructureD<2>* _blockStructure;  //* instead of &
public:
  /// virtual destructor for defined behaviour
  //~BlockF2D() override {};
  virtual BlockStructureD<2>& getBlockStructure(); // const;
  void setBlockStructure(BlockStructureD<2>* blockStructure);

  using GenericF<T,int>::operator();

  BlockF2D<T>& operator-(BlockF2D<T>& rhs);
  BlockF2D<T>& operator+(BlockF2D<T>& rhs);
  BlockF2D<T>& operator*(BlockF2D<T>& rhs);
  BlockF2D<T>& operator/(BlockF2D<T>& rhs);
};

/// BlockDataF2D can store data of any BlockFunctor2D
template <typename T, typename BaseType>
class BlockDataF2D : public BlockF2D<BaseType> {
protected:
  BlockDataF2D(int nx, int ny, int size=1);
  BlockData<2,T,BaseType>* _blockData;
  bool _owning;
public:
  /// Constructor
  BlockDataF2D(BlockData<2,T,BaseType>& blockData);
  ~BlockDataF2D();
  /// to store functor data, constuctor creates _blockData with functor data
  BlockDataF2D(BlockF2D<BaseType>& f);
  /// returns _blockData
  BlockData<2,T,BaseType>& getBlockData();
  /// access to _blockData via its get()
  bool operator() (BaseType output[], const int input[]) override;
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

/// functor to extract one component
template <typename T>
class BlockExtractComponentF2D : public BlockF2D<T> {
protected:
  BlockF2D<T>& _f;
  int _extractDim;
public:
  BlockExtractComponentF2D(BlockF2D<T>& f, int extractDim);
  int getExtractDim();
  bool operator() (T output[], const int input[]);
};

/// functor to extract one component inside an indicator
template <typename T>
class BlockExtractComponentIndicatorF2D : public BlockExtractComponentF2D<T> {
protected:
  BlockIndicatorF2D<T>& _indicatorF;
public:
  BlockExtractComponentIndicatorF2D(BlockF2D<T>& f, int extractDim,
                                    BlockIndicatorF2D<T>& indicatorF);
  bool operator() (T output[], const int input[]) override;
};

/// functor to extract data inside an indicator
template <typename T>
class BlockExtractIndicatorF2D : public BlockF2D<T> {
protected:
  BlockF2D<T>&          _f;
  BlockIndicatorF2D<T>& _indicatorF;
public:
  BlockExtractIndicatorF2D(BlockF2D<T>& f,
                           BlockIndicatorF2D<T>& indicatorF);
  bool operator() (T output[], const int input[]);
};

/// represents all functors that operate on a DESCRIPTOR in general, e.g. getVelocity(), getForce(), getPressure()
template <typename T, typename DESCRIPTOR>
class BlockLatticeF2D : public BlockF2D<T> {
protected:
  BlockLatticeF2D(BlockLattice<T,DESCRIPTOR>& blockLattice, int targetDim);
  BlockLattice<T,DESCRIPTOR>& _blockLattice;
public:
  /// Copy Constructor
  //BlockLatticeF2D(BlockLatticeF2D<T,DESCRIPTOR> const& rhs);
  /// Assignment Operator
  //BlockLatticeF2D<T,DESCRIPTOR>& operator=(BlockLatticeF2D<T,DESCRIPTOR> const& rhs);

  BlockLattice<T,DESCRIPTOR>& getBlock();
};


/// identity functor
template <typename T, typename DESCRIPTOR>
class BlockLatticeIdentity2D final : public BlockLatticeF2D<T,DESCRIPTOR> {
protected:
  BlockLatticeF2D<T,DESCRIPTOR>& _f;
public:
  BlockLatticeIdentity2D(BlockLatticeF2D<T,DESCRIPTOR>& f);
  bool operator() (T output[], const int input[]) override;
};


/// represents all functors that operate on a DESCRIPTOR with output in Phys, e.g. physVelocity(), physForce(), physPressure()
template <typename T, typename DESCRIPTOR>
class BlockLatticePhysF2D : public BlockLatticeF2D<T,DESCRIPTOR> {
protected:
  const UnitConverter<T,DESCRIPTOR>& _converter;
  BlockLatticePhysF2D(BlockLattice<T,DESCRIPTOR>& blockLattice,
                      const UnitConverter<T,DESCRIPTOR>& converter, int targetDim);
public:
  /// Copy Constructor
  //BlockLatticePhysF2D(BlockLatticePhysF2D<T,DESCRIPTOR> const& rhs);
  /// Assignment Operator
  //BlockLatticePhysF2D<T,DESCRIPTOR>& operator=(BlockLatticePhysF2D<T,DESCRIPTOR> const& rhs);
};

/// represents all thermal functors that operate on a DESCRIPTOR with output in Phys, e.g. physTemperature(), physHeatFlux()
template <typename T, typename DESCRIPTOR, typename TDESCRIPTOR>
class BlockLatticeThermalPhysF2D : public BlockLatticeF2D<T,TDESCRIPTOR> {
protected:
  BlockLatticeThermalPhysF2D(BlockLattice<T,TDESCRIPTOR>& blockLattice,
                             const ThermalUnitConverter<T,DESCRIPTOR,TDESCRIPTOR>& converter, int targetDim);
  const ThermalUnitConverter<T,DESCRIPTOR,TDESCRIPTOR>& _converter;
};


} // end namespace olb

#endif
