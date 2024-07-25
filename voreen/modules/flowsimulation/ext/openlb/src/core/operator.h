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

#ifndef CORE_OPERATOR_H
#define CORE_OPERATOR_H

namespace olb {

template <typename T, typename DESCRIPTOR, Platform PLATFORM> class ConcreteBlockLattice;
template <typename T, typename DESCRIPTOR, Platform PLATFORM> class ConcreteBlockMask;

/// Base of any block operator
template <typename T, typename DESCRIPTOR>
struct AbstractBlockO {
  virtual ~AbstractBlockO() = default;

  virtual std::type_index id() const = 0;
};

/// Base of block-wide operators such as post processors
template <typename T, typename DESCRIPTOR, Platform PLATFORM>
struct BlockO : public AbstractBlockO<T,DESCRIPTOR> {
  /// Set whether iCell is covered by the operator (optional)
  virtual void set(CellID iCell, bool state) = 0;
  /// Setup operator context
  virtual void setup(ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>& block) = 0;
  /// Apply operator on block
  virtual void apply(ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>& block) = 0;
};

/// Block-wide operator application scopes
/**
 * Declares how the actual OPERATOR::apply template wants to be called.
 **/
enum struct OperatorScope {
  /// Per-cell application, i.e. OPERATOR::apply is passed a CELL concept implementation
  PerCell,
  /// Per-block application, i.e. OPERATOR::apply is passed a ConcreteBlockLattice
  PerBlock,
  /// Per-cell application with parameters, i.e. OPERATOR::apply is passed a CELL concept implementation and parameters
  PerCellWithParameters,
};

/// Block application of concrete OPERATOR called using SCOPE on PLATFORM
template<typename T, typename DESCRIPTOR, Platform PLATFORM, typename OPERATOR, OperatorScope SCOPE> class ConcreteBlockO;

/// Base of collision operations performed by BlockDynamicsMap
template <typename T, typename DESCRIPTOR>
struct AbstractCollisionO : public AbstractBlockO<T,DESCRIPTOR> {
  virtual Dynamics<T,DESCRIPTOR>* getDynamics() = 0;

  /// Returns number of assigned cells
  /**
   * Used to determine the dominant dynamics to choose e.g. which
   * collision operator to vectorize or to prefer in GPU kernels.
   **/
  virtual std::size_t weight() const = 0;
  /// Set whether iCell is covered by the present collision step
  /**
   * \param iCell   Cell index
   * \param state   (De)activate for this dynamics / collision
   * \param overlap Cell index in overlap (set dynamics but do not collide)
   **/
  virtual void set(CellID iCell, bool state, bool overlap) = 0;
};


/// Collision dispatch strategy
enum struct CollisionDispatchStrategy {
  /// Apply dominant dynamics using mask and fallback to virtual dispatch for others
  Dominant,
  /// Apply all dynamics individually (async for Platform::GPU_CUDA)
  Individual
};

/// Collision operation on concrete blocks of PLATFORM
template <typename T, typename DESCRIPTOR, Platform PLATFORM>
struct BlockCollisionO : public AbstractCollisionO<T,DESCRIPTOR> {
  /// Setup collision on block
  virtual void setup(ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>& block) = 0;
  /// Apply collision on subdomain of block using strategy
  /**
   * Subdomain is currently assumed to be the non-overlap area of the block
   **/
  virtual void apply(ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>& block,
                     ConcreteBlockMask<T,DESCRIPTOR,PLATFORM>&    subdomain,
                     CollisionDispatchStrategy                    strategy) = 0;
};

/// Collision operation of concrete DYNAMICS on concrete block lattices of PLATFORM
template<typename T, typename DESCRIPTOR, Platform PLATFORM, typename DYNAMICS>
class ConcreteBlockCollisionO;


}

#endif
