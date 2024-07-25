/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006-2007 Jonas Latt,
 *                2015-2019 Mathias J. Krause, Adrian Kummerlaender
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
 * Definition of a LB cell -- header file.
 */
#ifndef CELL_H
#define CELL_H

#include <type_traits>

#include "util.h"
#include "meta.h"
#include "utilities/aliases.h"

#include "fieldArrayD.h"

namespace olb {

template<typename T> struct CellStatistic;
template<typename T> class LatticeStatistics;
template<typename T, typename DESCRIPTOR> struct Dynamics;
template<typename T, typename DESCRIPTOR> class BlockLattice;
template<typename T, typename DESCRIPTOR> class Cell;

/// Highest-level interface to read-only Cell data
/**
 * Only use where access to the specific ConcreteBlockLattice
 * resp. Dynamics information is not (yet) possible.
 **/
template<typename T, typename DESCRIPTOR>
class ConstCell {
private:
  ConstCell<T,DESCRIPTOR>& self() const;

protected:
  friend Cell<T,DESCRIPTOR>;

  /// Memory ID of currently represented cell
  std::size_t _iCell;
  /// Underlying BlockLattice
  BlockLattice<T,DESCRIPTOR>& _block;
  /// Pointers to population values of current cell
  /**
   * Workaround to increase performance of post-processors,
   * post-processing functors until they are ported to
   * work on concrete lattices again (as is already the case
   * for dynamics)
   *
   * Validity of neighborhood accesses via these pointers not
   * guaranteed (e.g. consider CPU_SIMD vs. CPU_SISD)
   **/
  Vector<T*,DESCRIPTOR::q> _populations;

public:
  using value_t = T;
  using descriptor_t = DESCRIPTOR;

  ConstCell(const BlockLattice<T,DESCRIPTOR>& block,
            std::size_t iCell);

  /// Return memory ID of the currently represented cell
  std::size_t getCellId() const;

  /// Jump to arbitrary cell memory ID
  /**
   * Caller is responsible that this is valid.
   **/
  void setCellId(std::size_t iCell)
  {
    _iCell = iCell;
    _populations = _block.getPopulationPointers(iCell);
  }

  bool operator==(ConstCell<T,DESCRIPTOR>& rhs) const;
  bool operator!=(ConstCell<T,DESCRIPTOR>& rhs) const;
  bool operator< (ConstCell<T,DESCRIPTOR>& rhs) const;
  bool operator<=(ConstCell<T,DESCRIPTOR>& rhs) const;

  ConstCell<T,DESCRIPTOR> neighbor(Vector<int,DESCRIPTOR::d> c) const;

  /// Read-only access to distribution functions.
  /**
   * \param iPop index of the accessed distribution function
   **/
  T operator[](unsigned iPop) const;

  /// Return read-only field accessor
  template <typename FIELD>
  auto getFieldPointer(FIELD field = FIELD()) const;
  /// Return copy of field component
  template <typename FIELD>
  auto getFieldComponent(unsigned iD) const;
  /// Return copy of descriptor-declared FIELD as a vector
  template <typename FIELD>
  auto getField(FIELD field = FIELD()) const;

  /// Get a pointer to the dynamics
  const Dynamics<T,DESCRIPTOR>* getDynamics() const;

  /// Copy FIELD content to given memory location
  template <typename FIELD>
  void computeField(T* data) const;
  /// Compute particle density on the cell
  T computeRho() const;
  /// Compute fluid velocity on the cell
  void computeU(T u[descriptors::d<DESCRIPTOR>()]) const;
  /// Compute fluid momentum (j = rho * u) on the cell
  void computeJ(T j[descriptors::d<DESCRIPTOR>()]) const;
  /// Compute components of the stress tensor on the cell
  void computeStress(T pi[util::TensorVal<DESCRIPTOR >::n]) const;
  /// Compute fluid velocity and particle density on the cell
  void computeRhoU(T& rho, T u[descriptors::d<DESCRIPTOR>()]) const;
  /// Compute equilibrium part of cell distribution
  void computeFeq(T fEq[descriptors::q<DESCRIPTOR>()]) const;
  /// Compute non-equilibrium part of cell distribution
  void computeFneq(T fNeq[descriptors::q<DESCRIPTOR>()]) const;
  /// Compute all momenta on the celll, up to second order
  void computeAllMomenta(
    T& rho, T u[descriptors::d<DESCRIPTOR>()],
    T pi[util::TensorVal<DESCRIPTOR >::n] ) const;

};

/// Highest-level interface to Cell data
template<typename T, typename DESCRIPTOR>
class Cell : public ConstCell<T,DESCRIPTOR> {
public:
  Cell(BlockLattice<T,DESCRIPTOR>& block,
       std::size_t iCell):
    ConstCell<T,DESCRIPTOR>(block, iCell)
  { }

  Cell<T,DESCRIPTOR> neighbor(Vector<int,DESCRIPTOR::d> c);

  /// Read-write access to distribution functions.
  /**
   * \param iPop index of the accessed distribution function
   **/
  T& operator[](unsigned iPop);

  /// Return field accessor
  template <typename FIELD>
  auto getFieldPointer(FIELD field = FIELD());

  /// Set value of FIELD from a vector
  template <typename FIELD>
  void setField(const FieldD<T,DESCRIPTOR,FIELD>& field);
  /// Set value of FIELD from a scalar
  template <typename FIELD>
  std::enable_if_t<(DESCRIPTOR::template size<FIELD>() == 1), void>
  setField(typename FIELD::template value_type<T> value);

  /// Get a pointer to the dynamics
  Dynamics<T,DESCRIPTOR>* getDynamics();

  /// Revert ("bounce-back") the distribution functions
  void revert();

  /// Apply LB collision to the cell according to local dynamics
  CellStatistic<T> collide();

  /// Set density on the cell
  void defineRho(T rho);
  /// Set fluid velocity on the cell
  void defineU(const T u[descriptors::d<DESCRIPTOR>()]);
  /// Define fluid velocity and particle density on the cell.
  void defineRhoU(T rho, const T u[descriptors::d<DESCRIPTOR>()]);
  /// Define particle populations through the dynamics object
  void definePopulations(const T* f_);
  /// Initialize all f values to their local equilibrium
  void iniEquilibrium(T rho, const T u[descriptors::d<DESCRIPTOR>()]);
  /// Initialize all f values with local equilibrium and non equilibrium part
  void iniRegularized(T rho, const T u[descriptors::d<DESCRIPTOR>()], const T pi[util::TensorVal<DESCRIPTOR >::n]);

};

}

#endif
