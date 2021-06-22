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

#include "blockStaticPopulationD2D.h"
#include "blockStaticFieldsD2D.h"
#include "blockDynamicFieldsD2D.h"

#include "blockStaticPopulationD3D.h"
#include "blockStaticFieldsD3D.h"
#include "blockDynamicFieldsD3D.h"

#include "blockDynamicsMap.h"

#include "utilities/meta.h"
#include "utilities/aliases.h"

namespace olb {

template<typename T, typename DESCRIPTOR> class Cell;

template<typename T, typename DESCRIPTOR>
class ConstCell {
private:
  ConstCell<T,DESCRIPTOR>& self() const;

protected:
  friend Cell<T,DESCRIPTOR>;

  std::size_t _iCell;

  BlockStaticPopulationD<T,DESCRIPTOR>& _staticPopulationD;
  BlockStaticFieldsD<T,DESCRIPTOR>& _staticFieldsD;

  BlockDynamicFieldsD<T,DESCRIPTOR>& _dynamicFieldsD;
  BlockDynamicsMap<T,DESCRIPTOR>&    _dynamicsMap;

public:
  ConstCell(const BlockStaticPopulationD<T,DESCRIPTOR>& staticPopulationD,
            const BlockStaticFieldsD<T,DESCRIPTOR>&     staticFieldsD,
            const BlockDynamicFieldsD<T,DESCRIPTOR>&    dynamicFieldsD,
            const BlockDynamicsMap<T,DESCRIPTOR>&       dynamicsMap,
            std::size_t iCell);

  /// Return memory ID of the currently represented cell
  std::size_t getCellId() const;

  /// Jump to arbitrary cell memory ID
  /**
   * Caller is responsible that this is valid.
   **/
  void setCellId(std::size_t iCell) {
    _iCell = iCell;
  }

  /// Jump to next cell in linearization sequence
  /**
   * Caller is responsible that this is valid.
   **/
  void advanceCellId() {
    ++_iCell;
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
  ConstFieldPtr<T,DESCRIPTOR,FIELD> getFieldPointer() const;
  /// Helper for accessing fields in a generic lambda expression
  template <typename FIELD_ID>
  ConstFieldPtr<T,DESCRIPTOR,typename FIELD_ID::type>
  getFieldPointer(FIELD_ID id) const;

  /// Return copy of descriptor-declared FIELD as a vector
  template <typename FIELD>
  std::enable_if_t<
    (DESCRIPTOR::template size<FIELD>() > 1),
    FieldD<T,DESCRIPTOR,FIELD>
  >
  getField() const;
  /// Return copy of descriptor-declared FIELD as a scalar
  template <typename FIELD>
  std::enable_if_t<(DESCRIPTOR::template size<FIELD>() == 1), typename FIELD::template value_type<T>>
  getField() const;

  /// Return pointer to dynamic FIELD of cell
  template <typename FIELD>
  ConstFieldPtr<T,DESCRIPTOR,FIELD>
  getDynamicFieldPointer() const;

  /// Return copy of dynamic FIELD as a scalar
  template <typename FIELD>
  std::enable_if_t<(DESCRIPTOR::template size<FIELD>() == 1), typename FIELD::template value_type<T>>
  getDynamicField() const;
  /// Return copy of dynamic FIELD as a vector
  template <typename FIELD>
  std::enable_if_t<
    (DESCRIPTOR::template size<FIELD>() > 1),
    FieldD<T,DESCRIPTOR,FIELD>
  >
  getDynamicField() const;

  /// Get a pointer to the dynamics
  const Dynamics<T,DESCRIPTOR>* getDynamics() const;

  /// Return serialized size of the currently represented cell
  /**
   * This value may vary depending on the number of dynamic fields.
   **/
  unsigned getSerializedSize() const;
  /// Serialize cell data (population and fields)
  void serialize(T* data) const;

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

template<typename T, typename DESCRIPTOR>
class Cell : public ConstCell<T,DESCRIPTOR> {
public:
  Cell(BlockStaticPopulationD<T,DESCRIPTOR>& staticPopulationD,
       BlockStaticFieldsD<T,DESCRIPTOR>&     staticFieldsD,
       BlockDynamicFieldsD<T,DESCRIPTOR>&    dynamicFieldsD,
       BlockDynamicsMap<T,DESCRIPTOR>&       dynamicsMap,
       std::size_t iCell):
    ConstCell<T,DESCRIPTOR>(staticPopulationD, staticFieldsD, dynamicFieldsD, dynamicsMap, iCell)
  { }

  /// Override all values with those of rhs
  Cell<T,DESCRIPTOR>& operator=(ConstCell<T,DESCRIPTOR>& rhs);

  Cell<T,DESCRIPTOR> neighbor(Vector<int,DESCRIPTOR::d> c);

  /// Read-write access to distribution functions.
  /**
   * \param iPop index of the accessed distribution function
   **/
  T& operator[](unsigned iPop);

  /// Zero-initialize memory of population and all cell fields
  void init();

  /// Return field accessor
  template <typename FIELD>
  FieldPtr<T,DESCRIPTOR,FIELD> getFieldPointer();
  /// Helper for accessing fields in a generic lambda expression
  template <typename FIELD_ID>
  FieldPtr<T,DESCRIPTOR,typename FIELD_ID::type>
  getFieldPointer(FIELD_ID id);

  /// Set value of FIELD from a vector
  template <typename FIELD>
  std::enable_if_t<(DESCRIPTOR::template size<FIELD>() > 1), void>
  setField(const FieldD<T,DESCRIPTOR,FIELD>& field);
  /// Set value of FIELD from a scalar
  template <typename FIELD>
  std::enable_if_t<(DESCRIPTOR::template size<FIELD>() == 1), void>
  setField(typename FIELD::template value_type<T> value);

  /// Return pointer to dynamic FIELD of cell
  template <typename FIELD>
  FieldPtr<T,DESCRIPTOR,FIELD> getDynamicFieldPointer();

  /// Set value of dynamic FIELD from a vector
  template <typename FIELD>
  std::enable_if_t<(DESCRIPTOR::template size<FIELD>() > 1), void>
  setDynamicField(const FieldD<T,DESCRIPTOR,FIELD>& field);
  /// Set value of dynamic FIELD from a scalar
  template <typename FIELD>
  std::enable_if_t<(DESCRIPTOR::template size<FIELD>() == 1), void>
  setDynamicField(typename FIELD::template value_type<T> value);

  /// Update cell (population and fields) using serialized data
  void unSerialize(const T* data);

  /// Set FIELD value from given memory location
  template <typename FIELD>
  void defineField(const T* data);
  /// Add to FIELD from given memory location
  /**
   * Similar to defineField(),but instead of replacing existing values
   * the data is added onto the existing values.
   **/
  template <typename FIELD>
  void addField(const T* data);
  /// Multiply FIELD with values at given memory location
  /**
   * Similar to defineField(), but instead of replacing existing values
   * the data is multiplied to the existing values.
   **/
  template <typename FIELD>
  void multiplyField(const T* data);

  /// Define or re-define dynamics of the cell.
  /**
   * \param dynamics_ a pointer to the dynamics object, whos memory management
   *                  falls under the responsibility of the user
   **/
  void defineDynamics(Dynamics<T,DESCRIPTOR>* dynamics_);
  /// Get a pointer to the dynamics
  Dynamics<T,DESCRIPTOR>* getDynamics();

  /// Revert ("bounce-back") the distribution functions
  void revert();

  /// Apply LB collision to the cell according to local dynamics
  void collide(LatticeStatistics<T>& statistics);

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


template<typename T, typename DESCRIPTOR>
struct WriteCellFunctional {
  virtual ~WriteCellFunctional() { };
  virtual void apply(Cell<T,DESCRIPTOR>& cell, int pos[descriptors::d<DESCRIPTOR>()]) const =0;
};


}

#endif
