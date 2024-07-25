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

#ifndef CORE_PLATFORM_CPU_CELL_H
#define CORE_PLATFORM_CPU_CELL_H

namespace olb {

/// Implementations of CPU specifics
namespace cpu {

template <typename T, typename DESCRIPTOR, Platform PLATFORM>
class Cell;

/// Virtual interface for dynamically-dispatched dynamics access on CPU targets
template <typename T, typename DESCRIPTOR, Platform PLATFORM>
struct Dynamics {
  virtual CellStatistic<T> collide(Cell<T,DESCRIPTOR,PLATFORM>& cell) = 0;

  virtual T    computeRho (Cell<T,DESCRIPTOR,PLATFORM>& cell              ) = 0;
  virtual void computeU   (Cell<T,DESCRIPTOR,PLATFORM>& cell,         T* u) = 0;
  virtual void computeJ   (Cell<T,DESCRIPTOR,PLATFORM>& cell,         T* j) = 0;
  virtual void computeRhoU(Cell<T,DESCRIPTOR,PLATFORM>& cell, T& rho, T* u) = 0;

  virtual void computeStress    (Cell<T,DESCRIPTOR,PLATFORM>& cell, T& rho, T* u, T* pi) = 0;
  virtual void computeAllMomenta(Cell<T,DESCRIPTOR,PLATFORM>& cell, T& rho, T* u, T* pi) = 0;

  virtual T computeEquilibrium(int iPop, T rho, T* u) = 0;

  virtual T getOmegaOrFallback(T fallback) = 0;

  void iniEquilibrium(Cell<T,DESCRIPTOR,PLATFORM>& cell, T rho, T* u) {
    for (unsigned iPop=0; iPop < DESCRIPTOR::q; ++iPop) {
      cell[iPop] = computeEquilibrium(iPop, rho, u);
    }
  };

};

/// CPU specific field mirroring BlockDynamicsMap
template <typename T, typename DESCRIPTOR, Platform PLATFORM>
struct DYNAMICS : public descriptors::TYPED_FIELD_BASE<Dynamics<T,DESCRIPTOR,PLATFORM>*,1> { };

/// Cell concept for concrete block lattices on CPU platforms
/**
 * Used for generic operators and non-legacy post processors
 **/
template <typename T, typename DESCRIPTOR, Platform PLATFORM>
class Cell {
private:
  ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>& _lattice;
  std::size_t _iCell;

public:
  using value_t = T;
  using descriptor_t = DESCRIPTOR;

  Cell(ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>& lattice, std::size_t iCell=0):
    _lattice(lattice),
    _iCell(iCell) { }

  CellID getCellId() const {
    return _iCell;
  }

  void setCellId(std::size_t iCell) {
    _iCell = iCell;
  }

  T& operator[](unsigned iPop) {
    return _lattice.template getField<descriptors::POPULATION>()[iPop][_iCell];
  }

  template <typename FIELD>
  auto getField() const {
    return _lattice.template getField<FIELD>().getRow(_iCell);
  }

  template <typename FIELD>
  void setField(const FieldD<T,DESCRIPTOR,FIELD>& v) {
    return _lattice.template getField<FIELD>().setRow(_iCell, v);
  }

  template <typename FIELD>
  auto getFieldPointer() {
    return _lattice.template getField<FIELD>().getRowPointer(_iCell);
  }

  template <typename FIELD>
  auto& getFieldComponent(unsigned iD) {
    return _lattice.template getField<FIELD>()[iD][_iCell];
  }

  Cell<T,DESCRIPTOR,PLATFORM> neighbor(LatticeR<DESCRIPTOR::d> offset) {
    return {_lattice, _iCell + _lattice.getNeighborDistance(offset)};
  }

  Dynamics<T,DESCRIPTOR,PLATFORM>& getDynamics() {
    return *getFieldComponent<DYNAMICS<T,DESCRIPTOR,PLATFORM>>(0);
  }

  T computeRho() {
    return getDynamics().computeRho(*this);
  }
  void computeU(T* u) {
    getDynamics().computeU(*this, u);
  }
  void computeJ(T* j) {
    getDynamics().computeJ(*this, j);
  }
  void computeRhoU(T& rho, T* u) {
    getDynamics().computeRhoU(*this, rho, u);
  }
  void computeStress(T* pi) {
    T rho, u[DESCRIPTOR::d] { };
    getDynamics().computeRhoU(*this, rho, u);
    getDynamics().computeStress(*this, rho, u, pi);
  }
  void computeAllMomenta(T& rho, T* u, T* pi) {
    getDynamics().computeAllMomenta(*this, rho, u, pi);
  }

  void iniEquilibrium(T rho, T* u) {
    getDynamics().iniEquilibrium(*this, rho, u);
  }

};

}

}

#endif
