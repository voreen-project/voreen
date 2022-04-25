/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006-2007 Jonas Latt
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
 * Definition of a LB cell -- generic implementation.
 */
#ifndef CELL_HH
#define CELL_HH

#include <algorithm>

#include "cell.h"
#include "util.h"
#include "dynamics/lbm.h"

#include "core/fieldArrayD.hh"

namespace olb {


template <typename T, typename DESCRIPTOR>
ConstCell<T,DESCRIPTOR>::ConstCell(const BlockLattice<T,DESCRIPTOR>& block,
                                   std::size_t iCell):
  _iCell(iCell),
  _block(const_cast<BlockLattice<T,DESCRIPTOR>&>(block)),
  _populations{_block.getPopulationPointers(iCell)}
{ }

template <typename T, typename DESCRIPTOR>
std::size_t ConstCell<T,DESCRIPTOR>::getCellId() const
{
  return _iCell;
}

template <typename T, typename DESCRIPTOR>
ConstCell<T,DESCRIPTOR>& ConstCell<T,DESCRIPTOR>::self() const
{
  return const_cast<ConstCell<T,DESCRIPTOR>&>(*this);
}

template <typename T, typename DESCRIPTOR>
template <typename FIELD>
auto ConstCell<T,DESCRIPTOR>::getFieldPointer(FIELD) const
{
  static_assert(descriptors::is_data_field<FIELD>::value,
                "FIELD must be structured data field");
  return _block.template getField<FIELD>().getPointer(_iCell);
}

template <typename T, typename DESCRIPTOR>
template <typename FIELD>
auto Cell<T,DESCRIPTOR>::getFieldPointer(FIELD)
{
  static_assert(descriptors::is_data_field<FIELD>::value,
                "FIELD must be structured data field");
  return this->_block.template getField<FIELD>().getPointer(this->_iCell);
}

template <typename T, typename DESCRIPTOR>
template <typename FIELD>
auto ConstCell<T,DESCRIPTOR>::getFieldComponent(unsigned iD) const
{
  static_assert(descriptors::is_data_field<FIELD>::value,
                "FIELD must be structured data field");
  return _block.template getField<FIELD>()[iD][_iCell];
}

template <typename T, typename DESCRIPTOR>
ConstCell<T,DESCRIPTOR> ConstCell<T,DESCRIPTOR>::neighbor(Vector<int,DESCRIPTOR::d> c) const
{
  return ConstCell<T,DESCRIPTOR>(_block, _iCell + _block.getNeighborDistance(c));
}

template <typename T, typename DESCRIPTOR>
Cell<T,DESCRIPTOR> Cell<T,DESCRIPTOR>::neighbor(Vector<int,DESCRIPTOR::d> c)
{
  return Cell<T,DESCRIPTOR>(this->_block, this->_iCell + this->_block.getNeighborDistance(c));
}

template <typename T, typename DESCRIPTOR>
T ConstCell<T,DESCRIPTOR>::operator[](unsigned iPop) const
{
  return *_populations[iPop];
}

template <typename T, typename DESCRIPTOR>
T& Cell<T,DESCRIPTOR>::operator[](unsigned iPop)
{
  return *this->_populations[iPop];
}

template <typename T, typename DESCRIPTOR>
bool ConstCell<T,DESCRIPTOR>::operator==(ConstCell<T,DESCRIPTOR>& rhs) const
{
  return getCellId() == rhs.getCellId();
}

template <typename T, typename DESCRIPTOR>
bool ConstCell<T,DESCRIPTOR>::operator!=(ConstCell<T,DESCRIPTOR>& rhs) const
{
  return getCellId() != rhs.getCellId();
}

template <typename T, typename DESCRIPTOR>
bool ConstCell<T,DESCRIPTOR>::operator<(ConstCell<T,DESCRIPTOR>& rhs) const
{
  return getCellId() < rhs.getCellId();
}

template <typename T, typename DESCRIPTOR>
bool ConstCell<T,DESCRIPTOR>::operator<=(ConstCell<T,DESCRIPTOR>& rhs) const
{
  return getCellId() <= rhs.getCellId();
}

template <typename T, typename DESCRIPTOR>
template <typename FIELD>
auto ConstCell<T,DESCRIPTOR>::getField(FIELD) const
{
  return _block.template getField<FIELD>().get(_iCell);
}

template <typename T, typename DESCRIPTOR>
template <typename FIELD>
void Cell<T,DESCRIPTOR>::setField(const FieldD<T,DESCRIPTOR,FIELD>& field)
{
  return this->_block.template getField<FIELD>().set(this->_iCell, field);
}

template <typename T, typename DESCRIPTOR>
template <typename FIELD>
std::enable_if_t<(DESCRIPTOR::template size<FIELD>() == 1), void>
Cell<T,DESCRIPTOR>::setField(typename FIELD::template value_type<T> value)
{
  return this->_block.template getField<FIELD>().set(this->_iCell,
                                                     FieldD<T,DESCRIPTOR,FIELD>(value));
}

template<typename T, typename DESCRIPTOR>
const Dynamics<T,DESCRIPTOR>* ConstCell<T,DESCRIPTOR>::getDynamics() const
{
  return _block.getDynamics(_iCell);
}

template<typename T, typename DESCRIPTOR>
Dynamics<T,DESCRIPTOR>* Cell<T,DESCRIPTOR>::getDynamics()
{
  return this->_block.getDynamics(this->_iCell);
}

template <typename T, typename DESCRIPTOR>
template <typename FIELD>
void ConstCell<T,DESCRIPTOR>::computeField(T* data) const
{
  auto field = getFieldPointer<FIELD>();
  for (long unsigned int i=0; i < DESCRIPTOR::template size<FIELD>(); ++i) {
    data[i] = field[i];
  }
}

template<typename T, typename DESCRIPTOR>
void Cell<T,DESCRIPTOR>::revert()
{
  for (int iPop=1; iPop <= descriptors::q<DESCRIPTOR>()/2; ++iPop) {
    std::swap(
      operator[](iPop),
      operator[](iPop+descriptors::q<DESCRIPTOR>()/2));
  }
}

template<typename T, typename DESCRIPTOR>
CellStatistic<T> Cell<T,DESCRIPTOR>::collide()
{
  return getDynamics()->collide(*this);
}

template<typename T, typename DESCRIPTOR>
T ConstCell<T,DESCRIPTOR>::computeRho() const
{
  return getDynamics()->computeRho(self());
}

template<typename T, typename DESCRIPTOR>
void ConstCell<T,DESCRIPTOR>::computeU(T u[descriptors::d<DESCRIPTOR>()]) const
{
  getDynamics()->computeU(self(), u);
}

template<typename T, typename DESCRIPTOR>
void ConstCell<T,DESCRIPTOR>::computeJ(T j[descriptors::d<DESCRIPTOR>()]) const
{
  getDynamics()->computeJ(self(), j);
}

template<typename T, typename DESCRIPTOR>
void ConstCell<T,DESCRIPTOR>::computeStress(T pi[util::TensorVal<DESCRIPTOR >::n]) const
{
  T rho, u[descriptors::d<DESCRIPTOR>()];
  getDynamics()->computeRhoU(self(), rho, u);
  getDynamics()->computeStress(self(), rho, u, pi);
}

template<typename T, typename DESCRIPTOR>
void ConstCell<T,DESCRIPTOR>::computeRhoU(T& rho, T u[descriptors::d<DESCRIPTOR>()]) const
{
  getDynamics()->computeRhoU(self(), rho, u);
}

template<typename T, typename DESCRIPTOR>
void ConstCell<T,DESCRIPTOR>::computeFeq(T fEq[descriptors::q<DESCRIPTOR>()]) const
{
  T rho{};
  Vector<T,descriptors::d<DESCRIPTOR>()> u;
  computeRhoU(rho, u.data());
  const T uSqr = norm_squared(u);
  for (int iPop=0; iPop < descriptors::q<DESCRIPTOR>(); ++iPop) {
    fEq[iPop] = lbm<DESCRIPTOR>::equilibrium(iPop, rho, u, uSqr);
  }
}

template <typename T, typename DESCRIPTOR>
void ConstCell<T,DESCRIPTOR>::computeFneq(T fNeq[descriptors::q<DESCRIPTOR>()]) const
{
  T rho{};
  T u[descriptors::d<DESCRIPTOR>()] { };
  computeRhoU(rho, u);
  lbm<DESCRIPTOR>::computeFneq(self(), fNeq, rho, u);
}

template<typename T, typename DESCRIPTOR>
void ConstCell<T,DESCRIPTOR>::computeAllMomenta(
  T& rho,
  T u[descriptors::d<DESCRIPTOR>()],
  T pi[util::TensorVal<DESCRIPTOR >::n]) const
{
  getDynamics()->computeAllMomenta(self(), rho, u, pi);
}

template<typename T, typename DESCRIPTOR>
void Cell<T,DESCRIPTOR>::defineRho(T rho)
{
  getDynamics()->defineRho(*this, rho);
}

template<typename T, typename DESCRIPTOR>
void Cell<T,DESCRIPTOR>::defineU(const T u[descriptors::d<DESCRIPTOR>()])
{
  getDynamics()->defineU(*this, u);
}

template<typename T, typename DESCRIPTOR>
void Cell<T,DESCRIPTOR>::defineRhoU(T rho, const T u[descriptors::d<DESCRIPTOR>()])
{
  getDynamics()->defineRhoU(*this, rho, u);
}

template<typename T, typename DESCRIPTOR>
void Cell<T,DESCRIPTOR>::definePopulations(const T* data)
{
  for (int iPop = 0; iPop < descriptors::q<DESCRIPTOR>(); ++iPop) {
    operator[](iPop) = data[iPop];
  }
}

template<typename T, typename DESCRIPTOR>
void Cell<T,DESCRIPTOR>::iniEquilibrium(T rho, const T u[descriptors::d<DESCRIPTOR>()])
{
  getDynamics()->iniEquilibrium(*this, rho, u);
}

template<typename T, typename DESCRIPTOR>
void Cell<T,DESCRIPTOR>::iniRegularized(
  T rho,
  const T u[descriptors::d<DESCRIPTOR>()],
  const T pi[util::TensorVal<DESCRIPTOR >::n])
{
  getDynamics()->iniRegularized(*this, rho, u, pi);
}


}

#endif
