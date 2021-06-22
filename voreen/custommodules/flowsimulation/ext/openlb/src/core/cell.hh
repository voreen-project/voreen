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
#include "dynamics/lbHelpers.h"
#include "core/blockStaticFieldsD2D.hh"
#include "core/blockStaticFieldsD3D.hh"

namespace olb {


template <typename T, typename DESCRIPTOR>
ConstCell<T,DESCRIPTOR>::ConstCell(
  const BlockStaticPopulationD<T,DESCRIPTOR>& staticPopulationD,
  const BlockStaticFieldsD<T,DESCRIPTOR>&     staticFieldsD,
  const BlockDynamicFieldsD<T,DESCRIPTOR>&    dynamicFieldsD,
  const BlockDynamicsMap<T,DESCRIPTOR>&       dynamicsMap,
  std::size_t iCell
):
  _iCell(iCell),
  _staticPopulationD(const_cast<BlockStaticPopulationD<T,DESCRIPTOR>&>(staticPopulationD)),
  _staticFieldsD(const_cast<BlockStaticFieldsD<T,DESCRIPTOR>&>(staticFieldsD)),
  _dynamicFieldsD(const_cast<BlockDynamicFieldsD<T,DESCRIPTOR>&>(dynamicFieldsD)),
  _dynamicsMap(const_cast<BlockDynamicsMap<T,DESCRIPTOR>&>(dynamicsMap))
{ }

template <typename T, typename DESCRIPTOR>
std::size_t ConstCell<T,DESCRIPTOR>::getCellId() const
{
  return _iCell;
}

template <typename T, typename DESCRIPTOR>
unsigned ConstCell<T,DESCRIPTOR>::getSerializedSize() const
{
  return DESCRIPTOR::q + DESCRIPTOR::size();
}

template <typename T, typename DESCRIPTOR>
ConstCell<T,DESCRIPTOR>& ConstCell<T,DESCRIPTOR>::self() const
{
  return const_cast<ConstCell<T,DESCRIPTOR>&>(*this);
}

template <typename T, typename DESCRIPTOR>
template <typename FIELD>
ConstFieldPtr<T,DESCRIPTOR,FIELD> ConstCell<T,DESCRIPTOR>::getFieldPointer() const
{
  return const_cast<const BlockStaticFieldsD<T,DESCRIPTOR>&>(_staticFieldsD).template getFieldPointer<FIELD>(_iCell);
}

template <typename T, typename DESCRIPTOR>
template <typename FIELD_ID>
ConstFieldPtr<T,DESCRIPTOR,typename FIELD_ID::type>
ConstCell<T,DESCRIPTOR>::getFieldPointer(FIELD_ID id) const
{
  return getFieldPointer<typename FIELD_ID::type>();
}

template <typename T, typename DESCRIPTOR>
template <typename FIELD>
FieldPtr<T,DESCRIPTOR,FIELD> Cell<T,DESCRIPTOR>::getFieldPointer()
{
  return this->_staticFieldsD.template getFieldPointer<FIELD>(this->_iCell);
}

template <typename T, typename DESCRIPTOR>
template <typename FIELD_ID>
FieldPtr<T,DESCRIPTOR,typename FIELD_ID::type>
Cell<T,DESCRIPTOR>::getFieldPointer(FIELD_ID id)
{
  return getFieldPointer<typename FIELD_ID::type>();
}

template <typename T, typename DESCRIPTOR>
ConstCell<T,DESCRIPTOR> ConstCell<T,DESCRIPTOR>::neighbor(Vector<int,DESCRIPTOR::d> c) const
{
  return ConstCell<T,DESCRIPTOR>(
    _staticPopulationD,
    _staticFieldsD,
    _dynamicFieldsD,
    _dynamicsMap,
    _iCell + _staticPopulationD.getNeighborDistance(c));
}

template <typename T, typename DESCRIPTOR>
Cell<T,DESCRIPTOR> Cell<T,DESCRIPTOR>::neighbor(Vector<int,DESCRIPTOR::d> c)
{
  return Cell<T,DESCRIPTOR>(
    this->_staticPopulationD,
    this->_staticFieldsD,
    this->_dynamicFieldsD,
    this->_dynamicsMap,
    this->_iCell + this->_staticPopulationD.getNeighborDistance(c));
}

template <typename T, typename DESCRIPTOR>
T ConstCell<T,DESCRIPTOR>::operator[](unsigned iPop) const
{
  OLB_PRECONDITION(iPop < descriptors::q<DESCRIPTOR>());
  return *_staticPopulationD.getPopulationPointer(iPop, _iCell);
}

template <typename T, typename DESCRIPTOR>
T& Cell<T,DESCRIPTOR>::operator[](unsigned iPop)
{
  OLB_PRECONDITION(iPop < descriptors::q<DESCRIPTOR>());
  return *this->_staticPopulationD.getPopulationPointer(iPop, this->_iCell);
}


template <typename T, typename DESCRIPTOR>
Cell<T,DESCRIPTOR>& Cell<T,DESCRIPTOR>::operator=(ConstCell<T,DESCRIPTOR>& rhs)
{
  for (unsigned iPop=0; iPop < DESCRIPTOR::template size<descriptors::POPULATION>(); ++iPop) {
    operator[](iPop) = rhs[iPop];
  }
  this->_staticFieldsD.forFieldsAt(this->_iCell, [&rhs](auto field, auto id) {
    field = rhs.getFieldPointer(id);
  });
  // TODO: reimplement dynamic field copy, potentially also using a forFieldsAt-like construct
  //this->getDynamicFields() = rhs.getDynamicFields();
  this->defineDynamics(const_cast<Dynamics<T,DESCRIPTOR>*>(rhs.getDynamics()));
  return *this;
}

template <typename T, typename DESCRIPTOR>
void Cell<T,DESCRIPTOR>::init()
{
  for (unsigned iPop=0; iPop < DESCRIPTOR::template size<descriptors::POPULATION>(); ++iPop) {
    operator[](iPop) = T();
  }
  this->_staticFieldsD.forFieldsAt(this->_iCell, [](auto field, auto id) {
    for (unsigned iDim=0; iDim < decltype(field)::d; ++iDim) {
      field[iDim] = T();
    }
  });
}

template <typename T, typename DESCRIPTOR>
void ConstCell<T,DESCRIPTOR>::serialize(T* data) const
{
  for (unsigned iPop = 0; iPop < DESCRIPTOR::template size<descriptors::POPULATION>(); ++iPop) {
    data[iPop] = operator[](iPop);
  }
  T* currData = data + DESCRIPTOR::template size<descriptors::POPULATION>();
  this->_staticFieldsD.forFieldsAt(this->_iCell, [&currData](auto field, auto id) {
    for (unsigned iDim=0; iDim < decltype(field)::d; ++iDim) {
      *(currData++) = field[iDim];
    }
  });
}

template <typename T, typename DESCRIPTOR>
void Cell<T,DESCRIPTOR>::unSerialize(const T* data)
{
  for (unsigned iPop = 0; iPop < DESCRIPTOR::template size<descriptors::POPULATION>(); ++iPop) {
    operator[](iPop) = data[iPop];
  }
  const T* currData = data + DESCRIPTOR::template size<descriptors::POPULATION>();
  this->_staticFieldsD.forFieldsAt(this->_iCell, [&currData](auto field, auto id) {
    for (unsigned iDim=0; iDim < decltype(field)::d; ++iDim) {
      field[iDim] = *(currData++);
    }
  });
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
std::enable_if_t<
  (DESCRIPTOR::template size<FIELD>() > 1),
  FieldD<T,DESCRIPTOR,FIELD>
>
ConstCell<T,DESCRIPTOR>::getField() const
{
  return _staticFieldsD.template getField<FIELD>(_iCell);
}

template <typename T, typename DESCRIPTOR>
template <typename FIELD>
std::enable_if_t<(DESCRIPTOR::template size<FIELD>() == 1), typename FIELD::template value_type<T>>
ConstCell<T,DESCRIPTOR>::getField() const
{
  return _staticFieldsD.template getField<FIELD>(_iCell)[0];
}

template <typename T, typename DESCRIPTOR>
template <typename FIELD>
std::enable_if_t<
  (DESCRIPTOR::template size<FIELD>() > 1),
  FieldD<T,DESCRIPTOR,FIELD>
>
ConstCell<T,DESCRIPTOR>::getDynamicField() const
{
  return _dynamicFieldsD.template getField<FIELD>(_iCell);
}

template <typename T, typename DESCRIPTOR>
template <typename FIELD>
std::enable_if_t<(DESCRIPTOR::template size<FIELD>() == 1), typename FIELD::template value_type<T>>
ConstCell<T,DESCRIPTOR>::getDynamicField() const
{
  return _dynamicFieldsD.template getField<FIELD>(_iCell)[0];
}

template <typename T, typename DESCRIPTOR>
template <typename FIELD>
std::enable_if_t<(DESCRIPTOR::template size<FIELD>() > 1), void>
Cell<T,DESCRIPTOR>::setField(const FieldD<T,DESCRIPTOR,FIELD>& field)
{
  this->_staticFieldsD.template setField<FIELD>(this->_iCell, field);
}

template <typename T, typename DESCRIPTOR>
template <typename FIELD>
std::enable_if_t<(DESCRIPTOR::template size<FIELD>() == 1), void>
Cell<T,DESCRIPTOR>::setField(typename FIELD::template value_type<T> value)
{
  this->_staticFieldsD.template setField<FIELD>(this->_iCell, FieldD<T,DESCRIPTOR,FIELD>(value));
}

template <typename T, typename DESCRIPTOR>
template <typename FIELD>
ConstFieldPtr<T,DESCRIPTOR,FIELD> ConstCell<T,DESCRIPTOR>::getDynamicFieldPointer() const
{
  return _dynamicFieldsD.template getFieldPointer<FIELD>(_iCell);
}

template <typename T, typename DESCRIPTOR>
template <typename FIELD>
FieldPtr<T,DESCRIPTOR,FIELD> Cell<T,DESCRIPTOR>::getDynamicFieldPointer()
{
  return this->_dynamicFieldsD.template getFieldPointer<FIELD>(this->getCellId());
}

template <typename T, typename DESCRIPTOR>
template <typename FIELD>
std::enable_if_t<(DESCRIPTOR::template size<FIELD>() > 1), void>
Cell<T,DESCRIPTOR>::setDynamicField(const FieldD<T,DESCRIPTOR,FIELD>& field)
{
  this->_dynamicFieldsD.template setField<FIELD>(this->_iCell, field);
}

template <typename T, typename DESCRIPTOR>
template <typename FIELD>
std::enable_if_t<(DESCRIPTOR::template size<FIELD>() == 1), void>
Cell<T,DESCRIPTOR>::setDynamicField(typename FIELD::template value_type<T> value)
{
  this->_dynamicFieldsD.template setField<FIELD>(this->_iCell, FieldD<T,DESCRIPTOR,FIELD>(value));
}

template<typename T, typename DESCRIPTOR>
void Cell<T,DESCRIPTOR>::defineDynamics(Dynamics<T,DESCRIPTOR>* dynamics)
{
  this->_dynamicsMap.set(this->getCellId(), dynamics);
}

template<typename T, typename DESCRIPTOR>
const Dynamics<T,DESCRIPTOR>* ConstCell<T,DESCRIPTOR>::getDynamics() const
{
  return &_dynamicsMap.get(this->getCellId());
}

template<typename T, typename DESCRIPTOR>
Dynamics<T,DESCRIPTOR>* Cell<T,DESCRIPTOR>::getDynamics()
{
  return &this->_dynamicsMap.get(this->getCellId());
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

template <typename T, typename DESCRIPTOR>
template <typename FIELD>
void Cell<T,DESCRIPTOR>::defineField(const T* data)
{
  this->_staticFieldsD.template setField<FIELD>(this->_iCell, Vector<T,DESCRIPTOR::template size<FIELD>()>(data));
}

template <typename T, typename DESCRIPTOR>
template <typename FIELD>
void Cell<T,DESCRIPTOR>::addField(const T* data)
{
  auto field = getFieldPointer<FIELD>();
  for (unsigned i=0; i < DESCRIPTOR::template size<FIELD>(); ++i) {
    field[i] += data[i];
  }
}

template <typename T, typename DESCRIPTOR>
template <typename FIELD>
void Cell<T,DESCRIPTOR>::multiplyField(const T* data)
{
  auto field = getFieldPointer<FIELD>();
  for (unsigned i=0; i < DESCRIPTOR::template size<FIELD>(); ++i) {
    field[i] *= data[i];
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
void Cell<T,DESCRIPTOR>::collide(LatticeStatistics<T>& statistics)
{
  getDynamics()->collide(*this, statistics);
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
    fEq[iPop] = lbHelpers<T,DESCRIPTOR>::equilibrium(iPop, rho, u.data(), uSqr);
  }
}

template <typename T, typename DESCRIPTOR>
void ConstCell<T,DESCRIPTOR>::computeFneq(T fNeq[descriptors::q<DESCRIPTOR>()]) const
{
  T rho{};
  T u[descriptors::d<DESCRIPTOR>()] { };
  computeRhoU(rho, u);
  lbHelpers<T,DESCRIPTOR>::computeFneq(self(), fNeq, rho, u);
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
