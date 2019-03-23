/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2008 Orestis Malaspinas, Andrea Parmigiani
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
 * A collection of dynamics classes (e.g. BGK) with which a Cell object
 * can be instantiated -- generic implementation.
 */
#ifndef ADVECTION_DIFFUSION_MOMENTA_HH
#define ADVECTION_DIFFUSION_MOMENTA_HH

#include <algorithm>
#include <limits>
#include "dynamics/dynamics.h"
#include "core/cell.h"
#include "lbHelpers.h"

namespace olb {

////////////////////// Class AdvectionDiffusionBulkMomenta //////////////////////////

template<typename T, template<typename U> class Lattice>
T AdvectionDiffusionBulkMomenta<T,Lattice>::computeRho(Cell<T,Lattice> const& cell) const
{
  return lbHelpers<T,Lattice>::computeRho(cell);
}

template<typename T, template<typename U> class Lattice>
void AdvectionDiffusionBulkMomenta<T,Lattice>::computeU(Cell<T,Lattice> const& cell, T u[Lattice<T>::d]) const
{
  const T *u_ = cell.getExternal(Lattice<T>::ExternalField::velocityBeginsAt);
  for (int iD = 0; iD < Lattice<T>::d; ++iD) {
    u[iD] = u_[iD];
  }
}

template<typename T, template<typename U> class Lattice>
void AdvectionDiffusionBulkMomenta<T,Lattice>::computeJ(Cell<T,Lattice> const& cell, T j[Lattice<T>::d]) const
{
}

template<typename T, template<typename U> class Lattice>
void AdvectionDiffusionBulkMomenta<T,Lattice>::computeStress (
  Cell<T,Lattice> const& cell,
  T rho, const T u[Lattice<T>::d],
  T pi[util::TensorVal<Lattice<T> >::n] ) const
{
}

template<typename T, template<typename U> class Lattice>
void AdvectionDiffusionBulkMomenta<T,Lattice>::computeRhoU (
  Cell<T,Lattice> const& cell,
  T& rho, T u[Lattice<T>::d] ) const
{
  rho = cell.computeRho();
  const T *u_ = cell.getExternal(Lattice<T>::ExternalField::velocityBeginsAt);
  for (int iD = 0; iD < Lattice<T>::d; ++iD) {
    u[iD] = u_[iD];
  }
}

template<typename T, template<typename U> class Lattice>
void AdvectionDiffusionBulkMomenta<T,Lattice>::computeAllMomenta (
  Cell<T,Lattice> const& cell,
  T& rho, T u[Lattice<T>::d],
  T pi[util::TensorVal<Lattice<T> >::n] ) const
{
  rho = cell.computeRho();
  const T *u_ = cell.getExternal(Lattice<T>::ExternalField::velocityBeginsAt);
  for (int iD = 0; iD < Lattice<T>::d; ++iD) {
    u[iD] = u_[iD];
  }
}

template<typename T, template<typename U> class Lattice>
void AdvectionDiffusionBulkMomenta<T,Lattice>::defineRho(Cell<T,Lattice>& cell, T rho)
{
  T *u = cell.getExternal(Lattice<T>::ExternalField::velocityBeginsAt);
  for (int iPop=0; iPop < Lattice<T>::q; ++iPop) {
    cell[iPop] = lbHelpers<T, Lattice>::equilibriumFirstOrder( iPop, rho, u );
  }
}

template<typename T, template<typename U> class Lattice>
void AdvectionDiffusionBulkMomenta<T,Lattice>::defineU (
  Cell<T,Lattice>& cell,
  const T u[Lattice<T>::d])
{
  T rho = cell.computeRho();
  T *u_ = cell.getExternal(Lattice<T>::ExternalField::velocityBeginsAt);
  for (int iD = 0; iD < Lattice<T>::d; ++iD) {
    u_[iD] = u[iD];
  }
  for (int iPop=0; iPop < Lattice<T>::q; ++iPop) {
    cell[iPop] = lbHelpers<T, Lattice>::equilibriumFirstOrder( iPop, rho, u );
  }
}

template<typename T, template<typename U> class Lattice>
void AdvectionDiffusionBulkMomenta<T,Lattice>::defineRhoU (
  Cell<T,Lattice>& cell,
  T rho, const T u[Lattice<T>::d])
{
  T *u_ = cell.getExternal(Lattice<T>::ExternalField::velocityBeginsAt);
  for (int iD = 0; iD < Lattice<T>::d; ++iD) {
    u_[iD] = u[iD];
  }
  for (int iPop = 0; iPop < Lattice<T>::q; ++iPop) {
    cell[iPop] = lbHelpers<T, Lattice>::equilibriumFirstOrder( iPop, rho, u );
  }
}

template<typename T, template<typename U> class Lattice>
void AdvectionDiffusionBulkMomenta<T,Lattice>::defineAllMomenta (
  Cell<T,Lattice>& cell,
  T rho, const T u[Lattice<T>::d],
  const T pi[util::TensorVal<Lattice<T> >::n] )
{
  T *u_ = cell.getExternal(Lattice<T>::ExternalField::velocityBeginsAt);
  for (int iD = 0; iD < Lattice<T>::d; ++iD) {
    u_[iD] = u[iD];
  }
  for (int iPop = 0; iPop < Lattice<T>::q; ++iPop) {
    cell[iPop] = lbHelpers<T, Lattice>::equilibriumFirstOrder( iPop, rho, u );
  }
}

/////////////// Singletons //////////////////////////////////

namespace instances {

template<typename T, template<typename U> class Lattice>
AdvectionDiffusionBulkMomenta<T,Lattice>& getAdvectionDiffusionBulkMomenta()
{
  static AdvectionDiffusionBulkMomenta<T,Lattice> bulkMomentaSingleton;
  return bulkMomentaSingleton;
}
}

}

#endif
