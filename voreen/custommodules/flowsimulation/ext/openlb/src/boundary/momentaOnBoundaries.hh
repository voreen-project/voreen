/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Jonas Latt
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
 * Implementation of boundary cell dynamics -- generic implementation.
 */
#ifndef MOMENTA_ON_BOUNDARIES_HH
#define MOMENTA_ON_BOUNDARIES_HH

#include <limits>
#include "momentaOnBoundaries.h"
#include "dynamics/lbHelpers.h"
#include "dynamics/firstOrderLbHelpers.h"

namespace olb {

////////////////////// Class EquilibriumBM //////////////////////

template<typename T, template<typename U> class Lattice>
EquilibriumBM<T,Lattice>::EquilibriumBM()
{
  _rho = (T)1;
  for (int iD=0; iD<Lattice<T>::d; ++iD) {
    _u[iD] = T();
  }
}

template<typename T, template<typename U> class Lattice>
EquilibriumBM<T,Lattice>::EquilibriumBM(T rho, const T u[Lattice<T>::d])
{
  _rho = rho;
  for (int iD=0; iD<Lattice<T>::d; ++iD) {
    _u[iD] = u[iD];
  }
}

template<typename T, template<typename U> class Lattice>
T EquilibriumBM<T,Lattice>::computeRho( Cell<T,Lattice> const& cell ) const
{
  return _rho;
}

template<typename T, template<typename U> class Lattice>
void EquilibriumBM<T,Lattice>::computeU(
  Cell<T,Lattice> const& cell, T u[Lattice<T>::d]) const
{
  for (int iD=0; iD<Lattice<T>::d; ++iD) {
    u[iD] = _u[iD];
  }
}

template<typename T, template<typename U> class Lattice>
void EquilibriumBM<T,Lattice>::computeJ (
  Cell<T,Lattice> const& cell, T j[Lattice<T>::d]) const
{
  for (int iD=0; iD<Lattice<T>::d; ++iD) {
    j[iD] = _u[iD]*_rho;
  }
}

template<typename T, template<typename U> class Lattice>
void EquilibriumBM<T,Lattice>::computeStress( Cell<T,Lattice> const& cell,
    T rho, const T u[Lattice<T>::d], T pi[util::TensorVal<Lattice<T> >::n] ) const
{
  for (int iPi=0; iPi<util::TensorVal<Lattice<T> >::n; ++iPi) {
    pi[iPi] = T();
  }
}


template<typename T, template<typename U> class Lattice>
void EquilibriumBM<T,Lattice>::defineRho( Cell<T,Lattice>& cell, T rho )
{
  _rho = rho;
}

template<typename T, template<typename U> class Lattice>
void EquilibriumBM<T,Lattice>::defineU(
  Cell<T,Lattice>& cell, const T u[Lattice<T>::d])
{
  for (int iD=0; iD<Lattice<T>::d; ++iD) {
    _u[iD] = u[iD];
  }
}

template<typename T, template<typename U> class Lattice>
void EquilibriumBM<T,Lattice>::defineAllMomenta( Cell<T,Lattice>& cell,
    T rho, const T u[Lattice<T>::d], const T pi[util::TensorVal<Lattice<T> >::n] )
{
  defineRho(cell, rho);
  defineU(cell, u);
}



////////////////////// Class VelocityBM //////////////////////

template<typename T, template<typename U> class Lattice, int direction, int orientation>
VelocityBM<T,Lattice,direction,orientation>::VelocityBM()
{
  for (int iD=0; iD<Lattice<T>::d; ++iD) {
    _u[iD] = T();
  }
}

/** It takes as argument the value of the velocity to be imposed on the
 * boundary.
 */
template<typename T, template<typename U> class Lattice, int direction, int orientation>
VelocityBM<T,Lattice,direction,orientation>::VelocityBM(const T u[Lattice<T>::d])
{
  for (int iD=0; iD<Lattice<T>::d; ++iD) {
    _u[iD] = u[iD];
  }
}

template<typename T, template<typename U> class Lattice, int direction, int orientation>
T VelocityBM<T,Lattice,direction,orientation>::computeRho( Cell<T,Lattice> const& cell ) const
{
  std::vector<int> const& onWallIndices
    = util::subIndex<Lattice<T>, direction, 0>();

  std::vector<int> const& normalIndices
    = util::subIndex<Lattice<T>, direction, orientation>();

  T rhoOnWall = T();
  for (unsigned fIndex=0; fIndex<onWallIndices.size(); ++fIndex) {
    rhoOnWall += cell[onWallIndices[fIndex]];
  }

  T rhoNormal = T();
  for (unsigned fIndex=0; fIndex<normalIndices.size(); ++fIndex) {
    rhoNormal += cell[normalIndices[fIndex]];
  }

  T rho =((T)2*rhoNormal+rhoOnWall+(T)1) /
         ((T)1+(T)orientation*this->_u[direction]);

  return rho;
}

template<typename T, template<typename U> class Lattice, int direction, int orientation>
void VelocityBM<T,Lattice,direction,orientation>::computeU (
  Cell<T,Lattice> const& cell, T u[Lattice<T>::d]) const
{
  computeU(u);
}

template<typename T, template<typename U> class Lattice, int direction, int orientation>
void VelocityBM<T,Lattice,direction,orientation>::computeJ (
  Cell<T,Lattice> const& cell, T j[Lattice<T>::d]) const
{
  T rho = computeRho(cell);
  for (int iD=0; iD<Lattice<T>::d; ++iD) {
    j[iD] = _u[iD]*rho;
  }
}

template<typename T, template<typename U> class Lattice, int direction, int orientation>
void VelocityBM<T,Lattice,direction,orientation>::computeU( T u[Lattice<T>::d] ) const
{
  for (int iD=0; iD<Lattice<T>::d; ++iD) {
    u[iD] = _u[iD];
  }
}

template<typename T, template<typename U> class Lattice, int direction, int orientation>
void VelocityBM<T,Lattice,direction,orientation>::defineRho(Cell<T,Lattice>& cell, T rho )
{
}

template<typename T, template<typename U> class Lattice, int direction, int orientation>
void VelocityBM<T,Lattice,direction,orientation>::defineU (
  Cell<T,Lattice>& cell, const T u[Lattice<T>::d])
{
  defineU(u);
}

template<typename T, template<typename U> class Lattice, int direction, int orientation>
void VelocityBM<T,Lattice,direction,orientation>::defineU(const T u[Lattice<T>::d])
{
  for (int iD=0; iD<Lattice<T>::d; ++iD) {
    _u[iD] = u[iD];
  }
}

template<typename T, template<typename U> class Lattice, int direction, int orientation>
void VelocityBM<T,Lattice,direction,orientation>::defineAllMomenta (
  Cell<T,Lattice>& cell, T rho, const T u[Lattice<T>::d], const T pi[util::TensorVal<Lattice<T> >::n] )
{
  this->defineRhoU(cell, rho, u);
}


////////////////////// Class PressureBM //////////////////////

template<typename T, template<typename U> class Lattice, int direction, int orientation>
PressureBM<T,Lattice,direction,orientation>::PressureBM()
{
  for (int iD=0; iD<Lattice<T>::d; ++iD) {
    _values[iD] = T();
  }
  _values[direction] = 1.;
}

/** It takes as argument the value of the tangential velocity
 * components, and the pressure, to be imposed on the boundary.
 */
template<typename T, template<typename U> class Lattice, int direction, int orientation>
PressureBM<T,Lattice,direction,orientation>::PressureBM(const T values[Lattice<T>::d])
{
  for (int iD=0; iD<Lattice<T>::d; ++iD) {
    _values[iD] = values[iD];
  }
}

template<typename T, template<typename U> class Lattice, int direction, int orientation>
T PressureBM<T,Lattice,direction,orientation>::computeRho( Cell<T,Lattice> const& cell ) const
{
  return computeRho();
}

template<typename T, template<typename U> class Lattice, int direction, int orientation>
T PressureBM<T,Lattice,direction,orientation>::computeRho() const
{
  return _values[direction];
}

template<typename T, template<typename U> class Lattice, int direction, int orientation>
void PressureBM<T,Lattice,direction,orientation>::computeU (
  Cell<T,Lattice> const& cell, T u[Lattice<T>::d]) const
{
  for (int iD=0; iD<Lattice<T>::d; ++iD) {
    u[iD] = _values[iD];
  }
  T rho = _values[direction];

  std::vector<int> const& onWallIndices
    = util::subIndex<Lattice<T>, direction, 0>();

  std::vector<int> const& normalIndices
    = util::subIndex<Lattice<T>, direction, orientation>();

  T rhoOnWall = T();
  for (unsigned fIndex=0; fIndex<onWallIndices.size(); ++fIndex) {
    rhoOnWall += cell[onWallIndices[fIndex]];
  }

  T rhoNormal = T();
  for (unsigned fIndex=0; fIndex<normalIndices.size(); ++fIndex) {
    rhoNormal += cell[normalIndices[fIndex]];
  }

  u[direction] = (T)orientation*( ((T)2*rhoNormal+rhoOnWall+(T)1 ) / rho-(T)1 );
}

template<typename T, template<typename U> class Lattice, int direction, int orientation>
void PressureBM<T,Lattice,direction,orientation>::computeJ (
  Cell<T,Lattice> const& cell, T j[Lattice<T>::d]) const
{
  computeU(cell, j);
  T rho = computeRho(cell);
  for (int iD=0; iD<Lattice<T>::d; ++iD) {
    j[iD] *= rho;
  }
}

template<typename T, template<typename U> class Lattice, int direction, int orientation>
void PressureBM<T,Lattice,direction,orientation>::defineRho(Cell<T,Lattice>& cell, T rho)
{
  defineRho(rho);
}

template<typename T, template<typename U> class Lattice, int direction, int orientation>
void PressureBM<T,Lattice,direction,orientation>::defineRho(T rho )
{
  _values[direction] = rho;
}

template<typename T, template<typename U> class Lattice, int direction, int orientation>
void PressureBM<T,Lattice,direction,orientation>::defineU(Cell<T,Lattice>& cell, const T u[Lattice<T>::d])
{
  for (int iD=0; iD<Lattice<T>::d; ++iD) {
    if (iD != direction) {
      _values[iD] = u[iD];
    }
  }
}

template<typename T, template<typename U> class Lattice, int direction, int orientation>
void PressureBM<T,Lattice,direction,orientation>::defineAllMomenta(
  Cell<T,Lattice>& cell, T rho, const T u[Lattice<T>::d], const T pi[util::TensorVal<Lattice<T> >::n] )
{
  this->defineRhoU(cell, rho, u);
}


////////  FreeStressBM //////////////////////////////////////////////

template<typename T, template<typename U> class Lattice>
void FreeStressBM<T,Lattice>::computeStress (
  Cell<T,Lattice> const& cell, T rho, const T u[Lattice<T>::d], T pi[util::TensorVal<Lattice<T> >::n] ) const
{
  lbHelpers<T,Lattice>::computeStress(cell, rho, u, pi);
}



////////////////////// Class RegularizedBM //////////////////////

template<typename T, template<typename U> class Lattice, int direction, int orientation>
void RegularizedBM<T,Lattice,direction,orientation>::computeStress (
  Cell<T,Lattice> const& cell, T rho, const T u[Lattice<T>::d],
  T pi[util::TensorVal<Lattice<T> >::n] ) const
{
  BoundaryHelpers<T,Lattice,direction,orientation>::computeStress(cell, rho, u, pi);
}

////////////////////// Class FixedVelocityBM //////////////////////////

template<typename T, template<typename U> class Lattice>
T FixedVelocityBM<T,Lattice>::computeRho(Cell<T,Lattice> const& cell) const
{
  return _basicMomenta.computeRho(cell);
}

template<typename T, template<typename U> class Lattice>
void FixedVelocityBM<T,Lattice>::computeU(Cell<T,Lattice> const& cell, T u[Lattice<T>::d]) const
{
  for (int iD=0; iD<Lattice<T>::d; ++iD) {
    u[iD] = _fixU[iD];
  }
}

template<typename T, template<typename U> class Lattice>
void FixedVelocityBM<T,Lattice>::computeJ(Cell<T,Lattice> const& cell, T j[Lattice<T>::d]) const
{
  T rho = computeRho(cell);
  for (int iD=0; iD<Lattice<T>::d; ++iD) {
    j[iD] = _fixU[iD]*rho;
  }
}

template<typename T, template<typename U> class Lattice>
void FixedVelocityBM<T,Lattice>::computeStress (
  Cell<T,Lattice> const& cell, T rho, const T u[Lattice<T>::d], T pi[util::TensorVal<Lattice<T> >::n] ) const
{
  _basicMomenta.computeStress(cell, rho, u, pi);
}

template<typename T, template<typename U> class Lattice>
void FixedVelocityBM<T,Lattice>::computeRhoU (
  Cell<T,Lattice> const& cell, T& rho, T u[Lattice<T>::d] ) const
{
  rho = computeRho(cell);
  computeU(cell,u);
}

template<typename T, template<typename U> class Lattice>
void FixedVelocityBM<T,Lattice>::computeAllMomenta (
  Cell<T,Lattice> const& cell, T& rho, T u[Lattice<T>::d], T pi[util::TensorVal<Lattice<T> >::n] ) const
{
  _basicMomenta.computeAllMomenta(cell, rho, u, pi);
  computeU(cell, u);
}

template<typename T, template<typename U> class Lattice>
void FixedVelocityBM<T,Lattice>::defineRho(Cell<T,Lattice>& cell, T rho)
{
  _basicMomenta.defineRho(cell, rho);
}

template<typename T, template<typename U> class Lattice>
void FixedVelocityBM<T,Lattice>::defineU(Cell<T,Lattice>& cell, const T u[Lattice<T>::d])
{
  for (int iD=0; iD<Lattice<T>::d; ++iD) {
    _fixU[iD] = u[iD];
  }
}

template<typename T, template<typename U> class Lattice>
void FixedVelocityBM<T,Lattice>::defineRhoU(Cell<T,Lattice>& cell, T rho, const T u[Lattice<T>::d])
{
  defineRho(cell,rho);
  defineU(cell,u);
}

template<typename T, template<typename U> class Lattice>
void FixedVelocityBM<T,Lattice>::defineAllMomenta( Cell<T,Lattice>& cell, T rho, const T u[Lattice<T>::d],
    const T pi[util::TensorVal<Lattice<T> >::n] )
{
  _basicMomenta.defineAllMomenta(cell, rho, u, pi);
  defineU(cell,u);
}



}  // namespace olb

#endif
