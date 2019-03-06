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
#ifndef MOMENTA_ON_BOUNDARIES_2D_HH
#define MOMENTA_ON_BOUNDARIES_2D_HH

#include "momentaOnBoundaries2D.h"
#include "dynamics/lbHelpers.h"
#include "dynamics/firstOrderLbHelpers.h"

namespace olb {

////////////////////// Class InnerCornerVelBM2D ///////////////

template<typename T, template<typename U> class Lattice,
         int normalX, int normalY>
InnerCornerVelBM2D<T,Lattice,normalX,normalY>::InnerCornerVelBM2D()
{ }

template<typename T, template<typename U> class Lattice,
         int normalX, int normalY>
InnerCornerVelBM2D<T,Lattice,normalX,normalY>::InnerCornerVelBM2D (
  const T u_[Lattice<T>::d])
  : xMomenta(u_), yMomenta(u_)
{ }

template<typename T, template<typename U> class Lattice,
         int normalX, int normalY>
T InnerCornerVelBM2D<T,Lattice,normalX,normalY>::computeRho (
  Cell<T,Lattice> const& cell ) const
{
  return (xMomenta.computeRho(cell) + yMomenta.computeRho(cell)) / (T)2;
}

template<typename T, template<typename U> class Lattice,
         int normalX, int normalY>
void InnerCornerVelBM2D<T,Lattice,normalX,normalY>::computeU (
  Cell<T,Lattice> const& cell,
  T u[Lattice<T>::d] ) const
{
  xMomenta.computeU(cell, u);
}

template<typename T, template<typename U> class Lattice,
         int normalX, int normalY>
void InnerCornerVelBM2D<T,Lattice,normalX,normalY>::computeJ (
  Cell<T,Lattice> const& cell,
  T j[Lattice<T>::d] ) const
{
  computeU(cell, j);
  T rho = computeRho(cell);
  for (int iD=0; iD<Lattice<T>::d; ++iD) {
    j[iD] *= rho;
  }
}

template<typename T, template<typename U> class Lattice,
         int normalX, int normalY>
void InnerCornerVelBM2D<T,Lattice,normalX,normalY>::computeU (
  T u[Lattice<T>::d] ) const
{
  xMomenta.computeU(u);
}

template<typename T, template<typename U> class Lattice,
         int normalX, int normalY>
void InnerCornerVelBM2D<T,Lattice,normalX,normalY>::defineRho (
  Cell<T,Lattice>& cell, T rho )
{ }

template<typename T, template<typename U> class Lattice,
         int normalX, int normalY>
void InnerCornerVelBM2D<T,Lattice,normalX,normalY>::defineU (
  Cell<T,Lattice>& cell,
  const T u[Lattice<T>::d] )
{
  xMomenta.defineU(cell, u);
  yMomenta.defineU(cell, u);
}

template<typename T, template<typename U> class Lattice,
         int normalX, int normalY>
void InnerCornerVelBM2D<T,Lattice,normalX,normalY>::defineU (
  const T u[Lattice<T>::d] )
{
  xMomenta.defineU(u);
  yMomenta.defineU(u);
}

template<typename T, template<typename U> class Lattice,
         int normalX, int normalY>
void InnerCornerVelBM2D<T,Lattice,normalX,normalY>::defineAllMomenta (
  Cell<T,Lattice>& cell,
  T rho, const T u[Lattice<T>::d],
  const T pi[util::TensorVal<Lattice<T> >::n] )
{
  xMomenta.defineU(u);
  yMomenta.defineU(u);
}

template<typename T, template<typename U> class Lattice,
         int normalX, int normalY>
void InnerCornerVelBM2D<T,Lattice,normalX,normalY>::computeStress (
  Cell<T,Lattice> const& cell,
  T rho, const T u[Lattice<T>::d],
  T pi[util::TensorVal<Lattice<T> >::n] ) const
{
  typedef lbHelpers<T,Lattice> lbH;
  Cell<T,Lattice> newCell(cell);
  int v[Lattice<T>::d] = { -normalX, -normalY };
  int unknownF  = util::findVelocity<Lattice<T> >(v);

  if (unknownF != Lattice<T>::q) {
    int oppositeF = util::opposite<Lattice<T> >(unknownF);

    T uSqr = util::normSqr<T,Lattice<T>::d>(u);

    newCell[unknownF] = newCell[oppositeF]
                        - lbH::equilibrium(oppositeF, rho, u, uSqr)
                        + lbH::equilibrium(unknownF, rho, u, uSqr);
  }

  lbH::computeStress(newCell, rho, u, pi);
}


}  // namespace olb

#endif
