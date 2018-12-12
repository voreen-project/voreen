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
 * Local boundary cell 3D dynamics -- generic implementation.
 */
#ifndef MOMENTA_ON_BOUNDARIES_3D_HH
#define MOMENTA_ON_BOUNDARIES_3D_HH

#include "momentaOnBoundaries3D.h"
#include "dynamics/lbHelpers.h"
#include "dynamics/firstOrderLbHelpers.h"

namespace olb {

////////////////////// Class InnerEdgeVelBM3D ///////////////

template<typename T, template<typename U> class Lattice,
         int plane, int normal1, int normal2>
InnerEdgeVelBM3D<T,Lattice,plane,normal1,normal2>::
InnerEdgeVelBM3D()
{ }

template<typename T, template<typename U> class Lattice,
         int plane, int normal1, int normal2>
InnerEdgeVelBM3D<T,Lattice,plane,normal1,normal2>::
InnerEdgeVelBM3D(const T u_[Lattice<T>::d])
  : momenta1(u_), momenta2(u_)
{ }

template<typename T, template<typename U> class Lattice,
         int plane, int normal1, int normal2>
T InnerEdgeVelBM3D<T,Lattice,plane,normal1,normal2>::computeRho (
  Cell<T,Lattice> const& cell ) const
{
  return (momenta1.computeRho(cell) + momenta2.computeRho(cell)) / (T)2;
}

template<typename T, template<typename U> class Lattice,
         int plane, int normal1, int normal2>
void InnerEdgeVelBM3D<T,Lattice,plane,normal1,normal2>::computeU (
  Cell<T,Lattice> const& cell,
  T u[Lattice<T>::d] ) const
{
  momenta1.computeU(cell, u);
}

template<typename T, template<typename U> class Lattice,
         int plane, int normal1, int normal2>
void InnerEdgeVelBM3D<T,Lattice,plane,normal1,normal2>::computeJ (
  Cell<T,Lattice> const& cell,
  T j[Lattice<T>::d] ) const
{
  momenta1.computeJ(cell, j);
}

template<typename T, template<typename U> class Lattice,
         int plane, int normal1, int normal2>
void InnerEdgeVelBM3D<T,Lattice,plane,normal1,normal2>::computeU (
  T u[Lattice<T>::d] ) const
{
  momenta1.computeU(u);
}

template<typename T, template<typename U> class Lattice,
         int plane, int normal1, int normal2>
void InnerEdgeVelBM3D<T,Lattice,plane,normal1,normal2>::defineRho (
  Cell<T,Lattice>& cell, T rho )
{ }

template<typename T, template<typename U> class Lattice,
         int plane, int normal1, int normal2>
void InnerEdgeVelBM3D<T,Lattice,plane,normal1,normal2>::defineU (
  Cell<T,Lattice>& cell,
  const T u[Lattice<T>::d] )
{
  momenta1.defineU(cell, u);
  momenta2.defineU(cell, u);
}

template<typename T, template<typename U> class Lattice,
         int plane, int normal1, int normal2>
void InnerEdgeVelBM3D<T,Lattice,plane,normal1,normal2>::defineU (
  const T u[Lattice<T>::d] )
{
  momenta1.defineU(u);
  momenta2.defineU(u);
}

template<typename T, template<typename U> class Lattice,
         int plane, int normal1, int normal2>
void InnerEdgeVelBM3D<T,Lattice,plane,normal1,normal2>::
defineAllMomenta (
  Cell<T,Lattice>& cell,
  T rho, const T u[Lattice<T>::d],
  const T pi[util::TensorVal<Lattice<T> >::n] )
{
  momenta1.defineU(u);
  momenta2.defineU(u);
}

template<typename T, template<typename U> class Lattice,
         int plane, int normal1, int normal2>
void InnerEdgeVelBM3D<T,Lattice,plane,normal1,normal2>::
computeStress (
  Cell<T,Lattice> const& cell,
  T rho, const T u[Lattice<T>::d],
  T pi[util::TensorVal<Lattice<T> >::n] ) const
{
  typedef lbHelpers<T,Lattice> lbH;

  T uSqr = util::normSqr<T,Lattice<T>::d>(u);

  Cell<T,Lattice> newCell(cell);
  for (int iPop=0; iPop<Lattice<T>::q; ++iPop) {
    if ( (Lattice<T>::c[iPop][direction1] == -normal1) &&
         (Lattice<T>::c[iPop][direction2] == -normal2) ) {
      int opp = util::opposite<Lattice<T> >(iPop);
      newCell[iPop] = newCell[opp]
                      - lbH::equilibrium(opp, rho, u, uSqr)
                      + lbH::equilibrium(iPop, rho, u, uSqr);
    }
  }
  lbH::computeStress(newCell, rho, u, pi);
}


////////////////////// Class InnerCornerVelBM3D ///////////////

template<typename T, template<typename U> class Lattice,
         int normalX, int normalY, int normalZ>
InnerCornerVelBM3D<T,Lattice,normalX,normalY,normalZ>::InnerCornerVelBM3D()
{ }

template<typename T, template<typename U> class Lattice,
         int normalX, int normalY, int normalZ>
InnerCornerVelBM3D<T,Lattice,normalX,normalY,normalZ>::InnerCornerVelBM3D (
  const T u_[Lattice<T>::d])
  : xMomenta(u_), yMomenta(u_), zMomenta(u_)
{ }

template<typename T, template<typename U> class Lattice,
         int normalX, int normalY, int normalZ>
T InnerCornerVelBM3D<T,Lattice,normalX,normalY,normalZ>::computeRho (
  Cell<T,Lattice> const& cell ) const
{
  return (xMomenta.computeRho(cell) +
          yMomenta.computeRho(cell) +
          zMomenta.computeRho(cell) ) / (T)3;
}

template<typename T, template<typename U> class Lattice,
         int normalX, int normalY, int normalZ>
void InnerCornerVelBM3D<T,Lattice,normalX,normalY,normalZ>::computeU (
  Cell<T,Lattice> const& cell,
  T u[Lattice<T>::d] ) const
{
  xMomenta.computeU(cell, u);
}

template<typename T, template<typename U> class Lattice,
         int normalX, int normalY, int normalZ>
void InnerCornerVelBM3D<T,Lattice,normalX,normalY,normalZ>::computeJ (
  Cell<T,Lattice> const& cell,
  T j[Lattice<T>::d] ) const
{
  xMomenta.computeJ(cell, j);
}

template<typename T, template<typename U> class Lattice,
         int normalX, int normalY, int normalZ>
void InnerCornerVelBM3D<T,Lattice,normalX,normalY,normalZ>::computeU (
  T u[Lattice<T>::d] ) const
{
  xMomenta.computeU(u);
}

template<typename T, template<typename U> class Lattice,
         int normalX, int normalY, int normalZ>
void InnerCornerVelBM3D<T,Lattice,normalX,normalY,normalZ>::defineRho (
  Cell<T,Lattice>& cell, T rho )
{ }

template<typename T, template<typename U> class Lattice,
         int normalX, int normalY, int normalZ>
void InnerCornerVelBM3D<T,Lattice,normalX,normalY,normalZ>::defineU (
  Cell<T,Lattice>& cell,
  const T u[Lattice<T>::d] )
{
  xMomenta.defineU(cell, u);
  yMomenta.defineU(cell, u);
  zMomenta.defineU(cell, u);
}

template<typename T, template<typename U> class Lattice,
         int normalX, int normalY, int normalZ>
void InnerCornerVelBM3D<T,Lattice,normalX,normalY,normalZ>::defineU (
  const T u[Lattice<T>::d] )
{
  xMomenta.defineU(u);
  yMomenta.defineU(u);
  zMomenta.defineU(u);
}

template<typename T, template<typename U> class Lattice,
         int normalX, int normalY, int normalZ>
void InnerCornerVelBM3D<T,Lattice,normalX,normalY,normalZ>::defineAllMomenta (
  Cell<T,Lattice>& cell,
  T rho, const T u[Lattice<T>::d],
  const T pi[util::TensorVal<Lattice<T> >::n] )
{
  xMomenta.defineU(u);
  yMomenta.defineU(u);
  zMomenta.defineU(u);
}

template<typename T, template<typename U> class Lattice,
         int normalX, int normalY, int normalZ>
void InnerCornerVelBM3D<T,Lattice,normalX,normalY,normalZ>::computeStress (
  Cell<T,Lattice> const& cell,
  T rho, const T u[Lattice<T>::d],
  T pi[util::TensorVal<Lattice<T> >::n] ) const
{
  typedef lbHelpers<T,Lattice> lbH;
  Cell<T,Lattice> newCell(cell);
  int v[Lattice<T>::d] = { -normalX, -normalY, -normalZ };
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
