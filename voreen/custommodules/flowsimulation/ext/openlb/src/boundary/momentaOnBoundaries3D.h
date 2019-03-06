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
 * Local boundary cell 3D dynamics -- header file.
 */
#ifndef MOMENTA_ON_BOUNDARIES_3D_H
#define MOMENTA_ON_BOUNDARIES_3D_H

#include "momentaOnBoundaries.h"

namespace olb {

template<typename T, template<typename U> class Lattice,
         int plane, int normal1, int normal2>
class InnerEdgeVelBM3D : public DirichletBoundaryMomenta<T,Lattice> {
public:
  enum { direction1 = (plane+1)%3, direction2 = (plane+2)%3 };
public:
  /// Default Constructor: initialization to zero
  InnerEdgeVelBM3D();
  /// Constructor with boundary initialization
  InnerEdgeVelBM3D(const T u_[Lattice<T>::d]);

  T computeRho(Cell<T,Lattice> const& cell) const override;
  void computeU (
    Cell<T,Lattice> const& cell,
    T u[Lattice<T>::d] ) const override;
  void computeJ (
    Cell<T,Lattice> const& cell,
    T j[Lattice<T>::d] ) const override;
  void computeU(T u[Lattice<T>::d]) const;
  void defineRho(Cell<T,Lattice>& cell, T rho) override ;
  void defineU(Cell<T,Lattice>& cell,
                       const T u[Lattice<T>::d]) override ;
  void defineU(const T u[Lattice<T>::d]);
  void defineAllMomenta (
    Cell<T,Lattice>& cell,
    T rho, const T u[Lattice<T>::d],
    const T pi[util::TensorVal<Lattice<T> >::n] ) override;
  /// Stress tensor
  void computeStress (
    Cell<T,Lattice> const& cell,
    T rho, const T u[Lattice<T>::d],
    T pi[util::TensorVal<Lattice<T> >::n] ) const override;
private:
  RegularizedVelocityBM<T,Lattice,direction1,normal1> momenta1;
  RegularizedVelocityBM<T,Lattice,direction2,normal2> momenta2;
};


template<typename T, template<typename U> class Lattice,
         int normalX, int normalY, int normalZ>
class InnerCornerVelBM3D : public DirichletBoundaryMomenta<T,Lattice> {
public:
  /// Default Constructor: initialization to zero
  InnerCornerVelBM3D();
  /// Constructor with boundary initialization
  InnerCornerVelBM3D(const T u_[Lattice<T>::d]);

  T computeRho(Cell<T,Lattice> const& cell) const override;
  void computeU (
    Cell<T,Lattice> const& cell,
    T u[Lattice<T>::d] ) const override;
  void computeJ (
    Cell<T,Lattice> const& cell,
    T j[Lattice<T>::d] ) const override;
  void computeU(T u[Lattice<T>::d]) const;
  void defineRho(Cell<T,Lattice>& cell, T rho) override ;
  void defineU(Cell<T,Lattice>& cell,
                       const T u[Lattice<T>::d]) override ;
  void defineU(const T u[Lattice<T>::d]);
  void defineAllMomenta (
    Cell<T,Lattice>& cell,
    T rho, const T u[Lattice<T>::d],
    const T pi[util::TensorVal<Lattice<T> >::n] ) override;
  /// Stress tensor
  void computeStress (
    Cell<T,Lattice> const& cell,
    T rho, const T u[Lattice<T>::d],
    T pi[util::TensorVal<Lattice<T> >::n] ) const override;
private:
  RegularizedVelocityBM<T,Lattice,0,normalX> xMomenta;
  RegularizedVelocityBM<T,Lattice,1,normalY> yMomenta;
  RegularizedVelocityBM<T,Lattice,2,normalZ> zMomenta;
};

}

#endif
