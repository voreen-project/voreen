/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012-2017 Lukas Baron, Tim Dornieden, Mathias J. Krause,
 *  Albert Mink, Benjamin Förster, Adrian Kummerländer
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

#ifndef SUPER_LATTICE_INTEGRAL_F_3D_H
#define SUPER_LATTICE_INTEGRAL_F_3D_H

#include <vector>

#include "functors/genericF.h"
#include "blockLatticeIntegralF3D.h"
#include "superBaseF3D.h"
#include "indicator/superIndicatorBaseF3D.h"
#include "indicator/indicatorBaseF3D.h"
#include "indicator/superIndicatorF3D.h"
#include "superLatticeLocalF3D.h"
#include "functors/analytical/interpolationF3D.h"
#include "functors/lattice/reductionF3D.h"
#include "integral/superIntegralF3D.h"
#include "core/superLattice3D.h"
#include "core/vector.h"
#include "io/ostreamManager.h"
#include "geometry/superGeometry3D.h"
#include "utilities/functorPtr.h"

/** Note: Throughout the whole source code directory genericFunctions, the
 *  template parameters for i/o dimensions are:
 *           F: S^m -> T^n  (S=source, T=target)
 */

namespace olb {

////////////////////////////////////////////////////////////////////////////////
//////if globIC is not on the local processor, the returned vector is empty/////
////////////////////////////////////////////////////////////////////////////////


/// SuperMin3D returns the min in each component of f on a indicated subset
template <typename T, typename W = T>
class SuperMin3D final : public SuperF3D<T,W> {
private:
  FunctorPtr<SuperF3D<T,W>>        _f;
  FunctorPtr<SuperIndicatorF3D<T>> _indicatorF;
public:
  /// Constructor for determining the minimum of f on a indicated subset
  /**
   * \param f          functor of which the minimum is to be determined
   * \param indicatorF indicator describing the subset on which to evaluate f
   **/
  SuperMin3D(FunctorPtr<SuperF3D<T,W>>&&        f,
             FunctorPtr<SuperIndicatorF3D<T>>&& indicatorF);
  /// Constructor for determining the minimum of f on a given material
  /**
   * \param f             functor of which the minimum is to be determined
   * \param superGeometry super geometry for constructing material indicator
   * \param material      number of the relevant material
   **/
  SuperMin3D(FunctorPtr<SuperF3D<T,W>>&& f,
             SuperGeometry3D<T>& superGeometry,
             const int material);

  bool operator() (W output[], const int input[]) override;
};


/// SuperMax3D returns the max in each component of f on a indicated subset
template <typename T, typename W = T>
class SuperMax3D final : public SuperF3D<T,W> {
private:
  FunctorPtr<SuperF3D<T,W>>        _f;
  FunctorPtr<SuperIndicatorF3D<T>> _indicatorF;
public:
  /// Constructor for determining the maximum of f on a indicated subset
  /**
   * \param f          functor of which the maximum is to be determined
   * \param indicatorF indicator describing the subset on which to evaluate f
   **/
  SuperMax3D(FunctorPtr<SuperF3D<T,W>>&&        f,
             FunctorPtr<SuperIndicatorF3D<T>>&& indicatorF);
  /// Constructor for determining the maximum of f on a given material
  /**
   * \param f             functor of which the maximum is to be determined
   * \param superGeometry super geometry for constructing material indicator
   * \param material      number of the relevant material
   **/
  SuperMax3D(FunctorPtr<SuperF3D<T,W>>&& f,
             SuperGeometry3D<T>& superGeometry,
             const int material);

  bool operator() (W output[], const int input[]) override;
};


/// SuperAverage3D returns the average in each component of f on a indicated subset
template <typename T, typename W = T>
class SuperAverage3D final : public SuperF3D<T,W> {
private:
  FunctorPtr<SuperF3D<T,W>>        _f;
  FunctorPtr<SuperIndicatorF3D<T>> _indicatorF;
public:
  /// Constructor for determining the average of f on a indicated subset
  /**
   * \param f          functor of which the average is to be determined
   * \param indicatorF indicator describing the subset on which to evaluate f
   **/
  SuperAverage3D(FunctorPtr<SuperF3D<T,W>>&&        f,
                 FunctorPtr<SuperIndicatorF3D<T>>&& indicatorF);
  /// Constructor for determining the average of f on a given material
  /**
   * \param f             functor of which the average is to be determined
   * \param superGeometry super geometry for constructing material indicator
   * \param material      number of the relevant material
   **/
  SuperAverage3D(FunctorPtr<SuperF3D<T,W>>&& f,
                 SuperGeometry3D<T>& superGeometry,
                 const int material);

  /// Global average operator
  /**
   * Note: While this functor exposes BlockAverage3D functors if possible, a call to
   * this function will not use them but calculate the global average by summing all
   * components and voxel counts.
   * Calling BlockAverage3D in this situation would unnecessarily complicate this as
   * we would have to weight the aggregated averages according to their share in the
   * global average.
   **/
  bool operator() (W output[], const int input[]) override;
};


/// sums over all cells of a SmoothIndicatorSphere
template <typename T, typename W = T>
class SuperSumIndicator3D final : public SuperF3D<T,W> {
private:
  SuperF3D<T,W>& _f;
  SuperGeometry3D<T>& _superGeometry;
  ParticleIndicatorF3D<T,T>& _indicator;
public:
  SuperSumIndicator3D(SuperF3D<T,W>& f, SuperGeometry3D<T>& superGeometry, ParticleIndicatorF3D<T,T>& indicator);
  bool operator() (T output[], const int input[]) override;
};


/// functor counts to get the discrete surface for a material no. in direction (1,0,0), (0,1,0), (0,0,1), (-1,0,0), (0,-1,0), (0,0,-1) and total surface, then it converts it into phys units
template <typename T>
class SuperGeometryFaces3D final : public GenericF<T,int> {
private:
  SuperGeometry3D<T>& _superGeometry;
  const int _material;
  T _latticeL;
  std::vector<BlockGeometryFaces3D<T>* > _blockGeometryFaces;
public:
  template<template<typename U> class DESCRIPTOR>
  SuperGeometryFaces3D(SuperGeometry3D<T>& superGeometry, const int material, const UnitConverter<T,DESCRIPTOR>& converter);
  SuperGeometryFaces3D(SuperGeometry3D<T>& superGeometry, const int material, T _latticeL);
  bool operator() (T output[], const int input[]) override;
  std::vector<T> operator() (std::vector<int> input)
  {
    std::vector<T> output(this->getTargetDim(), T());
    operator()(&output[0],&input[0]);
    return output;
  };
};


/// functor to get pointwise phys force acting on a boundary with a given material on local lattice
template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticePhysDrag3D final : public SuperLatticePhysF3D<T,DESCRIPTOR> {
private:
  SuperGeometry3D<T>& _superGeometry;
  const int _material;
  SuperGeometryFaces3D<T> _faces;
  SuperLatticePhysBoundaryForce3D<T,DESCRIPTOR> _pBoundForce;
  SuperSum3D<T,T> _sumF;
  T _factor;
public:
  SuperLatticePhysDrag3D(SuperLattice3D<T,DESCRIPTOR>& sLattice,
                         SuperGeometry3D<T>& superGeometry, const int material,
                         const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};

/// functor to get pointwise phys force acting on a boundary with a given indicator on local lattice
template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticePhysDragIndicator3D final : public SuperLatticePhysF3D<T,DESCRIPTOR> {
private:
  SuperGeometry3D<T>& _superGeometry;
  ParticleIndicatorSphere3D<T,T>& _indicator;
public:
  SuperLatticePhysDragIndicator3D(SuperLattice3D<T,DESCRIPTOR>& sLattice,
                                  SuperGeometry3D<T>& superGeometry,
                                  ParticleIndicatorSphere3D<T,T>& indicator,
                                  const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};


/**
 *  functor to get pointwise phys force acting on a boundary with a given material on local lattice
 *  see: Caiazzo, Junk: Boundary Forces in lattice Boltzmann: Analysis of MEA
 */
template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticePhysCorrDrag3D final : public SuperLatticePhysF3D<T,DESCRIPTOR> {
private:
  SuperGeometry3D<T>& _superGeometry;
  const int _material;
public:
  SuperLatticePhysCorrDrag3D(SuperLattice3D<T,DESCRIPTOR>& sLattice,
                             SuperGeometry3D<T>& superGeometry, const int material,
                             const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};


} // end namespace olb

#endif
