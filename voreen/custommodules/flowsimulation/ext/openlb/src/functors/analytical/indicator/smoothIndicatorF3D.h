/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014-2016 Cyril Masquelier, Mathias J. Krause, Albert Mink
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

#ifndef SMOOTH_INDICATOR_F_3D_H
#define SMOOTH_INDICATOR_F_3D_H

#include <vector>

#include "smoothIndicatorBaseF3D.h"
#include "io/xmlReader.h"

#include "core/blockData3D.h"
#include "core/unitConverter.h"
#include "functors/analytical/indicator/indicatorBaseF2D.h"
#include "functors/analytical/indicator/indicatorBaseF3D.h"

namespace olb {

/// implements a smooth particle cuboid in 3D with an _epsilon sector.
template <typename T, typename S, bool HLBM=false>
class SmoothIndicatorCuboid3D final: public SmoothIndicatorF3D<T, S, HLBM> {
private:
  S _xLength;
  S _yLength;
  S _zLength;
public:
  SmoothIndicatorCuboid3D(Vector<S,3> center, S xLength, S yLength, S zLength, S epsilon, Vector<S,3> theta = Vector<S,3> (0.,0.,0.), S density=0, Vector<S,3> vel = Vector<S,3> (0.,0.,0.));
  bool operator()(T output[],const S input[]) override;
};

/// implements a smooth particle ellipsoid in 3D with an _epsilon sector.
template <typename T, typename S, bool HLBM=false>
class SmoothIndicatorEllipsoid3D final: public SmoothIndicatorF3D<T, S, HLBM> {
private:
  S _xHalfAxis;
  S _yHalfAxis;
  S _zHalfAxis;
public:
  SmoothIndicatorEllipsoid3D(Vector<S,3> center, S xHalfAxis, S yHalfAxis, S zHalfAxis, S epsilon, Vector<S,3> theta = Vector<S,3> (0.,0.,0.), S density=0, Vector<S,3> vel = Vector<S,3> (0.,0.,0.));
  bool operator()(T output[],const S input[]) override;
};

/// implements a smooth particle super-ellipsoid in 3D. The epsilon sector is currently missing.
template <typename T, typename S, bool HLBM=false>
class SmoothIndicatorSuperEllipsoid3D final: public SmoothIndicatorF3D<T, S, HLBM> {
private:
  S _xHalfAxis;
  S _yHalfAxis;
  S _zHalfAxis;
  S _exp1;
  S _exp2;
public:
  SmoothIndicatorSuperEllipsoid3D(Vector<S,3> center, S xHalfAxis, S yHalfAxis, S zHalfAxis, S exponent1, S exponent2, S epsilon, Vector<S,3> theta = Vector<S,3> (0.,0.,0.), S density=0, Vector<S,3> vel = Vector<S,3> (0.,0.,0.));
  // this implements the beta function from the gamma function and will be deprecated when switching to c++17
  S beta(S arg1, S arg2);
  // calculates cartesian moments
  S moments(S p, S q, S r);
  bool operator()(T output[],const S input[]) override;
};

/// implements a smooth sphere in 3D with an _epsilon sector
template <typename T, typename S, bool HLBM=false>
class SmoothIndicatorSphere3D final: public SmoothIndicatorF3D<T, S, HLBM> {
private:
  S _radius;
public:
  SmoothIndicatorSphere3D(Vector<S, 3> center, S radius, S epsilon, S density=0, Vector<S,3> vel = Vector<S,3> (0.,0.,0.));
  bool operator()(T output[], const S input[]) override;
};

/// implements a smooth particle cylinder in 3D with an _epsilon sector.
template <typename T, typename S, bool HLBM=false>
class SmoothIndicatorCylinder3D final: public SmoothIndicatorF3D<T, S, HLBM> {
private:
  S _radius;
  S _length;
  void initIndicatorCylinder3D(Vector<S,3> normal, Vector<S,3> theta, S density, Vector<S,3> vel);
public:
  SmoothIndicatorCylinder3D(Vector<S,3> pointA, Vector<S,3> pointB, S radius, S epsilon, Vector<S,3> theta = Vector<S,3> (0.,0.,0.), S density=0, Vector<S,3> vel = Vector<S,3> (0.,0.,0.));
  SmoothIndicatorCylinder3D(Vector<S,3> center, Vector<S,3> normal, S radius, S length, S epsilon, Vector<S,3> theta = Vector<S,3> (0.,0.,0.), S density=0, Vector<S,3> vel = Vector<S,3> (0.,0.,0.));
  bool operator()(T output[], const S input[]) override;
};

/// implements a smooth particle cone in 3D with an _epsilon sector
template <typename T, typename S, bool HLBM=false>
class SmoothIndicatorCone3D : public SmoothIndicatorF3D<T, S, HLBM> {
private:
  S _length;
  S _radiusA;
  S _radiusB;
  void initIndicatorCone3D(Vector<S,3> normal, Vector<S,3> theta, S density, Vector<S,3> vel);
public:
  SmoothIndicatorCone3D(Vector<S,3> pointA, Vector<S,3> pointB,
                        S radiusA, S radiusB, S epsilon, Vector<S,3> theta = Vector<S,3> (0.,0.,0.), S density=0,
                        Vector<S,3> vel = Vector<S,3> (0.,0.,0.));
  SmoothIndicatorCone3D(Vector<S,3> center, Vector<S,3> normal, S lenght,
                        S radiusA, S radiusB, S epsilon, Vector<S,3> theta = Vector<S,3> (0.,0.,0.),
                        S density=0, Vector<S,3> vel = Vector<S,3> (0.,0.,0.));
  bool operator() (T output[], const S input[]) override;
};


//implements a custom shaped smooth particle //TODO: Check for consistency
//ALSO: adap .hh
template <typename T, typename S, typename DESCRIPTOR, bool HLBM=false>
class SmoothIndicatorCustom3D final: public SmoothIndicatorF3D<T, S, HLBM> {
private:
  /// Turn on/off additional output
  const bool _verbose;
  /// Turn on/off to use real boundary with eps = 0.5
  const bool _useRealBoudnary;
  /// Lattice spacing (in m) for the particle lattice (should be smaller or equal to the fluid lattice for best results)
  const T _latticeSpacing;
  /// Important parameter for the Gaussian point spread Function (standard deviations)
  const T _sigma;
  /// Local center
  std::vector<T> _center;
  /// Smoothed block data to store porosity
  BlockData3D<T, T> _blockData;

  void initRotationMatrix();
  void initBlockData(IndicatorF3D<T>& ind);
  void calcCenter();
  void calcMofi(T rhoP);
  void calcCircumRadius();

public:
  SmoothIndicatorCustom3D(T latticeSpacing,
                          IndicatorF3D<T>& ind, Vector<T,3> pos, T density, T epsilon,
                          Vector<T,3> theta, Vector<T,3> vel = Vector<T,3> (0.,0.,0.),
                          T sigma = 1., bool verbose=false, bool useRealBoundary=false);
  SmoothIndicatorCustom3D(UnitConverter<T,DESCRIPTOR> const& converter,
                          IndicatorF3D<T>& ind, Vector<T,3> pos, T density, T epsilon,
                          Vector<T,3> theta, Vector<T,3> vel = Vector<T,3> (0.,0.,0.),
                          T sigma = 1., bool verbose=false, bool useRealBoundary=false);
  Vector<T,3> getLocalCenter();
  Vector<T,3> getMofi();
  bool regardCell(BlockData3D<T,T>& blockData, int x, int y, int z);
  bool operator() (T output[], const S input[]);
};


}

#endif

