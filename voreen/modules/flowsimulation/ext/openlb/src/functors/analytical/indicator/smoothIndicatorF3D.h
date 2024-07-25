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
#include "utilities/functorPtr.h"

#include "core/blockData.h"
#include "core/unitConverter.h"
#include "functors/analytical/indicator/indicatorBaseF2D.h"
#include "functors/analytical/indicator/indicatorBaseF3D.h"
#include "functors/analytical/indicator/indicatorF3D.h"
#include "particles/functions/bodyMotionFunctions.h"
#include "sdf.h"

namespace olb {

/// implements a smooth particle cuboid in 3D with an _epsilon sector.
template <typename T, typename S, bool PARTICLE=false>
class SmoothIndicatorCuboid3D final: public SmoothIndicatorF3D<T, S, PARTICLE> {
private:
  IndicatorCuboid3D<S> _ind;

public:
  SmoothIndicatorCuboid3D(IndicatorCuboid3D<S>& ind, S epsilon, Vector<S,3> theta = Vector<S,3> (0.,0.,0.));
  SmoothIndicatorCuboid3D(S xLength, S yLength, S zLength, Vector<S,3> center, S epsilon, Vector<S,3> theta = Vector<S,3> (0.,0.,0.));
  Vector<S,3> surfaceNormal(const Vector<S,3>& pos, const S meshSize) override;
  const S signedDistance( const PhysR<S,3> input ) override;
  S getVolume( ) override;
  Vector<S,4> calcMofiAndMass(const S density) override;
};

/// implements a smooth particle ellipsoid in 3D with an _epsilon sector.
template <typename T, typename S, bool PARTICLE=false>
class SmoothIndicatorEllipsoid3D final: public SmoothIndicatorF3D<T, S, PARTICLE> {
private:
  IndicatorEllipsoid3D<S> _ind;

public:
  SmoothIndicatorEllipsoid3D(IndicatorEllipsoid3D<S>& ind, S epsilon, Vector<S,3> theta = Vector<S,3> (0.,0.,0.));
  SmoothIndicatorEllipsoid3D(Vector<S,3> center, Vector<S,3> radius, S epsilon, Vector<S,3> theta = Vector<S,3> (0.,0.,0.));
  Vector<S,3> surfaceNormal(const Vector<S,3>& pos, const S meshSize) override;
  const S signedDistance( const PhysR<S,3> input ) override;
  S getVolume( ) override;
  Vector<S,4> calcMofiAndMass(const S density) override;
  bool operator()(T output[],const S input[]) override;
};

/// implements a smooth particle super-ellipsoid in 3D. The epsilon sector is currently missing.
template <typename T, typename S, bool PARTICLE=false>
class SmoothIndicatorSuperEllipsoid3D final: public SmoothIndicatorF3D<T, S, PARTICLE> {
private:
  IndicatorSuperEllipsoid3D<S> _ind;
public:
  SmoothIndicatorSuperEllipsoid3D(IndicatorSuperEllipsoid3D<S>& ind, S epsilon, Vector<S,3> theta = Vector<S,3> (0.,0.,0.));
  SmoothIndicatorSuperEllipsoid3D(Vector<S,3> center, S xHalfAxis, S yHalfAxis, S zHalfAxis, S exponent1, S exponent2,
                                  S epsilon, Vector<S,3> theta = Vector<S,3> (0.,0.,0.));
  // this implements the beta function from the gamma function and will be deprecated when switching to c++17
  S beta(S arg1, S arg2);
  // calculates cartesian moments
  S moments(S p, S q, S r);
  Vector<S,3> surfaceNormal(const Vector<S,3>& pos, const S meshSize) override;
  const S signedDistance( const PhysR<S,3> input ) override;
  S getVolume( ) override;
  Vector<S,4> calcMofiAndMass(const S density) override;
  bool operator()(T output[],const S input[]) override;
};

/// implements a smooth sphere in 3D with an _epsilon sector
template <typename T, typename S, bool PARTICLE=false>
class SmoothIndicatorSphere3D final: public SmoothIndicatorF3D<T, S, PARTICLE> {
private:
  IndicatorSphere3D<S> _ind;
public:
  SmoothIndicatorSphere3D(IndicatorSphere3D<S>& ind, S epsilon);
  SmoothIndicatorSphere3D(Vector<S,3> center, S radius, S epsilon);
  Vector<S,3> surfaceNormal(const Vector<S,3>& pos, const S meshSize) override;
  const S signedDistance( const PhysR<S,3> input ) override;
  S getVolume( ) override;
  Vector<S,4> calcMofiAndMass(const S density) override;
};

/// implements a smooth particle cylinder in 3D with an _epsilon sector.
template <typename T, typename S, bool PARTICLE=false>
class SmoothIndicatorCylinder3D final: public SmoothIndicatorF3D<T, S, PARTICLE> {
private:
  IndicatorCylinder3D<S> _ind;
  void initIndicatorCylinder3D(Vector<S,3> theta, S length);

public:
  SmoothIndicatorCylinder3D(IndicatorCylinder3D<S>& ind, S epsilon, Vector<S,3> theta = Vector<S,3> (0.,0.,0.));
  SmoothIndicatorCylinder3D(Vector<S,3> center1, Vector<S,3> center2, S radius, S epsilon, Vector<S,3> theta = Vector<S,3> (0.,0.,0.));
  SmoothIndicatorCylinder3D(Vector<S,3> center, Vector<S,3> normal, S radius, S height, S epsilon, Vector<S,3> theta = Vector<S,3> (0.,0.,0.));
  Vector<S,3> surfaceNormal(const Vector<S,3>& pos, const S meshSize) override;
  const S signedDistance( const PhysR<S,3> input ) override;
  S getVolume( ) override;
  Vector<S,4> calcMofiAndMass(const S density) override;
};

/// implements a smooth particle cone in 3D with an _epsilon sector
template <typename T, typename S, bool PARTICLE=false>
class SmoothIndicatorCone3D : public SmoothIndicatorF3D<T, S, PARTICLE> {
private:
  Vector<S,3> _startPos;
  void initIndicatorCone3D(Vector<S,3> theta, S length);
  IndicatorCone3D<S> _ind;
public:
  SmoothIndicatorCone3D(IndicatorCone3D<S>& indPtr, S epsilon, Vector<S,3> theta = Vector<S,3> (0.,0.,0.));
  SmoothIndicatorCone3D(Vector<S,3> center1, Vector<S,3> center2, S radius1, S radius2, S epsilon, Vector<S,3> theta = Vector<S,3> (0.,0.,0.));
  Vector<S,3> surfaceNormal(const Vector<S,3>& pos, const S meshSize) override;
  const S signedDistance( const PhysR<S,3> input ) override;
  S getVolume( ) override;
  Vector<S,4> calcMofiAndMass(const S density) override;
  Vector<S,3> calcCenterOfMass() override;
};


//implements a custom shaped smooth particle //TODO: Check for consistency
//ALSO: adap .hh
template <typename T, typename S, bool PARTICLE=false>
class SmoothIndicatorCustom3D final: public SmoothIndicatorF3D<T, S, PARTICLE> {
private:
  std::shared_ptr<IndicatorF3D<T>> _indPtr;
  /// Lattice spacing (in m) for the particle lattice (should be smaller or equal to the fluid lattice for best results)
  const T _latticeSpacing;
  /// Local center
  PhysR<T,3> _center;
  /// Block data to store signed distance
  std::unique_ptr<BlockData<3,T,BaseType<T>>> _blockData;
  /// Cached particle volume (to avoid reiteration)
  T _volume;

  void initRotationMatrix();
  void initBlockData(IndicatorF3D<T>& ind);
  void calcCenter();
  void calcCircumRadius();

public:
  // TODO: Add specialized constructors (for PARTICLE = true and PARTICLE = false)
  SmoothIndicatorCustom3D(T latticeSpacing,
                          std::shared_ptr<IndicatorF3D<T>> indPtr, Vector<T,3> pos, T epsilon,
                          Vector<T,3> theta=Vector<T,3>(0.));

  Vector<T,3> getLocalCenter();
  S getVolume( ) override;
  Vector<T,4> calcMofiAndMass(T rhoP);
  Vector<S,3> surfaceNormal(const Vector<S,3>& pos, const S meshSize) override;
  Vector<S,3> surfaceNormal(const Vector<S,3>& pos, const S meshSize,
                            std::function<Vector<S,3>(const Vector<S,3>&)> transformPos) override;
  const S signedDistance( const PhysR<S,3> input ) override;
  bool regardCell(int input[3]);
};


}

#endif

