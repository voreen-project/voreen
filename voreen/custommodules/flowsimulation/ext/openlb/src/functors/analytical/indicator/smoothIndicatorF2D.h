/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2016 Thomas Henn, Cyril Masquelier, Jan Marquardt, Mathias J. Krause
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

#ifndef SMOOTH_INDICATOR_F_2D_H
#define SMOOTH_INDICATOR_F_2D_H

#include <vector>

#include "smoothIndicatorBaseF2D.h"
#include "io/xmlReader.h"
#include "utilities/functorPtr.h"

#include "core/blockData.h"
#include "core/unitConverter.h"
#include "smoothIndicatorBaseF3D.h"
#include "functors/analytical/indicator/indicatorBaseF2D.h"
#include "functors/analytical/indicator/indicatorBaseF3D.h"
#include "functors/analytical/indicator/indicatorF2D.h"
#include "particles/functions/bodyMotionFunctions.h"
#include "sdf.h"

namespace olb {


///////////////////////////SmoothIndicatorF/////////////////////////////////////

/** implements a smooth cuboid in 2D with an _epsilon sector.
 * \param mass    TODO
 * \param epsilon
 * \param theta   TODO
 *
 */
template <typename T, typename S, bool PARTICLE=false>
class SmoothIndicatorCuboid2D final : public SmoothIndicatorF2D<T,S,PARTICLE> {
private:
  IndicatorCuboid2D<S> _ind;
public:
  SmoothIndicatorCuboid2D(IndicatorCuboid2D<S>& ind, S epsilon, S theta=0);
  SmoothIndicatorCuboid2D(Vector<S,2>center, S xLength, S yLength, S epsilon, S theta=0);
  Vector<S,2> surfaceNormal(const Vector<S,2>& pos, const S meshSize) override;
  const S signedDistance( const Vector<S,2> input ) override;
  bool distance(S& distance, const Vector<S,2>& origin, const Vector<S,2>& direction, S precision, S pitch) override;
  S calcArea( ) override;
  Vector<S,2> calcMofiAndMass(const S density) override;
};

/// implements a smooth circle in 2D with an _epsilon sector
template <typename T, typename S, bool PARTICLE=false>
class SmoothIndicatorCircle2D final : public SmoothIndicatorF2D<T,S,PARTICLE> {
private:
  IndicatorCircle2D<S> _ind;
public:
  SmoothIndicatorCircle2D(IndicatorCircle2D<S>& indPtr, S epsilon);
  SmoothIndicatorCircle2D(Vector<S,2>center, S radius, S epsilon);
  Vector<S,2> surfaceNormal(const Vector<S,2>& pos, const S meshSize) override;
  const S signedDistance( const Vector<S,2> input ) override;
  bool distance(S& distance, const Vector<S,2>& origin, const Vector<S,2>& direction, S precision, S pitch) override;
  S calcArea( ) override;
  Vector<S,2> calcMofiAndMass(const S density) override;
};

/// implements a smooth triangle in 2D with an _epsilon sector
template <typename T, typename S, bool PARTICLE=false>
class SmoothIndicatorTriangle2D final : public SmoothIndicatorF2D<T,S,PARTICLE> {
private:
  IndicatorEquiTriangle2D<S> _ind;
public:
  SmoothIndicatorTriangle2D(IndicatorEquiTriangle2D<S>& indPtr, S epsilon, S theta=0);
  SmoothIndicatorTriangle2D(Vector<S,2>center, S radius, S epsilon, S theta=0);
  Vector<S,2> surfaceNormal(const Vector<S,2>& pos, const S meshSize) override;
  const S signedDistance( const Vector<S,2> input ) override;
  bool distance(S& distance, const Vector<S,2>& origin, const Vector<S,2>& direction, S precision, S pitch) override;
  S calcArea( ) override;
  Vector<S,2> calcMofiAndMass(const S density) override;
};

/*
/// implements a smooth triangle in 2D with an _epsilon sector
template <typename T, typename S, bool PARTICLE=false>
class SmoothIndicatorTriangle2D final : public SmoothIndicatorF2D<T,S,PARTICLE> {
private:
  /// corner points
  Vector<S, 2> _PointA, _PointB, _PointC;
  /// vectors connecting corner points (_ab: from a to b)
  Vector<S, 2> _ab, _bc, _ca;
  /// normal on _ab * _A  = _ab_d
  S _ab_d, _bc_d, _ca_d;
public:
  SmoothIndicatorTriangle2D(Vector<S,2> center, S radius, S epsilon, S theta=0, S density=0, Vector<S,2> vel = Vector<S,2> (0.,0.), S omega = 0);
  bool operator() (T output[], const S input[]) override;
};
*/

//implements a custom shaped smooth particle //TODO: Check for consistency
//ALSO: adap .hh
template <typename T, typename S, bool PARTICLE=false>
class SmoothIndicatorCustom2D : public SmoothIndicatorF2D<T, S, PARTICLE> {
private:
  std::shared_ptr<IndicatorF2D<T>> _indPtr;
  /// Lattice spacing (in m) for the particle lattice (should be smaller or equal to the fluid lattice for best results)
  const T _latticeSpacing;
  // _center is the local center
  Vector<T,2> _center;
  // _latticeCenter gives the center in local lattice coordinates
  Vector<int,2> _latticeCenter;
  /// Block data to store signed distance
  std::unique_ptr<BlockData<2,T,BaseType<T>>> _blockData;
  /// Cached particle volume (to avoid reiteration)
  T _area;

  void initData(IndicatorF2D<T>& ind);
  void initRotationMatrix();
  void initBlockData(IndicatorF2D<T>& ind);
  void calcCenter();
  void calcCircumRadius();

public:
  SmoothIndicatorCustom2D(T latticeSpacing,
                          std::shared_ptr<IndicatorF2D<T>> indPtr, Vector<T,2> pos, T epsilon,
                          T theta=0.);

  Vector<T,2> getLocalCenter();
  S calcArea( ) override;
  Vector<T,2> calcMofiAndMass(T rhoP) override;
  Vector<S,2> surfaceNormal(const Vector<S,2>& pos, const S meshSize) override;
  const S signedDistance( const Vector<S,2> input ) override;
  bool regardCell(int input[2]);
};

//Geng2019:
/// implements a smooth circle in 2D with an tangiant _epsilon sector
template <typename T, typename S, bool PARTICLE=false>
class SmoothIndicatorHTCircle2D final : public SmoothIndicatorF2D<T,S,PARTICLE> {
private:
  S _radius;
public:
  SmoothIndicatorHTCircle2D(Vector<S,2> center, S radius, S epsilon, S density=0, Vector<S,2> vel = Vector<S,2> (0.,0.), S omega = 0);
  bool operator() (T output[], const S input[]) override;
  Vector<S,2> calcMofiAndMass(const S density) override;
};

}

#endif

