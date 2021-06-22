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

#include "core/blockData2D.h"
#include "core/unitConverter.h"
#include "smoothIndicatorBaseF3D.h"
#include "functors/analytical/indicator/indicatorBaseF2D.h"
#include "functors/analytical/indicator/indicatorBaseF3D.h"

namespace olb {


///////////////////////////SmoothIndicatorF/////////////////////////////////////

/** implements a smooth cuboid in 2D with an _epsilon sector.
 * \param mass    TODO
 * \param epsilon
 * \param theta   TODO
 *
 */
template <typename T, typename S, bool HLBM=false>
class SmoothIndicatorCuboid2D final : public SmoothIndicatorF2D<T,S,HLBM> {
private:
  S _xLength;
  S _yLength;
public:
  SmoothIndicatorCuboid2D(Vector<S,2> center, S xLength, S yLength, S epsilon, S theta=0, S density=0, Vector<S,2> vel = Vector<S,2> (0.,0.));
  bool operator()(T output[],const S input[]) override;
};

/// implements a smooth circle in 2D with an _epsilon sector
template <typename T, typename S, bool HLBM=false>
class SmoothIndicatorCircle2D final : public SmoothIndicatorF2D<T,S,HLBM> {
private:
  S _radius;
public:
  SmoothIndicatorCircle2D(Vector<S,2> center, S radius, S epsilon, S density=0, Vector<S,2> vel = Vector<S,2> (0.,0.));
  bool operator() (T output[], const S input[]) override;
};


/// implements a smooth triangle in 2D with an _epsilon sector
template <typename T, typename S, bool HLBM=false>
class SmoothIndicatorTriangle2D final : public SmoothIndicatorF2D<T,S,HLBM> {
private:
  /// corner points
  Vector<S, 2> _PointA, _PointB, _PointC;
  /// vectors connecting corner points (_ab: from a to b)
  Vector<S, 2> _ab, _bc, _ca;
  /// normal on _ab * _A  = _ab_d
  S _ab_d, _bc_d, _ca_d;
public:
  SmoothIndicatorTriangle2D(Vector<S,2> center, S radius, S epsilon, S theta=0, S density=0, Vector<S,2> vel = Vector<S,2> (0.,0.));
  bool operator() (T output[], const S input[]) override;
};
// needs to be updated to current state
/*
/// implements a custom shaped smooth particle //TODO: Check for consistency
template <typename T, typename S, template<typename U> class DESCRIPTOR, bool HLBM=false>
class SmoothIndicatorCustom2D final : public SmoothIndicatorF2D<T,S,HLBM> {
private:
  // _center is the local center, _startPos the center at the start
  Vector<T,2> _center;
  // _latticeCenter gives the center in local lattice coordinates
  Vector<int,2> _latticeCenter;
  BlockData2D<T, T> _blockData;
  UnitConverter<T,DESCRIPTOR> const& _converter;
public:
  SmoothIndicatorCustom2D(UnitConverter<T,DESCRIPTOR> const& converter, IndicatorF3D<T>& ind, Vector<T,2> center, T epsilon, T slice, S theta=0, S density=0, Vector<S,2> vel = Vector<S,2> (0.,0.));
  bool operator() (T output[], const S input[]);
};
*/

//Geng2019:
/// implements a smooth circle in 2D with an tangiant _epsilon sector 
template <typename T, typename S, bool HLBM=false>
class SmoothIndicatorHTCircle2D final : public SmoothIndicatorF2D<T,S,HLBM> {
private:
  S _radius;
public:
  SmoothIndicatorHTCircle2D(Vector<S,2> center, S radius, S epsilon, S density=0, Vector<S,2> vel = Vector<S,2> (0.,0.));
  bool operator() (T output[], const S input[]) override;
};

}

#endif

