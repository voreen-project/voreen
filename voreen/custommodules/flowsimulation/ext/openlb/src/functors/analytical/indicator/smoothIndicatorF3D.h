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
#include "functors/lattice/indicator/indicatorBaseF2D.h"
#include "functors/lattice/indicator/indicatorBaseF3D.h"

namespace olb {


template<typename T, typename S> class SmoothIndicatorF3D;

///////////////////////////SmoothIndicatorF/////////////////////////////////////

/// implements a smooth sphere in 3D with an _epsilon sector
template<typename T, typename S>
class SmoothIndicatorSphere3D: public SmoothIndicatorF3D<T, S> {
private:
  S _innerRad;
  S _outerRad;
  S _epsilon;
public:
  SmoothIndicatorSphere3D(Vector<S, 3> center, S radius, S epsilon, S mass);
  SmoothIndicatorSphere3D(const SmoothIndicatorSphere3D<T, S>& rhs);
  bool operator()(T output[], const S input[]) override;
  Vector<S, 3>& getCenter() override;
  S getRadius() override;
  S getDiam() override;
};

/// implements a smooth cylinder in 3D with an _epsilon sector
template <typename T, typename S>
class SmoothIndicatorCylinder3D : public SmoothIndicatorF3D<T,S> {
private:
  Vector<S,3> _center1;
  Vector<S,3> _center2;
  Vector<S,3> _I;
  Vector<S,3> _J;
  Vector<S,3> _K;
  S _length;
  S _radius2;
  S _epsilon;
public:
  SmoothIndicatorCylinder3D(Vector<S,3> center1, Vector<S,3> center2,
                            S radius, S epsilon);
  bool operator() (T output[], const S input[]) override;
};

/// implements a smooth cone in 3D with an _epsilon sector
template <typename T, typename S>
class SmoothIndicatorCone3D : public SmoothIndicatorF3D<T,S> {
private:
  Vector<S,3> _center1;
  Vector<S,3> _center2;
  Vector<S,3> _I;
  Vector<S,3> _J;
  Vector<S,3> _K;
  S _length;
  S _radius1;
  S _radius2; /**< The 2nd radius is optional: if not defined, _center2 is the vertex of the cone */
  S _epsilon;
public:
  SmoothIndicatorCone3D(Vector<S,3> center1, Vector<S,3> center2,
                        S radius1, S radius2, S epsilon);
  bool operator() (T output[], const S input[]) override;
};



///////////////////////////ParticleIndicatorF/////////////////////////////////////

/// implements a smooth sphere in 3D with an _epsilon sector for particle simulations
template<typename T, typename S>
class ParticleIndicatorSphere3D: public ParticleIndicatorF3D<T, S> {
public:
  ParticleIndicatorSphere3D(Vector<S, 3> center, S radius, S epsilon, S mass);
  bool operator()(T output[], const S input[]) override;
};

/** implements a smooth particle cuboid in 3D with an _epsilon sector.
 * TODO construct by density
 * TODO rotation seems weird
 */
template <typename T, typename S>
class ParticleIndicatorCuboid3D : public ParticleIndicatorF3D<T, S> {
private:
  S _xLength;
  S _yLength;
  S _zLength;
public:
  ParticleIndicatorCuboid3D(Vector<S,3> center, S xLength, S yLength, S zLength, S mass, S epsilon, Vector<S,3> theta);
  bool operator()(T output[],const S x[]);
};

/** implements a smooth particle of shape given by in indicator (e.g. STL) in 3D with an _epsilon sector.
 * TODO construct by density
 * TODO check correctness of center and mofi
 */
template <typename T, typename S, template<typename U> class DESCRIPTOR>
class ParticleIndicatorCustom3D : public ParticleIndicatorF3D<T, S> {
private:
  // _center is the local center, _startPos the center at the start
  Vector<T,3> _center;
  // _latticeCenter gives the center in local lattice coordinates
  Vector<int,3> _latticeCenter;
  BlockData3D<T, T> _blockData;
  UnitConverter<T,DESCRIPTOR> const& _converter;

public:
  // for now epsilon has to be chosen to be latticeL
  ParticleIndicatorCustom3D(UnitConverter<T,DESCRIPTOR> const& converter, IndicatorF3D<T>& ind, Vector<T,3> center, T rhoP, T epsilon, Vector<T,3> theta);
  bool operator() (T output[], const S input[]);
};


}

#endif

