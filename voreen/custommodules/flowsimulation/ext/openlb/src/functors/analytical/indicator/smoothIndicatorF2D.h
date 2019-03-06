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
#include "functors/lattice/indicator/indicatorBaseF2D.h"
#include "functors/lattice/indicator/indicatorBaseF3D.h"

namespace olb {


///////////////////////////SmoothIndicatorF/////////////////////////////////////

/** implements a smooth cuboid in 2D with an _epsilon sector.
 * \param mass    TODO
 * \param epsilon
 * \param theta   TODO
 *
 * Must be final as getRadius in the constructor does not forward to derived class.
 */
template <typename T, typename S>
class SmoothIndicatorCuboid2D final : public SmoothIndicatorF2D<T,S> {
private:
  S _xLength;
  S _yLength;
public:
  SmoothIndicatorCuboid2D(Vector<S,2> center, S xLength, S yLength, S mass, S epsilon, S theta=0);
  bool operator()(T output[],const S x[]) override;
  /// \return radius of a circle at center, that contains the object
  S getRadius() override;
  S getDiam() override;
  Vector<S,2>& getMin() override;
  Vector<S,2>& getMax() override;
};


/// implements a smooth circle in 2D with an _epsilon sector
template <typename T, typename S>
class SmoothIndicatorCircle2D : public SmoothIndicatorF2D<T,S> {
public:
  SmoothIndicatorCircle2D(Vector<S,2> center, S radius, S mass, S epsilon);
  bool operator() (T output[], const S input[]) override;
  Vector<S,2>& getMin() override;
  Vector<S,2>& getMax() override;
};

/// implements a smooth triangle in 2D with an _epsilon sector
template <typename T, typename S>
class SmoothIndicatorTriangle2D : public SmoothIndicatorF2D<T, S> {
private:
  /// Eckpunkte des Dreiecks
  Vector<S, 2> _PointA, _PointB, _PointC;
  /// Verbindungsvektoren _ab von _A nach _B, etc.
  Vector<S, 2> _ab, _bc, _ca;
  /// normal on _ab * _A  = _ab_d
  S _ab_d, _bc_d, _ca_d;

public:
  SmoothIndicatorTriangle2D(Vector<S,2> center, S radius, S mass, S epsilon, S theta);
  bool operator() (T output[], const S input[]);
  Vector<S,2>& getMin();
  Vector<S,2>& getMax();

};


///////////////////////////ParticleIndicatorF/////////////////////////////////////

/// implements a smooth particle cuboid in 2D with an _epsilon sector.
template <typename T, typename S>
class ParticleIndicatorCuboid2D : public ParticleIndicatorF2D<T,S> {
private:
  S _xLength;
  S _yLength;
public:
  ParticleIndicatorCuboid2D(Vector<S,2> center, S xLength, S yLength, S density, S epsilon, S theta=0);
  bool operator()(T output[],const S x[]);
};

/// implements a smooth particle circle in 2D with an _epsilon sector.
template <typename T, typename S>
class ParticleIndicatorCircle2D : public ParticleIndicatorF2D<T,S> {
private:
  T _radius;
public:
  ParticleIndicatorCircle2D(Vector<S,2> center, S radius, S mass, S epsilon);
  bool operator() (T output[], const S input[]);
};

/** implements a smooth particle triangle in 2D with an _epsilon sector, constructed from circumradius
 * TODO generic constructor with angles
 */
template <typename T, typename S>
class ParticleIndicatorTriangle2D : public ParticleIndicatorF2D<T, S> {
private:
  /// Eckpunkte des Dreiecks
  Vector<S, 2> _PointA, _PointB, _PointC;
  /// Verbindungsvektoren _ab von _A nach _B, etc.
  Vector<S, 2> _ab, _bc, _ca;
  /// normal on _ab * _A  = _ab_d
  S _ab_d, _bc_d, _ca_d;

public:
  ParticleIndicatorTriangle2D(Vector<S,2> center, S radius, S density, S epsilon, S theta);
  bool operator() (T output[], const S input[]);
};

template <typename T, typename S, template<typename U> class DESCRIPTOR>
class ParticleIndicatorCustom2D : public ParticleIndicatorF2D<T, S> {
private:
  // _center is the local center, _startPos the center at the start
  Vector<T,2> _center;
  // _latticeCenter gives the center in local lattice coordinates
  Vector<int,2> _latticeCenter;
  BlockData2D<T, T> _blockData;
  UnitConverter<T,DESCRIPTOR> const& _converter;

public:
  ParticleIndicatorCustom2D(UnitConverter<T,DESCRIPTOR> const& converter, IndicatorF3D<T>& ind, Vector<T,2> center, T rhoP, T epsilon, T theta, T slice);
  bool operator() (T output[], const S input[]);
};


}

#endif

