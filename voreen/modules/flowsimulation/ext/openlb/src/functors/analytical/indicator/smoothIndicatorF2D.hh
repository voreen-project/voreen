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

#ifndef SMOOTH_INDICATOR_F_2D_HH
#define SMOOTH_INDICATOR_F_2D_HH

#include "utilities/omath.h"

#include "smoothIndicatorF2D.h"
#include "functors/lattice/reductionF3D.h"
#include "utilities/vectorHelpers.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define M_PI2 1.57079632679489661923

namespace olb {

template <typename T, typename S, bool PARTICLE>
SmoothIndicatorCuboid2D<T,S,PARTICLE>::SmoothIndicatorCuboid2D(IndicatorCuboid2D<S>& ind, S epsilon, S theta)
  : SmoothIndicatorCuboid2D<T,S,PARTICLE>(ind.getCenter(), ind.getxLength(), ind.getyLength(), epsilon, theta)
{ }

template <typename T, typename S, bool PARTICLE>
SmoothIndicatorCuboid2D<T,S,PARTICLE>::SmoothIndicatorCuboid2D(Vector<S,2>center, S xLength, S yLength, S epsilon, S theta)
  : _ind(xLength, yLength, center)
{
  this->_epsilon = epsilon;
  if constexpr (!PARTICLE) {
    this->_pos =  _ind.getCenter();
    this->_theta = theta * M_PI/180.;
  }

  this->_circumRadius = .5*(util::sqrt(util::pow( xLength, 2)+util::pow( yLength, 2))) + 0.5*epsilon;
  if constexpr (!PARTICLE) {
    this->_myMin = {
      this->_pos[0] - this->getCircumRadius(),
      this->_pos[1] - this->getCircumRadius()
    };
    this->_myMax = {
      this->_pos[0] + this->getCircumRadius(),
      this->_pos[1] + this->getCircumRadius()
    };
    this->init(this->_theta);
  }
}

template <typename T, typename S, bool PARTICLE>
S SmoothIndicatorCuboid2D<T, S, PARTICLE>::calcArea( )
{
  return _ind.getxLength()*_ind.getyLength();
}

template <typename T, typename S, bool PARTICLE>
Vector<S,2> SmoothIndicatorCuboid2D<T, S, PARTICLE>::calcMofiAndMass( const S density )
{
  const T xLength = _ind.getxLength();
  const T yLength = _ind.getyLength();
  const T mass = calcArea()*density;
  const T mofi = 1./12.*mass*(xLength*xLength+yLength*yLength);
  return Vector<S,2>(mofi, mass);
}

template <typename T, typename S, bool PARTICLE>
Vector<S,2> SmoothIndicatorCuboid2D<T, S, PARTICLE>::surfaceNormal( const Vector<S,2>& pos, const S meshSize )
{
  return _ind.surfaceNormal(pos, meshSize);
}

template <typename T, typename S, bool PARTICLE>
const S SmoothIndicatorCuboid2D<T, S, PARTICLE>::signedDistance( const Vector<S,2> input )
{
  Vector<S,2> p;
  if constexpr(!PARTICLE) {
    // rotation & translation
    p = body::motion::rotation<2,S,true>::execute(input, this->_rotMat, this->getPos());
  }
  else {
    p = input;
  }
  return _ind.signedDistance(p + _ind.getCenter());
}

template <typename T, typename S, bool PARTICLE>
bool SmoothIndicatorCuboid2D<T, S, PARTICLE>::distance(S& distance, const Vector<S,2>& origin, const Vector<S,2>& direction, S precision, S pitch)
{
  Vector<S,2> p;
  if constexpr(!PARTICLE) {
    // rotation & translation
    p = body::motion::rotation<2,S,true>::execute(origin, this->_rotMat, this->getPos());
  }
  else {
    p = origin;
  }
  return _ind.distance(distance, p, direction, precision, pitch);
}


template <typename T, typename S, bool PARTICLE>
SmoothIndicatorCircle2D<T,S,PARTICLE>::SmoothIndicatorCircle2D(IndicatorCircle2D<S>& ind, S epsilon)
  : SmoothIndicatorCircle2D(ind.getCenter(), ind.getRadius(), epsilon)
{ }

template <typename T, typename S, bool PARTICLE>
SmoothIndicatorCircle2D<T,S,PARTICLE>::SmoothIndicatorCircle2D(Vector<S,2>center, S radius, S epsilon)
  : _ind(center, radius)
{
  this->_epsilon = epsilon;
  if constexpr (!PARTICLE) {
    this->_pos = _ind.getCenter();
  }

  this->_circumRadius = radius + 0.5*epsilon;
  if constexpr (!PARTICLE) {
    this->_myMin = {this->_pos[0] - this->getCircumRadius(), this->_pos[1] - this->getCircumRadius()};
    this->_myMax = {this->_pos[0] + this->getCircumRadius(), this->_pos[1] + this->getCircumRadius()};
    this->init(0.);
  }
}

template <typename T, typename S, bool PARTICLE>
S SmoothIndicatorCircle2D<T, S, PARTICLE>::calcArea( )
{
  return M_PI*_ind.getRadius()*_ind.getRadius();
}

template <typename T, typename S, bool PARTICLE>
Vector<S,2> SmoothIndicatorCircle2D<T, S, PARTICLE>::calcMofiAndMass( const S density )
{
  const T radius = _ind.getRadius();
  const T mass = calcArea()*density;
  const T mofi = 0.5 * mass * radius * radius;
  return Vector<S,2>(mofi, mass);
}

template <typename T, typename S, bool PARTICLE>
Vector<S,2> SmoothIndicatorCircle2D<T, S, PARTICLE>::surfaceNormal( const Vector<S,2>& pos, const S meshSize )
{
  return _ind.surfaceNormal(pos, meshSize);
}

template <typename T, typename S, bool PARTICLE>
const S SmoothIndicatorCircle2D<T, S, PARTICLE>::signedDistance( const Vector<S,2> input )
{
  Vector<S,2> dist = input;
  if constexpr(!PARTICLE) {
    dist -= this->_pos;
  }
  return _ind.signedDistance(dist + _ind.getCenter());
}

template <typename T, typename S, bool PARTICLE>
bool SmoothIndicatorCircle2D<T, S, PARTICLE>::distance(S& distance, const Vector<S,2>& origin, const Vector<S,2>& direction, S precision, S pitch)
{
  Vector<S,2> dist = origin;
  if constexpr(!PARTICLE) {
    dist -= this->_pos;
  }
  return _ind.distance(distance, dist, direction, precision, pitch);
}


template <typename T, typename S, bool PARTICLE>
SmoothIndicatorTriangle2D<T,S,PARTICLE>::SmoothIndicatorTriangle2D(IndicatorEquiTriangle2D<S>& ind, S epsilon, S theta)
  : SmoothIndicatorTriangle2D(ind.getCenter(), ind.getRadius(), epsilon, theta)
{ }

template <typename T, typename S, bool PARTICLE>
SmoothIndicatorTriangle2D<T,S,PARTICLE>::SmoothIndicatorTriangle2D(Vector<S,2>center, S radius, S epsilon, S theta)
  : _ind(center, radius)
{
  this->_epsilon = epsilon;
  if constexpr (!PARTICLE) {
    this->_pos = _ind.getCenter();
    this->_theta = theta * M_PI/180.;
  }

  this->_circumRadius = radius + 0.5*epsilon;
  if constexpr (!PARTICLE) {
    this->_myMin = {center[0] - this->getCircumRadius(), center[1] - this->getCircumRadius()};
    this->_myMax = {center[0] + this->getCircumRadius(), center[1] + this->getCircumRadius()};
    this->init(this->_theta);
  }
}

template <typename T, typename S, bool PARTICLE>
S SmoothIndicatorTriangle2D<T, S, PARTICLE>::calcArea( )
{
  const T radius = _ind.getRadius();
  const T altitude = 1.5*radius;
  const T base = util::sqrt(3)*radius;
  return 0.5*base*altitude;
}

template <typename T, typename S, bool PARTICLE>
Vector<S,2> SmoothIndicatorTriangle2D<T, S, PARTICLE>::calcMofiAndMass( const S density )
{
  const T radius = _ind.getRadius();
  const T altitude = 1.5*radius;
  const T base = util::sqrt(3)*radius;
  const T mass = density*calcArea();
  const T mofi = mass*((altitude*altitude/18.)+(base*base/24.));
  return Vector<S,2>(mofi, mass);
}

template <typename T, typename S, bool PARTICLE>
Vector<S,2> SmoothIndicatorTriangle2D<T, S, PARTICLE>::surfaceNormal( const Vector<S,2>& pos, const S meshSize )
{
  return _ind.surfaceNormal(pos, meshSize);
}

template <typename T, typename S, bool PARTICLE>
const S SmoothIndicatorTriangle2D<T, S, PARTICLE>::signedDistance( const Vector<S,2> input )
{
  Vector<S,2> p;
  if constexpr(!PARTICLE) {
    // rotation & translation
    p = body::motion::rotation<2,S,true>::execute(input, this->_rotMat, this->getPos());
  }
  else {
    p = input;
  }
  return _ind.signedDistance(p + _ind.getCenter());
}

template <typename T, typename S, bool PARTICLE>
bool SmoothIndicatorTriangle2D<T, S, PARTICLE>::distance(S& distance, const Vector<S,2>& origin, const Vector<S,2>& direction, S precision, S pitch)
{
  Vector<S,2> p;
  if constexpr(!PARTICLE) {
    // rotation
    p = body::motion::rotation<2,S>::execute(origin, this->_rotMat, this->getPos());
    // translation
    p -= this->getPos();
  }
  else {
    p = origin;
  }
  return _ind.distance(distance, p, direction, precision, pitch);
}


//TODO: TO Be Repaired
//TODO: Check for consitency
template <typename T, typename S, bool PARTICLE>
SmoothIndicatorCustom2D<T,S,PARTICLE>::SmoothIndicatorCustom2D(T latticeSpacing,
    std::shared_ptr<IndicatorF2D<T>> indPtr,
    Vector<T,2> pos,
    T epsilon,
    T theta)
  :_indPtr(indPtr),
   _latticeSpacing(latticeSpacing)
{
  OstreamManager clout(std::cout,"createIndicatorCustom2D");
  this->_name = "custom2D";
  this->_epsilon = epsilon;
  if constexpr(!PARTICLE) {
    this->_pos = pos;         // global position of the local center
    this->_theta = theta * (M_PI/180.);
  }

  initData(*_indPtr);
}

template <typename T, typename S, bool PARTICLE>
void SmoothIndicatorCustom2D<T,S,PARTICLE>::initData(IndicatorF2D<T>& ind)
{
  initRotationMatrix();
  initBlockData(ind);

  // calculate mass and centerpoint for rotation
  calcCenter();
  // calculate min and max from circumRadius
  calcCircumRadius();
}

template <typename T, typename S, bool PARTICLE>
void SmoothIndicatorCustom2D<T,S,PARTICLE>::initRotationMatrix()
{
  // initialize rotation matrix
  if constexpr (!PARTICLE) {
    this->_rotMat[0] = util::cos(this->_theta);
    this->_rotMat[1] = util::sin(this->_theta);
    this->_rotMat[2] = -util::sin(this->_theta);
    this->_rotMat[3] = util::cos(this->_theta);
  }
}

template <typename T, typename S, bool PARTICLE>
void SmoothIndicatorCustom2D<T,S,PARTICLE>::initBlockData(IndicatorF2D<T>& ind)
{
  OstreamManager clout(std::cout,"createIndicatorCustom2D");

  // initialize temporary values
  int blockDataSize[3];
  int blockDataPadding[2];
  for (unsigned iD=0; iD<2; ++iD) {
    blockDataSize[iD] = util::ceil( (ind.getMax()[iD] - ind.getMin()[iD]) / _latticeSpacing );
    // Add a padding so that the distance can be cached in the vicinity of the geometry
    blockDataPadding[iD] = 2*util::ceil(0.05*blockDataSize[iD]);
    blockDataSize[iD] += blockDataPadding[iD];
  }
  blockDataSize[2] = 1;

  // create blockData containing signed distance information
  this->_blockData.reset(new BlockData<2,T,BaseType<T>>({blockDataSize, 0}));
  int iX[2];
  for (iX[0]=0; iX[0] < blockDataSize[0]; iX[0]++) {
    for (iX[1]=0; iX[1] < blockDataSize[1]; iX[1]++) {
      Vector<T,2> input;
      for (unsigned iD=0; iD<2; ++iD) {
        input[iD] = (iX[iD]-blockDataPadding[iD]/2)*_latticeSpacing+ind.getMin()[iD];
      }
      this->_blockData->get(iX) = ind.signedDistance(input);
    }
  }
}

template <typename T, typename S, bool PARTICLE>
void SmoothIndicatorCustom2D<T,S,PARTICLE>::calcCenter()
{
  // TODO check again for correctness of center due to smooth boundary and coordinate system
  unsigned nCells = 0;
  int input[2];
  this->_center = Vector<T,2>(0.);
  for (input[0] = 0; input[0] < this->_blockData->getNx(); ++input[0]) {
    for (input[1] = 0; input[1] < this->_blockData->getNy(); ++input[1]) {
      if (regardCell(input)) {
        // Always use real boundary as in other geometries too
        const unsigned short porosity = 1;
        this->_center[0] += porosity*this->_latticeSpacing*input[0];
        this->_center[1] += porosity*this->_latticeSpacing*input[1];
        nCells += porosity;
      }
    }
  }
  this->_center *= 1./nCells;
}

template <typename T, typename S, bool PARTICLE>
S SmoothIndicatorCustom2D<T,S,PARTICLE>::calcArea( )
{
  return _area;
}

/// calculates and returns mofi and mass of the particle (mofi is at index 0 and mass at index 1)
template <typename T, typename S, bool PARTICLE>
Vector<T,2> SmoothIndicatorCustom2D<T,S,PARTICLE>::calcMofiAndMass(T rhoP)
{
  // TODO - calculation
  T mofi = 0.;
  T mass;
  unsigned nCells = 0;
  int input[2];
  for (input[0] = 0; input[0] < this->_blockData->getNx(); ++input[0]) {
    const T dx = util::abs(input[0]*_latticeSpacing - this->_center[0]);
    for (input[1] = 0; input[1] < this->_blockData->getNy(); ++input[1]) {
      if (regardCell(input)) {
        const T dy = util::abs(input[1]*_latticeSpacing - this->_center[1]);
        mofi += (dx*dx+dy*dy);
        ++nCells;
      }
    }
  }
  _area = nCells * util::pow(_latticeSpacing, 2);
  mass = rhoP * _area;
  mofi += util::pow(_latticeSpacing, 4)/ 6.0;
  mofi *= mass/nCells;

  return Vector<T,2>(mofi,mass);
}

template <typename T, typename S, bool PARTICLE>
void SmoothIndicatorCustom2D<T,S,PARTICLE>::calcCircumRadius()
{
  Vector<T,2> min(std::numeric_limits<olb::BaseType<T>>::max());
  Vector<T,2> max(-std::numeric_limits<olb::BaseType<T>>::max());
  Vector<T,2> distance;
  T maxDistance = -std::numeric_limits<olb::BaseType<T>>::max();

  int input[2];
  for (input[0] = 0; input[0] < this->_blockData->getNx(); ++input[0]) {
    distance[0] = this->_center[0] - this->_latticeSpacing*input[0];
    for (input[1] = 0; input[1] < this->_blockData->getNy(); ++input[1]) {
      distance[1] = this->_center[1] - this->_latticeSpacing*input[1];
      if (regardCell(input)) {
        if constexpr (!PARTICLE) {
          for (unsigned iD=0; iD<2; ++iD) {
            min[iD] = util::min(distance[iD], min[iD]);
            max[iD] = util::max(distance[iD], max[iD]);
          }
        }
        maxDistance = util::max(norm(distance), maxDistance);
      }
    }
  }

  if constexpr (!PARTICLE) {
    min -= Vector<T,2>(-this->_epsilon);
    max += Vector<T,2>(-this->_epsilon);
    this->_myMin = this->_pos+min;
    this->_myMax = this->_pos+max;
  }

  this->_circumRadius = maxDistance + .5*this->_epsilon;
}

template <typename T, typename S, bool PARTICLE>
Vector<T,2> SmoothIndicatorCustom2D<T,S,PARTICLE>::getLocalCenter()
{
  return this->_center;
}

template <typename T, typename S, bool PARTICLE>
Vector<S,2> SmoothIndicatorCustom2D<T,S,PARTICLE>::surfaceNormal( const Vector<S,2>& pos, const S meshSize )
{
  return _indPtr->surfaceNormal(pos, meshSize);
}

template <typename T, typename S, bool PARTICLE>
const S SmoothIndicatorCustom2D<T,S,PARTICLE>::signedDistance( const Vector<S,2> input )
{
  // Translation
  T xDist = input[0];
  T yDist = input[1];
  if constexpr (!PARTICLE) {
    xDist -= this->getPos()[0];
    yDist -= this->getPos()[1];
  }

  // counter-clockwise rotation by _theta=-theta around (0/0) and movement from rotation center to local center
  int x,y;
  if constexpr(PARTICLE) {
    x = (xDist + this->_center[0]) / this->_latticeSpacing;
    y = (yDist + this->_center[1]) / this->_latticeSpacing;
  }
  else {
    x = (xDist*this->_rotMat[0] + yDist*this->_rotMat[2] + this->_center[0]) / this->_latticeSpacing;
    y = (xDist*this->_rotMat[1] + yDist*this->_rotMat[3] + this->_center[1]) / this->_latticeSpacing;
  }

  if (x >= 0 && x < this->_blockData->getNx() && y >= 0 && y < this->_blockData->getNy()) {
    return this->_blockData->get({x, y});
  }
  return 1;
}

template <typename T, typename S, bool PARTICLE>
bool SmoothIndicatorCustom2D<T,S,PARTICLE>::regardCell(int input[2])
{
  return this->_blockData->get(input) < std::numeric_limits<T>::epsilon();
}


//Geng2019:
template <typename T, typename S, bool PARTICLE>
SmoothIndicatorHTCircle2D<T,S,PARTICLE>::SmoothIndicatorHTCircle2D(Vector<S,2> center, S radius, S epsilon, S density, Vector<S,2> vel, S omega)
  : _radius(radius)
{
  this->_epsilon = epsilon;
  if constexpr (!PARTICLE) {
    this->_pos = center;
  }

  this->_circumRadius = radius + 0.5*epsilon;
  if constexpr (!PARTICLE) {
    this->_myMin = {center[0] - this->getCircumRadius(), center[1] - this->getCircumRadius()};
    this->_myMax = {center[0] + this->getCircumRadius(), center[1] + this->getCircumRadius()};
    this->init(0.);
  }

  T mass = M_PI*radius*radius*density;
  T mofi = 0.5 * mass * radius * radius;
}

template <typename T, typename S, bool PARTICLE>
Vector<S,2> SmoothIndicatorHTCircle2D<T, S, PARTICLE>::calcMofiAndMass( const S density )
{
  T mass = M_PI*_radius*_radius*density;
  T mofi = 0.5 * mass * _radius * _radius;
  return Vector<S,2>(mofi, mass);
}

// returns true if x is inside the sphere
template <typename T, typename S, bool PARTICLE>
bool SmoothIndicatorHTCircle2D<T,S,PARTICLE>::operator()(T output[], const S input[])
{
  double distToCenter2 = util::pow(this->getPos()[0]-input[0], 2) +
                         util::pow(this->getPos()[1]-input[1], 2);


  double d = util::sqrt(distToCenter2) - this->_radius;
  output[0] = T((1.-tanh(d/this->getEpsilon()))/2.);
  return true;
}


} // namespace olb

#endif
