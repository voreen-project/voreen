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

#include <cmath>

#include "smoothIndicatorF2D.h"
#include "functors/lattice/reductionF3D.h"
#include "utilities/vectorHelpers.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define M_PI2 1.57079632679489661923

namespace olb {

template <typename T, typename S, bool HLBM>
SmoothIndicatorCuboid2D<T,S,HLBM>::SmoothIndicatorCuboid2D(Vector<S,2> center, S xLength, S yLength, S epsilon, S theta, S density, Vector<S,2> vel)
  : _xLength(xLength),_yLength(yLength)
{
  this->_pos = center;
  this->_circumRadius = .5*(std::sqrt(std::pow(_xLength, 2)+std::pow(_yLength, 2))) + 0.5*epsilon;
  this->_myMin = {
    center[0] - this->getCircumRadius(),
    center[1] - this->getCircumRadius()
  };
  this->_myMax = {
    center[0] + this->getCircumRadius(),
    center[1] + this->getCircumRadius()
  };
  this->_epsilon = epsilon;
  this->_theta = theta * M_PI/180.;
  T mass = xLength*yLength*density;
  T mofi = 1./12.*mass*(xLength*xLength+yLength*yLength);
  this->init(theta, vel, mass, mofi);
}

template <typename T, typename S, bool HLBM>
bool SmoothIndicatorCuboid2D<T,S,HLBM>::operator()(T output[], const S input[])
{
  T xDist = input[0] - this->getPos()[0];
  T yDist = input[1] - this->getPos()[1];

  // counter-clockwise rotation by _theta=-theta around center
  T x = this->getPos()[0] + xDist*this->getRotationMatrix()[0] + yDist*this->getRotationMatrix()[2];
  T y = this->getPos()[1] + xDist*this->getRotationMatrix()[1] + yDist*this->getRotationMatrix()[3];

  xDist = fabs(x-this->getPos()[0]) - 0.5*(_xLength-this->getEpsilon());
  yDist = fabs(y-this->getPos()[1]) - 0.5*(_yLength-this->getEpsilon());

  if ( xDist <= 0 && yDist <= 0) {
    output[0] = 1.;
    return true;
  }
  if ( xDist >= this->getEpsilon() || yDist >= this->getEpsilon() ) {
    output[0] = 0.;
    return false;
  }
  // treating edges and corners as "rounded"
  if ( (xDist < this->getEpsilon() && xDist > 0) && (yDist < this->getEpsilon() && yDist > 0) ) {
    output[0] = T( (std::pow(cos(M_PI2*xDist/this->getEpsilon()), 2) *
                    std::pow(cos(M_PI2*yDist/this->getEpsilon()), 2)) );
    return true;
  }
  if ( xDist < this->getEpsilon() && xDist > 0 ) {
    output[0] = T( std::pow(cos(M_PI2*xDist/this->getEpsilon()), 2));
    return true;
  }
  if ( yDist < this->getEpsilon() && yDist > 0 ) {
    output[0] = T( std::pow(cos(M_PI2*yDist/this->getEpsilon()), 2));
    return true;
  }
  output[0] = 0.;
  return false;
}

template <typename T, typename S, bool HLBM>
SmoothIndicatorCircle2D<T,S,HLBM>::SmoothIndicatorCircle2D(Vector<S,2> center, S radius, S epsilon, S density, Vector<S,2> vel)
  : _radius(radius)
{
  this->_pos = center;
  this->_circumRadius = radius + 0.5*epsilon;
  this->_myMin = {center[0] - this->getCircumRadius(), center[1] - this->getCircumRadius()};
  this->_myMax = {center[0] + this->getCircumRadius(), center[1] + this->getCircumRadius()};
  this->_epsilon = epsilon;
  T mass = M_PI*radius*radius*density;
  T mofi = 0.5 * mass * radius * radius;
  this->init(0., vel, mass, mofi);
}

// returns true if x is inside the sphere
template <typename T, typename S, bool HLBM>
bool SmoothIndicatorCircle2D<T,S,HLBM>::operator()(T output[], const S input[])
{
  double distToCenter2 = std::pow(this->getPos()[0]-input[0], 2) +
                         std::pow(this->getPos()[1]-input[1], 2);
  if ( distToCenter2 >= std::pow(this->_radius + 0.5*this->getEpsilon(), 2)) {
    output[0] = 0.;
    return false;
  } else if ( distToCenter2 <= std::pow(this->_radius - 0.5*this->getEpsilon(), 2)) {
    output[0] = 1.;
    return true;
  } else {
    // d is between 0 and _epsilon
    double d = std::sqrt(distToCenter2) - this->_radius + 0.5*this->getEpsilon();
    output[0] = T(std::pow(cos(M_PI2*d/this->getEpsilon()), 2));
    return true;
  }
  return false;
}

template <typename T, typename S, bool HLBM>
SmoothIndicatorTriangle2D<T,S, HLBM>::SmoothIndicatorTriangle2D(Vector<S,2> center, S radius, S epsilon, S theta, S density, Vector<S,2> vel)
{
  this->_pos = center;
  this->_circumRadius = radius + 0.5*epsilon;
  this->_myMin = {center[0] - this->getCircumRadius(), center[1] - this->getCircumRadius()};
  this->_myMax = {center[0] + this->getCircumRadius(), center[1] + this->getCircumRadius()};
  this->_epsilon = epsilon;
  this->_theta = theta * M_PI/180.;

  T smallRad = radius * .5;    //sin(30)
  T halfEdge = radius * std::sqrt(3)/2.; // cos(30)
  T altitude = 1.5*radius;
  T base = std::sqrt(3)*radius;
  _PointA[0] = 0.;
  _PointA[1] = radius;
  _PointB[0] = - halfEdge;
  _PointB[1] = - smallRad;
  _PointC[0] = halfEdge;
  _PointC[1] = - smallRad;

  T invEps = 1./this->getEpsilon();

  _ab = _PointB - _PointA;
  _ab = normalize(_ab, invEps);
  _ab_d = _ab[1]*_PointA[0] - _ab[0]*_PointA[1];
  _bc = _PointC - _PointB;
  _bc = normalize(_bc, invEps);
  _bc_d = _bc[1]*_PointB[0] - _bc[0]*_PointB[1];
  _ca = _PointA - _PointC;
  _ca = normalize(_ca, invEps);
  _ca_d = _ca[1]*_PointC[0] - _ca[0]*_PointC[1];

  T mass = density*0.5*base*altitude;
  T mofi = mass*((altitude*altitude/18.)+(base*base/24.));
  this->init(theta, vel, mass, mofi);
}

// returns true if x is inside the sphere
// TODO: check if epsilon treatment is correct (currently epsilon not around but after the boundary)
template <typename T, typename S, bool HLBM>
bool SmoothIndicatorTriangle2D<T,S,HLBM>::operator()(T output[], const S input[])
{
  T xDist = input[0] - this->getPos()[0];
  T yDist = input[1] - this->getPos()[1];

  T x = xDist*this->getRotationMatrix()[0] + yDist*this->getRotationMatrix()[2];
  T y = xDist*this->getRotationMatrix()[1] + yDist*this->getRotationMatrix()[3];

  unsigned short area = 0;

  T dist_a = _bc[1]*x-_bc[0]*y - _bc_d;
  T dist_b = _ca[1]*x-_ca[0]*y - _ca_d;
  T dist_c = _ab[1]*x-_ab[0]*y - _ab_d;

  if (dist_c < 0) {
    area = (area | 100);
  }
  if (dist_a < 0) {
    area = (area | 10);
  }
  if (dist_b < 0) {
    area = (area | 1);
  }

  if (area == 111) {
    output[0] = 1.;
    return true;
  }

  if (area == 110 && dist_b < 1) {
    output[0] = T(std::pow(cos(M_PI2*dist_b), 2));
    return true;
  }

  if (area == 101 && dist_a < 1) {
    output[0] = T(std::pow(cos(M_PI2*dist_a), 2));
    return true;
  }

  if (area == 11 && dist_c < 1) {
    output[0] = T(std::pow(cos(M_PI2*dist_c), 2));
    return true;
  }

  if (area == 1 && dist_a < 1 && dist_c < 1) {
    output[0] = T(std::pow(cos(M_PI2*dist_a), 2)*std::pow(cos(M_PI2*dist_c), 2));
    return true;
  }

  if (area == 10 && dist_b < 1 && dist_c < 1) {
    output[0] = T(std::pow(cos(M_PI2*dist_b), 2)*std::pow(cos(M_PI2*dist_c), 2));
    return true;
  }

  if (area == 100 && dist_b < 1 && dist_a < 1) {
    output[0] = T(std::pow(cos(M_PI2*dist_b), 2)*std::pow(cos(M_PI2*dist_a), 2));
    return true;
  }

  output[0] = 0.;
  return false;
}
//needs to be updated to current state
/*
//TODO: Check for consitency
template <typename T, typename S, typename DESCRIPTOR, bool HLBM>
SmoothIndicatorCustom2D<T,S,DESCRIPTOR,HLBM>::SmoothIndicatorCustom2D(UnitConverter<T,DESCRIPTOR> const& converter,
    IndicatorF3D<T>& ind,
    Vector<T,2> center,
    T epsilon,
    T slice,
    S theta,
    S rhoP,
    Vector<S,2> vel)
  : _converter(converter)
{
  OstreamManager clout(std::cout,"createIndicatorCustom2D");
  this->_pos = center;
  this->_vel = vel;
  this->_epsilon = epsilon;
  this->_theta = theta;

  // initialize temporary values
  SmoothBlockIndicator3D<T,olb::descriptors::D3Q19<>> smoothBlock(ind, this->_epsilon);
  int _nX = smoothBlock.getBlockData().getNx();
  int _nY = smoothBlock.getBlockData().getNy();
  int tmpNcells = 0;
  if (slice<ind.getMin()[2] || slice>ind.getMax()[2]) {
    clout << "ERROR: Forbidden value, slice out of bounds. Value needs to be between " << ind.getMin()[2] << " and " << ind.getMax()[2] << std::endl;
    return;
  }

  // create smoothed blockData
  int tmpZ = int(slice/this->_converter.getConversionFactorLength());
  BlockData2D<T,BaseType> block_tmp(_nX, _nY);
  for (int iX=0; iX < _nX; iX++) {
    for (int iY=0; iY < _nY; iY++) {
      block_tmp.get(iX, iY) = smoothBlock.getBlockData().get(iX, iY, tmpZ);
      if (block_tmp.get(iX, iY) > 0) {
        tmpNcells++;
      }
    }
  }
  this->_blockData = block_tmp;
  T invNcells = 1./tmpNcells;

  // calculate mass and centerpoint for rotation
  this->_mass = rhoP * tmpNcells * std::pow(_converter.getConversionFactorLength(), 2);
  this->_center[0] = 0.0;
  this->_center[1] = 0.0;
  for (int iX= 0; iX < _nX; iX++) {
    for (int iY = 0; iY < _nY; iY++) {
      if (this->_blockData.get(iX,iY) > std::numeric_limits<T>::epsilon()) {
        this->_center[0] += _converter.getPhysLength(iX) * this->_blockData.get(iX, iY) * invNcells;
        this->_center[1] += _converter.getPhysLength(iY) * this->_blockData.get(iX, iY) * invNcells;
      }
    }
  }
  this->_latticeCenter[0] = _converter.numCells(this->_center[0]);  // TODO
  this->_latticeCenter[1] = _converter.numCells(this->_center[1]);

  // calculate moment of inertia
  this->_mofi = 0.;
  for (int iX = 0; iX < _nX; ++iX) {
    for (int iY = 0; iY < _nY; ++iY) {
      if (this->_blockData.get(iX,iY) > std::numeric_limits<T>::epsilon()) {
        T dx = std::abs(_converter.getPhysLength(iX) - this->_center[0]);
        T dy = std::abs(_converter.getPhysLength(iY) - this->_center[1]);
        this->_mofi += (dx*dx+dy*dy);
      }
    }
  }
  this->_mofi += pow(_converter.getPhysLength(1), 4)/ 6.0;
  this->_mofi *= this->_mass*invNcells;

  // calculate circumradius
  T distance = 0.;
  for (int iX = 0; iX < _nX; iX++) {
    T x = _converter.getPhysLength(iX);
    for (int iY = 0; iY < _nY; iY++) {
      T y = _converter.getPhysLength(iY);
      if (this->_blockData.get(iX,iY) > std::numeric_limits<T>::epsilon()) {
        T tmpDist = std::sqrt(std::pow(this->_center[0]-x,2)+std::pow(this->_center[1]-y,2));
        if (tmpDist > distance) {
          distance = tmpDist;
        }
      }
    }
  }
  this->_circumradius = 2.*this->_epsilon+distance;

  this->_rotMat[0] = std::cos(theta);
  this->_rotMat[1] = std::sin(theta);
  this->_rotMat[2] = -std::sin(theta);
  this->_rotMat[3] = std::cos(theta);
}

template <typename T, typename S, typename DESCRIPTOR, bool HLBM>
bool SmoothIndicatorCustom2D<T,S,DESCRIPTOR,HLBM>::operator() (T output[], const S input[])
{
  // Translation
  T xDist = input[0] - this->getPos()[0];
  T yDist = input[1] - this->getPos()[1];

  // counter-clockwise rotation by _theta=-theta around (0/0) and movement from rotation center to local center
  int x = this->_converter.numCells(xDist*this->_rotMat[0] + yDist*this->_rotMat[2]) + this->_latticeCenter[0];
  int y = this->_converter.numCells(xDist*this->_rotMat[1] + yDist*this->_rotMat[3]) + this->_latticeCenter[1];

  /// Checking if coordinates are inside the BlockData
  if (x >= 0 && x < _blockData.getNx() && y >= 0 && y < _blockData.getNy()) {
    if (this->_blockData.get(x, y) > std::numeric_limits<T>::epsilon()) {
      output[0] = T(this->_blockData.get(x, y));
      return true;
    }
  }
  output[0] = T(0);
  return false;
}
*/

//Geng2019:
template <typename T, typename S, bool HLBM>
SmoothIndicatorHTCircle2D<T,S,HLBM>::SmoothIndicatorHTCircle2D(Vector<S,2> center, S radius, S epsilon, S density, Vector<S,2> vel)
  : _radius(radius)
{
  this->_pos = center;
  this->_circumRadius = radius + 0.5*epsilon;
  this->_myMin = {center[0] - this->getCircumRadius(), center[1] - this->getCircumRadius()};
  this->_myMax = {center[0] + this->getCircumRadius(), center[1] + this->getCircumRadius()};
  this->_epsilon = epsilon;
  T mass = M_PI*radius*radius*density;
  T mofi = 0.5 * mass * radius * radius;
  this->init(0., vel, mass, mofi);
}

// returns true if x is inside the sphere
template <typename T, typename S, bool HLBM>
bool SmoothIndicatorHTCircle2D<T,S,HLBM>::operator()(T output[], const S input[])
{
  double distToCenter2 = std::pow(this->getPos()[0]-input[0], 2) +
                         std::pow(this->getPos()[1]-input[1], 2);
  
  
  double d = std::sqrt(distToCenter2) - this->_radius;
  output[0] = T((1.-tanh(d/this->getEpsilon()))/2.);
  return true;
}


} // namespace olb

#endif
