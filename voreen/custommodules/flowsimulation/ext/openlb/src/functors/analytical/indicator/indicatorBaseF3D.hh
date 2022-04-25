/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014-2016 Cyril Masquelier, Albert Mink, Mathias J. Krause, Benjamin Förster
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

#ifndef INDICATOR_BASE_F_3D_HH
#define INDICATOR_BASE_F_3D_HH


#include<cmath>
#include "indicatorBaseF3D.h"
#include "utilities/vectorHelpers.h"
#include "math.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define M_PI2 1.57079632679489661923

namespace olb {

template <typename S>
IndicatorF3D<S>::IndicatorF3D() : GenericF<bool,S>(1, 3)
{}

template <typename S>
Vector<S,3>& IndicatorF3D<S>::getMin()
{
  return _myMin;
}

template <typename S>
Vector<S,3>& IndicatorF3D<S>::getMax()
{
  return _myMax;
}

template <typename S>
bool IndicatorF3D<S>::distance(S& distance, const Vector<S,3>& origin, S precision, const Vector<S,3>& direction)
{
  return util::distance(distance, origin, direction, precision,
  [&](const Vector<S,3> input) {
    return this->signedDistance(input);
  },
  [&](const Vector<S,3>& pos) {
    return isInsideBox(pos);
  });
}

template <typename S>
bool IndicatorF3D<S>::distance(S& distance, const Vector<S,3>& origin, const Vector<S,3>& direction, S precision, S pitch)
{
  return util::distance(distance, origin, direction, precision, pitch,
  [&](bool output[1], const S input[3]) {
    return (*this)(output, input);
  },
  [&](const Vector<S,3>& pos) {
    return isInsideBox(pos);
  });
}

template <typename S>
bool IndicatorF3D<S>::distance(S& distance, const Vector<S,3>& origin, const Vector<S,3>& direction, int iC)
{
  S precision = .0001;
  S pitch = 0.5;

  return this->distance(distance, origin, direction, precision, pitch);
}

template <typename S>
bool IndicatorF3D<S>::distance(S& distance, const Vector<S,3>& origin)
{
  S input[3] = {origin[0],origin[1],origin[2]};
  return this->distance(distance, input);
}

template <typename S>
bool IndicatorF3D<S>::distance(S& distance, const S input[])
{
  distance = util::abs(signedDistance(input));
  return true;
}

template <typename S>
bool IndicatorF3D<S>::rotOnAxis(Vector<S,3>& vec_rot, const Vector<S,3>& vec, const Vector<S,3>& axis, S& theta)
{
  /// http://mathworld.wolfram.com/RodriguesRotationFormula.html https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula

  //Vector<S,3> axisN(axis*(1/(axis).norm())); // normalize rotation axis
  Vector<S,3> axisN(axis*(1/norm(axis))); // normalize rotation axis

  vec_rot = vec ;

  Vector<S,3> crossProd;

  crossProd[0] = axisN[1]*vec[2] - axisN[2]*vec[1];
  crossProd[1] = axisN[2]*vec[0] - axisN[0]*vec[2];
  crossProd[2] = axisN[0]*vec[1] - axisN[1]*vec[0];

  S dotProd = axisN[0]*vec[0] + axisN[1]*vec[1] + axisN[2]*vec[2];

  //v_rot = util::cos(theta)*vec + (crossProd)*util::sin(theta) + axisN*(dotProduct3D(axisN,vec))*(1 - util::cos(theta));
  vec_rot = util::cos(theta)*vec + (crossProd)*util::sin(theta) + axisN*(dotProd)*(1 - util::cos(theta));

  return true;

}

template <typename S>
bool IndicatorF3D<S>::normal(Vector<S,3>& normal, const Vector<S,3>& origin, const Vector<S,3>& direction, int iC)
{
  //OstreamManager clout(std::cout,"IndicatorF3D");
#ifdef OLB_DEBUG
  std::cout << "Calculating IndicatorF3D Normal" << std::endl;
#endif

  bool originValue;
  (*this)(&originValue, origin.data());
  Vector<S,3> currentPoint(origin);

  S precision = .0001;

  S dist;
  distance(dist, origin, direction, iC);

  S dirMag = norm(direction);
#ifdef OLB_DEBUG
  std::cout << "magnitude = " << dirMag << std::endl;
#endif

  //Vector<S,3> POS(origin + dist*direction*(1/const_cast<Vector<S,3>&> (direction).norm())); //Point on Surface
  //direction = direction*(1/const_cast<Vector<S,3>&> (direction).norm());
  Vector<S,3> directionN(direction*(1/dirMag));
  Vector<S,3> POS(origin + dist*directionN); //Point on Surface

  /// find perpendicular vector to direction
  Vector<S,3> directionPerp;
  if (    (util::nearZero(directionN[0]) && util::nearZero(directionN[1]) && !util::nearZero(directionN[2]))
          || (util::nearZero(directionN[0]) && !util::nearZero(directionN[1]) && util::nearZero(directionN[2]))    ) {
    directionPerp[0] = 1;
    directionPerp[1] = 0;
    directionPerp[2] = 0;
  }
  else if ( !util::nearZero(directionN[0]) && util::nearZero(directionN[1]) && util::nearZero(directionN[2]) ) {
    directionPerp[0] = 0;
    directionPerp[1] = 0;
    directionPerp[2] = 1;
  }
  else if ( ( !util::nearZero(directionN[0]) || !util::nearZero(directionN[1]) ) && !util::nearZero(directionN[2]) ) {
    directionPerp[0] = directionN[0];
    directionPerp[1] = directionN[1];
    directionPerp[2] = -(directionN[0] + directionN[1])/directionN[2];
  }
  else if ( ( !util::nearZero(directionN[0]) || !util::nearZero(directionN[1]) ) && util::nearZero(directionN[2]) ) {
    directionPerp[0] = directionN[0];
    directionPerp[1] = -(directionN[0] + directionN[2])/directionN[1];
    directionPerp[2] = directionN[2];
  }
  else {
    std::cout << "Error: unknown case for perpendicular check" << std::endl;
    return false;
  }

  Vector<S,3> directionPerpN(directionPerp*(dirMag/norm(directionPerp)));


  Vector<S,3> point1;
  Vector<S,3> point2;
  Vector<S,3> point3;

  bool currentValue;

  /// Loop 3 times to find three points on the surface to use for normal calc.
  /// orthogonal to direction vector 120 degree to each other


  for (int n: {
         0,120,240
       }) {
    S thetaMain = n*M_PI/180.;

    /// rotate directionPerpN through 3 angles {0,120,240}
    Vector<S,3> perp;
    rotOnAxis(perp, directionPerpN, directionN, thetaMain);
    Vector<S,3> perpPoint(POS + perp);

    //std::cout << "perp = [" << perp[0] << "," << perp[1]  << "," << perp[2] << "]" << std::endl;

    //S rotate = 90.;
    S rotate = 179.;
    S pitch = rotate/2.;

    Vector<S,3> vec(perp);

    Vector<S,3> rotAxis;

    //rotAxis = cross(perp,directionN);
    rotAxis[0] = perp[1]*directionN[2] - perp[2]*directionN[1];
    rotAxis[1] = perp[2]*directionN[0] - perp[0]*directionN[2];
    rotAxis[2] = perp[0]*directionN[1] - perp[1]*directionN[0];

    //normal = rotAxis;
    //return true;
    //std::cout << "rotAxis = [" << rotAxis[0] << "," << rotAxis[1]  << "," << rotAxis[2] << "]" << std::endl;

    /// Find 'positive' angle
    Vector<S,3> testPOS;
    S testAngle(45.*M_PI/180.);
    rotOnAxis(testPOS, perp, rotAxis, testAngle);
    Vector<S,3> testPoint( POS + testPOS);
    //std::cout << "testPOS = [" << testPOS[0] << "," << testPOS[1]  << "," << testPOS[2] << "]" << std::endl;
    //std::cout << "testPoint = [" << testPoint[0] << "," << testPoint[1]  << "," << testPoint[2] << "]" << std::endl;

    S distTestPoint = norm(testPoint - origin);
    S distPerpPoint = norm(perpPoint - origin);

    S mod = 0;
    if (distTestPoint < distPerpPoint) { // pos. angle rotates towards
      mod = -1;
    }
    else {
      mod = 1;
    }

    while (util::abs(pitch) >= precision) {

      S theta(pitch*M_PI/180);

      currentPoint = POS + vec;
      (*this)(&currentValue, currentPoint.data());

      S temp;
      if (currentValue == originValue) {
        temp = mod*theta;
        rotOnAxis(vec, vec, rotAxis, temp);
      }
      else {
        temp = -mod*theta;
        rotOnAxis(vec, vec, rotAxis, temp);
      }
      pitch /= 2.;



    }

    if (n == 0) {
      point1 = currentPoint;
    }
    else if (n == 120) {
      point2 = currentPoint;
    }
    else if (n == 240) {
      point3 = currentPoint;
    }
    else {
      std::cout << "Something broke" << std::endl;
      return false;
    }

  }

  /// Calculate Normal
  Vector<S,3> vec1 (point1 - point2);
  Vector<S,3> vec2 (point1 - point3);

  normal[0] = -(vec1[1]*vec2[2] - vec1[2]*vec2[1]);
  normal[1] = -(vec1[2]*vec2[0] - vec1[0]*vec2[2]);
  normal[2] = -(vec1[0]*vec2[1] - vec1[1]*vec2[0]);

  normalize(normal);


  //S dist;
  //Vector<S,3> dist;
  //distance(dist, origin, direction, iC);
  //normal = Vector<S,3>(dist);
  //normal = POS;
  //normal = directionPerpN;
  //normal = directionN;
  return true;

}

template <typename S>
bool IndicatorF3D<S>::isInsideBox(Vector<S,3> point)
{
  return point >= _myMin && point <= _myMax;
}

template <typename S>
S IndicatorF3D<S>::signedDistance(const Vector<S,3>& input)
{
  // TODO: Add fallback algorithm here
  assert(false);
  return std::numeric_limits<double>::quiet_NaN();
}

template <typename S>
Vector<S,3> IndicatorF3D<S>::surfaceNormal(const Vector<S,3>& pos, const S meshSize)
{
  return util::surfaceNormal(pos, meshSize, [&](const Vector<S,3>& pos) {
    return this->signedDistance(pos);
  });
}

template <typename S>
Vector<S,3> IndicatorF3D<S>::surfaceNormal(const Vector<S,3>& pos, const S meshSize,
    std::function<Vector<S,3>(const Vector<S,3>&)> transformPos)
{
  return this->surfaceNormal(transformPos(pos), meshSize);
}

template <typename S>
bool IndicatorF3D<S>::operator() (bool output[1], const S input[3])
{
  output[0] = this->signedDistance(input) <= 0;
  return output[0];
}

template <typename S>
IndicatorIdentity3D<S>::IndicatorIdentity3D(std::shared_ptr<IndicatorF3D<S>> f)
  : _f(f)
{
  this->_myMin = _f->getMin();
  this->_myMax = _f->getMax();
}

template <typename S>
bool IndicatorIdentity3D<S>::operator() (bool output[1], const S input[3])
{
  return (_f)->operator()(output, input);
}


} // namespace olb

#endif
