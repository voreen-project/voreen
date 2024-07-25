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
 *  but WITHOUS ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
*/


#ifndef INDICATOR_BASE_H
#define INDICATOR_BASE_H

#include "core/util.h"
#include "utilities/vectorHelpers.h"

namespace olb {

namespace util {

template <typename S, unsigned D, typename F1, typename F2>
bool distance(S& distance, const Vector<S,D>& origin, const Vector<S,D>& direction,
              S precision, S pitch, F1 isInside, F2 isInsideBoundingBox)
{
  // Check if point is on surface
  bool isInsideDistance;
  bool isInsideOppDistance;
  Vector<S,D> currentPoint = origin + 10*precision*direction;
  isInside(&isInsideDistance, currentPoint.data());
  currentPoint = origin - 10*precision*direction;
  isInside(&isInsideOppDistance, currentPoint.data());
  if (isInsideDistance == !isInsideOppDistance) {
    distance = 0;
    return true;
  }

  bool originValue;
  bool currentValue;

  isInside(&originValue, origin.data());

  // start at origin and move into given direction
  currentPoint = origin;
  currentValue = originValue;

  while (currentValue == originValue && isInsideBoundingBox(currentPoint)) {
    currentPoint += direction;
    // update currentValue until the first point on the other side (inside/outside) is found
    isInside(&currentValue, currentPoint.data());
  }

  // return false if no point was found in given direction
  if (!isInsideBoundingBox(currentPoint) && !originValue) {
    return false;
  }


  while (pitch >= precision) {
    if (!isInsideBoundingBox(currentPoint) && originValue) {
      currentPoint -= pitch * direction;
      pitch *= 0.5;
    }
    else {
      isInside(&currentValue, currentPoint.data());
      if (currentValue == originValue) {
        currentPoint += pitch * direction;
        pitch *= 0.5;
      }
      else {
        currentPoint -= pitch * direction;
        pitch *= 0.5;
      }
    }
  }

  distance = norm(currentPoint - origin);
  return true;
}

template <typename S, unsigned D, typename F1, typename F2>
bool distance(S& distance, const Vector<S,D>& origin, const Vector<S,D>& direction,
              S precision, F1 sdf, F2 isInsideBoundingBox)
{
  distance = S(0);
  S signedDistance = sdf(origin);
  short firstSign = util::sign(signedDistance);

  Vector<S,D> currPos = origin;
  Vector<S,D> dir = normalize(direction);

  while (util::fabs(signedDistance) > precision) {
    currPos += (signedDistance * firstSign) * dir;
    if (!isInsideBoundingBox(currPos)) {
      distance = 0;
      return false;
    }
    distance += signedDistance;
    signedDistance = sdf(currPos);
  }

  distance = util::fabs(distance);
  return true;
}

/** Using a bisect to find the unsigned distance (false if distance was not found, true if distance was found)
 *
 * \param distance              computed distance from origin to surface in direction \param direction is saved here
 * \param origin                the point we want to calculate the distance from
 * \param direction             the considered direction
 * \param pitch                 a first guess for the distance, must be > 0 (if origin is inside origin + normalize(direction) * pitch should be outside and if origin is outside origin + normalize(direction) * pitch should be inside)
 * \param precision             when the change of the found distance reaches the defined precision, the iteration stops
 * \param isInside              function that defines if a point is inside or outside a geometry (usually the ()-operator of an Indicator)
 * \param isInsideBoundingBox   function that defines if a point is inside the enclosing bounding box (maybe needs at least one extra layer between geometry surface and bounding box surface)
*/
template <typename S, unsigned D, typename F1, typename F2>
bool bisectDistance(S& distance, const Vector<S,D>& origin, const Vector<S,D>& direction,
                    S pitch, S precision, F1 isInside, F2 isInsideBoundingBox)
{
  distance = S(0);
  Vector<S,D> dir = normalize(direction);

  // Check if origin is on surface
  bool isInsideDistance = isInside(origin + 5*precision*dir);
  bool isInsideOppDistance = isInside(origin - 5*precision*dir);
  if (isInsideDistance == !isInsideOppDistance) {
    return true;
  }

  Vector<S,D> startPoint = origin;
  S fixedDistance = S(0);
  S oldDistance = distance;
  std::function<Vector<S,D>()> calcCurrPoint = [&startPoint, &pitch, &dir]() {
    return startPoint + pitch * dir;
  };
  Vector<S,D> currPoint = calcCurrPoint();

  bool originValue = isInside(startPoint);
  bool currentValue = isInside(currPoint);

  // if aforementioned requirement for pitch is not given, try to find a proper pitch
  while (originValue == currentValue) {
    pitch *= 1.5;
    currPoint = calcCurrPoint();
    currentValue = isInside(currPoint);

    // if both points are outside the bounding box return false indicating that no distance was found
    if (!currentValue && !isInsideBoundingBox(startPoint) && !isInsideBoundingBox(currPoint)) {
      return false;
    }
  }

  for (unsigned iD=0; iD<D; ++iD) {
    if (!std::isfinite(currPoint[iD])) {
      return false;
    }
  }

  // if we have one point inside and another outside, we can keep halving the line from origin to currPoint until we find the wanted precision
  distance = norm(startPoint - currPoint);
  while ( util::abs(distance - oldDistance) > precision ) {
    pitch *= 0.5;
    currPoint = calcCurrPoint();

    // set distances
    oldDistance = distance;
    distance = norm(startPoint - currPoint);

    // if both points are outside or inside run take result of same function with different start point (= current point)
    currentValue = isInside(currPoint);
    if (currentValue == originValue) {
      startPoint = currPoint;
      fixedDistance += distance;
      // convert distances to distance to new start point
      oldDistance -= distance;
      distance = S(0);
    }
  }
  distance += fixedDistance;

  return true;
}

template <typename S, unsigned D, typename F1>
Vector<S,D> surfaceNormal(const Vector<S,D>& pos, const S meshSize, F1 sdf)
{
  Vector<S,D> normal( S(0) );

  for (unsigned iD=0; iD<D; ++iD) {
    Vector<S,D> delta(S(0));
    delta[iD] = meshSize;
    normal[iD] = (sdf(pos+delta) - sdf(pos-delta)) / (2*meshSize);
  }

  return normal;
}

} //namespace util

} //namespace olb


#endif
