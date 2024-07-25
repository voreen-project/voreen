/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Nicolas Hafen, Mathias J. Krause
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

#ifndef DIMENSION_CONVERTER_H
#define DIMENSION_CONVERTER_H

#include "functors/analytical/indicator/smoothIndicatorBaseF2D.h"
#include "functors/analytical/indicator/smoothIndicatorBaseF3D.h"

namespace olb {

namespace utilities {

namespace dimensions {

/// Converts dimensions by deriving from given cartesian dimension D
// - differentiates between F2D and F3D indicators
// - differentiates between D4 and D9 Rotation Matrix
// - differentiates between D1 and D3 rotation related quantities
template <unsigned D>
struct convert;
template <>
struct convert<2> {
  // Rotational dimensions (rotational degrees of freedom)
  constexpr static const unsigned int rotation = 1;

  // Number of entries in rotation matrix
  constexpr static const unsigned int matrix = 4;

  // Number of direct neighbors of a cell
  constexpr static const unsigned short directNeighborsCount = 4;
  // Directions to direct neighbors
  constexpr static const short directNeighborDirections[4][2] = {
    {-1, 0}, {0,-1},
    {1, 0}, { 0, 1}
  };
  // Number of all neighboring cells
  constexpr static const unsigned short neighborsCount = 8;
  // Directions to neighbors
  constexpr static const short neighborDirections[8][2] = {
    {-1, 1}, {-1, 0}, {-1,-1}, { 0,-1},
    { 1,-1}, { 1, 0}, { 1, 1}, { 0, 1}
  };

  // Type used for surface representation
  template<typename T>
  using surfaceType = SmoothIndicatorF2D<T,T,true>;

  // Converting rotational vector to serial type for field access
  template<typename T>
  static constexpr T serialize_rotation( Vector<T,rotation> angle  )
  {
    return angle[0];
  }
};
template <>
struct convert<3> {
  // Rotational dimensions (rotational degrees of freedom)
  constexpr static const unsigned int rotation = 3;

  // Number of entries in rotation matrix
  constexpr static const unsigned int matrix = 9;

  // Number of direct neighbors of a cell
  constexpr static const unsigned short directNeighborsCount = 6;
  // Directions to direct neighbors
  constexpr static const short directNeighborDirections[6][3] = {
    {-1, 0, 0}, {0,-1, 0},
    { 0, 0,-1}, {1, 0, 0},
    { 0, 1, 0}, {0, 0, 1}
  };
  // Number of all neighboring cells
  constexpr static const unsigned short neighborsCount = 26;
  // Directions to neighbors
  constexpr static const short neighborDirections[26][3] = {
    {-1, 0, 0}, { 0,-1, 0}, { 0, 0,-1},
    {-1,-1, 0}, {-1, 1, 0}, {-1, 0,-1},
    {-1, 0, 1}, { 0,-1,-1}, { 0,-1, 1},
    {-1,-1,-1}, {-1,-1, 1}, {-1, 1,-1}, {-1, 1, 1},
    { 1, 0, 0}, { 0, 1, 0}, { 0, 0, 1},
    { 1, 1, 0}, { 1,-1, 0}, { 1, 0, 1},
    { 1, 0,-1}, { 0, 1, 1}, { 0, 1,-1},
    { 1, 1, 1}, { 1, 1,-1}, { 1,-1, 1}, { 1,-1,-1}
  };

  // Type used for surface representation
  template<typename T>
  using surfaceType = SmoothIndicatorF3D<T,T,true>;

  // Converting rotational vector to serial type for field access
  template<typename T>
  static constexpr Vector<T,rotation> serialize_rotation( Vector<T,rotation> angle  )
  {
    return angle;
  }
};


} //namespace dimensions

} //namespace utilities

} //namespace olb


#endif
