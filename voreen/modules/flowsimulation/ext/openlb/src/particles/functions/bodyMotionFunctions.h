/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Jan E. Marquardt, Mathias J. Krause
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


/* This file contains functions used for the calculation of body motions.
 *
*/

#ifndef BODY_MOTION_FUNCTIONS_H
#define BODY_MOTION_FUNCTIONS_H


#include <cassert>

namespace olb {

namespace body {

namespace motion {


////////////// Dimension sensitive Functions ////////////


template<unsigned D, typename T, bool OUTPUT_USES_ROTATION_CENTER_AS_ORIGIN=false>
struct rotation;

template<typename T, bool OUTPUT_USES_ROTATION_CENTER_AS_ORIGIN>
struct rotation<2,T,OUTPUT_USES_ROTATION_CENTER_AS_ORIGIN> {
  static constexpr Vector<T,2> execute( Vector<T,2> input, Vector<T,4> rotationMatrix, Vector<T,2> rotationCenter = Vector<T,2>(0.,0.) )
  {
    const Vector<T,2> dist = input - rotationCenter;
    const Vector<T,2> rotated = Vector<T,2>( dist[0]*rotationMatrix[0] + dist[1]*rotationMatrix[2],
                                dist[0]*rotationMatrix[1] + dist[1]*rotationMatrix[3] );
    if constexpr(!OUTPUT_USES_ROTATION_CENTER_AS_ORIGIN) {
      return rotationCenter + rotated;
    }
    else {
      return rotated;
    }
  }

  static constexpr Vector<T,2> invert( Vector<T,2> input, Vector<T,4> rotationMatrix, Vector<T,2> rotationCenter = Vector<T,2>(0.) )
  {
    std::cout << "Invert needs implementation." << std::endl;
    assert(false);
    return input;
  }
};

template<typename T, bool OUTPUT_USES_ROTATION_CENTER_AS_ORIGIN>
struct rotation<3,T,OUTPUT_USES_ROTATION_CENTER_AS_ORIGIN> {
  static constexpr Vector<T,3> execute( Vector<T,3> input, Vector<T,9> rotationMatrix, Vector<T,3> rotationCenter = Vector<T,3>(0.,0.,0.) )
  {
    const Vector<T,3> dist = input - rotationCenter;
    const Vector<T,3> rotated = Vector<T,3>( rotationMatrix[0]*dist[0] + rotationMatrix[3]*dist[1] + rotationMatrix[6]*dist[2],
                                rotationMatrix[1]*dist[0] + rotationMatrix[4]*dist[1] + rotationMatrix[7]*dist[2],
                                rotationMatrix[2]*dist[0] + rotationMatrix[5]*dist[1] + rotationMatrix[8]*dist[2] );
    if constexpr(!OUTPUT_USES_ROTATION_CENTER_AS_ORIGIN) {
      return rotationCenter + rotated;
    }
    else {
      return rotated;
    }
  }

  static constexpr Vector<T,3> invert( Vector<T,3> input, Vector<T,9> rotationMatrix, Vector<T,3> rotationCenter = Vector<T,3>(0.) )
  {
    std::cout << "Invert needs implementation." << std::endl;
    assert(false);
    return input;
  }
};

} //namespace motion

} //namespace body

} //namespace olb


#endif
