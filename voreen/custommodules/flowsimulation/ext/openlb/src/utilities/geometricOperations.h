/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Julius Jessberger
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

/** \file
 * Some arithmetic helper functions.
 */

#ifndef GEOMETRIC_OPERATIONS_H
#define GEOMETRIC_OPERATIONS_H

#include "core/vector.h"
#include "dimensionConverter.h"

namespace olb {

namespace util {

//Calculate local velocity
template <typename T, unsigned D>
Vector<T,D> calculateLocalVelocity( const Vector<T,D>& rotationCenter,
                                    const Vector<T,D>& velocity,
                                    const Vector<T,utilities::dimensions::convert<D>::rotation>& angularVelocity,
                                    const Vector<T,D>& position)
{
  if constexpr(D == 2) {
    // two dimensions: u = U + w x r = (Ux, Uy, 0) + (0,0,w) x (X,Y,0) = (Ux, Uy, 0) + (-w*Y, w*X, 0)
    return Vector<T,2>( velocity[0] - angularVelocity[0]*(position[1] - rotationCenter[1]),
                        velocity[1] + angularVelocity[0]*(position[0] - rotationCenter[0])
                      );
  }
  else {
    // three dimensions: u = U + w x r = (Ux, Uy, Uz) + (wx,wy,wz) x (X,Y,Z) = (Ux, Uy, Uz) + (wy*Z-wz*Y, wz*X-wx*Z, wx*Y-wy*X)
    return Vector<T,3>( velocity[0] + angularVelocity[1]*(position[2] - rotationCenter[2]) - angularVelocity[2]*(position[1] - rotationCenter[1]),
                        velocity[1] + angularVelocity[2]*(position[0] - rotationCenter[0]) - angularVelocity[0]*(position[2] - rotationCenter[2]),
                        velocity[2] + angularVelocity[0]*(position[1] - rotationCenter[1]) - angularVelocity[1]*(position[0] - rotationCenter[0])
                      );
  }
  __builtin_unreachable();
}

}

}

#endif
