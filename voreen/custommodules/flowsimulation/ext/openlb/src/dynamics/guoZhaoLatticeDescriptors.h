/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2016-2017 Davide Dapelo, Mathias J. Krause
 *  OpenLB e-mail contact: info@openlb.net
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
 * Groups all the 2D include files in the dynamics directory.
*/

/** \file
 * Descriptor for all types of 2D and 3D lattices for the Guo-Zhao
 * porous model. In principle, thanks
 * to the fact that the OpenLB code is generic, it is sufficient to
 * write a new descriptor when a new type of lattice is to be used.
 *  -- header file
 */

#ifndef GUOZHAO_LATTICE_DESCRIPTORS_H
#define GUOZHAO_LATTICE_DESCRIPTORS_H

#include <vector>
#include "core/olbDebug.h"

namespace olb {

/// Descriptors for the 2D and 3D lattices.
/** \warning Attention: The lattice directions must always be ordered in
 * such a way that c[i] = -c[i+(q-1)/2] for i=1..(q-1)/2, and c[0] = 0 must
 * be the rest velocity. Furthermore, the velocities c[i] for i=1..(q-1)/2
 * must verify
 *  - in 2D: (c[i][0]<0) || (c[i][0]==0 && c[i][1]<0)
 *  - in 3D: (c[i][0]<0) || (c[i][0]==0 && c[i][1]<0)
 *                       || (c[i][0]==0 && c[i][1]==0 && c[i][2]<0)
 * Otherwise some of the code will work erroneously, because the
 * aformentioned relations are taken as given to enable a few
 * optimizations.
*/
namespace descriptors {

struct GuoZhao2dDescriptor {
  static const int numScalars = 7;
  static const int numSpecies = 5;
  static const int forceBeginsAt = 0;
  static const int epsilonAt = 2;
  static const int KAt = 3;
  static const int nuAt = 4;
  static const int bodyForceBeginsAt = 5;
  static const int sizeOfForce   = 2;
};

struct GuoZhao2dDescriptorBase {
  typedef GuoZhao2dDescriptor ExternalField;
};

template <typename T> struct GuoZhaoD2Q9Descriptor
    : public D2Q9DescriptorBase<T>, public GuoZhao2dDescriptorBase {
};

}  // namespace descriptors

}  // namespace olb

#endif
