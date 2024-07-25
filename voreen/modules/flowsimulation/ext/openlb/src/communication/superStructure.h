/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014 Peter Weisbrod, Albert Mink, Mathias J. Krause
 *                2021 Adrian Kummerlaender
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

#ifndef SUPER_STRUCTURE_H
#define SUPER_STRUCTURE_H

#include "utilities/aliases.h"

#include "geometry/cuboidGeometry2D.h"
#include "geometry/cuboidGeometry3D.h"

#include "communication/loadBalancer.h"

namespace olb {


template<typename T, unsigned D>
class SuperStructure {
protected:
  /// The grid structure is stored here
  CuboidGeometry<T,D>& _cuboidGeometry;
  /// Distribution of the cuboids of the cuboid structure
  LoadBalancer<T>& _loadBalancer;
  /// Size of ghost cell layer (must be greater than 1 and
  /// greater_overlapBC, default =1)
  int _overlap;
  /// class specific output stream
  mutable OstreamManager clout;

public:
  using value_t = T;

  /// Virtual Destructor for inheritance
  virtual ~SuperStructure() {};
  /// Construction of a super structure
  SuperStructure(CuboidGeometry<T,D>& cuboidGeometry,
                 LoadBalancer<T>& loadBalancer, int overlap = 2);
  /// Default Constructor for empty SuperStructure
  SuperStructure(int overlap = 1);

  /// Read and write access to cuboid geometry
  CuboidGeometry<T,D>& getCuboidGeometry();
  /// Read only access to cuboid geometry
  CuboidGeometry<T,D> const& getCuboidGeometry() const;

  /// Read and write access to the overlap
  int getOverlap();
  /// Read only access to the overlap
  int getOverlap() const;

  /// Read and write access to the load balancer
  LoadBalancer<T>& getLoadBalancer();
  /// Read only access to the load balancer
  LoadBalancer<T> const& getLoadBalancer() const;

  virtual void communicate() { };

};


}

#endif
