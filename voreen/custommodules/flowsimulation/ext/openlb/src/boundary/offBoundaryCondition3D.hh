/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Jonas Kratzke, Mathias J. Krause
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
 * A helper for initialising 3D boundaries -- generic implementation.
 */
#ifndef OFF_BOUNDARY_CONDITION_3D_HH
#define OFF_BOUNDARY_CONDITION_3D_HH

#include "offBoundaryCondition3D.h"
#include "offBoundaryInstantiator3D.h"
#include "offBoundaryPostProcessors3D.h"

namespace olb {

/**
* Boundary Managers provide specific Boundary Processors by creating them
*/

template<typename T, template<typename U> class Lattice, class MixinDynamics>
class BouzidiBoundaryManager3D {
public:

  static PostProcessorGenerator3D<T,Lattice>*
  getOnePointZeroVelocityBoundaryProcessor(int x, int y, int z, int iPop, T dist);
  static PostProcessorGenerator3D<T,Lattice>*
  getTwoPointZeroVelocityBoundaryProcessor(int x, int y, int z, int iPop, T dist);
  static PostProcessorGenerator3D<T,Lattice>*
  getOnePointVelocityBoundaryProcessor(int x, int y, int z, int iPop, T dist);
  static PostProcessorGenerator3D<T,Lattice>*
  getTwoPointVelocityBoundaryProcessor(int x, int y, int z, int iPop, T dist);
  static Dynamics<T,Lattice>*
  getOffDynamics(T location[Lattice<T>::d]);
  static Dynamics<T,Lattice>*
  getOffDynamics(T location[Lattice<T>::d], T distances[Lattice<T>::q]);
};

////////// BouzidiBoundaryManager3D /////////////////////////////////////////

template<typename T, template<typename U> class Lattice, class MixinDynamics>
PostProcessorGenerator3D<T,Lattice>*
BouzidiBoundaryManager3D<T,Lattice,MixinDynamics>::
getOnePointZeroVelocityBoundaryProcessor(int x, int y, int z, int iPop, T dist)
{
  return new ZeroVelocityBounceBackPostProcessorGenerator3D
         <T, Lattice>(x, y, z, iPop, dist);
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
PostProcessorGenerator3D<T,Lattice>*
BouzidiBoundaryManager3D<T,Lattice,MixinDynamics>::
getTwoPointZeroVelocityBoundaryProcessor(int x, int y, int z, int iPop, T dist)
{
  return new ZeroVelocityBouzidiLinearPostProcessorGenerator3D
         <T, Lattice>(x, y, z, iPop, dist);
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
PostProcessorGenerator3D<T,Lattice>*
BouzidiBoundaryManager3D<T,Lattice,MixinDynamics>::
getOnePointVelocityBoundaryProcessor(int x, int y, int z, int iPop, T dist)
{
  return new VelocityBounceBackPostProcessorGenerator3D
         <T, Lattice>(x, y, z, iPop, dist);
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
PostProcessorGenerator3D<T,Lattice>*
BouzidiBoundaryManager3D<T,Lattice,MixinDynamics>::
getTwoPointVelocityBoundaryProcessor(int x, int y, int z, int iPop, T dist)
{
  return new VelocityBouzidiLinearPostProcessorGenerator3D
         <T, Lattice>(x, y, z, iPop, dist);
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
Dynamics<T,Lattice>*
BouzidiBoundaryManager3D<T,Lattice,MixinDynamics>::
getOffDynamics(T location[Lattice<T>::d])
{
  return new OffDynamics<T, Lattice>(location);
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
Dynamics<T,Lattice>*
BouzidiBoundaryManager3D<T,Lattice,MixinDynamics>::
getOffDynamics(T location[Lattice<T>::d], T distances[Lattice<T>::q])
{
  return new OffDynamics<T, Lattice>(location, distances);
}

////////// Factory functions //////////////////////////////////////////////////

template<typename T, template<typename U> class Lattice, typename MixinDynamics>
OffLatticeBoundaryCondition3D<T,Lattice>*
createBouzidiBoundaryCondition3D(BlockLatticeStructure3D<T,Lattice>& block)
{
  return new OffBoundaryConditionInstantiator3D <
         T, Lattice,
         BouzidiBoundaryManager3D<T,Lattice, MixinDynamics> > (block);
}

}  // namespace olb

#endif
