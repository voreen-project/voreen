/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012, 2016 Jonas Kratzke, Mathias J. Krause
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
 * A helper for initialising 2D boundaries -- generic implementation.
 */
#ifndef OFF_BOUNDARY_CONDITION_2D_HH
#define OFF_BOUNDARY_CONDITION_2D_HH

#include "offBoundaryCondition2D.h"
#include "offBoundaryInstantiator2D.h"
#include "offBoundaryPostProcessors2D.h"

namespace olb {

/**
* Boundary Managers provide specific Boundary Processors by creating them
*/

////////// BouzidiBoundaryManager2D /////////////////////////////////////////

template<typename T, template<typename U> class Lattice, class MixinDynamics>
class BouzidiBoundaryManager2D {
public:

  static PostProcessorGenerator2D<T,Lattice>*
  getOnePointZeroVelocityBoundaryProcessor(int iX, int iY, int iPop, T dist);
  static PostProcessorGenerator2D<T,Lattice>*
  getTwoPointZeroVelocityBoundaryProcessor(int iX, int iY, int iPop, T dist);
  static PostProcessorGenerator2D<T,Lattice>*
  getOnePointVelocityBoundaryProcessor(int iX, int iY, int iPop, T dist);
  static PostProcessorGenerator2D<T,Lattice>*
  getTwoPointVelocityBoundaryProcessor(int iX, int iY, int iPop, T dist);
  static Dynamics<T,Lattice>*
  getOffDynamics(T location[Lattice<T>::d]);
  static Dynamics<T,Lattice>*
  getOffDynamics(T location[Lattice<T>::d], T distances[Lattice<T>::q]);
};


template<typename T, template<typename U> class Lattice, class MixinDynamics>
PostProcessorGenerator2D<T,Lattice>*
BouzidiBoundaryManager2D<T,Lattice,MixinDynamics>::
getOnePointZeroVelocityBoundaryProcessor(int iX, int iY, int iPop, T dist)
{
  return new ZeroVelocityBounceBackPostProcessorGenerator2D
         <T, Lattice>(iX, iY, iPop, dist);
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
PostProcessorGenerator2D<T,Lattice>*
BouzidiBoundaryManager2D<T,Lattice,MixinDynamics>::
getTwoPointZeroVelocityBoundaryProcessor(int iX, int iY, int iPop, T dist)
{
  return new ZeroVelocityBouzidiLinearPostProcessorGenerator2D
         <T, Lattice>(iX, iY, iPop, dist);
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
PostProcessorGenerator2D<T,Lattice>*
BouzidiBoundaryManager2D<T,Lattice,MixinDynamics>::
getOnePointVelocityBoundaryProcessor(int iX, int iY, int iPop, T dist)
{
  return new VelocityBounceBackPostProcessorGenerator2D
         <T, Lattice>(iX, iY, iPop, dist);
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
PostProcessorGenerator2D<T,Lattice>*
BouzidiBoundaryManager2D<T,Lattice,MixinDynamics>::
getTwoPointVelocityBoundaryProcessor(int iX, int iY, int iPop, T dist)
{
  return new VelocityBouzidiLinearPostProcessorGenerator2D
         <T, Lattice>(iX, iY, iPop, dist);
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
Dynamics<T,Lattice>*
BouzidiBoundaryManager2D<T,Lattice,MixinDynamics>::
getOffDynamics(T location[Lattice<T>::d])
{
  return new OffDynamics<T, Lattice>(location);
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
Dynamics<T,Lattice>*
BouzidiBoundaryManager2D<T,Lattice,MixinDynamics>::
getOffDynamics(T location[Lattice<T>::d], T distances[Lattice<T>::q])
{
  return new OffDynamics<T, Lattice>(location, distances);
}

////////// BounceBackBoundaryManager2D /////////////////////////////////////////

template<typename T, template<typename U> class Lattice>
class BounceBackBoundaryManager2D {
public:

  static PostProcessorGenerator2D<T,Lattice>*
  getOnePointZeroVelocityBoundaryProcessor(int iX, int iY, int iPop, T dist);
  static PostProcessorGenerator2D<T,Lattice>*
  getTwoPointZeroVelocityBoundaryProcessor(int iX, int iY, int iPop, T dist);
  static PostProcessorGenerator2D<T,Lattice>*
  getOnePointVelocityBoundaryProcessor(int iX, int iY, int iPop, T dist);
  static PostProcessorGenerator2D<T,Lattice>*
  getTwoPointVelocityBoundaryProcessor(int iX, int iY, int iPop, T dist);
  static Dynamics<T,Lattice>*
  getOffDynamics(T location[Lattice<T>::d]);
  static Dynamics<T,Lattice>*
  getOffDynamics(T location[Lattice<T>::d], T distances[Lattice<T>::q]);
};

template<typename T, template<typename U> class Lattice>
PostProcessorGenerator2D<T,Lattice>*
BounceBackBoundaryManager2D<T,Lattice>::
getOnePointZeroVelocityBoundaryProcessor(int iX, int iY, int iPop, T dist)
{
  return new ZeroVelocityBounceBackPostProcessorGenerator2D
         <T, Lattice>(iX, iY, iPop, dist);
}

template<typename T, template<typename U> class Lattice>
PostProcessorGenerator2D<T,Lattice>*
BounceBackBoundaryManager2D<T,Lattice>::
getTwoPointZeroVelocityBoundaryProcessor(int iX, int iY, int iPop, T dist)
{
  return new ZeroVelocityBouzidiLinearPostProcessorGenerator2D
         <T, Lattice>(iX, iY, iPop, dist);
}

template<typename T, template<typename U> class Lattice>
PostProcessorGenerator2D<T,Lattice>*
BounceBackBoundaryManager2D<T,Lattice>::
getOnePointVelocityBoundaryProcessor(int iX, int iY, int iPop, T dist)
{
  return new VelocityBounceBackPostProcessorGenerator2D
         <T, Lattice>(iX, iY, iPop, dist);
}

template<typename T, template<typename U> class Lattice>
PostProcessorGenerator2D<T,Lattice>*
BounceBackBoundaryManager2D<T,Lattice>::
getTwoPointVelocityBoundaryProcessor(int iX, int iY, int iPop, T dist)
{
  return new VelocityBouzidiLinearPostProcessorGenerator2D
         <T, Lattice>(iX, iY, iPop, dist);
}

template<typename T, template<typename U> class Lattice>
Dynamics<T,Lattice>*
BounceBackBoundaryManager2D<T,Lattice>::
getOffDynamics(T location[Lattice<T>::d])
{
  return new OffDynamics<T, Lattice>(location);
}

template<typename T, template<typename U> class Lattice>
Dynamics<T,Lattice>*
BounceBackBoundaryManager2D<T,Lattice>::
getOffDynamics(T location[Lattice<T>::d], T distances[Lattice<T>::q])
{
  return new OffDynamics<T, Lattice>(location, distances);
}

////////// Factory functions //////////////////////////////////////////////////

template<typename T, template<typename U> class Lattice, typename MixinDynamics>
OffLatticeBoundaryCondition2D<T,Lattice>*
createBouzidiBoundaryCondition2D(BlockLatticeStructure2D<T,Lattice>& block)
{
  return new OffBoundaryConditionInstantiator2D <
         T, Lattice,
         BouzidiBoundaryManager2D<T,Lattice, MixinDynamics> > (block);
}


}  // namespace olb

#endif
