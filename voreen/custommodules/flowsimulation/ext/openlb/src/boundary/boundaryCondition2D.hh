/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Jonas Latt
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
#ifndef BOUNDARY_CONDITION_2D_HH
#define BOUNDARY_CONDITION_2D_HH

#include "boundaryCondition2D.h"
#include "boundaryInstantiator2D.h"


namespace olb {

template<typename T, template<typename U> class Lattice, class MixinDynamics>
class RegularizedBoundaryManager2D {
public:
  template<int direction, int orientation> static Momenta<T,Lattice>*
  getVelocityBoundaryMomenta();
  template<int direction, int orientation> static Dynamics<T,Lattice>*
  getVelocityBoundaryDynamics(T omega, Momenta<T,Lattice>& momenta);
  template<int direction, int orientation> static PostProcessorGenerator2D<T,Lattice>*
  getVelocityBoundaryProcessor(int x0, int x1, int y0, int y1);

  template<int direction, int orientation> static Momenta<T,Lattice>*
  getPressureBoundaryMomenta();
  template<int direction, int orientation> static Dynamics<T,Lattice>*
  getPressureBoundaryDynamics(T omega, Momenta<T,Lattice>& momenta);
  template<int direction, int orientation> static PostProcessorGenerator2D<T,Lattice>*
  getPressureBoundaryProcessor(int x0, int x1, int y0, int y1);

  template<int direction, int orientation> static PostProcessorGenerator2D<T,Lattice>*
  getConvectionBoundaryProcessor(int x0, int x1, int y0, int y1, T* uAv=NULL);

  template<int xNormal, int yNormal> static Momenta<T,Lattice>*
  getExternalVelocityCornerMomenta();
  template<int xNormal, int yNormal> static Dynamics<T,Lattice>*
  getExternalVelocityCornerDynamics(T omega, Momenta<T,Lattice>& momenta);
  template<int xNormal, int yNormal> static PostProcessorGenerator2D<T,Lattice>*
  getExternalVelocityCornerProcessor(int x, int y);

  template<int xNormal, int yNormal> static Momenta<T,Lattice>*
  getInternalVelocityCornerMomenta();
  template<int xNormal, int yNormal> static Dynamics<T,Lattice>*
  getInternalVelocityCornerDynamics(T omega, Momenta<T,Lattice>& momenta);
  template<int xNormal, int yNormal> static PostProcessorGenerator2D<T,Lattice>*
  getInternalVelocityCornerProcessor(int x, int y);
};

template<typename T, template<typename U> class Lattice, class MixinDynamics>
class InterpolationBoundaryManager2D {
public:
  template<int direction, int orientation> static Momenta<T,Lattice>*
  getVelocityBoundaryMomenta();
  template<int direction, int orientation> static Dynamics<T,Lattice>*
  getVelocityBoundaryDynamics(T omega, Momenta<T,Lattice>& momenta);
  template<int direction, int orientation> static PostProcessorGenerator2D<T,Lattice>*
  getVelocityBoundaryProcessor(int x0, int x1, int y0, int y1);

  template<int direction, int orientation> static Momenta<T,Lattice>*
  getPressureBoundaryMomenta();
  template<int direction, int orientation> static Dynamics<T,Lattice>*
  getPressureBoundaryDynamics(T omega, Momenta<T,Lattice>& momenta);
  template<int direction, int orientation> static PostProcessorGenerator2D<T,Lattice>*
  getPressureBoundaryProcessor(int x0, int x1, int y0, int y1);

  template<int direction, int orientation> static PostProcessorGenerator2D<T,Lattice>*
  getConvectionBoundaryProcessor(int x0, int x1, int y0, int y1, T* uAv=NULL);

  template<int xNormal, int yNormal> static Momenta<T,Lattice>*
  getExternalVelocityCornerMomenta();
  template<int xNormal, int yNormal> static Dynamics<T,Lattice>*
  getExternalVelocityCornerDynamics(T omega, Momenta<T,Lattice>& momenta);
  template<int xNormal, int yNormal> static PostProcessorGenerator2D<T,Lattice>*
  getExternalVelocityCornerProcessor(int x, int y);

  template<int xNormal, int yNormal> static Momenta<T,Lattice>*
  getInternalVelocityCornerMomenta();
  template<int xNormal, int yNormal> static Dynamics<T,Lattice>*
  getInternalVelocityCornerDynamics(T omega, Momenta<T,Lattice>& momenta);
  template<int xNormal, int yNormal> static PostProcessorGenerator2D<T,Lattice>*
  getInternalVelocityCornerProcessor(int x, int y);
};


////////// RegularizedBoundaryManager2D /////////////////////////////////////////

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int direction, int orientation>
Momenta<T,Lattice>*
RegularizedBoundaryManager2D<T,Lattice,MixinDynamics>::getVelocityBoundaryMomenta()
{
  return new RegularizedVelocityBM<T,Lattice, direction,orientation>;
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int direction, int orientation>
Dynamics<T,Lattice>* RegularizedBoundaryManager2D<T,Lattice,MixinDynamics>::
getVelocityBoundaryDynamics(T omega, Momenta<T,Lattice>& momenta)
{
  return new CombinedRLBdynamics<T,Lattice, MixinDynamics>(omega, momenta);
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int direction, int orientation>
PostProcessorGenerator2D<T,Lattice>*
RegularizedBoundaryManager2D<T,Lattice,MixinDynamics>::
getVelocityBoundaryProcessor(int x0, int x1, int y0, int y1)
{
  return nullptr;
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int direction, int orientation>
Momenta<T,Lattice>*
RegularizedBoundaryManager2D<T,Lattice,MixinDynamics>::getPressureBoundaryMomenta()
{
  return new RegularizedPressureBM<T,Lattice, direction,orientation>;
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int direction, int orientation>
Dynamics<T,Lattice>* RegularizedBoundaryManager2D<T,Lattice,MixinDynamics>::
getPressureBoundaryDynamics(T omega, Momenta<T,Lattice>& momenta)
{
  return new CombinedRLBdynamics<T,Lattice, MixinDynamics>(omega, momenta);
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int direction, int orientation>
PostProcessorGenerator2D<T,Lattice>*
RegularizedBoundaryManager2D<T,Lattice,MixinDynamics>::
getPressureBoundaryProcessor(int x0, int x1, int y0, int y1)
{
  return nullptr;
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int direction, int orientation>
PostProcessorGenerator2D<T,Lattice>*
RegularizedBoundaryManager2D<T,Lattice,MixinDynamics>::
getConvectionBoundaryProcessor(int x0, int x1, int y0, int y1, T* uAv)
{
  return nullptr;
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int xNormal, int yNormal>
Momenta<T,Lattice>*
RegularizedBoundaryManager2D<T,Lattice,MixinDynamics>::getExternalVelocityCornerMomenta()
{
  return new FixedVelocityBM<T,Lattice>;
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int xNormal, int yNormal>
Dynamics<T,Lattice>* RegularizedBoundaryManager2D<T,Lattice,MixinDynamics>::
getExternalVelocityCornerDynamics(T omega, Momenta<T,Lattice>& momenta)
{
  return new MixinDynamics(omega, momenta);
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int xNormal, int yNormal>
PostProcessorGenerator2D<T,Lattice>*
RegularizedBoundaryManager2D<T,Lattice,MixinDynamics>::
getExternalVelocityCornerProcessor(int x, int y)
{
  return new OuterVelocityCornerProcessorGenerator2D<T,Lattice, xNormal,yNormal> (x,y);
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int xNormal, int yNormal>
Momenta<T,Lattice>*
RegularizedBoundaryManager2D<T,Lattice,MixinDynamics>::getInternalVelocityCornerMomenta()
{
  return new InnerCornerVelBM2D<T,Lattice, xNormal,yNormal>;
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int xNormal, int yNormal>
Dynamics<T,Lattice>* RegularizedBoundaryManager2D<T,Lattice,MixinDynamics>::
getInternalVelocityCornerDynamics(T omega, Momenta<T,Lattice>& momenta)
{
  return new CombinedRLBdynamics<T,Lattice, MixinDynamics>(omega, momenta);
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int xNormal, int yNormal>
PostProcessorGenerator2D<T,Lattice>*
RegularizedBoundaryManager2D<T,Lattice,MixinDynamics>::getInternalVelocityCornerProcessor
(int x, int y)
{
  return nullptr;
}


////////// InterpolationBoundaryManager2D /////////////////////////////////////////

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int direction, int orientation>
Momenta<T,Lattice>*
InterpolationBoundaryManager2D<T,Lattice,MixinDynamics>::getVelocityBoundaryMomenta()
{
  return new BasicDirichletBM<T,Lattice,VelocityBM,direction,orientation>;
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int direction, int orientation>
Dynamics<T,Lattice>* InterpolationBoundaryManager2D<T,Lattice,MixinDynamics>::
getVelocityBoundaryDynamics(T omega, Momenta<T,Lattice>& momenta)
{
  return new MixinDynamics(omega, momenta);
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int direction, int orientation>
PostProcessorGenerator2D<T,Lattice>*
InterpolationBoundaryManager2D<T,Lattice,MixinDynamics>::
getVelocityBoundaryProcessor(int x0, int x1, int y0, int y1)
{
  return new StraightFdBoundaryProcessorGenerator2D
         < T,Lattice,direction,orientation  >  (x0, x1, y0, y1);
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int direction, int orientation>
Momenta<T,Lattice>*
InterpolationBoundaryManager2D<T,Lattice,MixinDynamics>::getPressureBoundaryMomenta()
{
  return new BasicDirichletBM<T,Lattice,PressureBM,direction,orientation>;
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int direction, int orientation>
Dynamics<T,Lattice>* InterpolationBoundaryManager2D<T,Lattice,MixinDynamics>::
getPressureBoundaryDynamics(T omega, Momenta<T,Lattice>& momenta)
{
  return new MixinDynamics(omega, momenta);
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int direction, int orientation>
PostProcessorGenerator2D<T,Lattice>*
InterpolationBoundaryManager2D<T,Lattice,MixinDynamics>::
getPressureBoundaryProcessor(int x0, int x1, int y0, int y1)
{
  return new StraightFdBoundaryProcessorGenerator2D
         < T,Lattice,direction,orientation >  (x0, x1, y0, y1);
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int direction, int orientation>
PostProcessorGenerator2D<T,Lattice>*
InterpolationBoundaryManager2D<T,Lattice,MixinDynamics>::
getConvectionBoundaryProcessor(int x0, int x1, int y0, int y1, T* uAv)
{
  return new StraightConvectionBoundaryProcessorGenerator2D
         < T,Lattice,direction,orientation >  (x0, x1, y0, y1, uAv);
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int xNormal, int yNormal>
Momenta<T,Lattice>*
InterpolationBoundaryManager2D<T,Lattice,MixinDynamics>::getExternalVelocityCornerMomenta()
{
  return new FixedVelocityBM<T,Lattice>;
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int xNormal, int yNormal>
Dynamics<T,Lattice>* InterpolationBoundaryManager2D<T,Lattice,MixinDynamics>::
getExternalVelocityCornerDynamics(T omega, Momenta<T,Lattice>& momenta)
{
  return new MixinDynamics(omega, momenta);
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int xNormal, int yNormal>
PostProcessorGenerator2D<T,Lattice>*
InterpolationBoundaryManager2D<T,Lattice,MixinDynamics>::getExternalVelocityCornerProcessor(int x, int y)
{
  return new OuterVelocityCornerProcessorGenerator2D<T,Lattice, xNormal,yNormal> (x,y);
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int xNormal, int yNormal>
Momenta<T,Lattice>*
InterpolationBoundaryManager2D<T,Lattice,MixinDynamics>::getInternalVelocityCornerMomenta()
{
  return new InnerCornerVelBM2D<T,Lattice, xNormal,yNormal>;
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int xNormal, int yNormal>
Dynamics<T,Lattice>* InterpolationBoundaryManager2D<T,Lattice,MixinDynamics>::
getInternalVelocityCornerDynamics(T omega, Momenta<T,Lattice>& momenta)
{
  return new CombinedRLBdynamics<T,Lattice, MixinDynamics>(omega, momenta);
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int xNormal, int yNormal>
PostProcessorGenerator2D<T,Lattice>*
InterpolationBoundaryManager2D<T,Lattice,MixinDynamics>::getInternalVelocityCornerProcessor (int x, int y)
{
  return nullptr;
}


////////// Factory functions //////////////////////////////////////////////////

template<typename T, template<typename U> class Lattice, typename MixinDynamics>
OnLatticeBoundaryCondition2D<T,Lattice>* createLocalBoundaryCondition2D(BlockLatticeStructure2D<T,Lattice>& block)
{
  return new BoundaryConditionInstantiator2D <
         T, Lattice,
         RegularizedBoundaryManager2D<T,Lattice, MixinDynamics> > (block);
}

template<typename T, template<typename U> class Lattice, typename MixinDynamics>
OnLatticeBoundaryCondition2D<T,Lattice>* createInterpBoundaryCondition2D(BlockLatticeStructure2D<T,Lattice>& block)
{
  return new BoundaryConditionInstantiator2D <
         T, Lattice,
         InterpolationBoundaryManager2D<T,Lattice, MixinDynamics> > (block);
}

}  // namespace olb

#endif
