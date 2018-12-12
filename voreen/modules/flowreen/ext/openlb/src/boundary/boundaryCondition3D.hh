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
 * A helper for initialising 3D boundaries -- generic implementation.
 */
#ifndef BOUNDARY_CONDITION_3D_HH
#define BOUNDARY_CONDITION_3D_HH

#include "boundaryCondition3D.h"
#include "boundaryInstantiator3D.h"
#include "momentaOnBoundaries3D.h"

namespace olb {

template<typename T, template<typename U> class Lattice, class MixinDynamics>
class RegularizedBoundaryManager3D {
public:
  template<int direction, int orientation> static Momenta<T,Lattice>*
  getVelocityBoundaryMomenta();
  template<int direction, int orientation> static Dynamics<T,Lattice>*
  getVelocityBoundaryDynamics(T omega, Momenta<T,Lattice>& momenta);
  template<int direction, int orientation> static PostProcessorGenerator3D<T,Lattice>*
  getVelocityBoundaryProcessor(int x0, int x1, int y0, int y1, int z0, int z1);

  template<int direction, int orientation> static Momenta<T,Lattice>*
  getPressureBoundaryMomenta();
  template<int direction, int orientation> static Dynamics<T,Lattice>*
  getPressureBoundaryDynamics(T omega, Momenta<T,Lattice>& momenta);
  template<int direction, int orientation> static PostProcessorGenerator3D<T,Lattice>*
  getPressureBoundaryProcessor(int x0, int x1, int y0, int y1, int z0, int z1);

  template<int direction, int orientation> static PostProcessorGenerator3D<T,Lattice>*
  getConvectionBoundaryProcessor(int x0, int x1, int y0, int y1, int z0, int z1, T* uAv=NULL);

  template<int plane, int normal1, int normal2> static Momenta<T,Lattice>*
  getExternalVelocityEdgeMomenta();
  template<int plane, int normal1, int normal2> static Dynamics<T,Lattice>*
  getExternalVelocityEdgeDynamics(T omega, Momenta<T,Lattice>& momenta);
  template<int plane, int normal1, int normal2> static PostProcessorGenerator3D<T,Lattice>*
  getExternalVelocityEdgeProcessor(int x0, int x1, int y0, int y1, int z0, int z1);

  template<int plane, int normal1, int normal2> static Momenta<T,Lattice>*
  getInternalVelocityEdgeMomenta();
  template<int plane, int normal1, int normal2> static Dynamics<T,Lattice>*
  getInternalVelocityEdgeDynamics(T omega, Momenta<T,Lattice>& momenta);
  template<int plane, int normal1, int normal2> static PostProcessorGenerator3D<T,Lattice>*
  getInternalVelocityEdgeProcessor(int x0, int x1, int y0, int y1, int z0, int z1);

  template<int xNormal, int yNormal, int zNormal> static Momenta<T,Lattice>*
  getExternalVelocityCornerMomenta();
  template<int xNormal, int yNormal, int zNormal> static Dynamics<T,Lattice>*
  getExternalVelocityCornerDynamics(T omega, Momenta<T,Lattice>& momenta);
  template<int xNormal, int yNormal, int zNormal> static PostProcessorGenerator3D<T,Lattice>*
  getExternalVelocityCornerProcessor(int x, int y, int z);

  template<int xNormal, int yNormal, int zNormal> static Momenta<T,Lattice>*
  getInternalVelocityCornerMomenta();
  template<int xNormal, int yNormal, int zNormal> static Dynamics<T,Lattice>*
  getInternalVelocityCornerDynamics(T omega, Momenta<T,Lattice>& momenta);
  template<int xNormal, int yNormal, int zNormal> static PostProcessorGenerator3D<T,Lattice>*
  getInternalVelocityCornerProcessor(int x, int y, int z);
};

template<typename T, template<typename U> class Lattice, class MixinDynamics>
class InterpolationBoundaryManager3D {
public:
  template<int direction, int orientation> static Momenta<T,Lattice>*
  getVelocityBoundaryMomenta();
  template<int direction, int orientation> static Dynamics<T,Lattice>*
  getVelocityBoundaryDynamics(T omega, Momenta<T,Lattice>& momenta);
  template<int direction, int orientation> static PostProcessorGenerator3D<T,Lattice>*
  getVelocityBoundaryProcessor(int x0, int x1, int y0, int y1, int z0, int z1);

  template<int direction, int orientation> static Momenta<T,Lattice>*
  getPressureBoundaryMomenta();
  template<int direction, int orientation> static Dynamics<T,Lattice>*
  getPressureBoundaryDynamics(T omega, Momenta<T,Lattice>& momenta);
  template<int direction, int orientation> static PostProcessorGenerator3D<T,Lattice>*
  getPressureBoundaryProcessor(int x0, int x1, int y0, int y1, int z0, int z1);

  template<int direction, int orientation> static PostProcessorGenerator3D<T,Lattice>*
  getConvectionBoundaryProcessor(int x0, int x1, int y0, int y1, int z0, int z1, T* uAv=NULL);

  template<int plane, int normal1, int normal2> static Momenta<T,Lattice>*
  getExternalVelocityEdgeMomenta();
  template<int plane, int normal1, int normal2> static Dynamics<T,Lattice>*
  getExternalVelocityEdgeDynamics(T omega, Momenta<T,Lattice>& momenta);
  template<int plane, int normal1, int normal2> static PostProcessorGenerator3D<T,Lattice>*
  getExternalVelocityEdgeProcessor(int x0, int x1, int y0, int y1, int z0, int z1);

  template<int plane, int normal1, int normal2> static Momenta<T,Lattice>*
  getInternalVelocityEdgeMomenta();
  template<int plane, int normal1, int normal2> static Dynamics<T,Lattice>*
  getInternalVelocityEdgeDynamics(T omega, Momenta<T,Lattice>& momenta);
  template<int plane, int normal1, int normal2> static PostProcessorGenerator3D<T,Lattice>*
  getInternalVelocityEdgeProcessor(int x0, int x1, int y0, int y1, int z0, int z1);

  template<int xNormal, int yNormal, int zNormal> static Momenta<T,Lattice>*
  getExternalVelocityCornerMomenta();
  template<int xNormal, int yNormal, int zNormal> static Dynamics<T,Lattice>*
  getExternalVelocityCornerDynamics(T omega, Momenta<T,Lattice>& momenta);
  template<int xNormal, int yNormal, int zNormal> static PostProcessorGenerator3D<T,Lattice>*
  getExternalVelocityCornerProcessor(int x, int y, int z);

  template<int xNormal, int yNormal, int zNormal> static Momenta<T,Lattice>*
  getInternalVelocityCornerMomenta();
  template<int xNormal, int yNormal, int zNormal> static Dynamics<T,Lattice>*
  getInternalVelocityCornerDynamics(T omega, Momenta<T,Lattice>& momenta);
  template<int xNormal, int yNormal, int zNormal> static PostProcessorGenerator3D<T,Lattice>*
  getInternalVelocityCornerProcessor(int x, int y, int z);
};


////////// RegularizedBoundaryManager3D /////////////////////////////////////////

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int direction, int orientation>
Momenta<T,Lattice>*
RegularizedBoundaryManager3D<T,Lattice,MixinDynamics>::getVelocityBoundaryMomenta()
{
  return new RegularizedVelocityBM<T,Lattice, direction,orientation>;
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int direction, int orientation>
Dynamics<T,Lattice>* RegularizedBoundaryManager3D<T,Lattice,MixinDynamics>::
getVelocityBoundaryDynamics(T omega, Momenta<T,Lattice>& momenta)
{
  return new CombinedRLBdynamics<T,Lattice, MixinDynamics>(omega, momenta);
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int direction, int orientation>
PostProcessorGenerator3D<T,Lattice>*
RegularizedBoundaryManager3D<T,Lattice,MixinDynamics>::
getVelocityBoundaryProcessor(int x0, int x1, int y0, int y1, int z0, int z1)
{
  return nullptr;
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int direction, int orientation>
Momenta<T,Lattice>*
RegularizedBoundaryManager3D<T,Lattice,MixinDynamics>::getPressureBoundaryMomenta()
{
  return new RegularizedPressureBM<T,Lattice, direction,orientation>;
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int direction, int orientation>
Dynamics<T,Lattice>* RegularizedBoundaryManager3D<T,Lattice,MixinDynamics>::
getPressureBoundaryDynamics(T omega, Momenta<T,Lattice>& momenta)
{
  return new CombinedRLBdynamics<T,Lattice, MixinDynamics>(omega, momenta);
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int direction, int orientation>
PostProcessorGenerator3D<T,Lattice>*
RegularizedBoundaryManager3D<T,Lattice,MixinDynamics>::
getPressureBoundaryProcessor(int x0, int x1, int y0, int y1, int z0, int z1)
{
  return nullptr;
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int direction, int orientation>
PostProcessorGenerator3D<T,Lattice>*
RegularizedBoundaryManager3D<T,Lattice,MixinDynamics>::
getConvectionBoundaryProcessor(int x0, int x1, int y0, int y1, int z0, int z1, T* uAv)
{
  return nullptr;
}


template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int plane, int normal1, int normal2>
Momenta<T,Lattice>*
RegularizedBoundaryManager3D<T,Lattice,MixinDynamics>::getExternalVelocityEdgeMomenta()
{
  return new FixedVelocityBM<T,Lattice>;
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int plane, int normal1, int normal2>
Dynamics<T,Lattice>*
RegularizedBoundaryManager3D<T,Lattice,MixinDynamics>::
getExternalVelocityEdgeDynamics(T omega, Momenta<T,Lattice>& momenta)
{
  return new MixinDynamics(omega, momenta);
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int plane, int normal1, int normal2>
PostProcessorGenerator3D<T,Lattice>*
RegularizedBoundaryManager3D<T,Lattice,MixinDynamics>::
getExternalVelocityEdgeProcessor(int x0, int x1, int y0, int y1, int z0, int z1)
{
  return new OuterVelocityEdgeProcessorGenerator3D<T,Lattice, plane,normal1,normal2>(x0,x1, y0,y1, z0,z1);
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int plane, int normal1, int normal2>
Momenta<T,Lattice>*
RegularizedBoundaryManager3D<T,Lattice,MixinDynamics>::getInternalVelocityEdgeMomenta()
{
  return new InnerEdgeVelBM3D<T,Lattice, plane,normal1,normal2>;
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int plane, int normal1, int normal2>
Dynamics<T,Lattice>*
RegularizedBoundaryManager3D<T,Lattice,MixinDynamics>::
getInternalVelocityEdgeDynamics(T omega, Momenta<T,Lattice>& momenta)
{
  return new CombinedRLBdynamics<T,Lattice, MixinDynamics>(omega, momenta);
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int plane, int normal1, int normal2>
PostProcessorGenerator3D<T,Lattice>*
RegularizedBoundaryManager3D<T,Lattice,MixinDynamics>::
getInternalVelocityEdgeProcessor(int x0, int x1, int y0, int y1, int z0, int z1)
{
  return nullptr;
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int xNormal, int yNormal, int zNormal>
Momenta<T,Lattice>*
RegularizedBoundaryManager3D<T,Lattice,MixinDynamics>::getExternalVelocityCornerMomenta()
{
  return new FixedVelocityBM<T,Lattice>;
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int xNormal, int yNormal, int zNormal>
Dynamics<T,Lattice>*
RegularizedBoundaryManager3D<T,Lattice,MixinDynamics>::
getExternalVelocityCornerDynamics(T omega, Momenta<T,Lattice>& momenta)
{
  return new MixinDynamics(omega, momenta);
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int xNormal, int yNormal, int zNormal>
PostProcessorGenerator3D<T,Lattice>*
RegularizedBoundaryManager3D<T,Lattice,MixinDynamics>::
getExternalVelocityCornerProcessor(int x, int y, int z)
{
  return new OuterVelocityCornerProcessorGenerator3D<T,Lattice, xNormal,yNormal,zNormal> (x,y,z);
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int xNormal, int yNormal, int zNormal>
Momenta<T,Lattice>*
RegularizedBoundaryManager3D<T,Lattice,MixinDynamics>::getInternalVelocityCornerMomenta()
{
  return new InnerCornerVelBM3D<T,Lattice, xNormal,yNormal,zNormal>;
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int xNormal, int yNormal, int zNormal>
Dynamics<T,Lattice>*
RegularizedBoundaryManager3D<T,Lattice,MixinDynamics>::
getInternalVelocityCornerDynamics(T omega, Momenta<T,Lattice>& momenta)
{
  return new CombinedRLBdynamics<T,Lattice, MixinDynamics>(omega, momenta);
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int xNormal, int yNormal, int zNormal>
PostProcessorGenerator3D<T,Lattice>*
RegularizedBoundaryManager3D<T,Lattice,MixinDynamics>::
getInternalVelocityCornerProcessor(int x, int y, int z)
{
  return nullptr;
}


////////// InterpolationBoundaryManager3D /////////////////////////////////////////

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int direction, int orientation>
Momenta<T,Lattice>*
InterpolationBoundaryManager3D<T,Lattice,MixinDynamics>::getVelocityBoundaryMomenta()
{
  return new BasicDirichletBM<T,Lattice,VelocityBM, direction,orientation>;
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int direction, int orientation>
Dynamics<T,Lattice>* InterpolationBoundaryManager3D<T,Lattice,MixinDynamics>::
getVelocityBoundaryDynamics(T omega, Momenta<T,Lattice>& momenta)
{
  return new MixinDynamics(omega, momenta);
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int direction, int orientation>
PostProcessorGenerator3D<T,Lattice>*
InterpolationBoundaryManager3D<T,Lattice,MixinDynamics>::
getVelocityBoundaryProcessor(int x0, int x1, int y0, int y1, int z0, int z1)
{
  return new PlaneFdBoundaryProcessorGenerator3D
         <T,Lattice, direction,orientation>(x0,x1, y0,y1, z0,z1);
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int direction, int orientation>
Momenta<T,Lattice>*
InterpolationBoundaryManager3D<T,Lattice,MixinDynamics>::getPressureBoundaryMomenta()
{
  return new BasicDirichletBM<T,Lattice,PressureBM, direction,orientation>;
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int direction, int orientation>
Dynamics<T,Lattice>* InterpolationBoundaryManager3D<T,Lattice,MixinDynamics>::
getPressureBoundaryDynamics(T omega, Momenta<T,Lattice>& momenta)
{
  return new MixinDynamics(omega, momenta);
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int direction, int orientation>
PostProcessorGenerator3D<T,Lattice>*
InterpolationBoundaryManager3D<T,Lattice,MixinDynamics>::
getPressureBoundaryProcessor(int x0, int x1, int y0, int y1, int z0, int z1)
{
  return new PlaneFdBoundaryProcessorGenerator3D
         <T,Lattice, direction,orientation>(x0, x1, y0, y1, z0, z1);
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int direction, int orientation>
PostProcessorGenerator3D<T,Lattice>*
InterpolationBoundaryManager3D<T,Lattice,MixinDynamics>::
getConvectionBoundaryProcessor(int x0, int x1, int y0, int y1, int z0, int z1, T* uAv)
{
  return new StraightConvectionBoundaryProcessorGenerator3D
         <T,Lattice,direction,orientation>(x0, x1, y0, y1, z0, z1, uAv);
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int plane, int normal1, int normal2>
Momenta<T,Lattice>*
InterpolationBoundaryManager3D<T,Lattice,MixinDynamics>::getExternalVelocityEdgeMomenta()
{
  return new FixedVelocityBM<T,Lattice>;
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int plane, int normal1, int normal2>
Dynamics<T,Lattice>*
InterpolationBoundaryManager3D<T,Lattice,MixinDynamics>::
getExternalVelocityEdgeDynamics(T omega, Momenta<T,Lattice>& momenta)
{
  return new MixinDynamics(omega, momenta);
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int plane, int normal1, int normal2>
PostProcessorGenerator3D<T,Lattice>*
InterpolationBoundaryManager3D<T,Lattice,MixinDynamics>::
getExternalVelocityEdgeProcessor(int x0, int x1, int y0, int y1, int z0, int z1)
{
  return new OuterVelocityEdgeProcessorGenerator3D<T,Lattice, plane,normal1,normal2>(x0,x1, y0,y1, z0,z1);
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int plane, int normal1, int normal2>
Momenta<T,Lattice>*
InterpolationBoundaryManager3D<T,Lattice,MixinDynamics>::getInternalVelocityEdgeMomenta()
{
  return new InnerEdgeVelBM3D<T,Lattice, plane,normal1,normal2>;
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int plane, int normal1, int normal2>
Dynamics<T,Lattice>*
InterpolationBoundaryManager3D<T,Lattice,MixinDynamics>::
getInternalVelocityEdgeDynamics(T omega, Momenta<T,Lattice>& momenta)
{
  return new CombinedRLBdynamics<T,Lattice, MixinDynamics>(omega, momenta);
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int plane, int normal1, int normal2>
PostProcessorGenerator3D<T,Lattice>*
InterpolationBoundaryManager3D<T,Lattice,MixinDynamics>::
getInternalVelocityEdgeProcessor(int x0, int x1, int y0, int y1, int z0, int z1)
{
  return nullptr;
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int xNormal, int yNormal, int zNormal>
Momenta<T,Lattice>*
InterpolationBoundaryManager3D<T,Lattice,MixinDynamics>::getExternalVelocityCornerMomenta()
{
  return new FixedVelocityBM<T,Lattice>;
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int xNormal, int yNormal, int zNormal>
Dynamics<T,Lattice>*
InterpolationBoundaryManager3D<T,Lattice,MixinDynamics>::
getExternalVelocityCornerDynamics(T omega, Momenta<T,Lattice>& momenta)
{
  return new MixinDynamics(omega, momenta);
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int xNormal, int yNormal, int zNormal>
PostProcessorGenerator3D<T,Lattice>*
InterpolationBoundaryManager3D<T,Lattice,MixinDynamics>::
getExternalVelocityCornerProcessor(int x, int y, int z)
{
  return new OuterVelocityCornerProcessorGenerator3D<T,Lattice, xNormal,yNormal,zNormal> (x,y,z);
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int xNormal, int yNormal, int zNormal>
Momenta<T,Lattice>*
InterpolationBoundaryManager3D<T,Lattice,MixinDynamics>::getInternalVelocityCornerMomenta()
{
  return new InnerCornerVelBM3D<T,Lattice, xNormal,yNormal,zNormal>;
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int xNormal, int yNormal, int zNormal>
Dynamics<T,Lattice>*
InterpolationBoundaryManager3D<T,Lattice,MixinDynamics>::
getInternalVelocityCornerDynamics(T omega, Momenta<T,Lattice>& momenta)
{
  return new CombinedRLBdynamics<T,Lattice, MixinDynamics>(omega, momenta);
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int xNormal, int yNormal, int zNormal>
PostProcessorGenerator3D<T,Lattice>*
InterpolationBoundaryManager3D<T,Lattice,MixinDynamics>::
getInternalVelocityCornerProcessor(int x, int y, int z)
{
  return nullptr;
}


////////// Factory functions //////////////////////////////////////////////////

template<typename T, template<typename U> class Lattice, typename MixinDynamics>
OnLatticeBoundaryCondition3D<T,Lattice>* createLocalBoundaryCondition3D(BlockLatticeStructure3D<T,Lattice>& block)
{
  return new BoundaryConditionInstantiator3D <
         T, Lattice,
         RegularizedBoundaryManager3D<T,Lattice, MixinDynamics> > (block);
}

template<typename T, template<typename U> class Lattice, typename MixinDynamics>
OnLatticeBoundaryCondition3D<T,Lattice>* createInterpBoundaryCondition3D(BlockLatticeStructure3D<T,Lattice>& block)
{
  return new BoundaryConditionInstantiator3D <
         T, Lattice,
         InterpolationBoundaryManager3D<T,Lattice, MixinDynamics> > (block);
}

}  // namespace olb

#endif
