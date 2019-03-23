/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2008 Orestis Malaspinas, Andrea Parmigiani
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
#ifndef ADVECTION_DIFFUSION_BOUNDARY_CONDITION_3D_HH
#define ADVECTION_DIFFUSION_BOUNDARY_CONDITION_3D_HH

#include "advectionDiffusionBoundaries.h"
#include "advectionDiffusionBoundaries.hh"
#include "advectionDiffusionBoundaryCondition3D.h"
#include "advectionDiffusionBoundaryInstantiator3D.h"
#include "advectionDiffusionBoundaryInstantiator3D.hh"


namespace olb {

template<typename T, template<typename U> class Lattice, class MixinDynamics>
class AdvectionDiffusionBoundaryManager3D {
public:
  template<int direction, int orientation>
  static Momenta<T,Lattice>* getTemperatureBoundaryMomenta();
  template<int direction, int orientation>
  static Dynamics<T,Lattice>* getTemperatureBoundaryDynamics(T omega, Momenta<T,Lattice>& momenta);
  template<int direction, int orientation>
  static PostProcessorGenerator3D<T,Lattice>* getTemperatureBoundaryProcessor(int x0, int x1, int y0, int y1, int z0, int z1);

  template<int plane, int normal1, int normal2>
  static Momenta<T,Lattice>* getTemperatureBoundaryEdgeMomenta();
  template<int plane, int normal1, int normal2>
  static Dynamics<T,Lattice>* getTemperatureBoundaryEdgeDynamics(T omega, Momenta<T,Lattice>& momenta);
  template<int plane, int normal1, int normal2>
  static PostProcessorGenerator3D<T,Lattice>* getTemperatureBoundaryEdgeProcessor(int x0, int x1, int y0, int y1, int z0, int z1);

  template<int normalX, int normalY, int normalZ>
  static Momenta<T,Lattice>* getTemperatureBoundaryCornerMomenta();
  template<int normalX, int normalY, int normalZ>
  static Dynamics<T,Lattice>* getTemperatureBoundaryCornerDynamics(T omega, Momenta<T,Lattice>& momenta);
  template<int normalX, int normalY, int normalZ>
  static PostProcessorGenerator3D<T,Lattice>* getTemperatureBoundaryCornerProcessor(int x, int y, int z);

};

//==================  Flat ================================
template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int direction, int orientation>
Momenta<T,Lattice>* AdvectionDiffusionBoundaryManager3D<T,Lattice,MixinDynamics>::
getTemperatureBoundaryMomenta()
{
  return new EquilibriumBM<T,Lattice>;
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int direction, int orientation>
Dynamics<T,Lattice>* AdvectionDiffusionBoundaryManager3D<T,Lattice,MixinDynamics>::
getTemperatureBoundaryDynamics(T omega, Momenta<T,Lattice>& momenta)
{
  return new AdvectionDiffusionBoundariesDynamics<T,Lattice,MixinDynamics,direction,orientation>(omega, momenta);
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int direction, int orientation>
PostProcessorGenerator3D<T,Lattice>* AdvectionDiffusionBoundaryManager3D<T,Lattice,MixinDynamics>::
getTemperatureBoundaryProcessor(int x0, int x1, int y0, int y1, int z0, int z1)
{
  return nullptr;
}

//==================  Edges ================================
template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int plane, int normal1, int normal2>
Momenta<T,Lattice>* AdvectionDiffusionBoundaryManager3D<T,Lattice,MixinDynamics>::
getTemperatureBoundaryEdgeMomenta()
{
  return new EquilibriumBM<T,Lattice>;
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int plane, int normal1, int normal2>
Dynamics<T,Lattice>* AdvectionDiffusionBoundaryManager3D<T,Lattice,MixinDynamics>::
getTemperatureBoundaryEdgeDynamics(T omega, Momenta<T,Lattice>& momenta)
{
  return new AdvectionDiffusionEdgesDynamics<T,Lattice,MixinDynamics,plane,normal1,normal2>(omega, momenta);
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int plane, int normal1, int normal2>
PostProcessorGenerator3D<T,Lattice>* AdvectionDiffusionBoundaryManager3D<T,Lattice,MixinDynamics>::
getTemperatureBoundaryEdgeProcessor(int x0, int x1, int y0, int y1, int z0, int z1)
{
  return nullptr;
}

//==================  Corners ================================
template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int xNormal, int yNormal, int zNormal>
Momenta<T,Lattice>* AdvectionDiffusionBoundaryManager3D<T,Lattice,MixinDynamics>::
getTemperatureBoundaryCornerMomenta()
{
  return new EquilibriumBM<T,Lattice>;
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int xNormal, int yNormal, int zNormal>
Dynamics<T,Lattice>* AdvectionDiffusionBoundaryManager3D<T,Lattice,MixinDynamics>::
getTemperatureBoundaryCornerDynamics(T omega, Momenta<T,Lattice>& momenta)
{
  return new AdvectionDiffusionCornerDynamics3D<T,Lattice,MixinDynamics,xNormal,yNormal,zNormal>(omega, momenta);
}



template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int xNormal, int yNormal, int zNormal>
PostProcessorGenerator3D<T,Lattice>* AdvectionDiffusionBoundaryManager3D<T,Lattice,MixinDynamics>::
getTemperatureBoundaryCornerProcessor(int x, int y, int z)
{
  return nullptr;
}


//================== Factory functions ================================
template<typename T, template<typename U> class Lattice, typename MixinDynamics>
OnLatticeAdvectionDiffusionBoundaryCondition3D<T,Lattice>*
createAdvectionDiffusionBoundaryCondition3D(BlockLatticeStructure3D<T,Lattice>& block)
{
  return new AdvectionDiffusionBoundaryConditionInstantiator3D<T, Lattice,
         AdvectionDiffusionBoundaryManager3D<T,Lattice, MixinDynamics> > (block);
}

template<typename T, template<typename U> class Lattice, typename MixinDynamics>
void createAdvectionDiffusionBoundaryCondition3D(sOnLatticeBoundaryCondition3D<T, Lattice>& sBC)
{
  int nC = sBC.getSuperLattice().getLoadBalancer().size();
  sBC.setOverlap(1);
  for (int iC = 0; iC < nC; iC++) {
    OnLatticeAdvectionDiffusionBoundaryCondition3D<T, Lattice>* blockBC =
      createAdvectionDiffusionBoundaryCondition3D<T,Lattice,MixinDynamics>(
        sBC.getSuperLattice().getExtendedBlockLattice(iC));
    sBC.getADblockBCs().push_back(blockBC);
  }
}


}  // namespace olb

#endif
