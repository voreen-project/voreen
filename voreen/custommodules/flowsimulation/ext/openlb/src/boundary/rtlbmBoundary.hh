/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2017 Albert Mink
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

#ifndef RTLBM_BOUNDARY_HH
#define RTLBM_BOUNDARY_HH

//#include "advectionDiffusionBoundaries.h"
//#include "advectionDiffusionLatticeDescriptors.h"
//#include "core/util.h"
//#include "advectionDiffusionLbHelpers.h"

#include "rtlbmBoundaryDynamics.h"
#include "advectionDiffusionBoundaries.h"   //TODO to be replaced
#include "advectionDiffusionBoundaryInstantiator3D.h"   //TODO to be replaced

namespace olb {


template<typename T, template<typename U> class Lattice>
class RtlbmBoundaryManager3D {
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

//================== Flat  ================================
template<typename T, template<typename U> class Lattice>
template<int direction, int orientation>
Momenta<T,Lattice>* RtlbmBoundaryManager3D<T,Lattice>::getTemperatureBoundaryMomenta()
{
  return new EquilibriumBM<T,Lattice>;
}

template<typename T, template<typename U> class Lattice>
template<int direction, int orientation>
Dynamics<T,Lattice>* RtlbmBoundaryManager3D<T,Lattice>::getTemperatureBoundaryDynamics(T omega, Momenta<T,Lattice>& momenta)
{
  return new RtlbmBoundaryDynamics<T,Lattice,direction,orientation>(omega, momenta);
}

template<typename T, template<typename U> class Lattice>
template<int direction, int orientation>
PostProcessorGenerator3D<T,Lattice>* RtlbmBoundaryManager3D<T,Lattice>::getTemperatureBoundaryProcessor(int x0, int x1, int y0, int y1, int z0, int z1)
{
  return nullptr;
}

//================== Edges ================================
template<typename T, template<typename U> class Lattice>
template<int plane, int normal1, int normal2>
Momenta<T,Lattice>* RtlbmBoundaryManager3D<T,Lattice>::getTemperatureBoundaryEdgeMomenta()
{
  return new EquilibriumBM<T,Lattice>;
}

template<typename T, template<typename U> class Lattice>
template<int plane, int normal1, int normal2>
Dynamics<T,Lattice>* RtlbmBoundaryManager3D<T,Lattice>::getTemperatureBoundaryEdgeDynamics(T omega, Momenta<T,Lattice>& momenta)
{
  return new AdvectionDiffusionEdgesDynamics<T,Lattice,BGKdynamics<T,Lattice>,plane,normal1,normal2>(omega, momenta);  // TODO AM mark as placeholder
}

template<typename T, template<typename U> class Lattice>
template<int plane, int normal1, int normal2>
PostProcessorGenerator3D<T,Lattice>* RtlbmBoundaryManager3D<T,Lattice>::getTemperatureBoundaryEdgeProcessor(int x0, int x1, int y0, int y1, int z0, int z1)
{
  return nullptr;
}

//================== Corners ================================
template<typename T, template<typename U> class Lattice>
template<int xNormal, int yNormal, int zNormal>
Momenta<T,Lattice>* RtlbmBoundaryManager3D<T,Lattice>::getTemperatureBoundaryCornerMomenta()
{
  return new EquilibriumBM<T,Lattice>;
}

template<typename T, template<typename U> class Lattice>
template<int xNormal, int yNormal, int zNormal>
Dynamics<T,Lattice>* RtlbmBoundaryManager3D<T,Lattice>::getTemperatureBoundaryCornerDynamics(T omega, Momenta<T,Lattice>& momenta)
{
  return new AdvectionDiffusionCornerDynamics3D<T,Lattice,BGKdynamics<T,Lattice>,xNormal,yNormal,zNormal>(omega, momenta); // TODO AM mark as placeholder
}

template<typename T, template<typename U> class Lattice>
template<int xNormal, int yNormal, int zNormal>
PostProcessorGenerator3D<T,Lattice>* RtlbmBoundaryManager3D<T,Lattice>::getTemperatureBoundaryCornerProcessor(int x, int y, int z)
{
  return nullptr;
}

//================== creator function ================================
// blockLattice creator
template<typename T, template<typename U> class Lattice>
OnLatticeAdvectionDiffusionBoundaryCondition3D<T,Lattice>* createRtlbmBoundaryCondition3D(BlockLatticeStructure3D<T,Lattice>& block)
{
  return new AdvectionDiffusionBoundaryConditionInstantiator3D<T, Lattice, RtlbmBoundaryManager3D<T,Lattice> > (block);
}

// superLattice creator
template<typename T, template<typename U> class Lattice>
void createRtlbmBoundaryCondition3D(sOnLatticeBoundaryCondition3D<T,Lattice>& sBC)
{
  sBC.setOverlap(1);
  for (int iC = 0; iC < sBC.getSuperLattice().getLoadBalancer().size(); iC++) {
    OnLatticeAdvectionDiffusionBoundaryCondition3D<T,Lattice>* blockBC =
      createRtlbmBoundaryCondition3D<T,Lattice>( sBC.getSuperLattice().getExtendedBlockLattice(iC) );
    sBC.getADblockBCs().push_back(blockBC);
  }
}


template<typename T, template<typename U> class Lattice>
class RtlbmDiffuseBoundaryManager3D {
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

//================== Flat  ================================
template<typename T, template<typename U> class Lattice>
template<int direction, int orientation>
Momenta<T,Lattice>* RtlbmDiffuseBoundaryManager3D<T,Lattice>::getTemperatureBoundaryMomenta()
{
  return new EquilibriumBM<T,Lattice>;
}

template<typename T, template<typename U> class Lattice>
template<int direction, int orientation>
Dynamics<T,Lattice>* RtlbmDiffuseBoundaryManager3D<T,Lattice>::getTemperatureBoundaryDynamics(T omega, Momenta<T,Lattice>& momenta)
{
  return new RtlbmDiffuseBoundaryDynamics<T,Lattice,direction,orientation>(omega, momenta);
}

template<typename T, template<typename U> class Lattice>
template<int direction, int orientation>
PostProcessorGenerator3D<T,Lattice>* RtlbmDiffuseBoundaryManager3D<T,Lattice>::getTemperatureBoundaryProcessor(int x0, int x1, int y0, int y1, int z0, int z1)
{
  return nullptr;
}

//================== Edges ================================
template<typename T, template<typename U> class Lattice>
template<int plane, int normal1, int normal2>
Momenta<T,Lattice>* RtlbmDiffuseBoundaryManager3D<T,Lattice>::getTemperatureBoundaryEdgeMomenta()
{
  return new EquilibriumBM<T,Lattice>;
}

template<typename T, template<typename U> class Lattice>
template<int plane, int normal1, int normal2>
Dynamics<T,Lattice>* RtlbmDiffuseBoundaryManager3D<T,Lattice>::getTemperatureBoundaryEdgeDynamics(T omega, Momenta<T,Lattice>& momenta)
{
  return new RtlbmDiffuseEdgeBoundaryDynamics<T,Lattice,plane,normal1,normal2>(omega, momenta);
}

template<typename T, template<typename U> class Lattice>
template<int plane, int normal1, int normal2>
PostProcessorGenerator3D<T,Lattice>* RtlbmDiffuseBoundaryManager3D<T,Lattice>::getTemperatureBoundaryEdgeProcessor(int x0, int x1, int y0, int y1, int z0, int z1)
{
  return nullptr;
}

//================== Corners ================================
template<typename T, template<typename U> class Lattice>
template<int xNormal, int yNormal, int zNormal>
Momenta<T,Lattice>* RtlbmDiffuseBoundaryManager3D<T,Lattice>::getTemperatureBoundaryCornerMomenta()
{
  return new EquilibriumBM<T,Lattice>;
}

template<typename T, template<typename U> class Lattice>
template<int xNormal, int yNormal, int zNormal>
Dynamics<T,Lattice>* RtlbmDiffuseBoundaryManager3D<T,Lattice>::getTemperatureBoundaryCornerDynamics(T omega, Momenta<T,Lattice>& momenta)
{
  return new RtlbmDiffuseCornerBoundaryDynamics<T,Lattice,xNormal,yNormal,zNormal>(omega, momenta);
}

template<typename T, template<typename U> class Lattice>
template<int xNormal, int yNormal, int zNormal>
PostProcessorGenerator3D<T,Lattice>* RtlbmDiffuseBoundaryManager3D<T,Lattice>::getTemperatureBoundaryCornerProcessor(int x, int y, int z)
{
  return nullptr;
}


//================== creator function ================================
// blockLattice creator
template<typename T, template<typename U> class Lattice>
OnLatticeAdvectionDiffusionBoundaryCondition3D<T,Lattice>* createRtlbmDiffuseBoundaryCondition3D(BlockLatticeStructure3D<T,Lattice>& block)
{
  return new AdvectionDiffusionBoundaryConditionInstantiator3D<T, Lattice, RtlbmDiffuseBoundaryManager3D<T,Lattice> > (block);   // TODO AM mark as placeholder
}

// superLattice creator
template<typename T, template<typename U> class Lattice>
void createRtlbmDiffuseBoundaryCondition3D(sOnLatticeBoundaryCondition3D<T,Lattice>& sBC)
{
  sBC.setOverlap(1);
  for (int iC = 0; iC < sBC.getSuperLattice().getLoadBalancer().size(); iC++) {
    OnLatticeAdvectionDiffusionBoundaryCondition3D<T,Lattice>* blockBC =
      createRtlbmDiffuseBoundaryCondition3D<T,Lattice>( sBC.getSuperLattice().getExtendedBlockLattice(iC) );
    sBC.getADblockBCs().push_back(blockBC);
  }
}

}  // namespace olb

#endif
