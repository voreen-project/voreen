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
 * A helper for initialising 2D boundaries -- generic implementation.
 */
#ifndef ADVECTION_DIFFUSION_BOUNDARY_CONDITION_2D_HH
#define ADVECTION_DIFFUSION_BOUNDARY_CONDITION_2D_HH

#include "advectionDiffusionBoundaries.h"
#include "advectionDiffusionBoundaries.hh"
#include "advectionDiffusionBoundaryCondition2D.h"
#include "advectionDiffusionBoundaryInstantiator2D.h"


namespace olb {
template<typename T, template<typename U> class Lattice, class MixinDynamics>
class AdvectionDiffusionBoundaryManager2D {
public:
  template<int direction, int orientation> static Momenta<T,Lattice>*
  getTemperatureBoundaryMomenta();
  template<int direction, int orientation> static Dynamics<T,Lattice>*
  getTemperatureBoundaryDynamics(T omega, Momenta<T,Lattice>& momenta);
  template<int direction, int orientation> static PostProcessorGenerator2D<T,Lattice>*
  getTemperatureBoundaryProcessor(int x0, int x1, int y0, int y1);

  template<int xNormal, int yNormal> static Momenta<T,Lattice>*
  getTemperatureCornerMomenta();
  template<int xNormal, int yNormal> static Dynamics<T,Lattice>*
  getTemperatureCornerDynamics(T omega, Momenta<T,Lattice>& momenta);
  template<int xNormal, int yNormal> static PostProcessorGenerator2D<T,Lattice>*
  getTemperatureCornerProcessor(int x, int y);

};


////////// AdvectionDiffusionBoundaryManager2D /////////////////////////////////////////

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int direction, int orientation>
Momenta<T,Lattice>*
AdvectionDiffusionBoundaryManager2D<T,Lattice,MixinDynamics>::getTemperatureBoundaryMomenta()
{
  return new EquilibriumBM<T,Lattice>;
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int direction, int orientation>
Dynamics<T,Lattice>* AdvectionDiffusionBoundaryManager2D<T,Lattice,MixinDynamics>::
getTemperatureBoundaryDynamics(T omega, Momenta<T,Lattice>& momenta)
{
  return new AdvectionDiffusionBoundariesDynamics<T,Lattice,MixinDynamics,direction,orientation>(omega, momenta);
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int direction, int orientation>
PostProcessorGenerator2D<T,Lattice>*
AdvectionDiffusionBoundaryManager2D<T,Lattice,MixinDynamics>::
getTemperatureBoundaryProcessor(int x0, int x1, int y0, int y1)
{
  return nullptr;
}

//==================  Corners ================================

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int xNormal, int yNormal>
Momenta<T,Lattice>*
AdvectionDiffusionBoundaryManager2D<T,Lattice,MixinDynamics>::getTemperatureCornerMomenta()
{
  return new EquilibriumBM<T,Lattice>;
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int xNormal, int yNormal>
Dynamics<T,Lattice>* AdvectionDiffusionBoundaryManager2D<T,Lattice,MixinDynamics>::
getTemperatureCornerDynamics(T omega, Momenta<T,Lattice>& momenta)
{
  return new AdvectionDiffusionCornerDynamics2D<T,Lattice,MixinDynamics,xNormal,yNormal>(omega, momenta);
}



template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int xNormal, int yNormal>
PostProcessorGenerator2D<T,Lattice>*
AdvectionDiffusionBoundaryManager2D<T,Lattice,MixinDynamics>::
getTemperatureCornerProcessor(int x, int y)
{
  return nullptr;
}

////////// Factory functions //////////////////////////////////////////////////

template<typename T, template<typename U> class Lattice, typename MixinDynamics>
OnLatticeAdvectionDiffusionBoundaryCondition2D<T,Lattice>*
createAdvectionDiffusionBoundaryCondition2D(BlockLatticeStructure2D<T,Lattice>& block)
{
  return new AdvectionDiffusionBoundaryConditionInstantiator2D<T, Lattice,
         AdvectionDiffusionBoundaryManager2D<T,Lattice, MixinDynamics> > (block);
}

template<typename T, template<typename U> class Lattice, typename MixinDynamics>
void createAdvectionDiffusionBoundaryCondition2D(sOnLatticeBoundaryCondition2D<T, Lattice>& sBC)
{
  int nC = sBC.getSuperLattice().getLoadBalancer().size();
  sBC.setOverlap(1);
  for (int iC = 0; iC < nC; iC++) {
    OnLatticeAdvectionDiffusionBoundaryCondition2D<T, Lattice>* blockBC =
      createAdvectionDiffusionBoundaryCondition2D<T,Lattice,MixinDynamics>(
        sBC.getSuperLattice().getExtendedBlockLattice(iC));
    sBC.getADblockBCs().push_back(blockBC);
  }
}


}  // namespace olb

#endif
