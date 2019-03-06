/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2008 Orestis Malaspinas, Andrea Parmigiani, Jonas Latt
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

#ifndef NAVIER_STOKES_ADVECTION_DIFFUSION_COUPLING_POST_PROCESSOR_2D_HH
#define NAVIER_STOKES_ADVECTION_DIFFUSION_COUPLING_POST_PROCESSOR_2D_HH

#include "latticeDescriptors.h"
#include "navierStokesAdvectionDiffusionCouplingPostProcessor2D.h"
#include "core/blockLattice2D.h"
#include "core/util.h"
#include "core/finiteDifference2D.h"

using namespace std;

namespace olb {

//=====================================================================================
//==============  NavierStokesAdvectionDiffusionCouplingPostProcessor2D ===============
//=====================================================================================

template<typename T, template<typename U> class Lattice>
NavierStokesAdvectionDiffusionCouplingPostProcessor2D<T,Lattice>::
NavierStokesAdvectionDiffusionCouplingPostProcessor2D(int x0_, int x1_, int y0_, int y1_,
    T gravity_, T T0_, T deltaTemp_, std::vector<T> dir_,
    std::vector<SpatiallyExtendedObject2D* > partners_)
  :  x0(x0_), x1(x1_), y0(y0_), y1(y1_),
     gravity(gravity_), T0(T0_), deltaTemp(deltaTemp_),
     dir(dir_), partners(partners_)
{
  // we normalize the direction of force vector
  T normDir = T();
  for (unsigned iD = 0; iD < dir.size(); ++iD) {
    normDir += dir[iD]*dir[iD];
  }
  normDir = sqrt(normDir);
  for (unsigned iD = 0; iD < dir.size(); ++iD) {
    dir[iD] /= normDir;
  }

  for (unsigned iD = 0; iD < dir.size(); ++iD) {
    forcePrefactor[iD] = gravity * dir[iD];
  }

  tPartner = dynamic_cast<BlockLattice2D<T,AdvectionDiffusionD2Q5Descriptor> *>(partners[0]);
}

template<typename T, template<typename U> class Lattice>
void NavierStokesAdvectionDiffusionCouplingPostProcessor2D<T,Lattice>::
processSubDomain(BlockLattice2D<T,Lattice>& blockLattice,
                 int x0_, int x1_, int y0_, int y1_)
{

  int newX0, newX1, newY0, newY1;
  if ( util::intersect (
         x0, x1, y0, y1,
         x0_, x1_, y0_, y1_,
         newX0, newX1, newY0, newY1 ) ) {

    for (int iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
        // Velocity coupling
        T *u = tPartner->get(iX,iY).getExternal(velOffset);
        blockLattice.get(iX,iY).computeU(u);

        // computation of the bousinessq force
        T *force = blockLattice.get(iX,iY).getExternal(forceOffset);
        T temperatureDifference = tPartner->get(iX,iY).computeRho() - T0;
        for (unsigned iD = 0; iD < L::d; ++iD) {
          force[iD] = forcePrefactor[iD] * temperatureDifference;
        }
      }
    }
  }
}

template<typename T, template<typename U> class Lattice>
void NavierStokesAdvectionDiffusionCouplingPostProcessor2D<T,Lattice>::
process(BlockLattice2D<T,Lattice>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1);
}

/// LatticeCouplingGenerator for advectionDiffusion coupling

template<typename T, template<typename U> class Lattice>
NavierStokesAdvectionDiffusionCouplingGenerator2D<T,Lattice>::
NavierStokesAdvectionDiffusionCouplingGenerator2D(int x0_, int x1_, int y0_, int y1_,
    T gravity_, T T0_, T deltaTemp_, std::vector<T> dir_)
  : LatticeCouplingGenerator2D<T,Lattice>(x0_, x1_, y0_, y1_),
    gravity(gravity_), T0(T0_), deltaTemp(deltaTemp_), dir(dir_)
{ }

template<typename T, template<typename U> class Lattice>
PostProcessor2D<T,Lattice>* NavierStokesAdvectionDiffusionCouplingGenerator2D<T,Lattice>::generate (
  std::vector<SpatiallyExtendedObject2D* > partners) const
{
  return new NavierStokesAdvectionDiffusionCouplingPostProcessor2D<T,Lattice>(
           this->x0,this->x1,this->y0,this->y1, gravity, T0, deltaTemp, dir,partners);
}

template<typename T, template<typename U> class Lattice>
LatticeCouplingGenerator2D<T,Lattice>* NavierStokesAdvectionDiffusionCouplingGenerator2D<T,Lattice>::clone() const
{
  return new NavierStokesAdvectionDiffusionCouplingGenerator2D<T,Lattice>(*this);
}

}  // namespace olb

#endif
