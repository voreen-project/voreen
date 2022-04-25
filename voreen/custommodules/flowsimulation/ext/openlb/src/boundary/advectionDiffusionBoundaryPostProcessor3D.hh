/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2016 Robin Trunk
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
#ifndef ADVECTION_DIFFUSION_BOUNDARY_POST_PROCESSOR_3D_HH
#define ADVECTION_DIFFUSION_BOUNDARY_POST_PROCESSOR_3D_HH

#include "advectionDiffusionBoundaryPostProcessor3D.h"
#include "core/util.h"
#include "dynamics/lbm.h"

namespace olb {

////////  ConvectionBoundaryProcessor3D ////////////////////////////////

template<typename T, typename DESCRIPTOR>
ConvectionBoundaryProcessor3D<T,DESCRIPTOR>::
ConvectionBoundaryProcessor3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
                              int discreteNormalX, int discreteNormalY, int discreteNormalZ)
  : x0(x0_), x1(x1_), y0(y0_), y1(y1_), z0(z0_), z1(z1_)
{
  OLB_PRECONDITION(x0==x1 || y0==y1 || z0==z1);

  this->getName() = "ConvectionBoundaryProcessor3D";

  interpolationPop[0] = 0;
  for (int iPop = 1; iPop < DESCRIPTOR::q; iPop++) {
    interpolationPop[iPop] = 0;
    // find incoming iPop from material 0
    if (descriptors::c<DESCRIPTOR>(iPop,0)*discreteNormalX + descriptors::c<DESCRIPTOR>(iPop,1)*discreteNormalY + descriptors::c<DESCRIPTOR>(iPop,2)*discreteNormalZ > 0) {
      // check for material number of neighbours has to be one level higher
      interpolationPop[iPop] = 1;
    }
  }
}

template<typename T, typename DESCRIPTOR>
void ConvectionBoundaryProcessor3D<T,DESCRIPTOR>::
processSubDomain(BlockLattice<T,DESCRIPTOR>& blockLattice, int x0_, int x1_, int y0_,
                 int y1_, int z0_, int z1_)
{
  int newX0, newX1, newY0, newY1, newZ0, newZ1;
  if ( util::intersect (
         x0, x1, y0, y1, z0, z1,
         x0_, x1_, y0_, y1_, z0_, z1_,
         newX0, newX1, newY0, newY1, newZ0, newZ1 ) ) {

    auto& pop = blockLattice.template getField<descriptors::POPULATION>();

#ifdef PARALLEL_MODE_OMP
    #pragma omp parallel for
#endif
    for (int iPop = 1; iPop < DESCRIPTOR::q; ++iPop) {
      if (interpolationPop[iPop] != 0) {
        for (int iX=newX0; iX<=newX1; ++iX) {
          for (int iY=newY0; iY<=newY1; ++iY) {
            for (int iZ=newZ0; iZ<=newZ1; ++iZ) {
              const auto c = descriptors::c<DESCRIPTOR>(iPop);
              if (blockLattice.isInside(iX+c[0],iY+c[0],iZ+c[0]) && blockLattice.isInside(iX+2*c[0],iY+2*c[0],iZ+2*c[0])) {
                //do reflection
                T v = 0.5 * (  pop[iPop][blockLattice.getCellId(iX+  c[0],iY+  c[1],iZ+  c[2])]
                             + pop[iPop][blockLattice.getCellId(iX+2*c[0],iY+2*c[1],iZ+2*c[2])]);
                pop[iPop][blockLattice.getCellId(iX,iY,iZ)] = v;
              }
            }
          }
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void ConvectionBoundaryProcessor3D<T,DESCRIPTOR>::
process(BlockLattice<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1, z0, z1);
}

////////  ZeroDistributionBoundaryProcessor3D ////////////////////////////////

template<typename T, typename DESCRIPTOR>
ZeroDistributionBoundaryProcessor3D<T,DESCRIPTOR>::
ZeroDistributionBoundaryProcessor3D(int x0_, int x1_, int y0_, int y1_, int z0_,
                                    int z1_, int discreteNormalX, int discreteNormalY, int discreteNormalZ)
  : x0(x0_), x1(x1_), y0(y0_), y1(y1_), z0(z0_), z1(z1_)
{
  this->getName() = "ZeroDistributionBoundaryProcessor3D";
  OLB_PRECONDITION(x0==x1 || y0==y1 || z0==z1);

  resetPop[0] = 0;
  for (int iPop = 1; iPop < DESCRIPTOR::q; iPop++) {
    resetPop[iPop] = 0;
    // find incoming iPop from material 0
    if (descriptors::c<DESCRIPTOR>(iPop,0)*discreteNormalX + descriptors::c<DESCRIPTOR>(iPop,1)*discreteNormalY + descriptors::c<DESCRIPTOR>(iPop,2)*discreteNormalZ > 0) {
      resetPop[iPop] = 1;
    }
  }
}

template<typename T, typename DESCRIPTOR>
void ZeroDistributionBoundaryProcessor3D<T,DESCRIPTOR>::
processSubDomain(BlockLattice<T,DESCRIPTOR>& blockLattice, int x0_, int x1_,
                 int y0_, int y1_, int z0_, int z1_)
{
  int newX0, newX1, newY0, newY1, newZ0, newZ1;
  if ( util::intersect (
         x0, x1, y0, y1, z0, z1,
         x0_, x1_, y0_, y1_, z0_, z1_,
         newX0, newX1, newY0, newY1, newZ0, newZ1 ) ) {

#ifdef PARALLEL_MODE_OMP
    #pragma omp parallel for
#endif
    for (int iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
        for (int iZ=newZ0; iZ<=newZ1; ++iZ) {
          for (int iPop = 1; iPop < DESCRIPTOR::q ; ++iPop) {
            if (resetPop[iPop]!=0) {
              blockLattice.get(iX,iY,iZ)[iPop] = -descriptors::t<T,DESCRIPTOR>(iPop);
            }
          }
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void ZeroDistributionBoundaryProcessor3D<T,DESCRIPTOR>::
process(BlockLattice<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1, z0, z1);
}

////////  ExtFieldBoundaryProcessor3D ////////////////////////////////

template<typename T, typename DESCRIPTOR, typename FIELD_A, typename FIELD_B>
ExtFieldBoundaryProcessor3D<T,DESCRIPTOR,FIELD_A,FIELD_B>::
ExtFieldBoundaryProcessor3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
                            int discreteNormalX_, int discreteNormalY_, int discreteNormalZ_)
  : x0(x0_), x1(x1_), y0(y0_), y1(y1_), z0(z0_), z1(z1_),
    discreteNormalX(discreteNormalX_), discreteNormalY(discreteNormalY_),
    discreteNormalZ(discreteNormalZ_)
{
  OLB_PRECONDITION(x0==x1 || y0==y1 || z0==z1);
  this->getName() = "ExtFieldBoundaryProcessor3D";
  tick = true;
}

template<typename T, typename DESCRIPTOR, typename FIELD_A, typename FIELD_B>
void ExtFieldBoundaryProcessor3D<T,DESCRIPTOR,FIELD_A,FIELD_B>::
processSubDomain(BlockLattice<T,DESCRIPTOR>& blockLattice,
                 int x0_, int x1_,
                 int y0_, int y1_,
                 int z0_, int z1_)
{
  int newX0, newX1, newY0, newY1, newZ0, newZ1;

  if ( util::intersect (
         x0, x1, y0, y1, z0, z1,
         x0_, x1_, y0_, y1_, z0_, z1_,
         newX0, newX1, newY0, newY1, newZ0, newZ1 ) ) {
#ifdef PARALLEL_MODE_OMP
    #pragma omp parallel for
#endif
    for (int iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
        for (int iZ=newZ0; iZ<=newZ1; ++iZ) {
          auto self = blockLattice.get(iX,iY,iZ);
          auto neighbor = self.neighbor({discreteNormalX, discreteNormalY, discreteNormalZ});
          if (tick) {
            self.template setField<FIELD_A>(neighbor.template getField<FIELD_B>());
          }
          else {
            self.template setField<FIELD_B>(neighbor.template getField<FIELD_A>());
          }
        }
      }
    }
  }
  tick = !tick;
}

template<typename T, typename DESCRIPTOR, typename FIELD_A, typename FIELD_B>
void ExtFieldBoundaryProcessor3D<T,DESCRIPTOR,FIELD_A,FIELD_B>::
process(BlockLattice<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1, z0, z1);
}

////////  ConvectionBoundaryProcessorGenerator3D ////////////////////////////////

template<typename T, typename DESCRIPTOR>
ConvectionBoundaryProcessorGenerator3D<T,DESCRIPTOR>::
ConvectionBoundaryProcessorGenerator3D(int x0_, int x1_, int y0_, int y1_, int z0_,
                                       int z1_, int discreteNormalX_, int discreteNormalY_, int discreteNormalZ_)
  : PostProcessorGenerator3D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_, z0_, z1_),
    discreteNormalX(discreteNormalX_), discreteNormalY(discreteNormalY_),
    discreteNormalZ(discreteNormalZ_)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor3D<T,DESCRIPTOR>*
ConvectionBoundaryProcessorGenerator3D<T,DESCRIPTOR>::generate() const
{
  return new ConvectionBoundaryProcessor3D<T,DESCRIPTOR>(this->x0, this->x1, this->y0,
         this->y1, this->z0, this->z1,
         discreteNormalX,
         discreteNormalY,
         discreteNormalZ);
}

template<typename T, typename DESCRIPTOR>
PostProcessorGenerator3D<T,DESCRIPTOR>*
ConvectionBoundaryProcessorGenerator3D<T,DESCRIPTOR>::clone() const
{
  return new ConvectionBoundaryProcessorGenerator3D<T,DESCRIPTOR>(this->x0, this->x1,
         this->y0, this->y1, this->z0, this->z1,
         discreteNormalX, discreteNormalY, discreteNormalZ);
}

////////  ZeroDistributionBoundaryProcessorGenerator3D ////////////////////////////////

template<typename T, typename DESCRIPTOR>
ZeroDistributionBoundaryProcessorGenerator3D<T,DESCRIPTOR>::
ZeroDistributionBoundaryProcessorGenerator3D(int x0_, int x1_, int y0_, int y1_,
    int z0_, int z1_, int discreteNormalX_, int discreteNormalY_, int discreteNormalZ_)
  : PostProcessorGenerator3D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_, z0_, z1_),
    discreteNormalX(discreteNormalX_), discreteNormalY(discreteNormalY_),
    discreteNormalZ(discreteNormalZ_)
{ }

template<typename T, typename DESCRIPTOR>
PostProcessor3D<T,DESCRIPTOR>*
ZeroDistributionBoundaryProcessorGenerator3D<T,DESCRIPTOR>::generate() const
{
  return new ZeroDistributionBoundaryProcessor3D<T,DESCRIPTOR>(this->x0, this->x1,
         this->y0, this->y1, this->z0, this->z1,
         discreteNormalX, discreteNormalY, discreteNormalZ);
}

template<typename T, typename DESCRIPTOR>
PostProcessorGenerator3D<T,DESCRIPTOR>*
ZeroDistributionBoundaryProcessorGenerator3D<T,DESCRIPTOR>::clone() const
{
  return new ZeroDistributionBoundaryProcessorGenerator3D<T,DESCRIPTOR>(this->x0,
         this->x1, this->y0, this->y1, this->z0, this->z1,
         discreteNormalX, discreteNormalY, discreteNormalZ);
}

////////  ExtFieldBoundaryProcessorGenerator3D ////////////////////////////////

template<typename T, typename DESCRIPTOR, typename FIELD_A, typename FIELD_B>
ExtFieldBoundaryProcessorGenerator3D<T,DESCRIPTOR,FIELD_A,FIELD_B>::
ExtFieldBoundaryProcessorGenerator3D(int x0_, int x1_, int y0_, int y1_, int z0_,
                                     int z1_, int discreteNormalX_, int discreteNormalY_, int discreteNormalZ_)
  : PostProcessorGenerator3D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_, z0_, z1_),
    discreteNormalX(discreteNormalX_),
    discreteNormalY(discreteNormalY_),
    discreteNormalZ(discreteNormalZ_)
{ }

template<typename T, typename DESCRIPTOR, typename FIELD_A, typename FIELD_B>
PostProcessor3D<T,DESCRIPTOR>*
ExtFieldBoundaryProcessorGenerator3D<T,DESCRIPTOR,FIELD_A,FIELD_B>::generate() const
{
  return new ExtFieldBoundaryProcessor3D<T,DESCRIPTOR,FIELD_A,FIELD_B>(this->x0, this->x1, this->y0,
         this->y1, this->z0, this->z1,
         discreteNormalX, discreteNormalY, discreteNormalZ);
}

template<typename T, typename DESCRIPTOR, typename FIELD_A, typename FIELD_B>
PostProcessorGenerator3D<T,DESCRIPTOR>*
ExtFieldBoundaryProcessorGenerator3D<T,DESCRIPTOR,FIELD_A,FIELD_B>::clone() const
{
  return new ExtFieldBoundaryProcessorGenerator3D<T,DESCRIPTOR,FIELD_A,FIELD_B>(this->x0, this->x1,
         this->y0, this->y1, this->z0, this->z1,
         discreteNormalX, discreteNormalY, discreteNormalZ);
}

}  // namespace olb


#endif
