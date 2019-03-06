/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007 Orestis Malaspinas, Jonas Latt
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

#ifndef EXTENDED_FINITE_DIFFERENCE_BOUNDARY_3D_HH
#define EXTENDED_FINITE_DIFFERENCE_BOUNDARY_3D_HH

#include "extendedFiniteDifferenceBoundary3D.h"
#include "core/finiteDifference3D.h"
#include "core/blockLattice3D.h"
#include "core/util.h"
#include "dynamics/lbHelpers.h"
#include "dynamics/firstOrderLbHelpers.h"
#include "boundaryInstantiator3D.h"
#include "momentaOnBoundaries3D.h"


namespace olb {

////////  ExtendedFdPlaneBoundaryPostProcessor3D ///////////////////////////////////

template<typename T, template<typename U> class Lattice, int direction, int orientation>
ExtendedFdPlaneBoundaryPostProcessor3D <T,Lattice,direction,orientation>::
ExtendedFdPlaneBoundaryPostProcessor3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
  : x0(x0_), x1(x1_), y0(y0_), y1(y1_), z0(z0_), z1(z1_)
{
  OLB_PRECONDITION(x0==x1 || y0==y1 || z0==z1);
}

template<typename T, template<typename U> class Lattice, int direction, int orientation>
void ExtendedFdPlaneBoundaryPostProcessor3D<T,Lattice,direction,orientation>::
processSubDomain(BlockLattice3D<T,Lattice>& blockLattice,
                 int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
{
  typedef Lattice<T> L;
  using namespace olb::util::tensorIndices3D;
  typedef lbHelpers<T,Lattice> lbH;
  enum {x,y,z};

  int newX0, newX1, newY0, newY1, newZ0, newZ1;
  if ( util::intersect (
         x0, x1, y0, y1, z0, z1,
         x0_, x1_, y0_, y1_, z0_, z1_,
         newX0, newX1, newY0, newY1, newZ0, newZ1 ) ) {

    for (int iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
        for (int iZ=newZ0; iZ<=newZ1; ++iZ) {
          Cell<T,Lattice>& cell = blockLattice.get(iX,iY,iZ);
          T rho, u[L::d];
          cell.computeRhoU(rho,u);
          T dx_U[Lattice<T>::d], dy_U[Lattice<T>::d], dz_U[Lattice<T>::d];
          interpolateGradients<0>(blockLattice, dx_U, iX, iY, iZ);
          interpolateGradients<1>(blockLattice, dy_U, iX, iY, iZ);
          interpolateGradients<2>(blockLattice, dz_U, iX, iY, iZ);

          T rhoGradU[L::d][L::d];
          rhoGradU[x][x] = rho *dx_U[x];
          rhoGradU[x][y] = rho *dx_U[y];
          rhoGradU[x][z] = rho *dx_U[z];
          rhoGradU[y][x] = rho *dy_U[x];
          rhoGradU[y][y] = rho *dy_U[y];
          rhoGradU[y][z] = rho *dy_U[z];
          rhoGradU[z][x] = rho *dz_U[x];
          rhoGradU[z][y] = rho *dz_U[y];
          rhoGradU[z][z] = rho *dz_U[z];

          T omega = blockLattice.getDynamics(iX, iY, iZ) -> getOmega();
          T sToPi = - (T)1 / Lattice<T>::invCs2 / omega;
          T pi[util::TensorVal<Lattice<T> >::n];

          pi[xx] = (T)2 * rhoGradU[x][x] * sToPi;
          pi[yy] = (T)2 * rhoGradU[y][y] * sToPi;
          pi[zz] = (T)2 * rhoGradU[z][z] * sToPi;
          pi[xy] = (rhoGradU[x][y] + rhoGradU[y][x]) * sToPi;
          pi[xz] = (rhoGradU[x][z] + rhoGradU[z][x]) * sToPi;
          pi[yz] = (rhoGradU[y][z] + rhoGradU[z][y]) * sToPi;

          // here ends the "regular" fdBoudaryCondition
          // implemented in OpenLB

          T uSqr = util::normSqr<T,Lattice<T>::d>(u);

          // first we compute the term
          // (c_{i\alpha} \nabla_\beta)(rho*u_\alpha*u_\beta)
          T dx_rho, dy_rho, dz_rho;
          interpolateGradients<0>(blockLattice, dx_rho, iX, iY, iZ);
          interpolateGradients<1>(blockLattice, dy_rho, iX, iY, iZ);
          interpolateGradients<2>(blockLattice, dz_rho, iX, iY, iZ);
          for (int iPop = 0; iPop < L::q; ++iPop) {
            T cGradRhoUU = T();
            for (int iAlpha=0; iAlpha < L::d; ++iAlpha) {
              cGradRhoUU += L::c[iPop][iAlpha] * (
                              dx_rho*u[iAlpha]*u[x] +
                              dx_U[iAlpha]*rho*u[x] +
                              dx_U[x]*rho*u[iAlpha] + //end of dx derivative
                              dy_rho*u[iAlpha]*u[y] +
                              dy_U[iAlpha]*rho*u[y] +
                              dy_U[y]*rho*u[iAlpha] +//end of dy derivative
                              dz_rho*u[iAlpha]*u[z] +
                              dz_U[iAlpha]*rho*u[z] +
                              dz_U[z]*rho*u[iAlpha]);
            }

            // then we compute the term
            // c_{i\gamma}\nabla_{\gamma}(\rho*u_\alpha * u_\beta)
            T cDivRhoUU[L::d][L::d]; //first step towards QcdivRhoUU
            for (int iAlpha = 0; iAlpha < L::d; ++iAlpha) {
              for (int iBeta = 0; iBeta < L::d; ++iBeta) {
                cDivRhoUU[iAlpha][iBeta] = L::c[iPop][x]*
                                           (dx_rho*u[iAlpha]*u[iBeta] +
                                            dx_U[iAlpha]*rho*u[iBeta] +
                                            dx_U[iBeta]*rho*u[iAlpha])
                                           +L::c[iPop][y]*
                                           (dy_rho*u[iAlpha]*u[iBeta] +
                                            dy_U[iAlpha]*rho*u[iBeta] +
                                            dy_U[iBeta]*rho*u[iAlpha])
                                           +L::c[iPop][z]*
                                           (dz_rho*u[iAlpha]*u[iBeta] +
                                            dz_U[iAlpha]*rho*u[iBeta] +
                                            dz_U[iBeta]*rho*u[iAlpha]);
              }
            }

            //Finally we can compute
            // Q_{i\alpha\beta}c_{i\gamma}\nabla_{\gamma}(\rho*u_\alpha * u_\beta)
            // and Q_{i\alpha\beta}\rho\nabla_{\alpha}u_\beta
            T qCdivRhoUU = T();
            T qRhoGradU = T();
            for (int iAlpha = 0; iAlpha < L::d; ++iAlpha) {
              for (int iBeta = 0; iBeta < L::d; ++iBeta) {
                int ci_ci = L::c[iPop][iAlpha]*L::c[iPop][iBeta];
                qCdivRhoUU  += ci_ci * cDivRhoUU[iAlpha][iBeta];
                qRhoGradU += ci_ci * rhoGradU[iAlpha][iBeta];
                if (iAlpha == iBeta) {
                  qCdivRhoUU -= cDivRhoUU[iAlpha][iBeta]/L::invCs2;
                  qRhoGradU  -= rhoGradU[iAlpha][iBeta]/L::invCs2;
                }
              }
            }

            // we then can reconstruct the value of the populations
            // according to the complete C-E expansion term
            cell[iPop] = lbH::equilibrium(iPop,rho,u,uSqr)
                         - L::t[iPop] * L::invCs2 / omega
                         * (qRhoGradU - cGradRhoUU + 0.5*L::invCs2*qCdivRhoUU);
          }
        }
      }
    }
  }
}

template<typename T, template<typename U> class Lattice, int direction, int orientation>
void ExtendedFdPlaneBoundaryPostProcessor3D<T,Lattice,direction,orientation>::
process(BlockLattice3D<T,Lattice>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1, z0, z1);
}

template<typename T, template<typename U> class Lattice, int direction, int orientation>
template<int deriveDirection>
void ExtendedFdPlaneBoundaryPostProcessor3D<T,Lattice,direction,orientation>::
interpolateGradients(BlockLattice3D<T,Lattice> const& blockLattice,
                     T velDeriv[Lattice<T>::d],
                     int iX, int iY, int iZ) const
{
  fd::DirectedGradients3D<T, Lattice, direction, orientation, deriveDirection, direction==deriveDirection>::
  interpolateVector(velDeriv, blockLattice, iX, iY, iZ);
}

template<typename T, template<typename U> class Lattice, int direction, int orientation>
template<int deriveDirection>
void ExtendedFdPlaneBoundaryPostProcessor3D<T,Lattice,direction,orientation>::
interpolateGradients(BlockLattice3D<T,Lattice> const& blockLattice,
                     T& rhoDeriv, int iX, int iY, int iZ) const
{
  fd::DirectedGradients3D<T, Lattice, direction, orientation, deriveDirection, direction==deriveDirection>::
  interpolateScalar(rhoDeriv, blockLattice, iX, iY, iZ);
}


////////  ExtendedFdPlaneBoundaryProcessorGenertor3D ///////////////////////////////

template<typename T, template<typename U> class Lattice, int direction, int orientation>
ExtendedFdPlaneBoundaryProcessorGenerator3D<T,Lattice,direction,orientation>::
ExtendedFdPlaneBoundaryProcessorGenerator3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
  : PostProcessorGenerator3D<T,Lattice>(x0_, x1_, y0_, y1_, z0_, z1_)
{ }

template<typename T, template<typename U> class Lattice, int direction, int orientation>
PostProcessor3D<T,Lattice>*
ExtendedFdPlaneBoundaryProcessorGenerator3D<T,Lattice,direction,orientation>::generate() const
{
  return new ExtendedFdPlaneBoundaryPostProcessor3D<T,Lattice,direction,orientation>
         (this->x0, this->x1, this->y0, this->y1, this->z0, this->z1);
}

template<typename T, template<typename U> class Lattice, int direction, int orientation>
PostProcessorGenerator3D<T,Lattice>*
ExtendedFdPlaneBoundaryProcessorGenerator3D<T,Lattice,direction,orientation>::clone() const
{
  return new ExtendedFdPlaneBoundaryProcessorGenerator3D<T,Lattice,direction,orientation>
         (this->x0, this->x1, this->y0, this->y1, this->z0, this->z1);
}


////////// Class ExtendedFdBoundaryManager3D /////////////////////////////////////////

template<typename T, template<typename U> class Lattice, class MixinDynamics>
class ExtendedFdBoundaryManager3D {
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
template<int direction, int orientation>
Momenta<T,Lattice>*
ExtendedFdBoundaryManager3D<T,Lattice,MixinDynamics>::getVelocityBoundaryMomenta()
{
  return new BasicDirichletBM<T,Lattice,VelocityBM, direction,orientation>;
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int direction, int orientation>
Dynamics<T,Lattice>* ExtendedFdBoundaryManager3D<T,Lattice,MixinDynamics>::
getVelocityBoundaryDynamics(T omega, Momenta<T,Lattice>& momenta)
{
  return new MixinDynamics(omega, momenta);
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int direction, int orientation>
PostProcessorGenerator3D<T,Lattice>*
ExtendedFdBoundaryManager3D<T,Lattice,MixinDynamics>::
getVelocityBoundaryProcessor(int x0, int x1, int y0, int y1, int z0, int z1)
{
  return new ExtendedFdPlaneBoundaryProcessorGenerator3D
         <T,Lattice, direction,orientation>(x0,x1, y0,y1, z0,z1);
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int direction, int orientation>
Momenta<T,Lattice>*
ExtendedFdBoundaryManager3D<T,Lattice,MixinDynamics>::getPressureBoundaryMomenta()
{
  return new BasicDirichletBM<T,Lattice,PressureBM, direction,orientation>;
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int direction, int orientation>
Dynamics<T,Lattice>* ExtendedFdBoundaryManager3D<T,Lattice,MixinDynamics>::
getPressureBoundaryDynamics(T omega, Momenta<T,Lattice>& momenta)
{
  return new MixinDynamics(omega, momenta);
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int direction, int orientation>
PostProcessorGenerator3D<T,Lattice>*
ExtendedFdBoundaryManager3D<T,Lattice,MixinDynamics>::
getPressureBoundaryProcessor(int x0, int x1, int y0, int y1, int z0, int z1)
{
  return new ExtendedFdPlaneBoundaryProcessorGenerator3D
         <T,Lattice, direction,orientation>(x0,x1, y0,y1, z0,z1);
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int direction, int orientation>
PostProcessorGenerator3D<T,Lattice>*
ExtendedFdBoundaryManager3D<T,Lattice,MixinDynamics>::
getConvectionBoundaryProcessor(int x0, int x1, int y0, int y1, int z0, int z1, T* uAv)
{
  return nullptr;
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int plane, int normal1, int normal2>
Momenta<T,Lattice>*
ExtendedFdBoundaryManager3D<T,Lattice,MixinDynamics>::getExternalVelocityEdgeMomenta()
{
  return new FixedVelocityBM<T,Lattice>;
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int plane, int normal1, int normal2>
Dynamics<T,Lattice>*
ExtendedFdBoundaryManager3D<T,Lattice,MixinDynamics>::
getExternalVelocityEdgeDynamics(T omega, Momenta<T,Lattice>& momenta)
{
  return new MixinDynamics(omega, momenta);
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int plane, int normal1, int normal2>
PostProcessorGenerator3D<T,Lattice>*
ExtendedFdBoundaryManager3D<T,Lattice,MixinDynamics>::
getExternalVelocityEdgeProcessor(int x0, int x1, int y0, int y1, int z0, int z1)
{
  return new OuterVelocityEdgeProcessorGenerator3D<T,Lattice, plane,normal1,normal2>(x0,x1, y0,y1, z0,z1);
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int plane, int normal1, int normal2>
Momenta<T,Lattice>*
ExtendedFdBoundaryManager3D<T,Lattice,MixinDynamics>::getInternalVelocityEdgeMomenta()
{
  return new InnerEdgeVelBM3D<T,Lattice, plane,normal1,normal2>;
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int plane, int normal1, int normal2>
Dynamics<T,Lattice>*
ExtendedFdBoundaryManager3D<T,Lattice,MixinDynamics>::
getInternalVelocityEdgeDynamics(T omega, Momenta<T,Lattice>& momenta)
{
  return new CombinedRLBdynamics<T,Lattice, MixinDynamics>(omega, momenta);
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int plane, int normal1, int normal2>
PostProcessorGenerator3D<T,Lattice>*
ExtendedFdBoundaryManager3D<T,Lattice,MixinDynamics>::
getInternalVelocityEdgeProcessor(int x0, int x1, int y0, int y1, int z0, int z1)
{
  return nullptr;
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int xNormal, int yNormal, int zNormal>
Momenta<T,Lattice>*
ExtendedFdBoundaryManager3D<T,Lattice,MixinDynamics>::getExternalVelocityCornerMomenta()
{
  return new FixedVelocityBM<T,Lattice>;
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int xNormal, int yNormal, int zNormal>
Dynamics<T,Lattice>*
ExtendedFdBoundaryManager3D<T,Lattice,MixinDynamics>::
getExternalVelocityCornerDynamics(T omega, Momenta<T,Lattice>& momenta)
{
  return new MixinDynamics(omega, momenta);
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int xNormal, int yNormal, int zNormal>
PostProcessorGenerator3D<T,Lattice>*
ExtendedFdBoundaryManager3D<T,Lattice,MixinDynamics>::
getExternalVelocityCornerProcessor(int x, int y, int z)
{
  return new OuterVelocityCornerProcessorGenerator3D<T,Lattice, xNormal,yNormal,zNormal> (x,y,z);
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int xNormal, int yNormal, int zNormal>
Momenta<T,Lattice>*
ExtendedFdBoundaryManager3D<T,Lattice,MixinDynamics>::getInternalVelocityCornerMomenta()
{
  return new InnerCornerVelBM3D<T,Lattice, xNormal,yNormal,zNormal>;
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int xNormal, int yNormal, int zNormal>
Dynamics<T,Lattice>*
ExtendedFdBoundaryManager3D<T,Lattice,MixinDynamics>::
getInternalVelocityCornerDynamics(T omega, Momenta<T,Lattice>& momenta)
{
  return new CombinedRLBdynamics<T,Lattice, MixinDynamics>(omega, momenta);
}

template<typename T, template<typename U> class Lattice, class MixinDynamics>
template<int xNormal, int yNormal, int zNormal>
PostProcessorGenerator3D<T,Lattice>*
ExtendedFdBoundaryManager3D<T,Lattice,MixinDynamics>::
getInternalVelocityCornerProcessor(int x, int y, int z)
{
  return nullptr;
}


////////// Factory function //////////////////////////////////////////////////

template<typename T, template<typename U> class Lattice, typename MixinDynamics>
OnLatticeBoundaryCondition3D<T,Lattice>*
createExtendedFdBoundaryCondition3D(BlockLatticeStructure3D<T,Lattice>& block)
{
  return new BoundaryConditionInstantiator3D <
         T, Lattice,
         ExtendedFdBoundaryManager3D<T,Lattice, MixinDynamics> > (block);
}


}  // namespace olb

#endif
