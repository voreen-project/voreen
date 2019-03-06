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

#ifndef BOUNDARY_POST_PROCESSORS_3D_HH
#define BOUNDARY_POST_PROCESSORS_3D_HH

#include "boundaryPostProcessors3D.h"
#include "core/finiteDifference3D.h"
#include "core/blockLattice3D.h"
#include "dynamics/firstOrderLbHelpers.h"
#include "core/util.h"
#include "utilities/vectorHelpers.h"

namespace olb {

////////  PlaneFdBoundaryProcessor3D ///////////////////////////////////

template<typename T, template<typename U> class Lattice, int direction, int orientation>
PlaneFdBoundaryProcessor3D<T,Lattice,direction,orientation>::
PlaneFdBoundaryProcessor3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
  : x0(x0_), x1(x1_), y0(y0_), y1(y1_), z0(z0_), z1(z1_)
{
  OLB_PRECONDITION(x0==x1 || y0==y1 || z0==z1);
}

template<typename T, template<typename U> class Lattice, int direction, int orientation>
void PlaneFdBoundaryProcessor3D<T,Lattice,direction,orientation>::
processSubDomain(BlockLattice3D<T,Lattice>& blockLattice,
                 int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
{
  using namespace olb::util::tensorIndices3D;

  int newX0, newX1, newY0, newY1, newZ0, newZ1;
  if ( util::intersect (
         x0, x1, y0, y1, z0, z1,
         x0_, x1_, y0_, y1_, z0_, z1_,
         newX0, newX1, newY0, newY1, newZ0, newZ1 ) ) {
    int iX;

#ifdef PARALLEL_MODE_OMP
    #pragma omp parallel for
#endif
    for (iX=newX0; iX<=newX1; ++iX) {
      T dx_u[Lattice<T>::d], dy_u[Lattice<T>::d], dz_u[Lattice<T>::d];
      for (int iY=newY0; iY<=newY1; ++iY) {
        for (int iZ=newZ0; iZ<=newZ1; ++iZ) {
          Cell<T,Lattice>& cell = blockLattice.get(iX,iY,iZ);
          Dynamics<T,Lattice>* dynamics = blockLattice.getDynamics(iX, iY, iZ);
          T rho, u[Lattice<T>::d];
          cell.computeRhoU(rho,u);

          interpolateGradients<0> ( blockLattice, dx_u, iX, iY, iZ );
          interpolateGradients<1> ( blockLattice, dy_u, iX, iY, iZ );
          interpolateGradients<2> ( blockLattice, dz_u, iX, iY, iZ );
          T dx_ux = dx_u[0];
          T dy_ux = dy_u[0];
          T dz_ux = dz_u[0];
          T dx_uy = dx_u[1];
          T dy_uy = dy_u[1];
          T dz_uy = dz_u[1];
          T dx_uz = dx_u[2];
          T dy_uz = dy_u[2];
          T dz_uz = dz_u[2];
          T omega = dynamics->getOmega();
          T sToPi = - rho / Lattice<T>::invCs2 / omega;
          T pi[util::TensorVal<Lattice<T> >::n];
          pi[xx] = (T)2 * dx_ux * sToPi;
          pi[yy] = (T)2 * dy_uy * sToPi;
          pi[zz] = (T)2 * dz_uz * sToPi;
          pi[xy] = (dx_uy + dy_ux) * sToPi;
          pi[xz] = (dx_uz + dz_ux) * sToPi;
          pi[yz] = (dy_uz + dz_uy) * sToPi;

          // Computation of the particle distribution functions
          // according to the regularized formula
          T uSqr = util::normSqr<T,Lattice<T>::d>(u);
          for (int iPop = 0; iPop < Lattice<T>::q; ++iPop)
            cell[iPop] = dynamics -> computeEquilibrium(iPop,rho,u,uSqr) +
                         firstOrderLbHelpers<T,Lattice>::fromPiToFneq(iPop, pi);
        }
      }
    }
  }
}

template<typename T, template<typename U> class Lattice, int direction, int orientation>
void PlaneFdBoundaryProcessor3D<T,Lattice,direction,orientation>::
process(BlockLattice3D<T,Lattice>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1, z0, z1);
}


template<typename T, template<typename U> class Lattice, int direction, int orientation>
template<int deriveDirection>
void PlaneFdBoundaryProcessor3D<T,Lattice,direction,orientation>::
interpolateGradients(BlockLattice3D<T,Lattice> const& blockLattice, T velDeriv[Lattice<T>::d],
                     int iX, int iY, int iZ) const
{
  fd::DirectedGradients3D<T, Lattice, direction, orientation, deriveDirection, direction==deriveDirection>::
  interpolateVector(velDeriv, blockLattice, iX, iY, iZ);
}


////////  PlaneFdBoundaryProcessorGenerator3D ///////////////////////////////

template<typename T, template<typename U> class Lattice, int direction, int orientation>
PlaneFdBoundaryProcessorGenerator3D<T,Lattice,direction,orientation>::
PlaneFdBoundaryProcessorGenerator3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
  : PostProcessorGenerator3D<T,Lattice>(x0_, x1_, y0_, y1_, z0_, z1_)
{ }

template<typename T, template<typename U> class Lattice, int direction, int orientation>
PostProcessor3D<T,Lattice>* PlaneFdBoundaryProcessorGenerator3D<T,Lattice,direction,orientation>::
generate() const
{
  return new PlaneFdBoundaryProcessor3D<T,Lattice, direction,orientation>
         ( this->x0, this->x1, this->y0, this->y1, this->z0, this->z1 );
}

template<typename T, template<typename U> class Lattice, int direction, int orientation>
PostProcessorGenerator3D<T,Lattice>*
PlaneFdBoundaryProcessorGenerator3D<T,Lattice,direction,orientation>::clone() const
{
  return new PlaneFdBoundaryProcessorGenerator3D<T,Lattice,direction,orientation>
         (this->x0, this->x1, this->y0, this->y1, this->z0, this->z1);
}


////////  StraightConvectionBoundaryProcessorGenerator3D ////////////////////////////////

template<typename T, template<typename U> class Lattice, int direction, int orientation>
StraightConvectionBoundaryProcessor3D<T,Lattice,direction,orientation>::
StraightConvectionBoundaryProcessor3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, T* uAv_)
  : x0(x0_), x1(x1_), y0(y0_), y1(y1_), z0(z0_), z1(z1_), uAv(uAv_)
{
  OLB_PRECONDITION(x0==x1 || y0==y1 || z0==z1);

  saveCell = new T*** [(size_t)(x1_-x0_+1)];
  for (int iX=0; iX<=x1_-x0_; ++iX) {
    saveCell[iX] = new T** [(size_t)(y1_-y0_+1)];
    for (int iY=0; iY<=y1_-y0_; ++iY) {
      saveCell[iX][iY] = new T* [(size_t)(z1_-z0_+1)];
      for (int iZ=0; iZ<=z1_-z0_; ++iZ) {
        saveCell[iX][iY][iZ] = new T [(size_t)(Lattice<T>::q)];
        for (int iPop=0; iPop<Lattice<T>::q; ++iPop) {
          // default set to -1 in order to avoid wrong results at first call
          saveCell[iX][iY][iZ][iPop] = T(-1);
        }
      }
    }
  }
}

template<typename T, template<typename U> class Lattice, int direction, int orientation>
StraightConvectionBoundaryProcessor3D<T,Lattice,direction,orientation>::
~StraightConvectionBoundaryProcessor3D()
{
  for (int iX=0; iX<=x1-x0; ++iX) {
    for (int iY=0; iY<=y1-y0; ++iY) {
      for (int iZ=0; iZ<=z1-z0; ++iZ) {
        delete [] saveCell[iX][iY][iZ];
      }
      delete [] saveCell[iX][iY];
    }
    delete [] saveCell[iX];
  }
  delete [] saveCell;
}

template<typename T, template<typename U> class Lattice, int direction,int orientation>
void StraightConvectionBoundaryProcessor3D<T,Lattice,direction,orientation>::
processSubDomain(BlockLattice3D<T,Lattice>& blockLattice, int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
{
  int newX0, newX1, newY0, newY1, newZ0, newZ1;
  if ( util::intersect (
         x0, x1, y0, y1, z0, z1,
         x0_, x1_, y0_, y1_, z0_, z1_,
         newX0, newX1, newY0, newY1, newZ0, newZ1 ) ) {

    int iX;
#ifdef PARALLEL_MODE_OMP
    #pragma omp parallel for
#endif
    for (iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
        for (int iZ=newZ0; iZ<=newZ1; ++iZ) {
          Cell<T,Lattice>& cell = blockLattice.get(iX,iY,iZ);
          for (int iPop = 0; iPop < Lattice<T>::q ; ++iPop) {
            if (Lattice<T>::c[iPop][direction]==-orientation) {
              // using default -1 avoids wrong first call
              if (!util::nearZero(1 + saveCell[iX-newX0][iY-newY0][iZ-newZ0][iPop]) ) {
                cell[iPop] = saveCell[iX-newX0][iY-newY0][iZ-newZ0][iPop];
              }
            }
          }

          T rho0, u0[3];
          T rho1, u1[3];
          T rho2, u2[3];
          if (direction==0) {
            blockLattice.get(iX,iY,iZ).computeRhoU(rho0,u0);
            blockLattice.get(iX-orientation,iY,iZ).computeRhoU(rho1,u1);
            blockLattice.get(iX-orientation*2,iY,iZ).computeRhoU(rho2,u2);
          } else if (direction==1) {
            blockLattice.get(iX,iY,iZ).computeRhoU(rho0,u0);
            blockLattice.get(iX,iY-orientation,iZ).computeRhoU(rho1,u1);
            blockLattice.get(iX,iY-orientation*2,iZ).computeRhoU(rho2,u2);
          } else {
            blockLattice.get(iX,iY,iZ).computeRhoU(rho0,u0);
            blockLattice.get(iX,iY,iZ-orientation).computeRhoU(rho1,u1);
            blockLattice.get(iX,iY,iZ-orientation*2).computeRhoU(rho2,u2);
          }

          // rho0 = T(1); rho1 = T(1); rho2 = T(1);

          T uDelta[3];
          T uAverage = rho0*u0[direction];
          if (uAv!=nullptr) {
            uAverage = *uAv * rho0;
          }
          uDelta[0]=-uAverage*0.5*(3*rho0*u0[0]-4*rho1*u1[0]+rho2*u2[0]);
          uDelta[1]=-uAverage*0.5*(3*rho0*u0[1]-4*rho1*u1[1]+rho2*u2[1]);
          uDelta[2]=-uAverage*0.5*(3*rho0*u0[2]-4*rho1*u1[2]+rho2*u2[2]);

          for (int iPop = 0; iPop < Lattice<T>::q ; ++iPop) {
            if (Lattice<T>::c[iPop][direction] == -orientation) {
              saveCell[iX-newX0][iY-newY0][iZ-newZ0][iPop] = cell[iPop] + Lattice<T>::invCs2*Lattice<T>::t[iPop]*(uDelta[0]*Lattice<T>::c[iPop][0]+uDelta[1]*Lattice<T>::c[iPop][1]+uDelta[2]*Lattice<T>::c[iPop][2]);
            }
          }
        }
      }
    }
  }
}

template<typename T, template<typename U> class Lattice, int direction,int orientation>
void StraightConvectionBoundaryProcessor3D<T,Lattice,direction,orientation>::
process(BlockLattice3D<T,Lattice>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1, z0, z1);
}


////////  StraightConvectionBoundaryProcessorGenerator2D ////////////////////////////////

template<typename T, template<typename U> class Lattice, int direction,int orientation>
StraightConvectionBoundaryProcessorGenerator3D<T,Lattice, direction,orientation>::
StraightConvectionBoundaryProcessorGenerator3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, T* uAv_)
  : PostProcessorGenerator3D<T,Lattice>(x0_, x1_, y0_, y1_, z0_, z1_), uAv(uAv_)
{ }

template<typename T, template<typename U> class Lattice, int direction,int orientation>
PostProcessor3D<T,Lattice>*
StraightConvectionBoundaryProcessorGenerator3D<T,Lattice,direction,orientation>::generate() const
{
  return new StraightConvectionBoundaryProcessor3D<T,Lattice,direction,orientation>
         ( this->x0, this->x1, this->y0, this->y1, this->z0, this->z1, uAv);
}

template<typename T, template<typename U> class Lattice, int direction,int orientation>
PostProcessorGenerator3D<T,Lattice>*
StraightConvectionBoundaryProcessorGenerator3D<T,Lattice,direction,orientation>::clone() const
{
  return new StraightConvectionBoundaryProcessorGenerator3D<T,Lattice,direction,orientation>
         (this->x0, this->x1, this->y0, this->y1, this->z0, this->z1, uAv);
}


////////  OuterVelocityEdgeProcessor3D ///////////////////////////////////

template<typename T, template<typename U> class Lattice, int plane, int normal1, int normal2>
OuterVelocityEdgeProcessor3D<T,Lattice, plane,normal1,normal2>::
OuterVelocityEdgeProcessor3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
  : x0(x0_), x1(x1_), y0(y0_), y1(y1_), z0(z0_), z1(z1_)
{
  OLB_PRECONDITION (
    (plane==2 && x0==x1 && y0==y1) ||
    (plane==1 && x0==x1 && z0==z1) ||
    (plane==0 && y0==y1 && z0==z1)     );

}

template<typename T, template<typename U> class Lattice, int plane, int normal1, int normal2>
void OuterVelocityEdgeProcessor3D<T,Lattice, plane,normal1,normal2>::
processSubDomain(BlockLattice3D<T,Lattice>& blockLattice,
                 int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
{
  using namespace olb::util::tensorIndices3D;

  int newX0, newX1, newY0, newY1, newZ0, newZ1;
  if ( util::intersect ( x0, x1, y0, y1, z0, z1,
                         x0_, x1_, y0_, y1_, z0_, z1_,
                         newX0, newX1, newY0, newY1, newZ0, newZ1 ) ) {
    int iX;

#ifdef PARALLEL_MODE_OMP
    #pragma omp parallel for
#endif
    for (iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
        for (int iZ=newZ0; iZ<=newZ1; ++iZ) {
          Cell<T,Lattice>& cell = blockLattice.get(iX,iY,iZ);
          Dynamics<T,Lattice>* dynamics = blockLattice.getDynamics(iX, iY, iZ);

          T rho10 = getNeighborRho(iX,iY,iZ,1,0, blockLattice);
          T rho01 = getNeighborRho(iX,iY,iZ,0,1, blockLattice);
          T rho20 = getNeighborRho(iX,iY,iZ,2,0, blockLattice);
          T rho02 = getNeighborRho(iX,iY,iZ,0,2, blockLattice);
          T rho = (T)2/(T)3*(rho01+rho10)-(T)1/(T)6*(rho02+rho20);

          T dA_uB_[3][3];
          interpolateGradients<plane,0>            ( blockLattice, dA_uB_[0], iX, iY, iZ );
          interpolateGradients<direction1,normal1> ( blockLattice, dA_uB_[1], iX, iY, iZ );
          interpolateGradients<direction2,normal2> ( blockLattice, dA_uB_[2], iX, iY, iZ );
          T dA_uB[3][3];
          for (int iBeta=0; iBeta<3; ++iBeta) {
            dA_uB[plane][iBeta]      = dA_uB_[0][iBeta];
            dA_uB[direction1][iBeta] = dA_uB_[1][iBeta];
            dA_uB[direction2][iBeta] = dA_uB_[2][iBeta];
          }
          T omega = dynamics -> getOmega();
          T sToPi = - rho / Lattice<T>::invCs2 / omega;
          T pi[util::TensorVal<Lattice<T> >::n];
          pi[xx] = (T)2 * dA_uB[0][0] * sToPi;
          pi[yy] = (T)2 * dA_uB[1][1] * sToPi;
          pi[zz] = (T)2 * dA_uB[2][2] * sToPi;
          pi[xy] = (dA_uB[0][1]+dA_uB[1][0]) * sToPi;
          pi[xz] = (dA_uB[0][2]+dA_uB[2][0]) * sToPi;
          pi[yz] = (dA_uB[1][2]+dA_uB[2][1]) * sToPi;

          // Computation of the particle distribution functions
          // according to the regularized formula
          T u[Lattice<T>::d];
          cell.computeU(u);
          T uSqr = util::normSqr<T,Lattice<T>::d>(u);

          for (int iPop = 0; iPop < Lattice<T>::q; ++iPop) {
            cell[iPop] = dynamics -> computeEquilibrium(iPop,rho,u,uSqr) +
                         firstOrderLbHelpers<T,Lattice>::fromPiToFneq(iPop, pi);
          }
        }
      }
    }
  }
}

template<typename T, template<typename U> class Lattice, int plane, int normal1, int normal2>
void OuterVelocityEdgeProcessor3D<T,Lattice, plane,normal1,normal2>::
process(BlockLattice3D<T,Lattice>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1, z0, z1);
}

template<typename T, template<typename U> class Lattice, int plane, int normal1, int normal2>
T OuterVelocityEdgeProcessor3D<T,Lattice, plane,normal1,normal2>::
getNeighborRho(int x, int y, int z, int step1, int step2, BlockLattice3D<T,Lattice> const& blockLattice)
{
  int coords[3] = {x, y, z};
  coords[direction1] += -normal1*step1;
  coords[direction2] += -normal2*step2;
  return blockLattice.get(coords[0], coords[1], coords[2]).computeRho();
}

template<typename T, template<typename U> class Lattice, int plane, int normal1, int normal2>
template<int deriveDirection, int orientation>
void OuterVelocityEdgeProcessor3D<T,Lattice, plane,normal1,normal2>::
interpolateGradients(BlockLattice3D<T,Lattice> const& blockLattice,
                     T velDeriv[Lattice<T>::d],
                     int iX, int iY, int iZ) const
{
  fd::DirectedGradients3D<T,Lattice,deriveDirection,orientation,deriveDirection,deriveDirection!=plane>::
  interpolateVector(velDeriv, blockLattice, iX, iY, iZ);
}

////////  OuterVelocityEdgeProcessorGenerator3D ///////////////////////////////

template<typename T, template<typename U> class Lattice, int plane, int normal1, int normal2>
OuterVelocityEdgeProcessorGenerator3D<T,Lattice, plane,normal1,normal2>::
OuterVelocityEdgeProcessorGenerator3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
  : PostProcessorGenerator3D<T,Lattice>(x0_, x1_, y0_, y1_, z0_, z1_)
{ }

template<typename T, template<typename U> class Lattice, int plane, int normal1, int normal2>
PostProcessor3D<T,Lattice>*
OuterVelocityEdgeProcessorGenerator3D<T,Lattice, plane,normal1,normal2>::
generate() const
{
  return new OuterVelocityEdgeProcessor3D < T,Lattice, plane,normal1,normal2 >
         ( this->x0, this->x1, this->y0, this->y1, this->z0, this->z1);
}

template<typename T, template<typename U> class Lattice, int plane, int normal1, int normal2>
PostProcessorGenerator3D<T,Lattice>*
OuterVelocityEdgeProcessorGenerator3D<T,Lattice, plane,normal1,normal2>::clone() const
{
  return new OuterVelocityEdgeProcessorGenerator3D<T,Lattice, plane,normal1,normal2 >
         (this->x0, this->x1, this->y0, this->y1, this->z0, this->z1);
}

/////////// OuterVelocityCornerProcessor3D /////////////////////////////////////

template<typename T, template<typename U> class Lattice, int xNormal, int yNormal, int zNormal>
OuterVelocityCornerProcessor3D<T, Lattice, xNormal, yNormal, zNormal>::
OuterVelocityCornerProcessor3D ( int x_, int y_, int z_ )
  : x(x_), y(y_), z(z_)
{ }

template<typename T, template<typename U> class Lattice, int xNormal, int yNormal, int zNormal>
void OuterVelocityCornerProcessor3D<T, Lattice, xNormal, yNormal, zNormal>::
process(BlockLattice3D<T,Lattice>& blockLattice)
{
  using namespace olb::util::tensorIndices3D;
  Cell<T,Lattice>& cell = blockLattice.get(x,y,z);
  Dynamics<T,Lattice>* dynamics = blockLattice.getDynamics(x, y, z);

  T rho100 = blockLattice.get(x - 1*xNormal, y - 0*yNormal, z - 0*zNormal).computeRho();
  T rho010 = blockLattice.get(x - 0*xNormal, y - 1*yNormal, z - 0*zNormal).computeRho();
  T rho001 = blockLattice.get(x - 0*xNormal, y - 0*yNormal, z - 1*zNormal).computeRho();
  T rho200 = blockLattice.get(x - 2*xNormal, y - 0*yNormal, z - 0*zNormal).computeRho();
  T rho020 = blockLattice.get(x - 0*xNormal, y - 2*yNormal, z - 0*zNormal).computeRho();
  T rho002 = blockLattice.get(x - 0*xNormal, y - 0*yNormal, z - 2*zNormal).computeRho();
  T rho = (T)4/(T)9 * (rho001 + rho010 + rho100) - (T)1/(T)9 * (rho002 + rho020 + rho200);

  T dx_u[Lattice<T>::d], dy_u[Lattice<T>::d], dz_u[Lattice<T>::d];
  fd::DirectedGradients3D<T, Lattice, 0, xNormal, 0, true>::interpolateVector(dx_u, blockLattice, x,y,z);
  fd::DirectedGradients3D<T, Lattice, 1, yNormal, 0, true>::interpolateVector(dy_u, blockLattice, x,y,z);
  fd::DirectedGradients3D<T, Lattice, 2, zNormal, 0, true>::interpolateVector(dz_u, blockLattice, x,y,z);

  T dx_ux = dx_u[0];
  T dy_ux = dy_u[0];
  T dz_ux = dz_u[0];
  T dx_uy = dx_u[1];
  T dy_uy = dy_u[1];
  T dz_uy = dz_u[1];
  T dx_uz = dx_u[2];
  T dy_uz = dy_u[2];
  T dz_uz = dz_u[2];
  T omega = dynamics -> getOmega();
  T sToPi = - rho / Lattice<T>::invCs2 / omega;
  T pi[util::TensorVal<Lattice<T> >::n];
  pi[xx] = (T)2 * dx_ux * sToPi;
  pi[yy] = (T)2 * dy_uy * sToPi;
  pi[zz] = (T)2 * dz_uz * sToPi;
  pi[xy] = (dx_uy + dy_ux) * sToPi;
  pi[xz] = (dx_uz + dz_ux) * sToPi;
  pi[yz] = (dy_uz + dz_uy) * sToPi;

  // Computation of the particle distribution functions
  // according to the regularized formula
  T u[Lattice<T>::d];
  cell.computeU(u);
  T uSqr = util::normSqr<T,Lattice<T>::d>(u);

  for (int iPop = 0; iPop < Lattice<T>::q; ++iPop) {
    cell[iPop] = dynamics -> computeEquilibrium(iPop,rho,u,uSqr) +
                 firstOrderLbHelpers<T,Lattice>::fromPiToFneq(iPop, pi);
  }

}

template<typename T, template<typename U> class Lattice, int xNormal, int yNormal, int zNormal>
void OuterVelocityCornerProcessor3D<T, Lattice, xNormal, yNormal, zNormal>::
processSubDomain(BlockLattice3D<T,Lattice>& blockLattice,
                 int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
{
  if (util::contained(x, y, z, x0_, x1_, y0_, y1_, z0_, z1_)) {
    process(blockLattice);
  }
}

////////  OuterVelocityCornerProcessorGenerator3D ///////////////////////////////

template<typename T, template<typename U> class Lattice, int xNormal, int yNormal, int zNormal>
OuterVelocityCornerProcessorGenerator3D<T,Lattice, xNormal,yNormal,zNormal>::
OuterVelocityCornerProcessorGenerator3D(int x_, int y_, int z_)
  : PostProcessorGenerator3D<T,Lattice>(x_,x_, y_,y_, z_,z_)
{ }

template<typename T, template<typename U> class Lattice, int xNormal, int yNormal, int zNormal>
PostProcessor3D<T,Lattice>*
OuterVelocityCornerProcessorGenerator3D<T,Lattice, xNormal,yNormal,zNormal>::
generate() const
{
  return new OuterVelocityCornerProcessor3D<T,Lattice, xNormal,yNormal,zNormal>
         ( this->x0, this->y0, this->z0 );
}

template<typename T, template<typename U> class Lattice, int xNormal, int yNormal, int zNormal>
PostProcessorGenerator3D<T,Lattice>*
OuterVelocityCornerProcessorGenerator3D<T,Lattice, xNormal,yNormal,zNormal>::clone() const
{
  return new OuterVelocityCornerProcessorGenerator3D<T,Lattice, xNormal, yNormal, zNormal>
         (this->x0, this->y0, this->z0);
}


////////  SlipBoundaryProcessor3D ////////////////////////////////

template<typename T, template<typename U> class Lattice>
SlipBoundaryProcessor3D<T,Lattice>::
SlipBoundaryProcessor3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, int discreteNormalX, int discreteNormalY, int discreteNormalZ)
  : x0(x0_), x1(x1_), y0(y0_), y1(y1_), z0(z0_), z1(z1_)
{
  OLB_PRECONDITION(x0==x1 || y0==y1 || z0==z1);
  reflectionPop[0] = 0;
  for (int iPop = 1; iPop < Lattice<T>::q; iPop++) {
    reflectionPop[iPop] = 0;
    // iPop are the directions which pointing into the fluid, discreteNormal is pointing outwarts
    if (Lattice<T>::c[iPop][0]*discreteNormalX + Lattice<T>::c[iPop][1]*discreteNormalY + Lattice<T>::c[iPop][2]*discreteNormalZ < 0) {
      // std::cout << "-----" <<s td::endl;
      int mirrorDirection0;
      int mirrorDirection1;
      int mirrorDirection2;
      int mult = 2 / (discreteNormalX*discreteNormalX + discreteNormalY*discreteNormalY + discreteNormalZ*discreteNormalZ);

      mirrorDirection0 = (Lattice<T>::c[iPop][0] - mult*(Lattice<T>::c[iPop][0]*discreteNormalX + Lattice<T>::c[iPop][1]*discreteNormalY + Lattice<T>::c[iPop][2]*discreteNormalZ)*discreteNormalX);
      mirrorDirection1 = (Lattice<T>::c[iPop][1] - mult*(Lattice<T>::c[iPop][0]*discreteNormalX + Lattice<T>::c[iPop][1]*discreteNormalY + Lattice<T>::c[iPop][2]*discreteNormalZ)*discreteNormalY);
      mirrorDirection2 = (Lattice<T>::c[iPop][2] - mult*(Lattice<T>::c[iPop][0]*discreteNormalX + Lattice<T>::c[iPop][1]*discreteNormalY + Lattice<T>::c[iPop][2]*discreteNormalZ)*discreteNormalZ);

      // bounce back for the case discreteNormalX = discreteNormalY = discreteNormalZ = 1, that is mult=0
      if (mult == 0) {
        mirrorDirection0 = -Lattice<T>::c[iPop][0];
        mirrorDirection1 = -Lattice<T>::c[iPop][1];
        mirrorDirection2 = -Lattice<T>::c[iPop][2];
      }

      // computes mirror jPop
      for (reflectionPop[iPop] = 1; reflectionPop[iPop] < Lattice<T>::q ; reflectionPop[iPop]++) {
        if (Lattice<T>::c[reflectionPop[iPop]][0]==mirrorDirection0 && Lattice<T>::c[reflectionPop[iPop]][1]==mirrorDirection1 && Lattice<T>::c[reflectionPop[iPop]][2]==mirrorDirection2) {
          break;
        }
      }
      //std::cout <<iPop << " to "<< jPop <<" for discreteNormal= "<< discreteNormalX << "/"<<discreteNormalY <<std::endl;
    }
  }
}

template<typename T, template<typename U> class Lattice>
void SlipBoundaryProcessor3D<T,Lattice>::
processSubDomain(BlockLattice3D<T,Lattice>& blockLattice, int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
{
  int newX0, newX1, newY0, newY1, newZ0, newZ1;
  if ( util::intersect (
         x0, x1, y0, y1, z0, z1,
         x0_, x1_, y0_, y1_, z0_, z1_,
         newX0, newX1, newY0, newY1, newZ0, newZ1 ) ) {

    int iX;
#ifdef PARALLEL_MODE_OMP
    #pragma omp parallel for
#endif
    for (iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
        for (int iZ=newZ0; iZ<=newZ1; ++iZ) {
          for (int iPop = 1; iPop < Lattice<T>::q ; ++iPop) {
            if (reflectionPop[iPop]!=0) {
              //do reflection
              blockLattice.get(iX,iY,iZ)[iPop] = blockLattice.get(iX,iY,iZ)[reflectionPop[iPop]];
            }
          }
        }
      }
    }
  }
}

template<typename T, template<typename U> class Lattice>
void SlipBoundaryProcessor3D<T,Lattice>::
process(BlockLattice3D<T,Lattice>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1, z0, z1);
}

////////  SlipBoundaryProcessorGenerator3D ////////////////////////////////

template<typename T, template<typename U> class Lattice>
SlipBoundaryProcessorGenerator3D<T,Lattice>::
SlipBoundaryProcessorGenerator3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, int discreteNormalX_, int discreteNormalY_, int discreteNormalZ_)
  : PostProcessorGenerator3D<T,Lattice>(x0_, x1_, y0_, y1_, z0_, z1_), discreteNormalX(discreteNormalX_), discreteNormalY(discreteNormalY_), discreteNormalZ(discreteNormalZ_)
{ }

template<typename T, template<typename U> class Lattice>
PostProcessor3D<T,Lattice>*
SlipBoundaryProcessorGenerator3D<T,Lattice>::generate() const
{
  return new SlipBoundaryProcessor3D<T,Lattice>
         ( this->x0, this->x1, this->y0, this->y1, this->z0, this->z1, discreteNormalX, discreteNormalY, discreteNormalZ);
}

template<typename T, template<typename U> class Lattice>
PostProcessorGenerator3D<T,Lattice>*
SlipBoundaryProcessorGenerator3D<T,Lattice>::clone() const
{
  return new SlipBoundaryProcessorGenerator3D<T,Lattice>
         (this->x0, this->x1, this->y0, this->y1, this->z0, this->z1, discreteNormalX, discreteNormalY, discreteNormalZ);
}


}  // namespace olb

#endif
