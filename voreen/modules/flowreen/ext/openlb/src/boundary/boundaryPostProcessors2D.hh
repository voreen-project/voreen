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

#ifndef FD_BOUNDARIES_2D_HH
#define FD_BOUNDARIES_2D_HH

#include "boundaryPostProcessors2D.h"
#include "core/finiteDifference2D.h"
#include "core/blockLattice2D.h"
#include "core/util.h"
#include "dynamics/lbHelpers.h"
#include "dynamics/firstOrderLbHelpers.h"

namespace olb {

///////////  StraightFdBoundaryProcessor2D ///////////////////////////////////

template<typename T, template<typename U> class Lattice, int direction, int orientation>
StraightFdBoundaryProcessor2D<T,Lattice,direction,orientation>::
StraightFdBoundaryProcessor2D(int x0_, int x1_, int y0_, int y1_)
  : x0(x0_), x1(x1_), y0(y0_), y1(y1_)
{
  OLB_PRECONDITION(x0==x1 || y0==y1);
}

template<typename T, template<typename U> class Lattice, int direction,int orientation>
void StraightFdBoundaryProcessor2D<T,Lattice,direction,orientation>::
processSubDomain(BlockLattice2D<T,Lattice>& blockLattice, int x0_, int x1_, int y0_, int y1_)
{
  using namespace olb::util::tensorIndices2D;

  int newX0, newX1, newY0, newY1;
  if ( util::intersect (
         x0, x1, y0, y1,
         x0_, x1_, y0_, y1_,
         newX0, newX1, newY0, newY1 ) ) {

    int iX;

#ifdef PARALLEL_MODE_OMP
    #pragma omp parallel for
#endif
    for (iX=newX0; iX<=newX1; ++iX) {
      T dx_u[Lattice<T>::d], dy_u[Lattice<T>::d];
      for (int iY=newY0; iY<=newY1; ++iY) {
        Cell<T,Lattice>& cell = blockLattice.get(iX,iY);
        Dynamics<T,Lattice>* dynamics = blockLattice.getDynamics(iX, iY);

        T rho, u[Lattice<T>::d];
        cell.computeRhoU(rho,u);

        interpolateGradients<0>(blockLattice, dx_u, iX, iY);
        interpolateGradients<1>(blockLattice, dy_u, iX, iY);
        T dx_ux = dx_u[0];
        T dy_ux = dy_u[0];
        T dx_uy = dx_u[1];
        T dy_uy = dy_u[1];
        T omega = dynamics->getOmega();
        T sToPi = - rho / Lattice<T>::invCs2 / omega;
        T pi[util::TensorVal<Lattice<T> >::n];
        pi[xx] = (T)2 * dx_ux * sToPi;
        pi[yy] = (T)2 * dy_uy * sToPi;
        pi[xy] = (dx_uy + dy_ux) * sToPi;

        // Computation of the particle distribution functions
        // according to the regularized formula

        T uSqr = util::normSqr<T,2>(u);
        for (int iPop = 0; iPop < Lattice<T>::q; ++iPop) {
          cell[iPop] = dynamics -> computeEquilibrium(iPop,rho,u,uSqr) +
                       firstOrderLbHelpers<T,Lattice>::fromPiToFneq(iPop, pi);
        }
      }
    }
  }
}

template<typename T, template<typename U> class Lattice, int direction,int orientation>
void StraightFdBoundaryProcessor2D<T,Lattice,direction,orientation>::
process(BlockLattice2D<T,Lattice>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1);
}

template<typename T, template<typename U> class Lattice, int direction,int orientation>
template<int deriveDirection>
void StraightFdBoundaryProcessor2D<T,Lattice,direction,orientation>::
interpolateGradients(BlockLattice2D<T,Lattice> const& blockLattice,
                     T velDeriv[Lattice<T>::d], int iX, int iY) const
{
  fd::DirectedGradients2D<T,Lattice,direction,orientation,direction==deriveDirection>::
  interpolateVector(velDeriv, blockLattice, iX, iY);
}

////////  StraightFdBoundaryProcessorGenerator2D ////////////////////////////////

template<typename T, template<typename U> class Lattice, int direction,int orientation>
StraightFdBoundaryProcessorGenerator2D<T,Lattice, direction,orientation>::
StraightFdBoundaryProcessorGenerator2D(int x0_, int x1_, int y0_, int y1_)
  : PostProcessorGenerator2D<T,Lattice>(x0_, x1_, y0_, y1_)
{ }

template<typename T, template<typename U> class Lattice, int direction,int orientation>
PostProcessor2D<T,Lattice>*
StraightFdBoundaryProcessorGenerator2D<T,Lattice,direction,orientation>::generate() const
{
  return new StraightFdBoundaryProcessor2D<T,Lattice,direction,orientation>
         ( this->x0, this->x1, this->y0, this->y1);
}

template<typename T, template<typename U> class Lattice, int direction,int orientation>
PostProcessorGenerator2D<T,Lattice>*
StraightFdBoundaryProcessorGenerator2D<T,Lattice,direction,orientation>::clone() const
{
  return new StraightFdBoundaryProcessorGenerator2D<T,Lattice,direction,orientation>
         (this->x0, this->x1, this->y0, this->y1);
}


////////  StraightConvectionBoundaryProcessor2D ////////////////////////////////

template<typename T, template<typename U> class Lattice, int direction, int orientation>
StraightConvectionBoundaryProcessor2D<T,Lattice,direction,orientation>::
StraightConvectionBoundaryProcessor2D(int x0_, int x1_, int y0_, int y1_, T* uAv_)
  : x0(x0_), x1(x1_), y0(y0_), y1(y1_), uAv(uAv_)
{
  OLB_PRECONDITION(x0==x1 || y0==y1);

  saveCell = new T** [(size_t)(x1_-x0_+1)];
  for (int iX=0; iX<=x1_-x0_; ++iX) {
    saveCell[iX] = new T* [(size_t)(y1_-y0_+1)];
    for (int iY=0; iY<=y1_-y0_; ++iY) {
      saveCell[iX][iY] = new T [(size_t)(Lattice<T>::q)];
      for (int iPop=0; iPop<Lattice<T>::q; ++iPop) {
        saveCell[iX][iY][iPop] = T();
      }
    }
  }
}

template<typename T, template<typename U> class Lattice, int direction, int orientation>
StraightConvectionBoundaryProcessor2D<T,Lattice,direction,orientation>::
~StraightConvectionBoundaryProcessor2D()
{
  for (int iX=0; iX<=x1-x0; ++iX) {
    for (int iY=0; iY<=y1-y0; ++iY) {
      delete [] saveCell[iX][iY];
    }
    delete [] saveCell[iX];
  }
  delete [] saveCell;
}

template<typename T, template<typename U> class Lattice, int direction,int orientation>
void StraightConvectionBoundaryProcessor2D<T,Lattice,direction,orientation>::
processSubDomain(BlockLattice2D<T,Lattice>& blockLattice, int x0_, int x1_, int y0_, int y1_)
{
  using namespace olb::util::tensorIndices2D;

  int newX0, newX1, newY0, newY1;
  if ( util::intersect (
         x0, x1, y0, y1,
         x0_, x1_, y0_, y1_,
         newX0, newX1, newY0, newY1 ) ) {

    int iX;
#ifdef PARALLEL_MODE_OMP
    #pragma omp parallel for
#endif
    for (iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
        Cell<T,Lattice>& cell = blockLattice.get(iX,iY);
        for (int iPop = 0; iPop < Lattice<T>::q ; ++iPop) {
          if (Lattice<T>::c[iPop][direction]==-orientation) {
            cell[iPop] = saveCell[iX-newX0][iY-newY0][iPop];
          }
        }

        T rho0, u0[2];
        T rho1, u1[2];
        T rho2, u2[2];
        if (direction==0) {
          blockLattice.get(iX,iY).computeRhoU(rho0,u0);
          blockLattice.get(iX-orientation,iY).computeRhoU(rho1,u1);
          blockLattice.get(iX-orientation*2,iY).computeRhoU(rho2,u2);
        } else {
          blockLattice.get(iX,iY).computeRhoU(rho0,u0);
          blockLattice.get(iX,iY-orientation).computeRhoU(rho1,u1);
          blockLattice.get(iX,iY-orientation*2).computeRhoU(rho2,u2);
        }

        rho0 = T(1);
        rho1 = T(1);
        rho2 = T(1);
        T uDelta[2];
        T uAverage = rho0*u0[direction];
        if (uAv!=nullptr) {
          uAverage = *uAv;
        }
        uDelta[0]=-uAverage*0.5*(3*rho0*u0[0]-4*rho1*u1[0]+rho2*u2[0]);
        uDelta[1]=-uAverage*0.5*(3*rho0*u0[1]-4*rho1*u1[1]+rho2*u2[1]);

        for (int iPop = 0; iPop < Lattice<T>::q ; ++iPop) {
          if (Lattice<T>::c[iPop][direction] == -orientation) {
            saveCell[iX-newX0][iY-newY0][iPop] = cell[iPop] + Lattice<T>::invCs2*Lattice<T>::t[iPop]*(uDelta[0]*Lattice<T>::c[iPop][0]+uDelta[1]*Lattice<T>::c[iPop][1]);
          }
        }
      }
    }
  }
}

template<typename T, template<typename U> class Lattice, int direction,int orientation>
void StraightConvectionBoundaryProcessor2D<T,Lattice,direction,orientation>::
process(BlockLattice2D<T,Lattice>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1);
}


////////  StraightConvectionBoundaryProcessorGenerator2D ////////////////////////////////

template<typename T, template<typename U> class Lattice, int direction,int orientation>
StraightConvectionBoundaryProcessorGenerator2D<T,Lattice, direction,orientation>::
StraightConvectionBoundaryProcessorGenerator2D(int x0_, int x1_, int y0_, int y1_, T* uAv_)
  : PostProcessorGenerator2D<T,Lattice>(x0_, x1_, y0_, y1_), uAv(uAv_)
{ }

template<typename T, template<typename U> class Lattice, int direction,int orientation>
PostProcessor2D<T,Lattice>*
StraightConvectionBoundaryProcessorGenerator2D<T,Lattice,direction,orientation>::generate() const
{
  return new StraightConvectionBoundaryProcessor2D<T,Lattice,direction,orientation>
         ( this->x0, this->x1, this->y0, this->y1, uAv);
}

template<typename T, template<typename U> class Lattice, int direction,int orientation>
PostProcessorGenerator2D<T,Lattice>*
StraightConvectionBoundaryProcessorGenerator2D<T,Lattice,direction,orientation>::clone() const
{
  return new StraightConvectionBoundaryProcessorGenerator2D<T,Lattice,direction,orientation>
         (this->x0, this->x1, this->y0, this->y1, uAv);
}


////////  SlipBoundaryProcessor2D ////////////////////////////////

template<typename T, template<typename U> class Lattice>
SlipBoundaryProcessor2D<T,Lattice>::
SlipBoundaryProcessor2D(int x0_, int x1_, int y0_, int y1_, int discreteNormalX, int discreteNormalY)
  : x0(x0_), x1(x1_), y0(y0_), y1(y1_)
{
  OLB_PRECONDITION(x0==x1 || y0==y1);
  reflectionPop[0] = 0;
  for (int iPop = 1; iPop < Lattice<T>::q; iPop++) {
    reflectionPop[iPop] = 0;
    // iPop are the directions which ointing into the fluid, discreteNormal is pointing outwarts
    if (Lattice<T>::c[iPop][0]*discreteNormalX + Lattice<T>::c[iPop][1]*discreteNormalY < 0) {
      // std::cout << "-----" <<s td::endl;
      int mirrorDirection0;
      int mirrorDirection1;
      int mult = 1;
      if (discreteNormalX*discreteNormalY==0) {
        mult = 2;
      }
      mirrorDirection0 = (Lattice<T>::c[iPop][0] - mult*(Lattice<T>::c[iPop][0]*discreteNormalX + Lattice<T>::c[iPop][1]*discreteNormalY )*discreteNormalX);
      mirrorDirection1 = (Lattice<T>::c[iPop][1] - mult*(Lattice<T>::c[iPop][0]*discreteNormalX + Lattice<T>::c[iPop][1]*discreteNormalY )*discreteNormalY);

      // computes mirror jPop
      for (reflectionPop[iPop] = 1; reflectionPop[iPop] < Lattice<T>::q ; reflectionPop[iPop]++) {
        if (Lattice<T>::c[reflectionPop[iPop]][0]==mirrorDirection0 && Lattice<T>::c[reflectionPop[iPop]][1]==mirrorDirection1 ) {
          break;
        }
      }
      //std::cout <<iPop << " to "<< jPop <<" for discreteNormal= "<< discreteNormalX << "/"<<discreteNormalY <<std::endl;
    }
  }
}

template<typename T, template<typename U> class Lattice>
void SlipBoundaryProcessor2D<T,Lattice>::
processSubDomain(BlockLattice2D<T,Lattice>& blockLattice, int x0_, int x1_, int y0_, int y1_)
{
  int newX0, newX1, newY0, newY1;
  if ( util::intersect (
         x0, x1, y0, y1,
         x0_, x1_, y0_, y1_,
         newX0, newX1, newY0, newY1 ) ) {

    int iX;
#ifdef PARALLEL_MODE_OMP
    #pragma omp parallel for
#endif
    for (iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
        for (int iPop = 1; iPop < Lattice<T>::q ; ++iPop) {
          if (reflectionPop[iPop]!=0) {
            //do reflection
            blockLattice.get(iX,iY)[iPop] = blockLattice.get(iX,iY)[reflectionPop[iPop]];
          }
        }
      }
    }
  }
}

template<typename T, template<typename U> class Lattice>
void SlipBoundaryProcessor2D<T,Lattice>::
process(BlockLattice2D<T,Lattice>& blockLattice)
{
  processSubDomain(blockLattice, x0, x1, y0, y1);
}

////////  SlipBoundaryProcessorGenerator2D ////////////////////////////////

template<typename T, template<typename U> class Lattice>
SlipBoundaryProcessorGenerator2D<T,Lattice>::
SlipBoundaryProcessorGenerator2D(int x0_, int x1_, int y0_, int y1_, int discreteNormalX_, int discreteNormalY_)
  : PostProcessorGenerator2D<T,Lattice>(x0_, x1_, y0_, y1_), discreteNormalX(discreteNormalX_), discreteNormalY(discreteNormalY_)
{ }

template<typename T, template<typename U> class Lattice>
PostProcessor2D<T,Lattice>*
SlipBoundaryProcessorGenerator2D<T,Lattice>::generate() const
{
  return new SlipBoundaryProcessor2D<T,Lattice>
         ( this->x0, this->x1, this->y0, this->y1, discreteNormalX, discreteNormalY);
}

template<typename T, template<typename U> class Lattice>
PostProcessorGenerator2D<T,Lattice>*
SlipBoundaryProcessorGenerator2D<T,Lattice>::clone() const
{
  return new SlipBoundaryProcessorGenerator2D<T,Lattice>
         (this->x0, this->x1, this->y0, this->y1, discreteNormalX, discreteNormalY);
}

/////////// OuterVelocityCornerProcessor2D /////////////////////////////////////

template<typename T, template<typename U> class Lattice, int xNormal,int yNormal>
OuterVelocityCornerProcessor2D<T, Lattice, xNormal, yNormal>::
OuterVelocityCornerProcessor2D(int x_, int y_)
  : x(x_), y(y_)
{ }

template<typename T, template<typename U> class Lattice, int xNormal,int yNormal>
void OuterVelocityCornerProcessor2D<T, Lattice, xNormal, yNormal>::
process(BlockLattice2D<T,Lattice>& blockLattice)
{
  using namespace olb::util::tensorIndices2D;

  T rho10 = blockLattice.get(x-1*xNormal, y-0*yNormal).computeRho();
  T rho01 = blockLattice.get(x-0*xNormal, y-1*yNormal).computeRho();

  T rho20 = blockLattice.get(x-2*xNormal, y-0*yNormal).computeRho();
  T rho02 = blockLattice.get(x-0*xNormal, y-2*yNormal).computeRho();

  T rho = (T)2/(T)3*(rho01+rho10) - (T)1/(T)6*(rho02+rho20);

  T dx_u[Lattice<T>::d], dy_u[Lattice<T>::d];
  fd::DirectedGradients2D<T, Lattice, 0, xNormal, true>::interpolateVector(dx_u, blockLattice, x,y);
  fd::DirectedGradients2D<T, Lattice, 1, yNormal, true>::interpolateVector(dy_u, blockLattice, x,y);
  T dx_ux = dx_u[0];
  T dy_ux = dy_u[0];
  T dx_uy = dx_u[1];
  T dy_uy = dy_u[1];

  Cell<T,Lattice>& cell = blockLattice.get(x,y);
  Dynamics<T,Lattice>* dynamics = blockLattice.getDynamics(x, y);
  T omega = dynamics -> getOmega();

  T sToPi = - rho / Lattice<T>::invCs2 / omega;
  T pi[util::TensorVal<Lattice<T> >::n];
  pi[xx] = (T)2 * dx_ux * sToPi;
  pi[yy] = (T)2 * dy_uy * sToPi;
  pi[xy] = (dx_uy + dy_ux) * sToPi;

  // Computation of the particle distribution functions
  // according to the regularized formula
  T u[Lattice<T>::d];
  blockLattice.get(x,y).computeU(u);

  T uSqr = util::normSqr<T,2>(u);
  for (int iPop = 0; iPop < Lattice<T>::q; ++iPop) {
    cell[iPop] =
      dynamics -> computeEquilibrium(iPop,rho,u,uSqr) +
      firstOrderLbHelpers<T,Lattice>::fromPiToFneq(iPop, pi);
  }
}

template<typename T, template<typename U> class Lattice, int xNormal,int yNormal>
void OuterVelocityCornerProcessor2D<T, Lattice, xNormal, yNormal>::
processSubDomain(BlockLattice2D<T,Lattice>& blockLattice,
                 int x0_, int x1_, int y0_, int y1_ )
{
  if (util::contained(x, y, x0_, x1_, y0_, y1_)) {
    process(blockLattice);
  }
}


////////  OuterVelocityCornerProcessorGenerator2D ////////////////////////////

template<typename T, template<typename U> class Lattice, int xNormal,int yNormal>
OuterVelocityCornerProcessorGenerator2D<T, Lattice, xNormal, yNormal>::
OuterVelocityCornerProcessorGenerator2D(int x_, int y_)
  : PostProcessorGenerator2D<T,Lattice>(x_, x_, y_, y_)
{ }

template<typename T, template<typename U> class Lattice, int xNormal,int yNormal>
PostProcessor2D<T,Lattice>*
OuterVelocityCornerProcessorGenerator2D<T, Lattice, xNormal, yNormal>::generate() const
{
  return new OuterVelocityCornerProcessor2D<T, Lattice, xNormal, yNormal>
         ( this->x0, this->y0);
}

template<typename T, template<typename U> class Lattice, int xNormal,int yNormal>
PostProcessorGenerator2D<T,Lattice>*
OuterVelocityCornerProcessorGenerator2D<T, Lattice, xNormal, yNormal>::
clone() const
{
  return new OuterVelocityCornerProcessorGenerator2D<T, Lattice, xNormal, yNormal>
         ( this->x0, this->y0);
}


}  // namespace olb

#endif
