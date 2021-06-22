/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006-2015 Jonas Latt, Mathias J. Krause
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
 * Template specializations for some computationally intensive LB
 * functions of the header file lbHelpers.h, for some D3Q19 grids.
 */

#ifndef LB_HELPERS_D3Q19_H
#define LB_HELPERS_D3Q19_H

namespace olb {

// Efficient specialization for D3Q19 lattice
template<typename T, typename... FIELDS>
struct lbDynamicsHelpers<T, descriptors::D3Q19<FIELDS...> > {
  using DESCRIPTOR = descriptors::D3Q19<FIELDS...>;

  static T equilibrium( int iPop, T rho, const T u[3], const T uSqr )
  {
    T c_u = descriptors::c<DESCRIPTOR>(iPop,0)*u[0] + descriptors::c<DESCRIPTOR>(iPop,1)*u[1] + descriptors::c<DESCRIPTOR>(iPop,2)*u[2];
    return rho * descriptors::t<T,DESCRIPTOR>(iPop) * ( 1. + 3.*c_u + 4.5*c_u*c_u - 1.5*uSqr ) - descriptors::t<T,DESCRIPTOR>(iPop);
  }

  static T equilibriumFirstOrder( int iPop, T rho, const T u[3] )
    {
      T c_u = descriptors::c<DESCRIPTOR>(iPop,0) * u[0] + descriptors::c<DESCRIPTOR>(iPop,1) * u[1] + descriptors::c<DESCRIPTOR>(iPop,2) * u[2];

      return rho * descriptors::t<T,DESCRIPTOR>(iPop) * ( ( T )1 + c_u * descriptors::invCs2<T,DESCRIPTOR>() ) - descriptors::t<T,DESCRIPTOR>(iPop);
    }

  static T incEquilibrium(int iPop, const T j[3], const T jSqr, const T pressure)
  {
    T c_j = descriptors::c<DESCRIPTOR>(iPop,0)*j[0] + descriptors::c<DESCRIPTOR>(iPop,1)*j[1] + descriptors::c<DESCRIPTOR>(iPop,2)*j[2];
    return descriptors::t<T,DESCRIPTOR>(iPop) * ( 3.*pressure + 3.*c_j + 4.5*c_j*c_j - 1.5*jSqr ) - descriptors::t<T,DESCRIPTOR>(iPop);
  }

  static void computeFneq (
    ConstCell<T,DESCRIPTOR>& cell, T fNeq[19], T rho, const T u[3] )
  {
    const T uSqr = u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
    for (int iPop=0; iPop < 19; ++iPop) {
      fNeq[iPop] = cell[iPop] - equilibrium(iPop, rho, u, uSqr);
    }
  }

 /// RLB advection diffusion collision step
  static T rlbCollision(Cell<T, descriptors::D3Q19<> >& cell, T rho, const T u[3], T omega )
  {
    const T uSqr = util::normSqr<T, DESCRIPTOR::d>( u );
    // First-order moment for the regularization
    T j1[DESCRIPTOR::d];
    for ( int iD = 0; iD < 3; ++iD ) {
      j1[iD] = T();
    }

    T fEq[DESCRIPTOR::q];
    for ( int iPop = 0; iPop < 19; ++iPop ) {
      fEq[iPop] = lbDynamicsHelpers<T, DESCRIPTOR>::equilibriumFirstOrder( iPop, rho, u );
      for ( int iD = 0; iD < DESCRIPTOR::d; ++iD ) {
        j1[iD] += descriptors::c<DESCRIPTOR>(iPop,iD) * ( cell[iPop] - fEq[iPop] );
      }
    }

    // Collision step
    for ( int iPop = 0; iPop < DESCRIPTOR::q; ++iPop ) {
      T fNeq = T();
      for ( int iD = 0; iD < DESCRIPTOR::d; ++iD ) {
        fNeq += descriptors::c<DESCRIPTOR>(iPop,iD) * j1[iD];
      }
      fNeq *= descriptors::t<T,DESCRIPTOR>(iPop) * descriptors::invCs2<T,DESCRIPTOR>();
      cell[iPop] = fEq[iPop] + ( (T)1 - omega ) * fNeq;
    }
    return uSqr;
  }

  static T bgkCollision(Cell<T,DESCRIPTOR>& cell, T const& rho, const T u[3], T const& omega)
  {
    T one_m_omega = (T)1 - omega;
    T t0_omega = descriptors::t<T,DESCRIPTOR>(0)*omega; // weight for i=0
    T t1_omega = descriptors::t<T,DESCRIPTOR>(1)*omega; // weight for i=1,2,3,10,11,12
    T t4_omega = descriptors::t<T,DESCRIPTOR>(4)*omega; // weight for i=4,5,6,7,8,9,13,14,15,16,17,18

    T uSqr     = u[0]*u[0] + u[1]*u[1] + u[2]*u[2]; // compute of usqr

    T u3x      = (T)3*u[0]; // compute of 3*ux
    T u3y      = (T)3*u[1]; // compute of 3*uy
    T u3z      = (T)3*u[2]; // compute of 3*uz

    T u3xSqr_  = .5*u3x*u3x; // compute 9/2*(ux)²
    T u3ySqr_  = .5*u3y*u3y; // compute 9/2*(uy)²
    T u3zSqr_  = .5*u3z*u3z;  // compute 9/2*(uz)²

    T u3xu3y_  = u3x*u3y; // compute 9(ux*uy)
    T u3xu3z_  = u3x*u3z; // compute 9(ux*uz)
    T u3yu3z_  = u3y*u3z; // compute 9(uy*uz)

    T C1 = (T)1 + (T)3*uSqr; // compute 1+3*((ux)²+(uy)²+(uz)²)
    T C2, C3;

    //**************************case i=0
    C3 = -u3xSqr_ - u3ySqr_ - u3zSqr_; // compute -9/2*((ux)²+(uy)²+(uz)²)

    // compute of feq0=rho*t0*(c1+c3)=rho*t0*(1-3/2*((ux)²+(uy)²+(uz)²))
    cell[0] *= one_m_omega;
    cell[0] += t0_omega*(rho*(C1 + C3) - (T)1);


    //**************************case i=1 and i=10
    C2 = -u3x; // compute -3*ux
    C3 = -u3ySqr_ - u3zSqr_; // compute -9/2*((uy)²+(uz)²)

    // compute of feq1=rho*t1*(c1+c2+c3)=rho*t0*(1-3*ux+9/2*(ux)²-3/2*((ux)²+(uy)²+(uz)²))
    cell[1]  *= one_m_omega;
    cell[1]  += t1_omega*(rho*(C1 + C2 + C3) - (T)1);

    // compute of feq10=rho*t1*(c1-c2+c3)=rho*t0*(1+3*ux+9/2*(ux)²-3/2*((ux)²+(uy)²+(uz)²))
    cell[10] *= one_m_omega;
    cell[10] += t1_omega*(rho*(C1 - C2 + C3) - (T)1);

    //************************case i=2 and i=11
    C2 = -u3y; // compute -3*uy
    C3 = -u3xSqr_ - u3zSqr_; // compute -9/2*((ux)²+(uz)²)

    // compute of feq2=rho*t1*(c1-c2+c3)=rho*t1*(1-3*uy+9/2*(uy)²-3/2*((ux)²+(uy)²+(uz)²))
    cell[2]  *= one_m_omega;
    cell[2]  += t1_omega*(rho*(C1 + C2 + C3) - (T)1);

    // compute of feq2=rho*t1*(c1+c2+c3)=rho*t1*(1+3*uy+9/2*(uy)²-3/2*((ux)²+(uy)²+(uz)²))
    cell[11] *= one_m_omega;
    cell[11] += t1_omega*(rho*(C1 - C2 + C3) - (T)1);

    //************************case i=3 and i=12
    C2 = -u3z; // compute -3*uz
    C3 = -u3xSqr_ - u3ySqr_; // compute -9/2*((ux)²+(uy)²)

    // compute of feq3=rho*t1*(c1+c2+c3)=rho*t1*(1-3*uz+9/2*(uz)²-3/2*((ux)²+(uy)²+(uz)²))
    cell[3]  *= one_m_omega;
    cell[3]  += t1_omega*(rho*(C1 + C2 + C3) - (T)1);

    // compute of feq12=rho*t1*(c1-c2+c3)=rho*t1*(1+3*uz+9/2*(uz)²-3/2*((ux)²+(uy)²+(uz)²))
    cell[12] *= one_m_omega;
    cell[12] += t1_omega*(rho*(C1 - C2 + C3) - (T)1);

    //************************case i=4 and i=13
    C2 = -u3x - u3y; // compute -3*(uz+uy)
    C3 = u3xu3y_ - u3zSqr_; // compute 9*(ux*uy)-9/2*(uz)²

    // compute of feq4=rho*t4*(c1+c2+c3)=rho*t4*(1-3*(ux+uy)+9/2*((ux)²+(uy)²+2*ux*uy)-3/2*((ux)²+(uy)²+(uz)²))
    cell[4]  *= one_m_omega;
    cell[4]  += t4_omega*(rho*(C1 + C2 + C3) - (T)1);

    // compute of feq13=rho*t4*(c1-c2+c3)=rho*t4*(1+3*(ux+uy)+9/2*((ux)²+(uy)²+2*ux*uy)-3/2*((ux)²+(uy)²+(uz)²))
    cell[13] *= one_m_omega;
    cell[13] += t4_omega*(rho*(C1 - C2 + C3) - (T)1);

    //************************case i=5 and i=14
    C2 = -u3x + u3y; // compute -3*(ux+uy)
    C3 = -u3xu3y_ - u3zSqr_; // compute -9*(ux*uy)-9/2*(uz)²

    // compute of feq5=rho*t4*(c1+c2+c3)=rho*t4*(1-3*(ux-uy)+9/2*((ux)²+(uy)²-2*ux*uy)-3/2*((ux)²+(uy)²+(uz)²))
    cell[5]  *= one_m_omega;
    cell[5]  += t4_omega*(rho*(C1 + C2 + C3) - (T)1);

    // compute of feq14=rho*t4*(c1-c2+c3)=rho*t4*(1+3*(ux-uy)+9/2*((ux)²+(uy)²-2*ux*uy)-3/2*((ux)²+(uy)²+(uz)²))
    cell[14] *= one_m_omega;
    cell[14] +=t4_omega*(rho*(C1 - C2 + C3) - (T)1);

    //************************case i=6 and i=15
    C2 = -u3x - u3z; // compute -3*(ux+uz)
    C3 = u3xu3z_ - u3ySqr_; // compute 9*(ux*uz)-9/2*(uy)²

    // compute of feq6=rho*t4*(c1+c2+c3)=rho*t4*(1-3*(ux+uz)+9/2*((ux)²+(uy)²+2*ux*uz)-3/2*((ux)²+(uy)²+(uz)²))
    cell[6]  *= one_m_omega;
    cell[6]  += t4_omega*(rho*(C1 + C2 + C3) - (T)1);

    // compute of feq15=rho*t4*(c1-c2+c3)=rho*t15*(1+3*(ux+uz)+9/2*((ux)²+(uy)²+2*ux*uz)-3/2*((ux)²+(uy)²+(uz)²))
    cell[15] *= one_m_omega;
    cell[15] +=t4_omega*(rho*(C1 - C2 + C3) - (T)1);

    //************************case i=7 and i=16
    C2 = -u3x + u3z; // compute -3*(ux-uz)
    C3 = -u3xu3z_ - u3ySqr_; // compute -9*(ux*uz)-9/2*(uy)²

    // compute of feq7=rho*t4*(c1+c2+c3)=rho*t4*(1-3*(ux-uz)+9/2*((ux)²+(uz)²-2*ux*uz)-3/2*((ux)²+(uy)²+(uz)²))
    cell[7]  *= one_m_omega;
    cell[7]  += t4_omega*(rho*(C1 + C2 + C3) - (T)1);

    // compute of feq16=rho*t4*(c1-c2+c3)=rho*t4*(1+3*(ux-uz)+9/2*((ux)²+(uz)²-2*ux*uz)-3/2*((ux)²+(uy)²+(uz)²))
    cell[16] *= one_m_omega;
    cell[16] += t4_omega*(rho*(C1 - C2 + C3) - (T)1);

    //************************case i=7 and i=16
    C2 = -u3y - u3z; // compute -3*(uy+uz)
    C3 = u3yu3z_ - u3xSqr_; // compute 9*(uy*uz)-9/2*(ux)²

    // compute of feq8=rho*t4*(c1+c2+c3)=rho*t4*(1-3*(uy+uz)+9/2*((uz)²+(uy)²+2*uy*uz)-3/2*((ux)²+(uy)²+(uz)²))
    cell[8]  *= one_m_omega;
    cell[8]  += t4_omega*(rho*(C1 + C2 + C3) - (T)1);

    // compute of feq17=rho*t4*(c1-c2+c3)=rho*t4*(1+3*(uy+uz)+9/2*((uz)²+(uy)²+2*uy*uz)-3/2*((ux)²+(uy)²+(uz)²))
    cell[17] *= one_m_omega;
    cell[17] += t4_omega*(rho*(C1 - C2 + C3) - (T)1);

    //************************case i=7 and i=16
    C2 = -u3y + u3z; // compute -3*(uy-uz)
    C3 = -u3yu3z_ - u3xSqr_; // compute -9*(uy*uz)-9/2*(ux)²

    // compute of feq9=rho*t4*(c1+c2+c3)=rho*t4*(1-3*(uy-uz)+9/2*((uy)²+(uz)²-2*uy*uz)-3/2*((ux)²+(uy)²+(uz)²))
    cell[9]  *= one_m_omega;
    cell[9]  += t4_omega*(rho*(C1 + C2 + C3) - (T)1);

    // compute of feq18=rho*t4*(c1+c2+c3)=rho*t4*(1+3*(uy-uz)+9/2*((uy)²+(uz)²-2*uy*uz)-3/2*((ux)²+(uy)²+(uz)²))
    cell[18] *= one_m_omega;
    cell[18] += t4_omega*(rho*(C1 - C2 + C3) - (T)1);

    return uSqr;
  }

  static T incBgkCollision(Cell<T,DESCRIPTOR>& cell, T pressure, const T j[3], T omega)
  {
    const T jSqr = util::normSqr<T,descriptors::D3Q19<>::d>(j);
    for (int iPop=0; iPop < descriptors::D3Q19<>::q; ++iPop) {
      cell[iPop] *= (T)1-omega;
      cell[iPop] += omega * lbDynamicsHelpers<T,descriptors::D3Q19<> >
                    ::incEquilibrium(iPop, j, jSqr, pressure);
    }
    return jSqr;
  }

  static T constRhoBgkCollision(Cell<T,DESCRIPTOR>& cell, T rho, const T u[3], T ratioRho, T omega)
  {
    const T uSqr = util::normSqr<T,descriptors::D3Q19<>::d>(u);
    for (int iPop=0; iPop < descriptors::D3Q19<>::q; ++iPop) {
      T feq = lbDynamicsHelpers<T,descriptors::D3Q19<> >::
              equilibrium(iPop, rho, u, uSqr );
      cell[iPop] =
        ratioRho*(feq+descriptors::t<T,DESCRIPTOR>(iPop))
        -descriptors::t<T,DESCRIPTOR>(iPop) +
        ((T)1-omega)*(cell[iPop]-feq);
    }
    return uSqr;
  }

  static void partial_rho ( ConstCell<T,DESCRIPTOR>& cell,
                            T& surfX_M1, T& surfX_0, T& surfX_P1,
                            T& surfY_M1, T& surfY_P1, T& surfZ_M1, T& surfZ_P1 )
  {
    surfX_M1 = cell[1] + cell[4] + cell[5] + cell[6] + cell[7];
    surfX_0  = cell[0] + cell[2] + cell[3] + cell[8] +
               cell[9] + cell[11] + cell[12] + cell[17] + cell[18];
    surfX_P1 = cell[10] + cell[13] + cell[14] + cell[15] + cell[16];

    surfY_M1 = cell[2] + cell[4] + cell[8] + cell[9] + cell[14];
    surfY_P1 = cell[5] + cell[11] + cell[13] + cell[17] + cell[18];

    surfZ_M1 = cell[3] + cell[6] + cell[8] + cell[16] + cell[18];
    surfZ_P1 = cell[7] + cell[9] + cell[12] + cell[15] + cell[17];

  }

  static void computeRhoU(ConstCell<T,DESCRIPTOR>& cell, T& rho, T u[3])
  {
    T surfX_M1, surfX_0, surfX_P1, surfY_M1, surfY_P1, surfZ_M1, surfZ_P1;

    surfX_M1 = cell[1] + cell[4] + cell[5] + cell[6] + cell[7];
    surfX_0  = cell[0] + cell[2] + cell[3] + cell[8] +
               cell[9] + cell[11] + cell[12] + cell[17] + cell[18];
    surfX_P1 = cell[10] + cell[13] + cell[14] + cell[15] + cell[16];

    surfY_M1 = cell[2] + cell[4] + cell[8] + cell[9] + cell[14];
    surfY_P1 = cell[5] + cell[11] + cell[13] + cell[17] + cell[18];

    surfZ_M1 = cell[3] + cell[6] + cell[8] + cell[16] + cell[18];
    surfZ_P1 = cell[7] + cell[9] + cell[12] + cell[15] + cell[17];

    rho = surfX_M1 + surfX_0 + surfX_P1 + (T)1;
    T invRho= 1./rho;

    u[0]  = ( surfX_P1 - surfX_M1 )*invRho;
    u[1]  = ( surfY_P1 - surfY_M1 )*invRho;
    u[2]  = ( surfZ_P1 - surfZ_M1 )*invRho;
  }

  static void computeRhoJ(ConstCell<T,DESCRIPTOR>& cell, T& rho, T j[3])
  {
    T surfX_M1, surfX_0, surfX_P1,
    surfY_M1, surfY_P1, surfZ_M1, surfZ_P1;

    partial_rho(cell, surfX_M1, surfX_0, surfX_P1,
                surfY_M1, surfY_P1, surfZ_M1, surfZ_P1);

    rho = surfX_M1 + surfX_0 + surfX_P1 + (T)1;

    j[0]  = ( surfX_P1 - surfX_M1 );
    j[1]  = ( surfY_P1 - surfY_M1 );
    j[2]  = ( surfZ_P1 - surfZ_M1 );
  }

  static void computeJ(ConstCell<T,DESCRIPTOR>& cell, T j[3])
  {
    T surfX_M1, surfX_P1, surfY_M1, surfY_P1, surfZ_M1, surfZ_P1;

    surfX_M1 = cell[1] + cell[4] + cell[5] + cell[6] + cell[7];
    surfX_P1 = cell[10] + cell[13] + cell[14] + cell[15] + cell[16];

    surfY_M1 = cell[2] + cell[4] + cell[8] + cell[9] + cell[14];
    surfY_P1 = cell[5] + cell[11] + cell[13] + cell[17] + cell[18];

    surfZ_M1 = cell[3] + cell[6] + cell[8] + cell[16] + cell[18];
    surfZ_P1 = cell[7] + cell[9] + cell[12] + cell[15] + cell[17];

    j[0]  = ( surfX_P1 - surfX_M1 );
    j[1]  = ( surfY_P1 - surfY_M1 );
    j[2]  = ( surfZ_P1 - surfZ_M1 );
  }

  static void computeStress(ConstCell<T,DESCRIPTOR>& cell, T rho, const T u[3], T pi[6])
  {
    // Workaround for Intel(r) compiler 9.1;
    // "using namespace util::tensorIndices3D" is not sufficient
    using util::tensorIndices3D::xx;
    using util::tensorIndices3D::yy;
    using util::tensorIndices3D::zz;
    using util::tensorIndices3D::xy;
    using util::tensorIndices3D::xz;
    using util::tensorIndices3D::yz;

    T surfX_M1, surfX_0, surfX_P1,
    surfY_M1, surfY_P1, surfZ_M1, surfZ_P1;

    partial_rho(cell, surfX_M1, surfX_0, surfX_P1,
                surfY_M1, surfY_P1, surfZ_M1, surfZ_P1);

    pi[xx] = surfX_P1+surfX_M1 - 1./descriptors::invCs2<T,DESCRIPTOR>()*(rho-(T)1) - rho*u[0]*u[0];
    pi[yy] = surfY_P1+surfY_M1 - 1./descriptors::invCs2<T,DESCRIPTOR>()*(rho-(T)1) - rho*u[1]*u[1];
    pi[zz] = surfZ_P1+surfZ_M1 - 1./descriptors::invCs2<T,DESCRIPTOR>()*(rho-(T)1) - rho*u[2]*u[2];

    pi[xy] = cell[4] - cell[5] + cell[13] - cell[14] - rho*u[0]*u[1];
    pi[xz] = cell[6] - cell[7] + cell[15] - cell[16] - rho*u[0]*u[2];
    pi[yz] = cell[8] - cell[9] + cell[17] - cell[18] - rho*u[1]*u[2];
  }

  static void computeAllMomenta(ConstCell<T,DESCRIPTOR>& cell, T& rho, T u[3], T pi[6])
  {
    // Workaround for Intel(r) compiler 9.1;
    // "using namespace util::tensorIndices3D" is not sufficient
    using util::tensorIndices3D::xx;
    using util::tensorIndices3D::yy;
    using util::tensorIndices3D::zz;
    using util::tensorIndices3D::xy;
    using util::tensorIndices3D::xz;
    using util::tensorIndices3D::yz;

    T surfX_M1, surfX_0, surfX_P1,
    surfY_M1, surfY_P1, surfZ_M1, surfZ_P1;

    partial_rho(cell, surfX_M1, surfX_0, surfX_P1,
                surfY_M1, surfY_P1, surfZ_M1, surfZ_P1);

    rho = surfX_M1 + surfX_0 + surfX_P1 + (T)1;
    T invRho = 1. / rho;

    u[0]  = ( surfX_P1 - surfX_M1 ) * invRho;
    u[1]  = ( surfY_P1 - surfY_M1 ) * invRho;
    u[2]  = ( surfZ_P1 - surfZ_M1 ) * invRho;

    pi[xx] = surfX_P1+surfX_M1 - 1./descriptors::invCs2<T,DESCRIPTOR>()*(rho-(T)1) - rho*u[0]*u[0];
    pi[yy] = surfY_P1+surfY_M1 - 1./descriptors::invCs2<T,DESCRIPTOR>()*(rho-(T)1) - rho*u[1]*u[1];
    pi[zz] = surfZ_P1+surfZ_M1 - 1./descriptors::invCs2<T,DESCRIPTOR>()*(rho-(T)1) - rho*u[2]*u[2];

    pi[xy] = cell[4] - cell[5] + cell[13] - cell[14] - rho*u[0]*u[1];
    pi[xz] = cell[6] - cell[7] + cell[15] - cell[16] - rho*u[0]*u[2];
    pi[yz] = cell[8] - cell[9] + cell[17] - cell[18] - rho*u[1]*u[2];
  }

  static T computeRho(ConstCell<T,DESCRIPTOR>& cell)
  {
    T rho = cell[0] + cell[1] + cell[2] + cell[3] + cell[4]
            + cell[5] + cell[6] + cell[7] + cell[8]
            + cell[9] + cell[10] + cell[11] + cell[12]
            + cell[13] + cell[14] + cell[15] + cell[16]
            + cell[17] + cell[18] + (T)1;
    return rho;
  }

  static void modifyVelocity(Cell<T,DESCRIPTOR>& cell, const T newU[3])
  {
    T rho, oldU[3];
    computeRhoU(cell, rho, oldU);
    const T oldUSqr = util::normSqr<T,3>(oldU);
    const T newUSqr = util::normSqr<T,3>(newU);
    for (int iPop=0; iPop<19; ++iPop) {
      cell[iPop] = cell[iPop]
                   - equilibrium(iPop, rho, oldU, oldUSqr)
                   + equilibrium(iPop, rho, newU, newUSqr);
    }
  }

};  //struct lbDynamicsHelpers<D3Q19<>>

// Efficient specialization for D3Q19 lattice with force
template<typename T>
struct lbExternalHelpers<T, descriptors::D3Q19<descriptors::FORCE>> {

  static void addExternalForce(
    Cell<T,descriptors::D3Q19<descriptors::FORCE>>& cell,
    const T u[descriptors::D3Q19<descriptors::FORCE>::d], T omega, T amplitude)
  {
    auto force = cell.template getFieldPointer<descriptors::FORCE>();
    T mu = amplitude*((T)1-omega/(T)2);

    cell[0]  += mu             *( force[0] * (-  u[0]                      ) +
                                  force[1] * (        -   u[1]             ) +
                                  force[2] * (                 -   u[2]    )   );
    cell[1]  += mu *(T)1/(T)6  *( force[0] * ( 2*u[0]                   - 1) +
                                  force[1] * (        -   u[1]             ) +
                                  force[2] * (                 -   u[2]    )   );
    cell[2]  += mu *(T)1/(T)6  *( force[0] * (-  u[0]                      ) +
                                  force[1] * (        + 2*u[1]          - 1) +
                                  force[2] * (                 -   u[2]    )   );
    cell[3]  += mu *(T)1/(T)6  *( force[0] * (-  u[0]                      ) +
                                  force[1] * (        -   u[1]             ) +
                                  force[2] * (                 + 2*u[2] - 1)   );
    cell[4]  += mu *(T)1/(T)12 *( force[0] * ( 2*u[0] + 3*u[1]          - 1) +
                                  force[1] * ( 3*u[0] + 2*u[1]          - 1) +
                                  force[2] * (                 -   u[2]    )   );
    cell[5]  += mu *(T)1/(T)12 *( force[0] * ( 2*u[0] - 3*u[1]          - 1) +
                                  force[1] * (-3*u[0] + 2*u[1]          + 1) +
                                  force[2] * (                 -   u[2]    )   );
    cell[6]  += mu *(T)1/(T)12 *( force[0] * ( 2*u[0]          + 3*u[2] - 1) +
                                  force[1] * (        -   u[1]             ) +
                                  force[2] * ( 3*u[0]          + 2*u[2] - 1)   );
    cell[7]  += mu *(T)1/(T)12 *( force[0] * ( 2*u[0]          - 3*u[2] - 1) +
                                  force[1] * (        -   u[1]             ) +
                                  force[2] * (-3*u[0]          + 2*u[2] + 1)   );
    cell[8]  += mu *(T)1/(T)12 *( force[0] * (-  u[0]                      ) +
                                  force[1] * (          2*u[1] + 3*u[2] - 1) +
                                  force[2] * (          3*u[1] + 2*u[2] - 1)   );
    cell[9]  += mu *(T)1/(T)12 *( force[0] * (-  u[0]                      ) +
                                  force[1] * (          2*u[1] - 3*u[2] - 1) +
                                  force[2] * (        - 3*u[1] + 2*u[2] + 1)   );
    cell[10] += mu *(T)1/(T)6  *( force[0] * ( 2*u[0]                   + 1) +
                                  force[1] * (        -   u[1]             ) +
                                  force[2] * (                 -   u[2]    )   );
    cell[11] += mu *(T)1/(T)6  *( force[0] * (-  u[0]                      ) +
                                  force[1] * (        + 2*u[1]          + 1) +
                                  force[2] * (                 -   u[2]    )   );
    cell[12] += mu *(T)1/(T)6  *( force[0] * (-  u[0]                      ) +
                                  force[1] * (        -   u[1]             ) +
                                  force[2] * (                 + 2*u[2] + 1)   );
    cell[13] += mu *(T)1/(T)12 *( force[0] * ( 2*u[0] + 3*u[1]          + 1) +
                                  force[1] * ( 3*u[0] + 2*u[1]          + 1) +
                                  force[2] * (                 -   u[2]    )   );
    cell[14] += mu *(T)1/(T)12 *( force[0] * ( 2*u[0] - 3*u[1]          + 1) +
                                  force[1] * (-3*u[0] + 2*u[1]          - 1) +
                                  force[2] * (                 -   u[2]    )   );
    cell[15] += mu *(T)1/(T)12 *( force[0] * ( 2*u[0]          + 3*u[2] + 1) +
                                  force[1] * (        -   u[1]             ) +
                                  force[2] * ( 3*u[0]          + 2*u[2] + 1)   );
    cell[16] += mu *(T)1/(T)12 *( force[0] * ( 2*u[0]          - 3*u[2] + 1) +
                                  force[1] * (        -   u[1]             ) +
                                  force[2] * (-3*u[0]          + 2*u[2] - 1)   );
    cell[17] += mu *(T)1/(T)12 *( force[0] * (-  u[0]                      ) +
                                  force[1] * (          2*u[1] + 3*u[2] + 1) +
                                  force[2] * (          3*u[1] + 2*u[2] + 1)   );
    cell[18] += mu *(T)1/(T)12 *( force[0] * (-  u[0]                      ) +
                                  force[1] * (          2*u[1] - 3*u[2] + 1) +
                                  force[2] * (        - 3*u[1] + 2*u[2] - 1)   );
  }
};

// Efficient specialization for D3Q19 lattice and for forced D3Q19 lattice
//   (operations applying to the whole lattice)

template<typename T>
struct lbLatticeHelpers<T, descriptors::D3Q19<>> {

  static void swapAndStreamCell (
    BlockLattice3D<T,descriptors::D3Q19<>>& lattice,
    int iX, int iY, int iZ, int jX, int jY, int jZ, int iPop, T& fTmp )
  {
    auto iCell = lattice.get(iX,iY,iZ);
    auto jCell = lattice.get(jX,jY,jZ);
    fTmp          = iCell[iPop];
    iCell[iPop]   = iCell[iPop+9];
    iCell[iPop+9] = jCell[iPop];
    jCell[iPop]   = fTmp;
  }

  static void swapAndStream3D(BlockLattice3D<T,descriptors::D3Q19<>>& lattice,
                              int iX, int iY, int iZ)
  {
    T fTmp;
    swapAndStreamCell(lattice, iX, iY, iZ, iX-1, iY,   iZ,   1, fTmp);
    swapAndStreamCell(lattice, iX, iY, iZ, iX,   iY-1, iZ,   2, fTmp);
    swapAndStreamCell(lattice, iX, iY, iZ, iX,   iY  , iZ-1, 3, fTmp);
    swapAndStreamCell(lattice, iX, iY, iZ, iX-1, iY-1, iZ,   4, fTmp);
    swapAndStreamCell(lattice, iX, iY, iZ, iX-1, iY+1, iZ,   5, fTmp);
    swapAndStreamCell(lattice, iX, iY, iZ, iX-1, iY  , iZ-1, 6, fTmp);
    swapAndStreamCell(lattice, iX, iY, iZ, iX-1, iY  , iZ+1, 7, fTmp);
    swapAndStreamCell(lattice, iX, iY, iZ, iX  , iY-1, iZ-1, 8, fTmp);
    swapAndStreamCell(lattice, iX, iY, iZ, iX  , iY-1, iZ+1, 9, fTmp);
  }

};

template<typename T>
struct lbLatticeHelpers<T, descriptors::D3Q19<descriptors::FORCE>> {

  static void swapAndStreamCell (
    BlockLattice3D<T,descriptors::D3Q19<descriptors::FORCE>>& lattice,
    int iX, int iY, int iZ, int jX, int jY, int jZ, int iPop, T& fTmp )
  {
    auto iCell = lattice.get(iX,iY,iZ);
    auto jCell = lattice.get(jX,jY,jZ);
    fTmp          = iCell[iPop];
    iCell[iPop]   = iCell[iPop+9];
    iCell[iPop+9] = jCell[iPop];
    jCell[iPop]   = fTmp;
  }

  static void swapAndStream3D(BlockLattice3D<T,descriptors::D3Q19<descriptors::FORCE>>& lattice,
                              int iX, int iY, int iZ)
  {
    T fTmp;
    swapAndStreamCell(lattice, iX, iY, iZ, iX-1, iY,   iZ,   1, fTmp);
    swapAndStreamCell(lattice, iX, iY, iZ, iX,   iY-1, iZ,   2, fTmp);
    swapAndStreamCell(lattice, iX, iY, iZ, iX,   iY  , iZ-1, 3, fTmp);
    swapAndStreamCell(lattice, iX, iY, iZ, iX-1, iY-1, iZ,   4, fTmp);
    swapAndStreamCell(lattice, iX, iY, iZ, iX-1, iY+1, iZ,   5, fTmp);
    swapAndStreamCell(lattice, iX, iY, iZ, iX-1, iY  , iZ-1, 6, fTmp);
    swapAndStreamCell(lattice, iX, iY, iZ, iX-1, iY  , iZ+1, 7, fTmp);
    swapAndStreamCell(lattice, iX, iY, iZ, iX  , iY-1, iZ-1, 8, fTmp);
    swapAndStreamCell(lattice, iX, iY, iZ, iX  , iY-1, iZ+1, 9, fTmp);
  }

};

}  // namespace olb

#endif
