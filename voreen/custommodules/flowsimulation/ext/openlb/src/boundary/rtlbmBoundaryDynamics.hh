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

#ifndef RTLBM_BOUNDARY_DYNAMICS_HH
#define RTLBM_BOUNDARY_DYNAMICS_HH


namespace olb {



// For flat Walls
template<typename T, template<typename U> class Lattice, int direction, int orientation>
RtlbmBoundaryDynamics<T,Lattice,direction,orientation>::RtlbmBoundaryDynamics( T omega_, Momenta<T,Lattice>& momenta_)
  : BasicDynamics<T,Lattice>(momenta_)
{
}

template<typename T, template<typename U> class Lattice, int direction, int orientation>
T RtlbmBoundaryDynamics<T,Lattice,direction,orientation>::computeEquilibrium(int iPop, T rho, const T u[Lattice<T>::d], T uSqr) const
{
  return lbHelpers<T,Lattice>::equilibriumFirstOrder( iPop, rho, u );
}

template<typename T, template<typename U> class Lattice, int direction, int orientation>
void RtlbmBoundaryDynamics<T,Lattice,direction,orientation>::collide(Cell<T,Lattice>& cell,LatticeStatistics<T>& statistics)
{
  typedef Lattice<T> L;
  T dirichletTemperature = this->_momenta.computeRho(cell);

  for( int iPop = 0; iPop < L::q; ++iPop ) {
    cell[iPop] = - L::t[iPop];
  }

  std::vector<int> const missingDiagonal = util::subIndexOutgoing<L,direction,orientation>();
  for ( int i : missingDiagonal ) {
    // compute norm of c_iPopMissing
    // is direction axis parallel
    if ( util::normSqr<int,L::d>(L::c[i]) == 1 ) {
      cell[i] = dirichletTemperature - L::t[i];
    }
  }
}

template<typename T, template<typename U> class Lattice, int direction, int orientation>
void RtlbmBoundaryDynamics<T,Lattice,direction,orientation>::staticCollide( Cell<T,Lattice>& cell, const T u[Lattice<T>::d], LatticeStatistics<T>& statistics)
{
  assert(false);
}

template<typename T, template<typename U> class Lattice, int direction, int orientation>
T RtlbmBoundaryDynamics<T,Lattice,direction,orientation>::getOmega() const
{
  return T(-1);
}

template<typename T, template<typename U> class Lattice, int direction, int orientation>
void RtlbmBoundaryDynamics<T,Lattice,direction,orientation>::setOmega(T omega_)
{
}



// for flat diffuse walls
template<typename T, template<typename U> class Lattice, int direction, int orientation>
RtlbmDiffuseBoundaryDynamics<T,Lattice,direction,orientation>::RtlbmDiffuseBoundaryDynamics( T omega_, Momenta<T,Lattice>& momenta_)
  : BasicDynamics<T,Lattice>(momenta_)
{
}

template<typename T, template<typename U> class Lattice, int direction, int orientation>
T RtlbmDiffuseBoundaryDynamics<T,Lattice,direction,orientation>::computeEquilibrium(int iPop, T rho, const T u[Lattice<T>::d], T uSqr) const
{
  return lbHelpers<T,Lattice>::equilibriumFirstOrder( iPop, rho, u );
}

template<typename T, template<typename U> class Lattice, int direction, int orientation>
void RtlbmDiffuseBoundaryDynamics<T,Lattice,direction,orientation>::collide(Cell<T,Lattice>& cell,LatticeStatistics<T>& statistics)
{
    // For direction i \in I_in define
    // cell_i = w_i * dirichlet/sumWeights - w_i
    // For direction i \in I_out defube
    // cell_i = - w_i
    // This construction yields
    // sum_{i=0}^{q-1} cell_i == dirichlet - 1


  // TODO AM, goal more consistent code reading/writting
  // for int i 0 < L::q; if (i \in missing_iPop) then else

  typedef Lattice<T> L;
  // shift all: cell_i = f_i - weight_i
  for ( int iPop = 0; iPop < L::q; ++iPop ) {
    cell[iPop] = - L::t[iPop];
  }

  std::vector<int> const missing_iPop = util::subIndexOutgoing<L,direction,orientation>();
  double sumWeights = 0;
  for ( int i : missing_iPop ) {
    sumWeights += L::t[i];
  }

  T dirichletTemperature = this->_momenta.computeRho(cell);
  for ( int i : missing_iPop ) {
    cell[i] = L::t[i]*dirichletTemperature/sumWeights - L::t[i];
  }
}

template<typename T, template<typename U> class Lattice, int direction, int orientation>
void RtlbmDiffuseBoundaryDynamics<T,Lattice,direction,orientation>::staticCollide( Cell<T,Lattice>& cell, const T u[Lattice<T>::d], LatticeStatistics<T>& statistics)
{
  assert(false);
}

template<typename T, template<typename U> class Lattice, int direction, int orientation>
T RtlbmDiffuseBoundaryDynamics<T,Lattice,direction,orientation>::getOmega() const
{
  return T(-1);
}

template<typename T, template<typename U> class Lattice, int direction, int orientation>
void RtlbmDiffuseBoundaryDynamics<T,Lattice,direction,orientation>::setOmega(T omega_)
{
}



// for edge diffuse walls
template<typename T, template<typename U> class Lattice, int plane, int normal1, int normal2>
RtlbmDiffuseEdgeBoundaryDynamics<T,Lattice,plane,normal1,normal2>::RtlbmDiffuseEdgeBoundaryDynamics( T omega_, Momenta<T,Lattice>& momenta_)
  : BasicDynamics<T,Lattice>(momenta_)
{
}

template<typename T, template<typename U> class Lattice, int plane, int normal1, int normal2>
T RtlbmDiffuseEdgeBoundaryDynamics<T,Lattice,plane,normal1,normal2>::computeEquilibrium(int iPop, T rho, const T u[Lattice<T>::d], T uSqr) const
{
  return lbHelpers<T,Lattice>::equilibriumFirstOrder( iPop, rho, u );
}

template<typename T, template<typename U> class Lattice, int plane, int normal1, int normal2>
void RtlbmDiffuseEdgeBoundaryDynamics<T,Lattice,plane,normal1,normal2>::collide(Cell<T,Lattice>& cell,LatticeStatistics<T>& statistics)
{
    // For direction i \in I_in define
    // cell_i = w_i * dirichlet/sumWeights - w_i
    // For direction i \in I_out defube
    // cell_i = - w_i
    // This construction yields
    // sum_{i=0}^{q-1} cell_i == dirichlet - 1

  typedef Lattice<T> L;

  // shift all: cell_i = f_i - weight_i
  for ( int iPop = 0; iPop < L::q; ++iPop ) {
    cell[iPop] = - L::t[iPop];
  }

  std::vector<int> missing_iPop = util::subIndexOutgoing3DonEdges<L,plane,normal1,normal2>();
  double sumWeights = 0;
  for ( int i : missing_iPop ) {
    sumWeights += L::t[i];
  }

  T dirichletTemperature = this->_momenta.computeRho(cell);
  for ( int i : missing_iPop ) {
    cell[i] = L::t[i]*dirichletTemperature/sumWeights - L::t[i];
  }
}

template<typename T, template<typename U> class Lattice, int plane, int normal1, int normal2>
void RtlbmDiffuseEdgeBoundaryDynamics<T,Lattice,plane,normal1,normal2>::staticCollide( Cell<T,Lattice>& cell, const T u[Lattice<T>::d], LatticeStatistics<T>& statistics)
{
  assert(false);
}

template<typename T, template<typename U> class Lattice, int plane, int normal1, int normal2>
T RtlbmDiffuseEdgeBoundaryDynamics<T,Lattice,plane,normal1,normal2>::getOmega() const
{
  return T(-1);
}

template<typename T, template<typename U> class Lattice, int plane, int normal1, int normal2>
void RtlbmDiffuseEdgeBoundaryDynamics<T,Lattice,plane,normal1,normal2>::setOmega(T omega_)
{
}



// for corner diffuse walls
template<typename T, template<typename U> class Lattice, int xNormal, int yNormal, int zNormal>
RtlbmDiffuseCornerBoundaryDynamics<T,Lattice,xNormal,yNormal,zNormal>::RtlbmDiffuseCornerBoundaryDynamics( T omega_, Momenta<T,Lattice>& momenta_)
  : BasicDynamics<T,Lattice>(momenta_)
{
}

template<typename T, template<typename U> class Lattice, int xNormal, int yNormal, int zNormal>
T RtlbmDiffuseCornerBoundaryDynamics<T,Lattice,xNormal,yNormal,zNormal>::computeEquilibrium(int iPop, T rho, const T u[Lattice<T>::d], T uSqr) const
{
  return lbHelpers<T,Lattice>::equilibriumFirstOrder( iPop, rho, u );
}

template<typename T, template<typename U> class Lattice, int xNormal, int yNormal, int zNormal>
void RtlbmDiffuseCornerBoundaryDynamics<T,Lattice,xNormal,yNormal,zNormal>::collide(Cell<T,Lattice>& cell,LatticeStatistics<T>& statistics)
{
    // For direction i \in I_in define
    // cell_i = w_i * dirichlet/sumWeights - w_i
    // For direction i \in I_out defube
    // cell_i = - w_i
    // This construction yields
    // sum_{i=0}^{q-1} cell_i == dirichlet - 1

  typedef Lattice<T> L;

  // shift all: cell_i = f_i - weight_i
  for ( int iPop = 0; iPop < L::q; ++iPop ) {
    cell[iPop] = - L::t[iPop];
  }

  std::vector<int> const missing_iPop = util::subIndexOutgoing3DonCorners<L,xNormal,yNormal,zNormal>();
  double sumWeights = 0;
  for ( int i : missing_iPop ) {
    sumWeights += L::t[i];
  }

  T dirichletTemperature = this->_momenta.computeRho(cell);
  for ( int i : missing_iPop ) {
    cell[i] = L::t[i]*dirichletTemperature/sumWeights - L::t[i];
  }
}

template<typename T, template<typename U> class Lattice, int xNormal, int yNormal, int zNormal>
void RtlbmDiffuseCornerBoundaryDynamics<T,Lattice,xNormal,yNormal,zNormal>::staticCollide( Cell<T,Lattice>& cell, const T u[Lattice<T>::d], LatticeStatistics<T>& statistics)
{
  assert(false);
}

template<typename T, template<typename U> class Lattice, int xNormal, int yNormal, int zNormal>
T RtlbmDiffuseCornerBoundaryDynamics<T,Lattice,xNormal,yNormal,zNormal>::getOmega() const
{
  return T(-1);
}

template<typename T, template<typename U> class Lattice, int xNormal, int yNormal, int zNormal>
void RtlbmDiffuseCornerBoundaryDynamics<T,Lattice,xNormal,yNormal,zNormal>::setOmega(T omega_)
{
}

}  // namespace olb


#endif
