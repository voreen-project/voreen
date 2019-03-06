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

#ifndef ADVECTION_DIFFUSION_BOUNDARIES_HH
#define ADVECTION_DIFFUSION_BOUNDARIES_HH

#include "advectionDiffusionBoundaries.h"
#include "dynamics/latticeDescriptors.h"
#include "core/util.h"
#include "dynamics/lbHelpers.h"

namespace olb {

using namespace descriptors;

//==================================================================================================
//==================== For regularized Advection Diffusion Boundary Condition ======================
//============================================================================================


// For flat Walls

template<typename T, template<typename U> class Lattice, typename Dynamics, int direction, int orientation>
AdvectionDiffusionBoundariesDynamics<T,Lattice,Dynamics,direction,orientation>::
AdvectionDiffusionBoundariesDynamics( T omega_, Momenta<T,Lattice>& momenta_)
  : BasicDynamics<T,Lattice>(momenta_), boundaryDynamics(omega_, momenta_)
{
}

template<typename T, template<typename U> class Lattice, typename Dynamics, int direction, int orientation>
T AdvectionDiffusionBoundariesDynamics<T,Lattice,Dynamics,direction,orientation>::
computeEquilibrium(int iPop, T rho, const T u[Lattice<T>::d], T uSqr) const
{
  return lbHelpers<T, Lattice>::equilibriumFirstOrder( iPop, rho, u );
}


template<typename T, template<typename U> class Lattice, typename Dynamics, int direction, int orientation>
void AdvectionDiffusionBoundariesDynamics<T,Lattice,Dynamics,direction,orientation>::
collide(Cell<T,Lattice>& cell,LatticeStatistics<T>& statistics)
{
  typedef Lattice<T> L;
  T dirichletTemperature = this->_momenta.computeRho(cell);

  int missingNormal = 0;
  std::vector<int> missingDiagonal = util::subIndexOutgoing<L,direction,orientation>();
  std::vector<int> knownIndexes = util::remainingIndexes<L>(missingDiagonal);
  // here I know all missing and non missing f_i
  for (unsigned iPop = 0; iPop < missingDiagonal.size(); ++iPop) {
    int numOfNonNullComp = 0;
    for (int iDim = 0; iDim < L::d; ++iDim) {
      numOfNonNullComp += abs(L::c[missingDiagonal[iPop]][iDim]);
    }

    if (numOfNonNullComp == 1) {
      missingNormal = missingDiagonal[iPop];
      missingDiagonal.erase(missingDiagonal.begin()+iPop);
      break;
    }
  }

  T sum = T();
  for (unsigned iPop = 0; iPop < knownIndexes.size(); ++iPop) {
    sum += cell[knownIndexes[iPop]];
  }
  cell[missingNormal] = dirichletTemperature - sum -(T)1;

  // Once all the f_i are known, I can call the collision for the default=Regularized Model.
  boundaryDynamics.collide(cell, statistics);
}

template<typename T, template<typename U> class Lattice, typename Dynamics, int direction, int orientation>
void AdvectionDiffusionBoundariesDynamics<T,Lattice,Dynamics,direction,orientation>::
staticCollide( Cell<T,Lattice>& cell, const T u[Lattice<T>::d], LatticeStatistics<T>& statistics)
{
  assert(false);
}

template<typename T, template<typename U> class Lattice, typename Dynamics, int direction, int orientation>
T AdvectionDiffusionBoundariesDynamics<T,Lattice,Dynamics,direction,orientation>::
getOmega() const
{
  return boundaryDynamics.getOmega();
}

template<typename T, template<typename U> class Lattice, typename Dynamics, int direction, int orientation>
void AdvectionDiffusionBoundariesDynamics<T,Lattice,Dynamics,direction,orientation>::
setOmega(T omega_)
{
  boundaryDynamics.setOmega(omega_);
}

//=================================================================
// For 2D Corners with regularized Dynamic ==============================================
//=================================================================
template<typename T, template<typename U> class Lattice, typename Dynamics, int xNormal, int yNormal>
AdvectionDiffusionCornerDynamics2D<T,Lattice,Dynamics,xNormal,yNormal>::AdvectionDiffusionCornerDynamics2D(
  T omega_, Momenta<T,Lattice>& momenta_)
  : BasicDynamics<T,Lattice>(momenta_),
    boundaryDynamics(omega_, momenta_)
{
}

template<typename T, template<typename U> class Lattice, typename Dynamics,  int xNormal, int yNormal>
T AdvectionDiffusionCornerDynamics2D<T,Lattice,Dynamics,xNormal,yNormal>::computeEquilibrium(int iPop, T rho, const T u[Lattice<T>::d], T uSqr) const
{
  return lbHelpers<T, Lattice>::equilibriumFirstOrder( iPop, rho, u );
}


template<typename T, template<typename U> class Lattice, typename Dynamics,  int xNormal, int yNormal>
void AdvectionDiffusionCornerDynamics2D<T,Lattice,Dynamics,xNormal,yNormal>::collide(Cell<T,Lattice>& cell,LatticeStatistics<T>& statistics)
{
  typedef Lattice<T> L;
  typedef lbHelpers<T,Lattice> lbH;

  T temperature = this->_momenta.computeRho(cell);
  T* u = cell.getExternal(Lattice<T>::ExternalField::velocityBeginsAt);
  // I need to get Missing information on the corners !!!!
  std::vector<int> unknownIndexes = util::subIndexOutgoing2DonCorners<L,xNormal,yNormal>();
  // here I know all missing and non missing f_i


  // The collision procedure for D2Q5 and D3Q7 lattice is the same ...
  // Given the rule f_i_neq = -f_opposite(i)_neq
  // I have the right number of equations for the number of unknowns using these lattices

  for (unsigned iPop = 0; iPop < unknownIndexes.size(); ++iPop) {
    cell[unknownIndexes[iPop]] = lbH::equilibriumFirstOrder(unknownIndexes[iPop], temperature, u)
                                 -(cell[util::opposite<L>(unknownIndexes[iPop])]
                                 - lbH::equilibriumFirstOrder(util::opposite<L>(unknownIndexes[iPop]), temperature, u) ) ;
  }

  // Once all the f_i are known, I can call the collision for the Regularized Model.
  boundaryDynamics.collide(cell, statistics);

}

template<typename T, template<typename U> class Lattice, typename Dynamics, int xNormal, int yNormal>
void AdvectionDiffusionCornerDynamics2D<T,Lattice,Dynamics,xNormal,yNormal>::staticCollide (
  Cell<T,Lattice>& cell,
  const T u[Lattice<T>::d],
  LatticeStatistics<T>& statistics )
{
  assert(false);
}

template<typename T, template<typename U> class Lattice, typename Dynamics, int xNormal, int yNormal>
T AdvectionDiffusionCornerDynamics2D<T,Lattice,Dynamics,xNormal,yNormal>::getOmega() const
{
  return boundaryDynamics.getOmega();
}

template<typename T, template<typename U> class Lattice, typename Dynamics, int xNormal, int yNormal>
void AdvectionDiffusionCornerDynamics2D<T,Lattice,Dynamics,xNormal,yNormal>::setOmega(T omega_)
{
  boundaryDynamics.setOmega(omega_);
}



//=================================================================
// For 3D Corners with regularized Dynamic ==============================================
//=================================================================
template<typename T, template<typename U> class Lattice, typename Dynamics, int xNormal, int yNormal, int zNormal>
AdvectionDiffusionCornerDynamics3D<T,Lattice,Dynamics,xNormal,yNormal,zNormal>::AdvectionDiffusionCornerDynamics3D(
  T omega_, Momenta<T,Lattice>& momenta_)
  : BasicDynamics<T,Lattice>(momenta_),
    boundaryDynamics(omega_, momenta_)
{
}

template<typename T, template<typename U> class Lattice, typename Dynamics,  int xNormal, int yNormal, int zNormal>
T AdvectionDiffusionCornerDynamics3D<T,Lattice,Dynamics,xNormal,yNormal,zNormal>::computeEquilibrium(int iPop, T rho, const T u[Lattice<T>::d], T uSqr) const
{
  return lbHelpers<T, Lattice>::equilibriumFirstOrder( iPop, rho, u );
}


template<typename T, template<typename U> class Lattice, typename Dynamics,  int xNormal, int yNormal, int zNormal>
void AdvectionDiffusionCornerDynamics3D<T,Lattice,Dynamics,xNormal,yNormal,zNormal>::collide(Cell<T,Lattice>& cell,LatticeStatistics<T>& statistics)
{
  typedef Lattice<T> L;
  typedef lbHelpers<T,Lattice> lbH;

  T temperature = this->_momenta.computeRho(cell);
  T* u = cell.getExternal(Lattice<T>::ExternalField::velocityBeginsAt);
  // I need to get Missing information on the corners !!!!
  std::vector<int> unknownIndexes = util::subIndexOutgoing3DonCorners<L,xNormal,yNormal,zNormal>();
  // here I know all missing and non missing f_i


  // The collision procedure for D2Q5 and D3Q7 lattice is the same ...
  // Given the rule f_i_neq = -f_opposite(i)_neq
  // I have the right number of equations for the number of unknowns using these lattices

  for (unsigned iPop = 0; iPop < unknownIndexes.size(); ++iPop) {
    cell[unknownIndexes[iPop]] = lbH::equilibriumFirstOrder(unknownIndexes[iPop], temperature, u)
                                 -(cell[util::opposite<L>(unknownIndexes[iPop])]
                                   - lbH::equilibriumFirstOrder(util::opposite<L>(unknownIndexes[iPop]), temperature, u) ) ;
  }

  // Once all the f_i are known, I can call the collision for the Regularized Model.
  boundaryDynamics.collide(cell, statistics);

}

template<typename T, template<typename U> class Lattice, typename Dynamics, int xNormal, int yNormal, int zNormal>
void AdvectionDiffusionCornerDynamics3D<T,Lattice,Dynamics,xNormal,yNormal,zNormal>::staticCollide (
  Cell<T,Lattice>& cell,
  const T u[Lattice<T>::d],
  LatticeStatistics<T>& statistics )
{
  assert(false);
}

template<typename T, template<typename U> class Lattice, typename Dynamics, int xNormal, int yNormal, int zNormal>
T AdvectionDiffusionCornerDynamics3D<T,Lattice,Dynamics,xNormal,yNormal,zNormal>::getOmega() const
{
  return boundaryDynamics.getOmega();
}

template<typename T, template<typename U> class Lattice, typename Dynamics, int xNormal, int yNormal, int zNormal>
void AdvectionDiffusionCornerDynamics3D<T,Lattice,Dynamics,xNormal,yNormal,zNormal>::setOmega(T omega_)
{
  boundaryDynamics.setOmega(omega_);
}

//=================================================================
// For 3D Edges with regularized Dynamic ==============================================
//=================================================================
template<typename T, template<typename U> class Lattice, typename Dynamics, int plane, int normal1, int normal2>
AdvectionDiffusionEdgesDynamics<T,Lattice,Dynamics,plane,normal1, normal2>::AdvectionDiffusionEdgesDynamics(
  T omega_, Momenta<T,Lattice>& momenta_)
  : BasicDynamics<T,Lattice>(momenta_),
    boundaryDynamics(omega_, momenta_)
{
}

template<typename T, template<typename U> class Lattice, typename Dynamics, int plane, int normal1, int normal2>
T AdvectionDiffusionEdgesDynamics<T,Lattice,Dynamics,plane,normal1, normal2>::computeEquilibrium(int iPop, T rho, const T u[Lattice<T>::d], T uSqr) const
{
  return lbHelpers<T, Lattice>::equilibriumFirstOrder( iPop, rho, u );
}


template<typename T, template<typename U> class Lattice, typename Dynamics,  int plane, int normal1, int normal2>
void AdvectionDiffusionEdgesDynamics<T,Lattice,Dynamics,plane,normal1, normal2>::collide(Cell<T,Lattice>& cell,LatticeStatistics<T>& statistics)
{
  typedef Lattice<T> L;
  typedef lbHelpers<T,Lattice> lbH;

  T temperature = this->_momenta.computeRho(cell);
  T* u = cell.getExternal(Lattice<T>::ExternalField::velocityBeginsAt);
  // I need to get Missing information on the corners !!!!
  std::vector<int> unknownIndexes = util::subIndexOutgoing3DonEdges<L,plane,normal1, normal2>();
  // here I know all missing and non missing f_i


  // The collision procedure for D2Q5 and D3Q7 lattice is the same ...
  // Given the rule f_i_neq = -f_opposite(i)_neq
  // I have the right number of equations for the number of unknowns using these lattices

  for (unsigned iPop = 0; iPop < unknownIndexes.size(); ++iPop) {
    cell[unknownIndexes[iPop]] = lbH::equilibriumFirstOrder(unknownIndexes[iPop], temperature, u)
                                 -(cell[util::opposite<L>(unknownIndexes[iPop])]
                                 - lbH::equilibriumFirstOrder(util::opposite<L>(unknownIndexes[iPop]), temperature, u) ) ;
  }

  // Once all the f_i are known, I can call the collision for the Regularized Model.
  boundaryDynamics.collide(cell, statistics);

}

template<typename T, template<typename U> class Lattice, typename Dynamics, int plane, int normal1, int normal2>
void AdvectionDiffusionEdgesDynamics<T,Lattice,Dynamics,plane,normal1, normal2>::staticCollide (
  Cell<T,Lattice>& cell,
  const T u[Lattice<T>::d],
  LatticeStatistics<T>& statistics )
{
  assert(false);
}

template<typename T, template<typename U> class Lattice, typename Dynamics, int plane, int normal1, int normal2>
T AdvectionDiffusionEdgesDynamics<T,Lattice,Dynamics,plane,normal1, normal2>::getOmega() const
{
  return boundaryDynamics.getOmega();
}

template<typename T, template<typename U> class Lattice, typename Dynamics, int plane, int normal1, int normal2>
void AdvectionDiffusionEdgesDynamics<T,Lattice,Dynamics,plane,normal1, normal2>::setOmega(T omega_)
{
  boundaryDynamics.setOmega(omega_);
}




}  // namespace olb




#endif
