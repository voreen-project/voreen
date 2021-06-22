/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006-2008 Jonas Latt
 *                2008-2020 Mathias Krause
 *                2020 Adrian Kummerlaender
 *  OMP parallel code by Mathias Krause, Copyright (C) 2007
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
 * The dynamics of a 3D block lattice -- generic implementation.
 */
#ifndef BLOCK_LATTICE_3D_HH
#define BLOCK_LATTICE_3D_HH

#include <algorithm>

#include "util.h"
#include "blockLattice3D.h"
#include "functors/lattice/indicator/blockIndicatorF3D.h"
#include "dynamics/dynamics.h"
#include "dynamics/lbHelpers.h"
#include "communication/loadBalancer.h"
#include "communication/blockLoadBalancer.h"
#include "communication/ompManager.h"

namespace olb {

////////////////////// Class BlockLattice3D /////////////////////////

/** \param nx_ lattice width (first index)
 *  \param ny_ lattice height (second index)
 *  \param nz_ lattice depth (third index)
 */
template<typename T, typename DESCRIPTOR>
BlockLattice3D<T,DESCRIPTOR>::BlockLattice3D(int nX, int nY, int nZ)
  : BlockLatticeStructure3D<T,DESCRIPTOR>(nX, nY, nZ),
    _staticPopulationD(nX, nY, nZ),
    _staticFieldsD(nX, nY, nZ),
    _dynamicFieldsD(nX, nY, nZ),
    _dynamicsMap(this->getNcells())
{
  resetPostProcessors();

#ifdef PARALLEL_MODE_OMP
  statistics = new LatticeStatistics<T>* [3*omp.get_size()];
  #pragma omp parallel
  {
    statistics[omp.get_rank() + omp.get_size()]
      = new LatticeStatistics<T>;
    statistics[omp.get_rank()] = new LatticeStatistics<T>;
    statistics[omp.get_rank() + 2*omp.get_size()]
      = new LatticeStatistics<T>;
  }
#else
  statistics = new LatticeStatistics<T>;
  statistics->initialize();
#endif
}

/** During destruction, the memory for the lattice and the contained
 * cells is released. However, the dynamics objects pointed to by
 * the cells must be deleted manually by the user.
 */
template<typename T, typename DESCRIPTOR>
BlockLattice3D<T,DESCRIPTOR>::~BlockLattice3D()
{
  clearPostProcessors();
  clearLatticeCouplings();
#ifdef PARALLEL_MODE_OMP
  #pragma omp parallel
  {
    delete statistics[omp.get_rank()];
  }
  delete statistics;
#else
  delete statistics;
#endif
}

template<typename T, typename DESCRIPTOR>
void BlockLattice3D<T,DESCRIPTOR>::initialize()
{
  std::stable_sort(_postProcessors.begin(),
                   _postProcessors.end(),
                   [](PostProcessor3D<T,DESCRIPTOR>* lhs, PostProcessor3D<T,DESCRIPTOR>* rhs) -> bool {
                     return lhs->getPriority() <= rhs->getPriority();
                   });
  postProcess();
}

template<typename T, typename DESCRIPTOR>
Dynamics<T,DESCRIPTOR>* BlockLattice3D<T,DESCRIPTOR>::getDynamics (
  int iX, int iY, int iZ)
{
  return get(iX, iY, iZ).getDynamics();
}

template<typename T, typename DESCRIPTOR>
void BlockLattice3D<T,DESCRIPTOR>::defineDynamics (
  int iX, int iY, int iZ, Dynamics<T,DESCRIPTOR>* dynamics )
{
  OLB_PRECONDITION(iX>=0 && iX<this->_nx);
  OLB_PRECONDITION(iY>=0 && iY<this->_ny);
  OLB_PRECONDITION(iZ>=0 && iZ<this->_nz);

  get(iX, iY, iZ).defineDynamics(dynamics);
}

/** The dynamics object is not duplicated: all cells of the rectangular
 * domain point to the same dynamics.
 *
 * The dynamics object is not owned by the BlockLattice3D object, its
 * memory management must be taken care of by the user.
 */
template<typename T, typename DESCRIPTOR>
void BlockLattice3D<T,DESCRIPTOR>::defineDynamics (
  int x0, int x1, int y0, int y1, int z0, int z1,
  Dynamics<T,DESCRIPTOR>* dynamics )
{
  OLB_PRECONDITION(x0>=0 && x1<this->_nx);
  OLB_PRECONDITION(x1>=x0);
  OLB_PRECONDITION(y0>=0 && y1<this->_ny);
  OLB_PRECONDITION(y1>=y0);
  OLB_PRECONDITION(z0>=0 && z1<this->_nz);
  OLB_PRECONDITION(z1>=z0);

  for (int iX = x0; iX <= x1; ++iX) {
    for (int iY = y0; iY <= y1; ++iY) {
      for (int iZ = z0; iZ <= z1; ++iZ) {
        get(iX, iY, iZ).defineDynamics(dynamics);
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLattice3D<T,DESCRIPTOR>::defineDynamics (
  BlockIndicatorF3D<T>& indicator, Dynamics<T,DESCRIPTOR>* dynamics)
{
  for (int iX = 0; iX < this->_nx; ++iX) {
    for (int iY = 0; iY < this->_ny; ++iY) {
      for (int iZ = 0; iZ < this->_nz; ++iZ) {
        if (indicator(iX, iY, iZ)) {
          get(iX, iY, iZ).defineDynamics(dynamics);
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLattice3D<T,DESCRIPTOR>::defineDynamics (
  BlockGeometryStructure3D<T>& blockGeometry, int material, Dynamics<T,DESCRIPTOR>* dynamics)
{
  BlockIndicatorMaterial3D<T> indicator(blockGeometry, std::vector<int>(1, material));
  defineDynamics(indicator, dynamics);
}

template<typename T, typename DESCRIPTOR>
void BlockLattice3D<T,DESCRIPTOR>::collide (
  int x0, int x1, int y0, int y1, int z0, int z1)
{
  OLB_PRECONDITION(x0>=0 && x1<this->_nx);
  OLB_PRECONDITION(x1>=x0);
  OLB_PRECONDITION(y0>=0 && y1<this->_ny);
  OLB_PRECONDITION(y1>=y0);
  OLB_PRECONDITION(z0>=0 && z1<this->_nz);
  OLB_PRECONDITION(z1>=z0);

  int iX;
  auto cell = get(0,0,0);
#ifdef PARALLEL_MODE_OMP
  #pragma omp parallel for schedule(dynamic,1) firstprivate(cell)
#endif
  for (iX=x0; iX<=x1; ++iX) {
    for (int iY=y0; iY<=y1; ++iY) {
      cell.setCellId(this->getCellId(iX,iY,z0));
      for (int iZ=z0; iZ<=z1; ++iZ) {
        cell.collide(getStatistics());
        cell.advanceCellId();
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLattice3D<T,DESCRIPTOR>::stream()
{
  _staticPopulationD.shift();
}

template<typename T, typename DESCRIPTOR>
void BlockLattice3D<T,DESCRIPTOR>::collideAndStream (
  int x0, int x1, int y0, int y1, int z0, int z1)
{
  collide(x0, x1, y0, y1, z0, z1);
  stream();
  postProcess();
}

template<typename T, typename DESCRIPTOR>
void BlockLattice3D<T,DESCRIPTOR>::collide()
{
  collide(0, this->_nx-1, 0, this->_ny-1, 0, this->_nz-1);
}

template<typename T, typename DESCRIPTOR>
void BlockLattice3D<T,DESCRIPTOR>::collideAndStream()
{
  collideAndStream(0, this->_nx-1, 0, this->_ny-1, 0, this->_nz-1);
}

template<typename T, typename DESCRIPTOR>
T BlockLattice3D<T,DESCRIPTOR>::computeAverageDensity (
  int x0, int x1, int y0, int y1, int z0, int z1) const
{
  T sumRho = T();
  for (int iX=x0; iX<=x1; ++iX) {
    for (int iY=y0; iY<=y1; ++iY) {
      for (int iZ=z0; iZ<=z1; ++iZ) {
        T rho, u[DESCRIPTOR::d];
        get(iX,iY,iZ).computeRhoU(rho, u);
        sumRho += rho;
      }
    }
  }
  return sumRho / (T)(x1-x0+1) / (T)(y1-y0+1) / (T)(z1-z0+1);
}

template<typename T, typename DESCRIPTOR>
T BlockLattice3D<T,DESCRIPTOR>::computeAverageDensity() const
{
  return computeAverageDensity(0, this->_nx-1, 0, this->_ny-1, 0, this->_nz-1);
}

template<typename T, typename DESCRIPTOR>
void BlockLattice3D<T,DESCRIPTOR>::computeStress(int iX, int iY, int iZ,
    T pi[util::TensorVal<DESCRIPTOR >::n])
{
  get(iX,iY,iZ).computeStress(pi);
}

template<typename T, typename DESCRIPTOR>
void BlockLattice3D<T,DESCRIPTOR>::stripeOffDensityOffset (
  int x0, int x1, int y0, int y1, int z0, int z1, T offset )
{
  for (int iX=x0; iX<=x1; ++iX) {
    for (int iY=y0; iY<=y1; ++iY) {
      for (int iZ=z0; iZ<=z1; ++iZ) {
        //if (offset<-42000.) {
        //T rho = get(iX,iY,iZ).computeRho();
        // if (rho<0) {
        //for (int iPop=0; iPop<DESCRIPTOR::q; ++iPop) {
        //  if (get(iX,iY,iZ)[iPop] + descriptors::t<T,DESCRIPTOR>(iPop) < T() ) {
        //    get(iX,iY,iZ)[iPop] = -descriptors::t<T,DESCRIPTOR>(iPop)+0.0000001;
        //  }
        //  else if(rho>1.)
        // get(iX,iY,iZ)[iPop] -= descriptors::t<T,DESCRIPTOR>(iPop) * (rho-1.);
        //}
        //}
        //}
        //else {
        // only stripe off if rho stays positive
        //if (get(iX,iY,iZ).computeRho()>offset) {
        for (int iPop=0; iPop<DESCRIPTOR::q; ++iPop) {
          get(iX,iY,iZ)[iPop] -= descriptors::t<T,DESCRIPTOR>(iPop) * offset;
        }
        // }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLattice3D<T,DESCRIPTOR>::stripeOffDensityOffset(T offset)
{
  stripeOffDensityOffset(0, this->_nx-1, 0, this->_ny-1, 0, this->_nz-1, offset);
}

template<typename T, typename DESCRIPTOR>
void BlockLattice3D<T,DESCRIPTOR>::forAll (
  int x0, int x1, int y0, int y1, int z0, int z1,
  WriteCellFunctional<T,DESCRIPTOR> const& application )
{
  for (int iX=x0; iX<=x1; ++iX) {
    for (int iY=y0; iY<=y1; ++iY) {
      for (int iZ=z0; iZ<=z1; ++iZ) {
        int pos[] = {iX, iY, iZ};
        auto cell = get(iX,iY,iZ);
        application.apply(cell, pos);
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLattice3D<T,DESCRIPTOR>::forAll(WriteCellFunctional<T,DESCRIPTOR> const& application)
{
  forAll(0, this->_nx-1, 0, this->_ny-1, 0, this->_nz-1, application);
}


template<typename T, typename DESCRIPTOR>
void BlockLattice3D<T,DESCRIPTOR>::addPostProcessor (
  PostProcessorGenerator3D<T,DESCRIPTOR> const& ppGen )
{
  _postProcessors.push_back(ppGen.generate());
}

template<typename T, typename DESCRIPTOR>
void BlockLattice3D<T,DESCRIPTOR>::resetPostProcessors()
{
  clearPostProcessors();
  StatPPGenerator3D<T,DESCRIPTOR> statPPGenerator;
  addPostProcessor(statPPGenerator);
}

template<typename T, typename DESCRIPTOR>
void BlockLattice3D<T,DESCRIPTOR>::clearPostProcessors()
{
  for (PostProcessor3D<T,DESCRIPTOR>* postProcessor : _postProcessors) {
    delete postProcessor;
  }
  _postProcessors.clear();
}

template<typename T, typename DESCRIPTOR>
void BlockLattice3D<T,DESCRIPTOR>::postProcess()
{
  for (PostProcessor3D<T,DESCRIPTOR>* postProcessor : _postProcessors) {
    postProcessor->process(*this);
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLattice3D<T,DESCRIPTOR>::postProcess (
  int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
{
  for (PostProcessor3D<T,DESCRIPTOR>* postProcessor : _postProcessors) {
    postProcessor->processSubDomain(*this, x0_, x1_, y0_, y1_, z0_, z1_);
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLattice3D<T,DESCRIPTOR>::addLatticeCoupling (
  LatticeCouplingGenerator3D<T,DESCRIPTOR> const& lcGen,
  std::vector<SpatiallyExtendedObject3D*> partners )
{
  _latticeCouplings.push_back(lcGen.generate(partners));
}

template<typename T, typename DESCRIPTOR>
void BlockLattice3D<T,DESCRIPTOR>::executeCoupling()
{
  for (PostProcessor3D<T,DESCRIPTOR>* coupling : _latticeCouplings) {
    coupling->process(*this);
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLattice3D<T,DESCRIPTOR>::executeCoupling (
  int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
{
  for (PostProcessor3D<T,DESCRIPTOR>* coupling : _latticeCouplings) {
    coupling->processSubDomain(*this, x0_, x1_, y0_, y1_, z0_, z1_);
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLattice3D<T,DESCRIPTOR>::clearLatticeCouplings()
{
  for (PostProcessor3D<T,DESCRIPTOR>* coupling : _latticeCouplings) {
    delete coupling;
  }
  _latticeCouplings.clear();
}

template<typename T, typename DESCRIPTOR>
LatticeStatistics<T>& BlockLattice3D<T,DESCRIPTOR>::getStatistics()
{
#ifdef PARALLEL_MODE_OMP
  return *statistics[omp.get_rank()];
#else
  return *statistics;
#endif
}

template<typename T, typename DESCRIPTOR>
LatticeStatistics<T> const&
BlockLattice3D<T,DESCRIPTOR>::getStatistics() const
{
#ifdef PARALLEL_MODE_OMP
  return *statistics[omp.get_rank()];
#else
  return *statistics;
#endif
}


template<typename T, typename DESCRIPTOR>
void BlockLattice3D<T,DESCRIPTOR>::makePeriodic()
{
  static const int vicinity = descriptors::vicinity<DESCRIPTOR>();
  int maxX = this->getNx()-1;
  int maxY = this->getNy()-1;
  int maxZ = this->getNz()-1;
  periodicSurface(0,             vicinity-1,    0,    maxY,              0,       maxZ);
  periodicSurface(maxX-vicinity+1,     maxX,    0,    maxY,              0,       maxZ);
  periodicSurface(vicinity,      maxX-vicinity, 0,    vicinity-1,        0,       maxZ);
  periodicSurface(vicinity,      maxX-vicinity, maxY-vicinity+1, maxY,   0,       maxZ);
  periodicSurface(vicinity,      maxX-vicinity, vicinity, maxY-vicinity, 0, vicinity-1);
  periodicSurface(vicinity,      maxX-vicinity, vicinity, maxY-vicinity, maxZ-vicinity+1, maxZ);
}

template<typename T, typename DESCRIPTOR>
void BlockLattice3D<T,DESCRIPTOR>::periodicSurface (
  int x0, int x1, int y0, int y1, int z0, int z1)
{
  for (int iX=x0; iX<=x1; ++iX) {
    for (int iY=y0; iY<=y1; ++iY) {
      for (int iZ=z0; iZ<=z1; ++iZ) {
        for (int iPop=1; iPop<=DESCRIPTOR::q/2; ++iPop) {
          int nextX = iX + descriptors::c<DESCRIPTOR>(iPop,0);
          int nextY = iY + descriptors::c<DESCRIPTOR>(iPop,1);
          int nextZ = iZ + descriptors::c<DESCRIPTOR>(iPop,2);
          if ( nextX<0 || nextX>=this->getNx() ||
               nextY<0 || nextY>=this->getNy() ||
               nextZ<0 || nextZ>=this->getNz() ) {
            nextX = (nextX+this->getNx())%this->getNx();
            nextY = (nextY+this->getNy())%this->getNy();
            nextZ = (nextZ+this->getNz())%this->getNz();
            std::swap (
              get(iX,iY,iZ)         [iPop+DESCRIPTOR::q/2],
              get(nextX,nextY,nextZ)[iPop] );
          }
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
std::size_t BlockLattice3D<T,DESCRIPTOR>::getNblock() const
{
  return 3
       + _staticPopulationD.getNblock()
       + _staticFieldsD.getNblock()
       + _dynamicFieldsD.getNblock();
}


template<typename T, typename DESCRIPTOR>
std::size_t BlockLattice3D<T,DESCRIPTOR>::getSerializableSize() const
{
  return 3 * sizeof(int)
       + _staticPopulationD.getSerializableSize()
       + _staticFieldsD.getSerializableSize()
       + _dynamicFieldsD.getSerializableSize();
}


template<typename T, typename DESCRIPTOR>
bool* BlockLattice3D<T,DESCRIPTOR>::getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode)
{
  std::size_t currentBlock = 0;
  bool* dataPtr = nullptr;

  registerVar                     (iBlock, sizeBlock, currentBlock, dataPtr, this->_nx);
  registerVar                     (iBlock, sizeBlock, currentBlock, dataPtr, this->_ny);
  registerVar                     (iBlock, sizeBlock, currentBlock, dataPtr, this->_nz);
  registerSerializableOfConstSize (iBlock, sizeBlock, currentBlock, dataPtr, _staticPopulationD, loadingMode);
  registerSerializableOfConstSize (iBlock, sizeBlock, currentBlock, dataPtr, _staticFieldsD, loadingMode);
  registerSerializableOfConstSize (iBlock, sizeBlock, currentBlock, dataPtr, _dynamicFieldsD, loadingMode);

  return dataPtr;
}


}  // namespace olb

#endif
