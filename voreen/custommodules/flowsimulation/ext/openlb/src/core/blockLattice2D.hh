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
 * The dynamics of a 2D block lattice -- generic implementation.
 */
#ifndef BLOCK_LATTICE_2D_HH
#define BLOCK_LATTICE_2D_HH

#include <algorithm>

#include "blockLattice2D.h"
#include "dynamics/dynamics.h"
#include "dynamics/lbHelpers.h"
#include "util.h"
#include "communication/loadBalancer.h"
#include "communication/blockLoadBalancer.h"
#include "functors/lattice/indicator/blockIndicatorF2D.h"
#include "communication/ompManager.h"

namespace olb {

////////////////////// Class BlockLattice2D /////////////////////////

template<typename T, typename DESCRIPTOR>
BlockLattice2D<T,DESCRIPTOR>::BlockLattice2D(int nX, int nY)
  : BlockLatticeStructure2D<T,DESCRIPTOR>(nX, nY),
    _staticPopulationD(nX, nY),
    _staticFieldsD(nX, nY),
    _dynamicFieldsD(nX, nY),
    _dynamicsMap(this->getNcells())
{
  resetPostProcessors();

#ifdef PARALLEL_MODE_OMP
  _statistics = new LatticeStatistics<T>* [3*omp.get_size()];
  #pragma omp parallel
  {
    _statistics[omp.get_rank() + omp.get_size()]
      = new LatticeStatistics<T>;
    _statistics[omp.get_rank()] = new LatticeStatistics<T>;
    _statistics[omp.get_rank() + 2*omp.get_size()]
      = new LatticeStatistics<T>;
  }
#else
  _statistics = new LatticeStatistics<T>;
#endif
}

/** During destruction, the memory for the lattice and the contained
 * cells is released. However, the dynamics objects pointed to by
 * the cells must be deleted manually by the user.
 */
template<typename T, typename DESCRIPTOR>
BlockLattice2D<T,DESCRIPTOR>::~BlockLattice2D()
{
  clearPostProcessors();
  clearLatticeCouplings();
#ifdef PARALLEL_MODE_OMP
  #pragma omp parallel
  {
    delete _statistics[omp.get_rank()];
  }
  delete _statistics;
#else
  delete _statistics;
#endif
}

template<typename T, typename DESCRIPTOR>
void BlockLattice2D<T,DESCRIPTOR>::initialize()
{
  std::stable_sort(_postProcessors.begin(),
                   _postProcessors.end(),
                   [](PostProcessor2D<T,DESCRIPTOR>* lhs, PostProcessor2D<T,DESCRIPTOR>* rhs) -> bool {
                     return lhs->getPriority() <= rhs->getPriority();
                   });
  postProcess();
}

template<typename T, typename DESCRIPTOR>
Dynamics<T,DESCRIPTOR>* BlockLattice2D<T,DESCRIPTOR>::getDynamics(int iX, int iY)
{
  return &_dynamicsMap.get(this->getCellId(iX,iY));
}

/** The dynamics object is not duplicated: all cells of the rectangular
 * domain point to the same dynamics.
 *
 * The dynamics object is not owned by the BlockLattice2D object, its
 * memory management is in charge of the user.
 */
template<typename T, typename DESCRIPTOR>
void BlockLattice2D<T,DESCRIPTOR>::defineDynamics (
  int x0, int x1, int y0, int y1, Dynamics<T,DESCRIPTOR>* dynamics )
{
  OLB_PRECONDITION(x0>=0 && x1<this->_nx);
  OLB_PRECONDITION(x1>=x0);
  OLB_PRECONDITION(y0>=0 && y1<this->_ny);
  OLB_PRECONDITION(y1>=y0);

  for (int iX=x0; iX<=x1; ++iX) {
    for (int iY=y0; iY<=y1; ++iY) {
      defineDynamics(iX, iY, dynamics);
    }
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLattice2D<T,DESCRIPTOR>::defineDynamics (
  int iX, int iY, Dynamics<T,DESCRIPTOR>* dynamics )
{
  OLB_PRECONDITION(iX>=0 && iX<this->_nx);
  OLB_PRECONDITION(iY>=0 && iY<this->_ny);

  _dynamicsMap.set(this->getCellId(iX,iY), dynamics);
}

template<typename T, typename DESCRIPTOR>
void BlockLattice2D<T,DESCRIPTOR>::defineDynamics (
  BlockIndicatorF2D<T>& indicator, Dynamics<T,DESCRIPTOR>* dynamics)
{
  int latticeR[2];
  for (latticeR[0] = 0; latticeR[0] < this->_nx; ++latticeR[0]) {
    for (latticeR[1] = 0; latticeR[1] < this->_ny; ++latticeR[1]) {
      if (indicator(latticeR)) {
        get(latticeR).defineDynamics(dynamics);
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLattice2D<T,DESCRIPTOR>::defineDynamics(
  BlockGeometryStructure2D<T>& blockGeometry, int material, Dynamics<T,DESCRIPTOR>* dynamics)
{
  BlockIndicatorMaterial2D<T> indicator(blockGeometry, 1);
  defineDynamics(indicator, dynamics);
}

template<typename T, typename DESCRIPTOR>
void BlockLattice2D<T,DESCRIPTOR>::collide(int x0, int x1, int y0, int y1)
{
  OLB_PRECONDITION(x0>=0 && x1<this->_nx);
  OLB_PRECONDITION(x1>=x0);
  OLB_PRECONDITION(y0>=0 && y1<this->_ny);
  OLB_PRECONDITION(y1>=y0);

  auto cell = get(0,0);
#ifdef PARALLEL_MODE_OMP
  #pragma omp parallel for firstprivate(cell)
#endif
  for (int iX=x0; iX <= x1; ++iX) {
    cell.setCellId(this->getCellId(iX,y0));
    for (int iY=y0; iY <= y1; ++iY) {
      cell.collide(getStatistics());
      cell.advanceCellId();
    }
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLattice2D<T,DESCRIPTOR>::stream()
{
  _staticPopulationD.shift();
}

template<typename T, typename DESCRIPTOR>
void BlockLattice2D<T,DESCRIPTOR>::collideAndStream(int x0, int x1, int y0, int y1)
{
  collide(x0, x1, y0, y1);
  stream();
  postProcess();
}

template<typename T, typename DESCRIPTOR>
void BlockLattice2D<T,DESCRIPTOR>::collide()
{
#ifdef PARALLEL_MODE_OMP
  #pragma omp parallel for
#endif
  for (std::size_t iCell=0; iCell < this->getNcells(); ++iCell) {
    auto cell = get(iCell);
    cell.collide(getStatistics());
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLattice2D<T,DESCRIPTOR>::collideAndStream()
{
  collideAndStream(0, this->_nx-1, 0, this->_ny-1);
}

template<typename T, typename DESCRIPTOR>
T BlockLattice2D<T,DESCRIPTOR>::computeAverageDensity ( int x0, int x1, int y0, int y1) const
{
  T sumRho = T();
  for (int iX=x0; iX<=x1; ++iX) {
    for (int iY=y0; iY<=y1; ++iY) {
      T rho, u[DESCRIPTOR::d];
      get(iX,iY).computeRhoU(rho, u);
      sumRho += rho;
    }
  }
  return sumRho / (T)(x1-x0+1) / (T)(y1-y0+1);
}

template<typename T, typename DESCRIPTOR>
T BlockLattice2D<T,DESCRIPTOR>::computeAverageDensity() const
{
  return computeAverageDensity(0, this->_nx-1, 0, this->_ny-1);
}

template<typename T, typename DESCRIPTOR>
void BlockLattice2D<T,DESCRIPTOR>::computeStress(int iX, int iY, T pi[util::TensorVal<DESCRIPTOR>::n])
{
  get(iX,iY).computeStress(pi);
}

template<typename T, typename DESCRIPTOR>
void BlockLattice2D<T,DESCRIPTOR>::stripeOffDensityOffset ( int x0, int x1, int y0, int y1, T offset )
{
  for (int iX=x0; iX<=x1; ++iX) {
    for (int iY=y0; iY<=y1; ++iY) {
      for (int iPop=0; iPop<DESCRIPTOR::q; ++iPop) {
        get(iX,iY)[iPop] -= descriptors::t<T,DESCRIPTOR>(iPop) * offset;
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLattice2D<T,DESCRIPTOR>::stripeOffDensityOffset(T offset)
{
  stripeOffDensityOffset(0, this->_nx-1, 0, this->_ny-1, offset);
}

template<typename T, typename DESCRIPTOR>
void BlockLattice2D<T,DESCRIPTOR>::forAll (
  int x0, int x1, int y0, int y1, WriteCellFunctional<T,DESCRIPTOR> const& application )
{
  for (int iX=x0; iX<=x1; ++iX) {
    for (int iY=y0; iY<=y1; ++iY) {
      int pos[] = {iX, iY};
      auto cell = get(iX,iY);
      application.apply( cell, pos );
    }
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLattice2D<T,DESCRIPTOR>::forAll(WriteCellFunctional<T,DESCRIPTOR> const& application)
{
  forAll(0, this->_nx-1, 0, this->_ny-1, application);
}

template<typename T, typename DESCRIPTOR>
void BlockLattice2D<T,DESCRIPTOR>::addPostProcessor (
  PostProcessorGenerator2D<T,DESCRIPTOR> const& ppGen )
{
  _postProcessors.push_back(ppGen.generate());
}

template<typename T, typename DESCRIPTOR>
void BlockLattice2D<T,DESCRIPTOR>::resetPostProcessors()
{
  clearPostProcessors();
  StatPPGenerator2D<T,DESCRIPTOR> statPPGenerator;
  addPostProcessor(statPPGenerator);
}

template<typename T, typename DESCRIPTOR>
void BlockLattice2D<T,DESCRIPTOR>::clearPostProcessors()
{
  for (PostProcessor2D<T,DESCRIPTOR>* postProcessor : _postProcessors) {
    delete postProcessor;
  }
  _postProcessors.clear();
}

template<typename T, typename DESCRIPTOR>
void BlockLattice2D<T,DESCRIPTOR>::postProcess()
{
  for (PostProcessor2D<T,DESCRIPTOR>* postProcessor : _postProcessors) {
    postProcessor->process(*this);
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLattice2D<T,DESCRIPTOR>::postProcess(int x0_, int x1_, int y0_, int y1_)
{
  for (PostProcessor2D<T,DESCRIPTOR>* postProcessor : _postProcessors) {
    postProcessor->processSubDomain(*this, x0_, x1_, y0_, y1_);
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLattice2D<T,DESCRIPTOR>::addLatticeCoupling (
  LatticeCouplingGenerator2D<T,DESCRIPTOR> const& lcGen,
  std::vector<SpatiallyExtendedObject2D*> partners )
{
  _latticeCouplings.push_back(lcGen.generate(partners));
}

template<typename T, typename DESCRIPTOR>
void BlockLattice2D<T,DESCRIPTOR>::executeCoupling()
{
  for (PostProcessor2D<T,DESCRIPTOR>* coupling : _latticeCouplings) {
    coupling->process(*this);
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLattice2D<T,DESCRIPTOR>::executeCoupling(int x0_, int x1_, int y0_, int y1_)
{
  for (PostProcessor2D<T,DESCRIPTOR>* coupling : _latticeCouplings) {
    coupling->processSubDomain(*this, x0_, x1_, y0_, y1_);
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLattice2D<T,DESCRIPTOR>::clearLatticeCouplings()
{
  for (PostProcessor2D<T,DESCRIPTOR>* coupling : _latticeCouplings) {
    delete coupling;
  }
  _latticeCouplings.clear();
}

template<typename T, typename DESCRIPTOR>
LatticeStatistics<T>& BlockLattice2D<T,DESCRIPTOR>::getStatistics()
{
#ifdef PARALLEL_MODE_OMP
  return *_statistics[omp.get_rank()];
#else
  return *_statistics;
#endif
}

template<typename T, typename DESCRIPTOR>
LatticeStatistics<T> const& BlockLattice2D<T,DESCRIPTOR>::getStatistics() const
{
#ifdef PARALLEL_MODE_OMP
  return *_statistics[omp.get_rank()];
#else
  return *_statistics;
#endif
}

template<typename T, typename DESCRIPTOR>
std::size_t BlockLattice2D<T,DESCRIPTOR>::getNblock() const
{
  return 2
       + _staticPopulationD.getNblock()
       + _staticFieldsD.getNblock()
       + _dynamicFieldsD.getNblock();
}

template<typename T, typename DESCRIPTOR>
std::size_t BlockLattice2D<T,DESCRIPTOR>::getSerializableSize() const
{
  return 2 * sizeof(int)
       + _staticPopulationD.getSerializableSize()
       + _staticFieldsD.getSerializableSize()
       + _dynamicFieldsD.getSerializableSize();
}

template<typename T, typename DESCRIPTOR>
bool* BlockLattice2D<T,DESCRIPTOR>::getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode)
{
  std::size_t currentBlock = 0;
  bool* dataPtr = nullptr;

  registerVar                     (iBlock, sizeBlock, currentBlock, dataPtr, this->_nx);
  registerVar                     (iBlock, sizeBlock, currentBlock, dataPtr, this->_ny);
  registerSerializableOfConstSize (iBlock, sizeBlock, currentBlock, dataPtr, _staticPopulationD, loadingMode);
  registerSerializableOfConstSize (iBlock, sizeBlock, currentBlock, dataPtr, _staticFieldsD, loadingMode);
  registerSerializableOfConstSize (iBlock, sizeBlock, currentBlock, dataPtr, _dynamicFieldsD, loadingMode);

  return dataPtr;
}

template<typename T, typename DESCRIPTOR>
void BlockLattice2D<T,DESCRIPTOR>::periodicEdge(int x0, int x1, int y0, int y1)
{
  for (int iX=x0; iX<=x1; ++iX) {
    for (int iY=y0; iY<=y1; ++iY) {
      auto iCell = get(iX,iY);
      for (int iPop=1; iPop<=DESCRIPTOR::q/2; ++iPop) {
        int nextX = iX + descriptors::c<DESCRIPTOR>(iPop,0);
        int nextY = iY + descriptors::c<DESCRIPTOR>(iPop,1);
        if ( nextX<0 || nextX>=this->_nx ||
             nextY<0 || nextY>=this->_ny ) {
          nextX = (nextX+this->_nx)%this->_nx;
          nextY = (nextY+this->_ny)%this->_ny;
          std::swap (
            iCell[iPop+DESCRIPTOR::q/2],
            get(nextX,nextY)[iPop] );
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
void BlockLattice2D<T,DESCRIPTOR>::makePeriodic()
{
  static const int vicinity = descriptors::vicinity<DESCRIPTOR>();
  int maxX = this->_nx-1;
  int maxY = this->_ny-1;
  periodicEdge(0,vicinity-1, 0,maxY);
  periodicEdge(maxX-vicinity+1,maxX, 0,maxY);
  periodicEdge(vicinity,maxX-vicinity, 0,vicinity-1);
  periodicEdge(vicinity,maxX-vicinity, maxY-vicinity+1,maxY);
}

}  // namespace olb

#endif
