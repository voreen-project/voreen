/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007, 2008 Jonas Latt
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
 * Dynamics for a generic 3D block lattice view -- generic implementation.
 */
#ifndef BLOCK_LATTICE_VIEW_3D_HH
#define BLOCK_LATTICE_VIEW_3D_HH

#include "blockLatticeView3D.h"
#include "cell.h"

namespace olb {

////////////////////// Class BlockLatticeView3D /////////////////////////

template<typename T, template<typename U> class Lattice>
BlockLatticeView3D<T,Lattice>::BlockLatticeView3D(BlockLatticeStructure3D<T,Lattice>& originalLattice_)
  : BlockLatticeStructure3D<T,Lattice>(originalLattice_.getNx(),
                                       originalLattice_.getNy(),
                                       originalLattice_.getNz()),
    originalLattice(&originalLattice_), x0(0), y0(0), z0(0)
{ }

template<typename T, template<typename U> class Lattice>
BlockLatticeView3D<T,Lattice>::BlockLatticeView3D (
  BlockLatticeStructure3D<T,Lattice>& originalLattice_,
  int x0_, int x1_, int y0_, int y1_, int z0_, int z1_ )
  : BlockLatticeStructure3D<T,Lattice>(x1_-x0_+1, y1_-y0_+1, z1_-z0_+1),
    originalLattice(&originalLattice_), x0(x0_), y0(y0_), z0(z0_)
{
  OLB_PRECONDITION(x0 < originalLattice->getNx());
  //  OLB_PRECONDITION(nx >= 1);
  OLB_PRECONDITION(y0 < originalLattice->getNy());
  //  OLB_PRECONDITION(ny >= 1);
  OLB_PRECONDITION(z0 < originalLattice->getNz());
  //  OLB_PRECONDITION(nz >= 1);
}

template<typename T, template<typename U> class Lattice>
BlockLatticeView3D<T,Lattice>::~BlockLatticeView3D()
{
}

template<typename T, template<typename U> class Lattice>
BlockLatticeView3D<T,Lattice>::BlockLatticeView3D(BlockLatticeView3D<T,Lattice> const& rhs)
  : BlockLatticeStructure3D<T,Lattice>(rhs._nx,rhs._ny,rhs._nz), originalLattice(rhs.originalLattice),
    x0(rhs.x0), y0(rhs.y0), z0(rhs.z0)
{ }

template<typename T, template<typename U> class Lattice>
BlockLatticeView3D<T,Lattice>& BlockLatticeView3D<T,Lattice>::operator= (
  BlockLatticeView3D<T,Lattice> const& rhs )
{
  BlockLatticeView3D<T,Lattice> tmp(rhs);
  swap(tmp);
  return *this;
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeView3D<T,Lattice>::swap (
  BlockLatticeView3D<T,Lattice>& rhs)
{
  std::swap(this->_nx, rhs._nx);
  std::swap(this->_ny, rhs._ny);
  std::swap(this->_nz, rhs._nz);
  std::swap(x0, rhs.x0);
  std::swap(y0, rhs.y0);
  std::swap(z0, rhs.z0);
  std::swap(originalLattice, rhs.originalLattice);
}

template<typename T, template<typename U> class Lattice>
Cell<T,Lattice>& BlockLatticeView3D<T,Lattice>::get(int iX, int iY, int iZ)
{
  OLB_PRECONDITION(iX<originalLattice->getNx());
  OLB_PRECONDITION(iY<originalLattice->getNy());
  OLB_PRECONDITION(iZ<originalLattice->getNz());
  return originalLattice->get(iX+x0, iY+y0, iZ+z0);
}

template<typename T, template<typename U> class Lattice>
Cell<T,Lattice> const& BlockLatticeView3D<T,Lattice>::get (
  int iX, int iY, int iZ ) const
{
  OLB_PRECONDITION(iX<originalLattice->getNx());
  OLB_PRECONDITION(iY<originalLattice->getNy());
  OLB_PRECONDITION(iZ<originalLattice->getNz());
  return originalLattice->get(iX+x0, iY+y0, iZ+z0);
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeView3D<T,Lattice>::initialize()
{
  originalLattice->initialize();
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeView3D<T,Lattice>::defineDynamics (
  int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
  Dynamics<T,Lattice>* dynamics )
{
  originalLattice->defineDynamics( x0_+x0, x1_+x0,
                                   y0_+y0, y1_+y0,
                                   z0_+z0, z1_+z0,
                                   dynamics );

}

template<typename T, template<typename U> class Lattice>
void BlockLatticeView3D<T,Lattice>::defineDynamics (
  int iX, int iY, int iZ, Dynamics<T,Lattice>* dynamics )
{
  originalLattice->defineDynamics( iX+x0, iY+y0, iZ+z0, dynamics );
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeView3D<T,Lattice>::defineDynamics (
  BlockGeometryStructure3D<T>& blockGeometry, int material, Dynamics<T,Lattice>* dynamics)
{
  for (int iX=0; iX<this->_nx; ++iX) {
    for (int iY=0; iY<this->_ny; ++iY) {
      for (int iZ=0; iZ<this->_nz; ++iZ) {
        if (blockGeometry.getMaterial(iX, iY, iZ)==material) {
          defineDynamics(iX,iY,iZ, dynamics);
        }
      }
    }
  }
}

template<typename T, template<typename U> class Lattice>
Dynamics<T,Lattice>* BlockLatticeView3D<T,Lattice>::getDynamics (
        int iX, int iY, int iZ)
{
    return originalLattice->getDynamics(iX+x0, iY+y0, iZ+z0);
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeView3D<T,Lattice>::collide (
  int x0_, int x1_, int y0_, int y1_, int z0_, int z1_ )
{
  originalLattice->collide( x0_+x0, x1_+x0,
                            y0_+y0, y1_+y0,
                            z0_+z0, z1_+z0 );

}

template<typename T, template<typename U> class Lattice>
void BlockLatticeView3D<T,Lattice>::collide()
{
  originalLattice->collide( x0, x0+this->_nx-1, y0, y0+this->_ny-1, z0, z0+this->_nz-1);
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeView3D<T,Lattice>::staticCollide (
  int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
  BlockData3D<T,T> const& u)
{
  originalLattice->staticCollide( x0_+x0, x1_+x0,
                                  y0_+y0, y1_+y0,
                                  z0_+z0, z1_+z0, u);

}

template<typename T, template<typename U> class Lattice>
void BlockLatticeView3D<T,Lattice>::staticCollide (
  BlockData3D<T,T> const& u )
{
  originalLattice->staticCollide( x0, x0+this->_nx-1, y0, y0+this->_ny-1,
                                  z0, z0+this->_nz-1, u);
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeView3D<T,Lattice>::stream(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
{
  originalLattice->stream( x0_+x0, x1_+x0,
                           y0_+y0, y1_+y0,
                           z0_+z0, z1_+z0 );
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeView3D<T,Lattice>::stream(bool periodic)
{
  OLB_PRECONDITION(!periodic);
  originalLattice->stream( x0, x0+this->_nx-1, y0, y0+this->_ny-1, z0, z0+this->_nz-1);
  postProcess();
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeView3D<T,Lattice>::collideAndStream (
  int x0_, int x1_, int y0_, int y1_, int z0_, int z1_ )
{
  originalLattice->collideAndStream( x0_+x0, x1_+x0,
                                     y0_+y0, y1_+y0,
                                     z0_+z0, z1_+z0 );
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeView3D<T,Lattice>::collideAndStream(bool periodic)
{
  OLB_PRECONDITION(!periodic);
  originalLattice->collideAndStream(x0, x0+this->_nx-1, y0, y0+this->_ny-1, z0, z0+this->_nz-1);
  postProcess();
}

template<typename T, template<typename U> class Lattice>
T BlockLatticeView3D<T,Lattice>::computeAverageDensity (
  int x0_, int x1_, int y0_, int y1_, int z0_, int z1_ ) const
{
  return originalLattice->computeAverageDensity( x0_+x0, x1_+x0,
         y0_+y0, y1_+y0,
         z0_+z0, z1_+z0 );
}

template<typename T, template<typename U> class Lattice>
T BlockLatticeView3D<T,Lattice>::computeAverageDensity() const
{
  return originalLattice->computeAverageDensity (
           x0, x0+this->_nx-1, y0, y0+this->_ny-1, z0, z0+this->_nz-1 );
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeView3D<T,Lattice>::computeStress(int iX, int iY, int iZ,
        T pi[util::TensorVal<Lattice<T> >::n])
{
    originalLattice->computeStress(iX + x0, iY + y0, iZ + z0, pi);
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeView3D<T,Lattice>::stripeOffDensityOffset (
  int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, T offset )
{
  originalLattice->stripeOffDensityOffset( x0_+x0, x1_+x0,
      y0_+y0, y1_+y0,
      z0_+z0, z1_+z0, offset);
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeView3D<T,Lattice>::stripeOffDensityOffset(T offset)
{
  originalLattice->stripeOffDensityOffset (
    x0, x0+this->_nx-1, y0, y0+this->_ny-1, z0, z0+this->_nz-1, offset);
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeView3D<T,Lattice>::forAll (
  int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
  WriteCellFunctional<T,Lattice> const& application )
{
  originalLattice->forAll( x0_+x0, x1_+x0,
                           y0_+y0, y1_+y0,
                           z0_+z0, z1_+z0, application);
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeView3D<T,Lattice>::forAll(WriteCellFunctional<T,Lattice> const& application)
{
  originalLattice->forAll (
    x0, x0+this->_nx-1, y0, y0+this->_ny-1, z0, z0+this->_nz-1, application );
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeView3D<T,Lattice>::addPostProcessor (
  PostProcessorGenerator3D<T,Lattice> const& ppGen )
{
  PostProcessorGenerator3D<T,Lattice>* shiftedGen = ppGen.clone();
  shiftedGen->shift(x0,y0,z0);
  originalLattice->addPostProcessor(*shiftedGen);
  delete shiftedGen;
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeView3D<T,Lattice>::resetPostProcessors()
{
  originalLattice->resetPostProcessors();
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeView3D<T,Lattice>::postProcess()
{
  originalLattice -> postProcess(x0, x0+this->_nx-1, y0, y0+this->_ny-1, z0, z0+this->_nz-1);
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeView3D<T,Lattice>::postProcess(
  int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
{
  originalLattice -> postProcess( x0_+x0, x1_+x0, y0_+y0, y1_+y0, z0_+z0, z1_+z0 );
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeView3D<T,Lattice>::addLatticeCoupling (
  LatticeCouplingGenerator3D<T,Lattice> const& lcGen,
  std::vector<SpatiallyExtendedObject3D*> partners )
{
  LatticeCouplingGenerator3D<T,Lattice>* shiftedGen = lcGen.clone();
  shiftedGen->shift(x0,y0,z0);
  originalLattice->addLatticeCoupling(*shiftedGen, partners);
  delete shiftedGen;
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeView3D<T,Lattice>::executeCoupling()
{
  originalLattice -> executeCoupling(x0, x0+this->_nx-1, y0, y0+this->_ny-1, z0, z0+this->_nz-1);
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeView3D<T,Lattice>::executeCoupling(
  int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
{
  originalLattice -> executeCoupling( x0_+x0, x1_+x0, y0_+y0, y1_+y0, z0_+z0, z1_+z0 );
}

template<typename T, template<typename U> class Lattice>
void BlockLatticeView3D<T,Lattice>::
subscribeReductions(Reductor<T>& reductor)
{
  originalLattice -> subscribeReductions(reductor);
}

template<typename T, template<typename U> class Lattice>
LatticeStatistics<T>& BlockLatticeView3D<T,Lattice>::
getStatistics()
{
  return originalLattice->getStatistics();
}

template<typename T, template<typename U> class Lattice>
LatticeStatistics<T> const& BlockLatticeView3D<T,Lattice>::getStatistics() const
{
  return originalLattice->getStatistics();
}

}  // namespace olb

#endif
