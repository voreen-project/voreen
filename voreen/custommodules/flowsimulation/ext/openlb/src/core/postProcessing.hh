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

#ifndef POST_PROCESSING_HH
#define POST_PROCESSING_HH

#include "blockLattice2D.h"
#include "blockLattice3D.h"
//#include <cmath>
//#include <numeric>
//#include <limits>
//#include "util.h"

namespace olb {

////////////////////// Class PostProcessorGenerator2D /////////////////

template<typename T, template<typename U> class Lattice>
PostProcessorGenerator2D<T,Lattice>::PostProcessorGenerator2D (
  int x0_, int x1_, int y0_, int y1_)
  : x0(x0_), x1(x1_), y0(y0_), y1(y1_)
{ }

template<typename T, template<typename U> class Lattice>
void PostProcessorGenerator2D<T,Lattice>::shift(int deltaX, int deltaY)
{
  x0 += deltaX;
  x1 += deltaX;
  y0 += deltaY;
  y1 += deltaY;
}

template<typename T, template<typename U> class Lattice>
bool PostProcessorGenerator2D<T,Lattice>::
extract(int x0_, int x1_, int y0_, int y1_)
{
  int newX0, newX1, newY0, newY1;
  if ( util::intersect (
         x0, x1, y0, y1,
         x0_, x1_, y0_, y1_,
         newX0, newX1, newY0, newY1 ) ) {
    x0 = newX0;
    x1 = newX1;
    y0 = newY0;
    y1 = newY1;
    return true;
  } else {
    return false;
  }
}


////////////////////// Class LatticeCouplingGenerator2D /////////////////

template<typename T, template<typename U> class Lattice>
LatticeCouplingGenerator2D<T,Lattice>::LatticeCouplingGenerator2D (
  int x0_, int x1_, int y0_, int y1_)
  : x0(x0_), x1(x1_), y0(y0_), y1(y1_)
{ }

template<typename T, template<typename U> class Lattice>
void LatticeCouplingGenerator2D<T,Lattice>::shift(int deltaX, int deltaY)
{
  x0 += deltaX;
  x1 += deltaX;
  y0 += deltaY;
  y1 += deltaY;
}

template<typename T, template<typename U> class Lattice>
bool LatticeCouplingGenerator2D<T,Lattice>::extract(int x0_, int x1_, int y0_, int y1_)
{
  int newX0, newX1, newY0, newY1;
  if ( util::intersect (
         x0, x1, y0, y1,
         x0_, x1_, y0_, y1_,
         newX0, newX1, newY0, newY1 ) ) {
    x0 = newX0;
    x1 = newX1;
    y0 = newY0;
    y1 = newY1;
    return true;
  } else {
    return false;
  }
}


////////////////////// Class StatisticsPostProcessor2D //////////////

template<typename T, template<typename U> class Lattice>
StatisticsPostProcessor2D<T,Lattice>::StatisticsPostProcessor2D()
{ }

#ifndef PARALLEL_MODE_OMP
template<typename T, template<typename U> class Lattice>
void StatisticsPostProcessor2D<T,Lattice>::process (
  BlockLattice2D<T,Lattice>& blockLattice )
{
  blockLattice.getStatistics().reset();
}
#endif


#ifdef PARALLEL_MODE_OMP
template<typename T, template<typename U> class Lattice>
void StatisticsPostProcessor2D<T,Lattice>::process (
  BlockLattice2D<T,Lattice>& blockLattice )
{
  #pragma omp parallel
  blockLattice.getStatistics().reset();


  int numCells     = 0;
  T avRho    = T();
  T avEnergy = T();
  T maxU     = T();

  #pragma omp parallel
  {
    #pragma omp critical
    {
      numCells       += blockLattice.getStatistics().getNumCells();
      avRho          += blockLattice.getStatistics().getAverageRho()
      *blockLattice.getStatistics().getNumCells();
      avEnergy       += blockLattice.getStatistics().getAverageEnergy()
      *blockLattice.getStatistics().getNumCells();
      if (maxU<blockLattice.getStatistics().getMaxU() )
      {
        maxU        = blockLattice.getStatistics().getMaxU();
      }
    }
  }
  if (numCells==0) {
    // avoid division by zero
    avRho = T();
    avEnergy = T();
    maxU = T();
    numCells = 0;
  } else {
    avRho    = avRho / numCells;
    avEnergy = avEnergy / numCells;
  }
  #pragma omp parallel
  blockLattice.getStatistics().reset(avRho,avEnergy, maxU, numCells);
}
#endif

template<typename T, template<typename U> class Lattice>
void StatisticsPostProcessor2D<T,Lattice>::
subscribeReductions(BlockLattice2D<T,Lattice>& blockLattice, Reductor<T>* reductor)
{
  std::vector<T>& averageVect = blockLattice.getStatistics().getAverageVect();
  for (size_t iVect=0; iVect<averageVect.size(); ++iVect) {
    reductor->subscribeAverage(blockLattice.getStatistics().getNumCells(), averageVect[iVect]);
  }
  std::vector<T>& sumVect = blockLattice.getStatistics().getSumVect();
  for (size_t iVect=0; iVect<sumVect.size(); ++iVect) {
    reductor->subscribeSum(sumVect[iVect]);
  }
  std::vector<T>& minVect = blockLattice.getStatistics().getMinVect();
  for (size_t iVect=0; iVect<minVect.size(); ++iVect) {
    reductor->subscribeMin(minVect[iVect]);
  }
  std::vector<T>& maxVect = blockLattice.getStatistics().getMaxVect();
  for (size_t iVect=0; iVect<maxVect.size(); ++iVect) {
    reductor->subscribeMax(maxVect[iVect]);
  }
}


////////////////////// Class StatPPGenerator2D //////////////

template<typename T, template<typename U> class Lattice>
StatPPGenerator2D<T,Lattice>::StatPPGenerator2D()
  : PostProcessorGenerator2D<T,Lattice>(-1,-1,-1,-1)
{ }

template<typename T, template<typename U> class Lattice>
PostProcessor2D<T,Lattice>* StatPPGenerator2D<T,Lattice>::generate() const
{
  return new StatisticsPostProcessor2D<T,Lattice>;
}


template<typename T, template<typename U> class Lattice>
PostProcessorGenerator2D<T,Lattice>*
StatPPGenerator2D<T,Lattice>::clone() const
{
  return new StatPPGenerator2D;
}


////////////////////// Class PostProcessorGenerator3D /////////////////

template<typename T, template<typename U> class Lattice>
PostProcessorGenerator3D<T,Lattice>::PostProcessorGenerator3D (
  int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
  : x0(x0_), x1(x1_), y0(y0_), y1(y1_), z0(z0_), z1(z1_)
{ }

template<typename T, template<typename U> class Lattice>
void PostProcessorGenerator3D<T,Lattice>::shift (
  int deltaX, int deltaY, int deltaZ )
{
  x0 += deltaX;
  x1 += deltaX;
  y0 += deltaY;
  y1 += deltaY;
  z0 += deltaZ;
  z1 += deltaZ;
}

template<typename T, template<typename U> class Lattice>
bool PostProcessorGenerator3D<T,Lattice>::
extract(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
{
  int newX0, newX1, newY0, newY1, newZ0, newZ1;
  if ( util::intersect (
         x0, x1, y0, y1, z0, z1,
         x0_, x1_, y0_, y1_, z0_, z1_,
         newX0, newX1, newY0, newY1, newZ0, newZ1 ) ) {
    x0 = newX0;
    x1 = newX1;
    y0 = newY0;
    y1 = newY1;
    z0 = newZ0;
    z1 = newZ1;
    return true;
  } else {
    return false;
  }
}

////////////////////// Class LatticeCouplingGenerator3D /////////////////

template<typename T, template<typename U> class Lattice>
LatticeCouplingGenerator3D<T,Lattice>::LatticeCouplingGenerator3D (
  int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
  : x0(x0_), x1(x1_), y0(y0_), y1(y1_), z0(z0_), z1(z1_)
{ }

template<typename T, template<typename U> class Lattice>
void LatticeCouplingGenerator3D<T,Lattice>::shift (
  int deltaX, int deltaY, int deltaZ, int iC_)
{
  x0 += deltaX;
  x1 += deltaX;
  y0 += deltaY;
  y1 += deltaY;
  z0 += deltaZ;
  z1 += deltaZ;
  iC = iC_;
}

template<typename T, template<typename U> class Lattice>
bool LatticeCouplingGenerator3D<T,Lattice>::
extract(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
{
  int newX0, newX1, newY0, newY1, newZ0, newZ1;
  if ( util::intersect (
         x0, x1, y0, y1, z0, z1,
         x0_, x1_, y0_, y1_, z0_, z1_,
         newX0, newX1, newY0, newY1, newZ0, newZ1 ) ) {
    x0 = newX0;
    x1 = newX1;
    y0 = newY0;
    y1 = newY1;
    z0 = newZ0;
    z1 = newZ1;
    return true;
  } else {
    return false;
  }
}


////////////////////// Class StatisticsPostProcessor3D //////////////

template<typename T, template<typename U> class Lattice>
StatisticsPostProcessor3D<T,Lattice>::StatisticsPostProcessor3D()
{ }

#ifndef PARALLEL_MODE_OMP
template<typename T, template<typename U> class Lattice>
void StatisticsPostProcessor3D<T,Lattice>::process (
  BlockLattice3D<T,Lattice>& blockLattice )
{
  blockLattice.getStatistics().reset();
}
#endif
#ifdef PARALLEL_MODE_OMP
template<typename T, template<typename U> class Lattice>
void StatisticsPostProcessor3D<T,Lattice>::process (
  BlockLattice3D<T,Lattice>& blockLattice )
{
  #pragma omp parallel
  blockLattice.getStatistics().reset();


  int numCells     = 0;
  T avRho    = T();
  T avEnergy = T();
  T maxU     = T();

  #pragma omp parallel
  {
    #pragma omp critical
    {
      numCells       += blockLattice.getStatistics().getNumCells();
      avRho          += blockLattice.getStatistics().getAverageRho()
      *blockLattice.getStatistics().getNumCells();
      avEnergy       += blockLattice.getStatistics().getAverageEnergy()
      *blockLattice.getStatistics().getNumCells();
      if (maxU<blockLattice.getStatistics().getMaxU() )
      {
        maxU        = blockLattice.getStatistics().getMaxU();
      }
    }
  }
  if (numCells==0) {
    // avoid division by zero
    avRho = T();
    avEnergy = T();
    maxU = T();
    numCells = 0;
  } else {
    avRho    = avRho / numCells;
    avEnergy = avEnergy / numCells;
  }
  #pragma omp parallel
  blockLattice.getStatistics().reset(avRho,avEnergy, maxU, numCells);
}
#endif


template<typename T, template<typename U> class Lattice>
void StatisticsPostProcessor3D<T,Lattice>::
subscribeReductions(BlockLattice3D<T,Lattice>& blockLattice, Reductor<T>* reductor)
{
  std::vector<T>& averageVect = blockLattice.getStatistics().getAverageVect();
  for (size_t iVect=0; iVect<averageVect.size(); ++iVect) {
    reductor->subscribeAverage(blockLattice.getStatistics().getNumCells(), averageVect[iVect]);
  }
  std::vector<T>& sumVect = blockLattice.getStatistics().getSumVect();
  for (size_t iVect=0; iVect<sumVect.size(); ++iVect) {
    reductor->subscribeSum(sumVect[iVect]);
  }
  std::vector<T>& minVect = blockLattice.getStatistics().getMinVect();
  for (size_t iVect=0; iVect<minVect.size(); ++iVect) {
    reductor->subscribeMin(minVect[iVect]);
  }
  std::vector<T>& maxVect = blockLattice.getStatistics().getMaxVect();
  for (size_t iVect=0; iVect<maxVect.size(); ++iVect) {
    reductor->subscribeMax(maxVect[iVect]);
  }
}


////////////////////// Class StatPPGenerator3D //////////////

template<typename T, template<typename U> class Lattice>
StatPPGenerator3D<T,Lattice>::StatPPGenerator3D()
  : PostProcessorGenerator3D<T,Lattice>(-1,-1,-1,-1,-1,-1)
{ }

template<typename T, template<typename U> class Lattice>
PostProcessor3D<T,Lattice>* StatPPGenerator3D<T,Lattice>::generate() const
{
  return new StatisticsPostProcessor3D<T,Lattice>;
}


template<typename T, template<typename U> class Lattice>
PostProcessorGenerator3D<T,Lattice>* StatPPGenerator3D<T,Lattice>::clone() const
{
  return new StatPPGenerator3D;
}




}  // namespace olb

#endif
