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


namespace olb {

////////////////////// Class PostProcessor2D /////////////////

template <typename T, typename DESCRIPTOR>
std::string& PostProcessor2D<T,DESCRIPTOR>::getName()
{
  return _name;
}

template <typename T, typename DESCRIPTOR>
std::string const& PostProcessor2D<T,DESCRIPTOR>::getName() const
{
  return _name;
}

template <typename T, typename DESCRIPTOR>
int PostProcessor2D<T,DESCRIPTOR>::getPriority() const
{
  return _priority;
}

////////////////////// Class PostProcessorGenerator2D /////////////////

template<typename T, typename DESCRIPTOR>
PostProcessorGenerator2D<T,DESCRIPTOR>::PostProcessorGenerator2D (
  int x0_, int x1_, int y0_, int y1_)
  : x0(x0_), x1(x1_), y0(y0_), y1(y1_)
{ }

template<typename T, typename DESCRIPTOR>
void PostProcessorGenerator2D<T,DESCRIPTOR>::shift(int deltaX, int deltaY)
{
  x0 += deltaX;
  x1 += deltaX;
  y0 += deltaY;
  y1 += deltaY;
}

template<typename T, typename DESCRIPTOR>
void PostProcessorGenerator2D<T,DESCRIPTOR>::shift(LatticeR<2> delta)
{
  x0 += delta[0];
  x1 += delta[0];
  y0 += delta[1];
  y1 += delta[1];
}

template<typename T, typename DESCRIPTOR>
bool PostProcessorGenerator2D<T,DESCRIPTOR>::
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
  }
  else {
    return false;
  }
}

template<typename T, typename DESCRIPTOR>
bool PostProcessorGenerator2D<T,DESCRIPTOR>::extract(LatticeR<2> lower, LatticeR<2> upper)
{
  return extract(lower[0], upper[0], lower[1], upper[1]);
}

template<typename T, typename DESCRIPTOR>
void PostProcessorGenerator2D<T,DESCRIPTOR>::reset(LatticeR<2> lower, LatticeR<2> upper)
{
  x0 = lower[0];
  x1 = upper[0];
  y0 = lower[1];
  y1 = upper[1];
}

template<typename T, typename DESCRIPTOR>
void PostProcessorGenerator2D<T,DESCRIPTOR>::reset(int x0_, int x1_, int y0_, int y1_)
{
  x0 = x0_;
  x1 = x1_;
  y0 = y0_;
  y1 = y1_;
}


////////////////////// Class LatticeCouplingGenerator2D /////////////////

template<typename T, typename DESCRIPTOR>
LatticeCouplingGenerator2D<T,DESCRIPTOR>::LatticeCouplingGenerator2D (
  int x0_, int x1_, int y0_, int y1_)
  : x0(x0_), x1(x1_), y0(y0_), y1(y1_)
{ }

template<typename T, typename DESCRIPTOR>
void LatticeCouplingGenerator2D<T,DESCRIPTOR>::shift(LatticeR<2> delta)
{
  x0 += delta[0];
  x1 += delta[0];
  y0 += delta[1];
  y1 += delta[1];
}

template<typename T, typename DESCRIPTOR>
void LatticeCouplingGenerator2D<T,DESCRIPTOR>::shift(int deltaX, int deltaY)
{
  x0 += deltaX;
  x1 += deltaX;
  y0 += deltaY;
  y1 += deltaY;
}

template<typename T, typename DESCRIPTOR>
bool LatticeCouplingGenerator2D<T,DESCRIPTOR>::extract(int x0_, int x1_, int y0_, int y1_)
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
  }
  else {
    return false;
  }
}

template<typename T, typename DESCRIPTOR>
bool LatticeCouplingGenerator2D<T,DESCRIPTOR>::extract(LatticeR<2> lower, LatticeR<2> upper)
{
  return extract(lower[0], upper[0], lower[1], upper[1]);
}

template<typename T, typename DESCRIPTOR>
void LatticeCouplingGenerator2D<T,DESCRIPTOR>::reset(int x0_, int x1_, int y0_, int y1_)
{
  x0 = x0_;
  x1 = x1_;
  y0 = y0_;
  y1 = y1_;
}

template<typename T, typename DESCRIPTOR>
void LatticeCouplingGenerator2D<T,DESCRIPTOR>::reset(LatticeR<2> lower, LatticeR<2> upper)
{
  reset(lower[0], upper[0], lower[1], upper[1]);
}



////////////////////// Class PostProcessor3D /////////////////

template <typename T, typename DESCRIPTOR>
std::string& PostProcessor3D<T,DESCRIPTOR>::getName()
{
  return _name;
}

template <typename T, typename DESCRIPTOR>
std::string const& PostProcessor3D<T,DESCRIPTOR>::getName() const
{
  return _name;
}

template <typename T, typename DESCRIPTOR>
int PostProcessor3D<T,DESCRIPTOR>::getPriority() const
{
  return _priority;
}

////////////////////// Class PostProcessorGenerator3D /////////////////

template<typename T, typename DESCRIPTOR>
PostProcessorGenerator3D<T,DESCRIPTOR>::PostProcessorGenerator3D (
  int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
  : x0(x0_), x1(x1_), y0(y0_), y1(y1_), z0(z0_), z1(z1_)
{ }

template<typename T, typename DESCRIPTOR>
void PostProcessorGenerator3D<T,DESCRIPTOR>::shift (
  int deltaX, int deltaY, int deltaZ, int iC_ )
{
  x0 += deltaX;
  x1 += deltaX;
  y0 += deltaY;
  y1 += deltaY;
  z0 += deltaZ;
  z1 += deltaZ;
  iC = iC_;
}

template<typename T, typename DESCRIPTOR>
void PostProcessorGenerator3D<T,DESCRIPTOR>::shift(
  LatticeR<3> delta, int iC_)
{
  x0 += delta[0];
  x1 += delta[0];
  y0 += delta[1];
  y1 += delta[1];
  z0 += delta[2];
  z1 += delta[2];
  iC = iC_;
}

template<typename T, typename DESCRIPTOR>
bool PostProcessorGenerator3D<T,DESCRIPTOR>::
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
  }
  else {
    return false;
  }
}

template<typename T, typename DESCRIPTOR>
bool PostProcessorGenerator3D<T,DESCRIPTOR>::extract(LatticeR<3> lower, LatticeR<3> upper)
{
  return extract(lower[0], upper[0], lower[1], upper[1], lower[2], upper[2]);
}

template<typename T, typename DESCRIPTOR>
void PostProcessorGenerator3D<T,DESCRIPTOR>::reset(LatticeR<3> lower, LatticeR<3> upper)
{
  x0 = lower[0];
  x1 = upper[0];
  y0 = lower[1];
  y1 = upper[1];
  z0 = lower[2];
  z1 = upper[2];
}

template<typename T, typename DESCRIPTOR>
void PostProcessorGenerator3D<T,DESCRIPTOR>::
reset(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
{
  x0 = x0_;
  x1 = x1_;
  y0 = y0_;
  y1 = y1_;
  z0 = z0_;
  z1 = z1_;
}

////////////////////// Class LatticeCouplingGenerator3D /////////////////

template<typename T, typename DESCRIPTOR>
LatticeCouplingGenerator3D<T,DESCRIPTOR>::LatticeCouplingGenerator3D (
  int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
  : x0(x0_), x1(x1_), y0(y0_), y1(y1_), z0(z0_), z1(z1_), iC(-1)
{ }

template<typename T, typename DESCRIPTOR>
void LatticeCouplingGenerator3D<T,DESCRIPTOR>::shift(LatticeR<3> delta, int iC_)
{
  x0 += delta[0];
  x1 += delta[0];
  y0 += delta[1];
  y1 += delta[1];
  z0 += delta[2];
  z1 += delta[2];
  iC = iC_;
}

template<typename T, typename DESCRIPTOR>
void LatticeCouplingGenerator3D<T,DESCRIPTOR>::shift (
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

template<typename T, typename DESCRIPTOR>
bool LatticeCouplingGenerator3D<T,DESCRIPTOR>::
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
  }
  else {
    return false;
  }
}

template<typename T, typename DESCRIPTOR>
bool LatticeCouplingGenerator3D<T,DESCRIPTOR>::extract(LatticeR<3> lower, LatticeR<3> upper)
{
  return extract(lower[0], upper[0], lower[1], upper[1], lower[2], upper[2]);
}

template<typename T, typename DESCRIPTOR>
void LatticeCouplingGenerator3D<T,DESCRIPTOR>::reset(LatticeR<3> lower, LatticeR<3> upper)
{
  x0 = lower[0];
  x1 = upper[0];
  y0 = lower[1];
  y1 = upper[1];
  z0 = lower[2];
  z1 = upper[2];
}

template<typename T, typename DESCRIPTOR>
void LatticeCouplingGenerator3D<T,DESCRIPTOR>::
reset(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_)
{
  x0 = x0_;
  x1 = x1_;
  y0 = y0_;
  y1 = y1_;
  z0 = z0_;
  z1 = z1_;
}

////////////////////// Class StatisticsPostProcessor3D //////////////

template <typename BLOCK>
void StatisticsPostProcessor::type<BLOCK>::setup(BLOCK& blockLattice)
{ }

template <typename BLOCK>
void StatisticsPostProcessor::type<BLOCK>::apply(BLOCK& blockLattice)
{
  if (!blockLattice.statisticsEnabled()) {
    return;
  }

  #ifndef PARALLEL_MODE_OMP
  blockLattice.getStatistics().reset();
  #endif
  #ifdef PARALLEL_MODE_OMP
  using V = typename BLOCK::value_t;
  #pragma omp parallel
  blockLattice.getStatistics().reset();

  int numCells     = 0;
  V avRho    = V();
  V avEnergy = V();
  V maxU     = V();

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
    avRho = V();
    avEnergy = V();
    maxU = V();
    numCells = 0;
  }
  else {
    avRho    = avRho / numCells;
    avEnergy = avEnergy / numCells;
  }
  #pragma omp parallel
  blockLattice.getStatistics().reset(avRho, avEnergy, maxU, numCells);
  #endif
}

}  // namespace olb

#endif
