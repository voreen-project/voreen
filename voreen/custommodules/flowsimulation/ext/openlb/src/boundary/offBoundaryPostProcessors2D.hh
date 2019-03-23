/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Jonas Kratzke, Mathias J. Krause
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

#ifndef OFF_BOUNDARY_POST_PROCESSORS_2D_HH
#define OFF_BOUNDARY_POST_PROCESSORS_2D_HH

#include "offBoundaryPostProcessors2D.h"
#include "core/blockLattice2D.h"
#include "core/util.h"
#include "core/cell.h"

namespace olb {

/////////// LinearBouzidiPostProcessor2D /////////////////////////////////////

/* Bouzidi Interpolation scheme of first order
 *
 * fluid nodes               wall  solid node
 * --o-------<-o->-----<-o->--|----x----
 *            xB         x        xN
 * directions: --> iPop
 *             <-- opp
 *
*/

template<typename T, template<typename U> class Lattice>
ZeroVelocityBouzidiLinearPostProcessor2D<T,Lattice>::
ZeroVelocityBouzidiLinearPostProcessor2D(int x_, int y_, int iPop_, T dist_)
  : x(x_), y(y_), iPop(iPop_), dist(dist_)
{
  if (dist < 0 || dist > 1)
    std::cout << "WARNING: Bogus distance at (" << x << "," << y << "): "
              << dist << std::endl;
  typedef Lattice<T> L;
  const int* c = L::c[iPop];
  opp = util::opposite<L>(iPop);
  xN = x + c[0];
  yN = y + c[1];

  if (dist >= 0.5) {
    xB = x - c[0];
    yB = y - c[1];
    q = 1/(2*dist);
    iPop2 = opp;
  } else {
    xB = x;
    yB = y;
    q = 2*dist;
    iPop2 = iPop;
  }
  /*
    std::cout << "ZeroVelocityLinear (" << x << "," << y << "," <<
      "), iPop: " << iPop << ", nP: (" << xN << "," << yN << "," <<
      "), opp: " << opp << ", bP: (" << xB << "," << yB << "," <<
      "), dist: " << dist << ", q: " << q << std::endl;
  */
}

template<typename T, template<typename U> class Lattice>
void ZeroVelocityBouzidiLinearPostProcessor2D<T,Lattice>::
processSubDomain(BlockLattice2D<T,Lattice>& blockLattice, int x0_, int x1_, int y0_, int y1_)
{
  if (util::contained(x, y, x0_, x1_, y0_, y1_) ) {
    process(blockLattice);
  }
}

template<typename T, template<typename U> class Lattice>
void ZeroVelocityBouzidiLinearPostProcessor2D<T,Lattice>::
process(BlockLattice2D<T,Lattice>& blockLattice)
{
  blockLattice.get(x, y)[opp] = q*blockLattice.get(xN, yN)[iPop] +
                                (1-q)*blockLattice.get(xB, yB)[iPop2];
}

template<typename T, template<typename U> class Lattice>
VelocityBouzidiLinearPostProcessor2D<T,Lattice>::
VelocityBouzidiLinearPostProcessor2D(int x_, int y_, int iPop_, T dist_)
  : x(x_), y(y_), iPop(iPop_), dist(dist_)
{
  if (dist < 0 || dist > 1)
    std::cout << "WARNING: Bogus distance at (" << x << "," << y << "," << "): "
              << dist << std::endl;
  typedef Lattice<T> L;
  const int* c = L::c[iPop];
  opp = util::opposite<L>(iPop);
  xN = x + c[0];
  yN = y + c[1];

  if (dist >= 0.5) {
    xB = x - c[0];
    yB = y - c[1];
    q = 1/(2*dist);
    ufrac = q;
    iPop2 = opp;
  } else {
    xB = x;
    yB = y;
    q = 2*dist;
    iPop2 = iPop;
    ufrac = 1;
  }
  /*
    std::cout << "VelocityLinear (" << x << "," << y << "," <<
      "), iPop: " << iPop << ", nP: (" << xN << "," << yN << "," <<
      "), opp: " << opp << ", bP: (" << xB << "," << yB << "," <<
      "), dist: " << dist << ", q: " << q << std::endl;
  */
}

template<typename T, template<typename U> class Lattice>
void VelocityBouzidiLinearPostProcessor2D<T,Lattice>::
processSubDomain(BlockLattice2D<T,Lattice>& blockLattice, int x0_, int x1_, int y0_, int y1_)
{
  if (util::contained(x, y, x0_, x1_, y0_, y1_) ) {
    process(blockLattice);
  }
}

template<typename T, template<typename U> class Lattice>
void VelocityBouzidiLinearPostProcessor2D<T,Lattice>::
process(BlockLattice2D<T,Lattice>& blockLattice)
{
  Dynamics<T,Lattice>* dynamics = blockLattice.getDynamics(xN, yN);
  T u = ufrac*dynamics->getVelocityCoefficient(iPop);
  dynamics->defineRho( blockLattice.get(xN, yN), blockLattice.get(x, y).computeRho() );
  T j = u;// * blockLattice.get(x, y).computeRho();
  blockLattice.get(x, y)[opp] = q*blockLattice.get(xN, yN)[iPop] +
                                (1-q)*blockLattice.get(xB, yB)[iPop2] + j;
}


//////// CornerBouzidiPostProcessor2D ///////////////////

template<typename T, template<typename U> class Lattice>
ZeroVelocityBounceBackPostProcessor2D<T,Lattice>::
ZeroVelocityBounceBackPostProcessor2D(int x_, int y_, int iPop_, T dist_)
  : x(x_), y(y_), iPop(iPop_), dist(dist_)
{
  if (dist < 0 || dist > 1)
    std::cout << "WARNING: Bogus distance at (" << x << "," << y << "," << "): "
              << dist << std::endl;
  typedef Lattice<T> L;
  const int* c = L::c[iPop];
  opp = util::opposite<L>(iPop);
  xN = x + c[0];
  yN = y + c[1];
  /*
    std::cout << "Corner (" << x << "," << y << "," <<
        "), iPop: " << iPop << ", nP: (" << xN << "," << yN << "," <<
        "), dist: " << dist << std::endl;
  */
}

template<typename T, template<typename U> class Lattice>
void ZeroVelocityBounceBackPostProcessor2D<T,Lattice>::
processSubDomain(BlockLattice2D<T,Lattice>& blockLattice, int x0_, int x1_, int y0_, int y1_)
{
  if (util::contained(x, y, x0_, x1_, y0_, y1_) ) {
    process(blockLattice);
  }
}

template<typename T, template<typename U> class Lattice>
void ZeroVelocityBounceBackPostProcessor2D<T,Lattice>::
process(BlockLattice2D<T,Lattice>& blockLattice)
{
  blockLattice.get(x, y)[opp] = blockLattice.get(xN, yN)[iPop];
}


template<typename T, template<typename U> class Lattice>
VelocityBounceBackPostProcessor2D<T,Lattice>::
VelocityBounceBackPostProcessor2D(int x_, int y_, int iPop_, T dist_)
  : x(x_), y(y_), iPop(iPop_), dist(dist_)
{
  if (dist < 0 || dist > 1)
    std::cout << "WARNING: Bogus distance at (" << x << "," << y << "," << "): "
              << dist << std::endl;
  typedef Lattice<T> L;
  const int* c = L::c[iPop];
  opp = util::opposite<L>(iPop);
  xN = x + c[0];
  yN = y + c[1];

  /*
    std::cout << "Corner (" << x << "," << y << "," <<
        "), iPop: " << iPop << ", nP: (" << xN << "," << yN << "," <<
        "), dist: " << dist << std::endl;
  */
}

template<typename T, template<typename U> class Lattice>
void VelocityBounceBackPostProcessor2D<T,Lattice>::
processSubDomain(BlockLattice2D<T,Lattice>& blockLattice, int x0_, int x1_, int y0_, int y1_)
{
  if (util::contained(x, y, x0_, x1_, y0_, y1_) ) {
    process(blockLattice);
  }
}

template<typename T, template<typename U> class Lattice>
void VelocityBounceBackPostProcessor2D<T,Lattice>::
process(BlockLattice2D<T,Lattice>& blockLattice)
{
  Dynamics<T,Lattice>* dynamics = blockLattice.getDynamics(xN, yN);
  T u = dynamics->getVelocityCoefficient(iPop);
  dynamics->defineRho( blockLattice.get(xN, yN), blockLattice.get(x, y).computeRho() );
  T j = u;//*blockLattice.get(x, y).computeRho();
  blockLattice.get(x, y)[opp] = blockLattice.get(xN, yN)[iPop] + j;
}


template<typename T, template<typename U> class Lattice>
AntiBounceBackPostProcessor2D<T,Lattice>::
AntiBounceBackPostProcessor2D(int x_, int y_, int iPop_)
  : x(x_), y(y_), iPop(iPop_)
{
  typedef Lattice<T> L;
  const int* c = L::c[iPop];
  opp = util::opposite<L>(iPop);
  xN = x + c[0];
  yN = y + c[1];

  /*
    std::cout << "Corner (" << x << "," << y << "," <<
        "), iPop: " << iPop << ", nP: (" << xN << "," << yN << "," <<
        "), dist: " << dist << std::endl;
  */
}

template<typename T, template<typename U> class Lattice>
void AntiBounceBackPostProcessor2D<T,Lattice>::
processSubDomain(BlockLattice2D<T,Lattice>& blockLattice, int x0_, int x1_, int y0_, int y1_)
{
  if (util::contained(x, y, x0_, x1_, y0_, y1_) ) {
    process(blockLattice);
  }
}

template<typename T, template<typename U> class Lattice>
void AntiBounceBackPostProcessor2D<T,Lattice>::
process(BlockLattice2D<T,Lattice>& blockLattice)
{
  /*Dynamics<T,Lattice>* dynamics = blockLattice.getDynamics(xN, yN);
  T u = dynamics->getVelocityCoefficient(iPop);
  dynamics->defineRho( blockLattice.get(xN, yN), blockLattice.get(x, y).computeRho() );*/
  //T j = u;//*blockLattice.get(x, y).computeRho();
  if (Lattice<T>::c[iPop][1]==0) {
    blockLattice.get(x, y)[opp] = -blockLattice.get(xN, yN)[iPop];  // + j;
  }
  //std::cout << "here" << std::endl;
}


template<typename T, template<typename U> class Lattice>
BoundaryStreamPostProcessor2D<T,Lattice>::
BoundaryStreamPostProcessor2D(int x_, int y_, const bool streamDirection[Lattice<T>::q])
  : x(x_), y(y_)
{
  for (int iPop = 0; iPop < Lattice<T>::q ; ++iPop) {
    this->_streamDirections[iPop] = streamDirection[iPop];
  }
}

template<typename T, template<typename U> class Lattice>
void BoundaryStreamPostProcessor2D<T,Lattice>::
processSubDomain(BlockLattice2D<T,Lattice>& blockLattice, int x0_, int x1_, int y0_, int y1_)
{
  if (util::contained(x, y, x0_, x1_, y0_, y1_) ) {
    process(blockLattice);
  }
}

template<typename T, template<typename U> class Lattice>
void BoundaryStreamPostProcessor2D<T,Lattice>::
process(BlockLattice2D<T,Lattice>& blockLattice)
{
  for (int iPop = 1; iPop < Lattice<T>::q ; ++iPop) {
    if (_streamDirections[iPop]) {
      blockLattice.get(x + Lattice<T>::c[iPop][0], y + Lattice<T>::c[iPop][1])[iPop] = blockLattice.get(x, y)[iPop];
    }
  }
}


////////  LinearBouzidiBoundaryPostProcessorGenerator ////////////////////////////////

template<typename T, template<typename U> class Lattice>
ZeroVelocityBouzidiLinearPostProcessorGenerator2D<T,Lattice>::
ZeroVelocityBouzidiLinearPostProcessorGenerator2D(int x_, int y_, int iPop_, T dist_)
  : PostProcessorGenerator2D<T,Lattice>(x_, x_, y_, y_),
    x(x_), y(y_), iPop(iPop_), dist(dist_)
{ }

template<typename T, template<typename U> class Lattice>
PostProcessor2D<T,Lattice>*
ZeroVelocityBouzidiLinearPostProcessorGenerator2D<T,Lattice>::generate() const
{
  return new ZeroVelocityBouzidiLinearPostProcessor2D<T,Lattice>
         ( this->x, this->y, this->iPop, this->dist);
}

template<typename T, template<typename U> class Lattice>
PostProcessorGenerator2D<T,Lattice>*
ZeroVelocityBouzidiLinearPostProcessorGenerator2D<T,Lattice>::clone() const
{
  return new ZeroVelocityBouzidiLinearPostProcessorGenerator2D<T,Lattice>
         (this->x, this->y, this->iPop, this->dist);
}

template<typename T, template<typename U> class Lattice>
VelocityBouzidiLinearPostProcessorGenerator2D<T,Lattice>::
VelocityBouzidiLinearPostProcessorGenerator2D(int x_, int y_, int iPop_, T dist_)
  : PostProcessorGenerator2D<T,Lattice>(x_, x_, y_, y_),
    x(x_), y(y_), iPop(iPop_), dist(dist_)
{ }

template<typename T, template<typename U> class Lattice>
PostProcessor2D<T,Lattice>*
VelocityBouzidiLinearPostProcessorGenerator2D<T,Lattice>::generate() const
{
  return new VelocityBouzidiLinearPostProcessor2D<T,Lattice>
         ( this->x, this->y, this->iPop, this->dist);
}

template<typename T, template<typename U> class Lattice>
PostProcessorGenerator2D<T,Lattice>*
VelocityBouzidiLinearPostProcessorGenerator2D<T,Lattice>::clone() const
{
  return new VelocityBouzidiLinearPostProcessorGenerator2D<T,Lattice>
         (this->x, this->y, this->iPop, this->dist);
}

/////////// CornerBouzidiBoundaryPostProcessorGenerator /////////////////////////////////////

template<typename T, template<typename U> class Lattice>
ZeroVelocityBounceBackPostProcessorGenerator2D<T,Lattice>::
ZeroVelocityBounceBackPostProcessorGenerator2D(int x_, int y_, int iPop_, T dist_)
  : PostProcessorGenerator2D<T,Lattice>(x_, x_, y_, y_),
    x(x_), y(y_), iPop(iPop_), dist(dist_)
{ }

template<typename T, template<typename U> class Lattice>
PostProcessor2D<T,Lattice>*
ZeroVelocityBounceBackPostProcessorGenerator2D<T,Lattice>::generate() const
{
  return new ZeroVelocityBounceBackPostProcessor2D<T,Lattice>
         ( this->x, this->y, this->iPop, this->dist);
}

template<typename T, template<typename U> class Lattice>
PostProcessorGenerator2D<T,Lattice>*
ZeroVelocityBounceBackPostProcessorGenerator2D<T,Lattice>::clone() const
{
  return new ZeroVelocityBounceBackPostProcessorGenerator2D<T,Lattice>
         (this->x, this->y, this->iPop, this->dist);
}

template<typename T, template<typename U> class Lattice>
VelocityBounceBackPostProcessorGenerator2D<T,Lattice>::
VelocityBounceBackPostProcessorGenerator2D(int x_, int y_, int iPop_, T dist_)
  : PostProcessorGenerator2D<T,Lattice>(x_, x_, y_, y_),
    x(x_), y(y_), iPop(iPop_), dist(dist_)
{ }

template<typename T, template<typename U> class Lattice>
PostProcessor2D<T,Lattice>*
VelocityBounceBackPostProcessorGenerator2D<T,Lattice>::generate() const
{
  return new VelocityBounceBackPostProcessor2D<T,Lattice>
         ( this->x, this->y, this->iPop, this->dist);
}

template<typename T, template<typename U> class Lattice>
PostProcessorGenerator2D<T,Lattice>*
VelocityBounceBackPostProcessorGenerator2D<T,Lattice>::clone() const
{
  return new VelocityBounceBackPostProcessorGenerator2D<T,Lattice>
         (this->x, this->y, this->iPop, this->dist);
}


template<typename T, template<typename U> class Lattice>
AntiBounceBackPostProcessorGenerator2D<T,Lattice>::
AntiBounceBackPostProcessorGenerator2D(int x_, int y_, int iPop_)
  : PostProcessorGenerator2D<T,Lattice>(x_, x_, y_, y_),
    x(x_), y(y_), iPop(iPop_)
{ }

template<typename T, template<typename U> class Lattice>
PostProcessor2D<T,Lattice>*
AntiBounceBackPostProcessorGenerator2D<T,Lattice>::generate() const
{
  return new AntiBounceBackPostProcessor2D<T,Lattice>
         ( this->x, this->y, this->iPop);
}

template<typename T, template<typename U> class Lattice>
PostProcessorGenerator2D<T,Lattice>*
AntiBounceBackPostProcessorGenerator2D<T,Lattice>::clone() const
{
  return new AntiBounceBackPostProcessorGenerator2D<T,Lattice>
         (this->x, this->y, this->iPop);
}

template<typename T, template<typename U> class Lattice>
BoundaryStreamPostProcessorGenerator2D<T,Lattice>::
BoundaryStreamPostProcessorGenerator2D(int x_, int y_, const bool streamDirections[Lattice<T>::q])
  : PostProcessorGenerator2D<T,Lattice>(x_, x_, y_, y_),
    x(x_), y(y_)
{
  for (int iPop = 0; iPop < Lattice<T>::q ; ++iPop) {
    this->_streamDirections[iPop] = streamDirections[iPop];
  }
}

template<typename T, template<typename U> class Lattice>
PostProcessor2D<T,Lattice>*
BoundaryStreamPostProcessorGenerator2D<T,Lattice>::generate() const
{
  return new BoundaryStreamPostProcessor2D<T,Lattice>
         ( this->x, this->y, this->_streamDirections);
}

template<typename T, template<typename U> class Lattice>
PostProcessorGenerator2D<T,Lattice>*
BoundaryStreamPostProcessorGenerator2D<T,Lattice>::clone() const
{
  return new BoundaryStreamPostProcessorGenerator2D<T,Lattice>
         (this->x, this->y, this->_streamDirections);
}
}  // namespace olb

#endif
