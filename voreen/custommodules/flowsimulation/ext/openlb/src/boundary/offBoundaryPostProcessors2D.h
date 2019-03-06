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

#ifndef OFF_BOUNDARY_POST_PROCESSORS_2D_H
#define OFF_BOUNDARY_POST_PROCESSORS_2D_H

#include "core/postProcessing.h"
#include "core/blockLattice2D.h"

namespace olb {

/**
* This class computes the Linear Bouzidi BC
*/

template<typename T, template<typename U> class Lattice>
class ZeroVelocityBouzidiLinearPostProcessor2D : public LocalPostProcessor2D<T,Lattice> {
public:
  ZeroVelocityBouzidiLinearPostProcessor2D(int x_, int y_, int iPop_, T dist_);
  int extent() const override
  {
    return 1;
  }
  int extent(int whichDirection) const override
  {
    return 1;
  }
  void process(BlockLattice2D<T,Lattice>& blockLattice) override;
  void processSubDomain(BlockLattice2D<T,Lattice>& blockLattice,
                                int x0_, int x1_, int y0_, int y1_ ) override;
private:
  int x, y;
  int xN, yN, xB, yB;
  int iPop, opp, iPop2;
  T q, dist;
};

template<typename T, template<typename U> class Lattice>
class ZeroVelocityBounceBackPostProcessor2D : public LocalPostProcessor2D<T,Lattice> {
public:
  ZeroVelocityBounceBackPostProcessor2D(int x_, int y_, int iPop_, T dist_);
  int extent() const override
  {
    return 1;
  }
  int extent(int whichDirection) const override
  {
    return 1;
  }
  void process(BlockLattice2D<T,Lattice>& blockLattice) override;
  void processSubDomain(BlockLattice2D<T,Lattice>& blockLattice,
                                int x0_, int x1_, int y0_, int y1_ ) override;
private:
  int x, y;
  int xN, yN;
  int iPop, opp;
  T dist;
};

template<typename T, template<typename U> class Lattice>
class VelocityBouzidiLinearPostProcessor2D : public LocalPostProcessor2D<T,Lattice> {
public:
  VelocityBouzidiLinearPostProcessor2D(int x_, int y_, int iPop_, T dist_);
  int extent() const override
  {
    return 1;
  }
  int extent(int whichDirection) const override
  {
    return 1;
  }
  void process(BlockLattice2D<T,Lattice>& blockLattice) override;
  void processSubDomain(BlockLattice2D<T,Lattice>& blockLattice,
                                int x0_, int x1_, int y0_, int y1_ ) override;
private:
  int x, y;
  int xN, yN, xB, yB;
  int iPop, opp, iPop2;
  T q, dist;
  T ufrac;
};

template<typename T, template<typename U> class Lattice>
class VelocityBounceBackPostProcessor2D : public LocalPostProcessor2D<T,Lattice> {
public:
  VelocityBounceBackPostProcessor2D(int x_, int y_, int iPop_, T dist_);
  int extent() const override
  {
    return 1;
  }
  int extent(int whichDirection) const override
  {
    return 1;
  }
  void process(BlockLattice2D<T,Lattice>& blockLattice) override;
  void processSubDomain(BlockLattice2D<T,Lattice>& blockLattice,
                                int x0_, int x1_, int y0_, int y1_ ) override;
private:
  int x, y;
  int xN, yN;
  int iPop, opp;
  T dist;
};

template<typename T, template<typename U> class Lattice>
class AntiBounceBackPostProcessor2D : public LocalPostProcessor2D<T,Lattice> {
public:
  AntiBounceBackPostProcessor2D(int x_, int y_, int iPop_);
  int extent() const override
  {
    return 1;
  }
  int extent(int whichDirection) const override
  {
    return 1;
  }
  void process(BlockLattice2D<T,Lattice>& blockLattice) override;
  void processSubDomain(BlockLattice2D<T,Lattice>& blockLattice,
                                int x0_, int x1_, int y0_, int y1_ ) override;
private:
  int x, y;
  int xN, yN;
  int iPop, opp;
};

template<typename T, template<typename U> class Lattice>
class BoundaryStreamPostProcessor2D : public LocalPostProcessor2D<T,Lattice> {
public:
  BoundaryStreamPostProcessor2D(int x_, int y_, const bool streamDirections[Lattice<T>::q]);
  int extent() const override
  {
    return 1;
  }
  int extent(int whichDirection) const override
  {
    return 1;
  }
  void process(BlockLattice2D<T,Lattice>& blockLattice) override;
  void processSubDomain(BlockLattice2D<T,Lattice>& blockLattice,
                                int x0_, int x1_, int y0_, int y1_ ) override;
private:
  int x, y;
  bool _streamDirections[Lattice<T>::q];
};

/**
* Linear Bouzidi BC Generator
*/

template<typename T, template<typename U> class Lattice>
class ZeroVelocityBouzidiLinearPostProcessorGenerator2D : public PostProcessorGenerator2D<T,Lattice> {
public:
  ZeroVelocityBouzidiLinearPostProcessorGenerator2D(int x_, int y_, int iPop_, T dist_);
  PostProcessor2D<T,Lattice>* generate() const override;
  PostProcessorGenerator2D<T,Lattice>*  clone() const override;
private:
  int x, y;
  int iPop;
  T dist;
};

template<typename T, template<typename U> class Lattice>
class ZeroVelocityBounceBackPostProcessorGenerator2D : public PostProcessorGenerator2D<T,Lattice> {
public:
  ZeroVelocityBounceBackPostProcessorGenerator2D(int x_, int y_, int iPop_, T dist_);
  PostProcessor2D<T,Lattice>* generate() const override;
  PostProcessorGenerator2D<T,Lattice>*  clone() const override;
private:
  int x, y;
  int iPop;
  T dist;
};

template<typename T, template<typename U> class Lattice>
class VelocityBouzidiLinearPostProcessorGenerator2D : public PostProcessorGenerator2D<T,Lattice> {
public:
  VelocityBouzidiLinearPostProcessorGenerator2D(int x_, int y_, int iPop_, T dist_);
  PostProcessor2D<T,Lattice>* generate() const override;
  PostProcessorGenerator2D<T,Lattice>*  clone() const override;
private:
  int x, y;
  int iPop;
  T dist;
};

template<typename T, template<typename U> class Lattice>
class VelocityBounceBackPostProcessorGenerator2D : public PostProcessorGenerator2D<T,Lattice> {
public:
  VelocityBounceBackPostProcessorGenerator2D(int x_, int y_, int iPop_, T dist_);
  PostProcessor2D<T,Lattice>* generate() const override;
  PostProcessorGenerator2D<T,Lattice>*  clone() const override;
private:
  int x, y;
  int iPop;
  T dist;
};

template<typename T, template<typename U> class Lattice>
class AntiBounceBackPostProcessorGenerator2D : public PostProcessorGenerator2D<T,Lattice> {
public:
  AntiBounceBackPostProcessorGenerator2D(int x_, int y_, int iPop_);
  PostProcessor2D<T,Lattice>* generate() const override;
  PostProcessorGenerator2D<T,Lattice>*  clone() const override;
private:
  int x, y;
  int iPop;
};

template<typename T, template<typename U> class Lattice>
class BoundaryStreamPostProcessorGenerator2D : public PostProcessorGenerator2D<T,Lattice> {
public:
  BoundaryStreamPostProcessorGenerator2D(int x_, int y_, const bool _streamDirections[Lattice<T>::q]);
  PostProcessor2D<T,Lattice>* generate() const override;
  PostProcessorGenerator2D<T,Lattice>*  clone() const override;
private:
  int x, y;
  bool _streamDirections[Lattice<T>::q];
};

}

#endif
