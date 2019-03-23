/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Jonas Latt
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  Generic collision, which modifies the particle distribution
 *  functions, implemented by Orestis Malaspinas, 2007
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

#ifndef BOUNDARY_POST_PROCESSORS_3D_H
#define BOUNDARY_POST_PROCESSORS_3D_H

#include "core/postProcessing.h"
#include "momentaOnBoundaries.h"
#include "core/blockLattice3D.h"

namespace olb {

/**
* This class computes the skordos BC
* on a plane wall in 3D but with a limited number of terms added to the
* equilibrium distributions (i.e. only the Q_i : Pi term)
*/
template<typename T, template<typename U> class Lattice, int direction, int orientation>
class PlaneFdBoundaryProcessor3D : public LocalPostProcessor3D<T,Lattice> {
public:
  PlaneFdBoundaryProcessor3D (int x0_, int x1_, int y0_, int y1_, int z0_, int z1_);
  int extent() const override
  {
    return 1;
  }
  int extent(int whichDirection) const override
  {
    return 1;
  }
  void process(BlockLattice3D<T,Lattice>& blockLattice) override;
  void processSubDomain(BlockLattice3D<T,Lattice>& blockLattice,
                                int x0_, int x1_, int y0_, int y1_, int z0_, int z1_ ) override;
private:
  template<int deriveDirection>
  void interpolateGradients (
    BlockLattice3D<T,Lattice> const& blockLattice,
    T velDeriv[Lattice<T>::d], int iX, int iY, int iZ ) const;
private:
  int x0, x1, y0, y1, z0, z1;
};

template<typename T, template<typename U> class Lattice, int direction, int orientation>
class PlaneFdBoundaryProcessorGenerator3D : public PostProcessorGenerator3D<T,Lattice> {
public:
  PlaneFdBoundaryProcessorGenerator3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_);
  PostProcessor3D<T,Lattice>* generate() const override;
  PostProcessorGenerator3D<T,Lattice>*  clone() const override;
};


/**
* This class computes a convection BC on a flat wall in 2D
*/
template<typename T, template<typename U> class Lattice, int direction, int orientation>
class StraightConvectionBoundaryProcessor3D : public LocalPostProcessor3D<T,Lattice> {
public:
  StraightConvectionBoundaryProcessor3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, T* uAv_ = NULL);
  ~StraightConvectionBoundaryProcessor3D() override;
  int extent() const override
  {
    return 1;
  }
  int extent(int whichDirection) const override
  {
    return 1;
  }
  void process(BlockLattice3D<T,Lattice>& blockLattice) override;
  void processSubDomain ( BlockLattice3D<T,Lattice>& blockLattice,
                                  int x0_, int x1_, int y0_, int y1_ , int z0_, int z1_) override;
private:
  int x0, x1, y0, y1, z0, z1;
  T**** saveCell;
  T* uAv;
};

template<typename T, template<typename U> class Lattice, int direction, int orientation>
class StraightConvectionBoundaryProcessorGenerator3D : public PostProcessorGenerator3D<T,Lattice> {
public:
  StraightConvectionBoundaryProcessorGenerator3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, T* uAv_ = NULL);
  PostProcessor3D<T,Lattice>* generate() const override;
  PostProcessorGenerator3D<T,Lattice>*  clone() const override;
private:
  T* uAv;
};

/**
* This class computes the skordos BC
* on a convex edge wall in 3D but with a limited number of terms added to the
* equilibrium distributions (i.e. only the Q_i : Pi term)
*/
template<typename T, template<typename U> class Lattice,
         int plane, int normal1, int normal2>
class OuterVelocityEdgeProcessor3D : public LocalPostProcessor3D<T,Lattice> {
public:
  enum { direction1 = (plane+1)%3, direction2 = (plane+2)%3 };
public:
  OuterVelocityEdgeProcessor3D (
    int x0_, int x1_, int y0_, int y1_, int z0_, int z1_ );
  int extent() const override
  {
    return 2;
  }
  int extent(int whichDirection) const override
  {
    return 2;
  }
  void process(BlockLattice3D<T,Lattice>& blockLattice) override;
  void processSubDomain(BlockLattice3D<T,Lattice>& blockLattice,
                                int x0_, int x1_, int y0_, int y1_,
                                int z0_, int z1_ ) override;
private:
  T getNeighborRho(int x, int y, int z, int step1, int step2,
                   BlockLattice3D<T,Lattice> const& blockLattice);
  template<int deriveDirection, int orientation>
  void interpolateGradients (
    BlockLattice3D<T,Lattice> const& blockLattice,
    T velDeriv[Lattice<T>::d], int iX, int iY, int iZ ) const;
private:
  int x0, x1, y0, y1, z0, z1;
};

template<typename T, template<typename U> class Lattice,
         int plane, int normal1, int normal2>
class OuterVelocityEdgeProcessorGenerator3D
  : public PostProcessorGenerator3D<T,Lattice> {
public:
  OuterVelocityEdgeProcessorGenerator3D(int x0_, int x1_, int y0_, int y1_,
                                        int z0_, int z1_);
  PostProcessor3D<T,Lattice>* generate() const override;
  PostProcessorGenerator3D<T,Lattice>*  clone() const override;
};


template<typename T, template<typename U> class Lattice,
         int xNormal, int yNormal, int zNormal>
class OuterVelocityCornerProcessor3D : public LocalPostProcessor3D<T,Lattice> {
public:
  OuterVelocityCornerProcessor3D(int x_, int y_, int z_);
  int extent() const override
  {
    return 2;
  }
  int extent(int whichDirection) const override
  {
    return 2;
  }
  void process(BlockLattice3D<T,Lattice>& blockLattice) override;
  void processSubDomain(BlockLattice3D<T,Lattice>& blockLattice,
                                int x0_, int x1_, int y0_, int y1_,
                                int z0_, int z1_ ) override;
private:
  int x,y,z;
};

template<typename T, template<typename U> class Lattice,
         int xNormal, int yNormal, int zNormal>
class OuterVelocityCornerProcessorGenerator3D
  : public PostProcessorGenerator3D<T,Lattice> {
public:
  OuterVelocityCornerProcessorGenerator3D(int x_, int y_, int z_);
  PostProcessor3D<T,Lattice>* generate() const override;
  PostProcessorGenerator3D<T,Lattice>*  clone() const override;
};

/**
* This class computes a slip BC in 3D
*/

template<typename T, template<typename U> class Lattice>
class SlipBoundaryProcessor3D : public LocalPostProcessor3D<T,Lattice> {
public:
  SlipBoundaryProcessor3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, int discreteNormalX_, int discreteNormalY_, int discreteNormalZ_);
  int extent() const override
  {
    return 0;
  }
  int extent(int whichDirection) const override
  {
    return 0;
  }
  void process(BlockLattice3D<T,Lattice>& blockLattice) override;
  void processSubDomain ( BlockLattice3D<T,Lattice>& blockLattice,
                                  int x0_, int x1_, int y0_, int y1_, int z0_, int z1_ ) override;
private:
  int reflectionPop[Lattice<T>::q];
  int x0, x1, y0, y1, z0, z1;
};


template<typename T, template<typename U> class Lattice>
class SlipBoundaryProcessorGenerator3D : public PostProcessorGenerator3D<T,Lattice> {
public:
  SlipBoundaryProcessorGenerator3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, int discreteNormalX_, int discreteNormalY_, int discreteNormalZ_);
  PostProcessor3D<T,Lattice>* generate() const override;
  PostProcessorGenerator3D<T,Lattice>*  clone() const override;
private:
  int discreteNormalX;
  int discreteNormalY;
  int discreteNormalZ;
};
}

#endif
