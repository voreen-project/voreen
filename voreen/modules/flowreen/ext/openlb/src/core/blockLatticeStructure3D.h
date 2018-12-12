/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006-2008 Jonas Latt
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
 * Dynamics for a generic 3D block structure -- header file.
 */
#ifndef BLOCK_LATTICE_STRUCTURE_3D_H
#define BLOCK_LATTICE_STRUCTURE_3D_H

#include <vector>
#include "cell.h"
#include "blockData3D.h"
#include "blockStructure3D.h"
#include "postProcessing.h"
#include "serializer.h"
#include "spatiallyExtendedObject3D.h"
#include "geometry/blockGeometryStructure3D.h"
#include "latticeStatistics.h"
#include "functors/analytical/analyticalF.h"
#include "functors/lattice/indicator/indicatorBaseF3D.h"


namespace olb {
template<typename T, typename S> class AnalyticalF3D;

template<typename T, template<typename U> class Lattice> struct Dynamics;
template<typename T, template<typename U> class Lattice> class Cell;
template<typename T, template<typename U> class Lattice> struct WriteCellFunctional;
template<typename T> class BlockIndicatorF3D;
template<typename T> class IndicatorSphere3D;
template<typename T> class BlockGeometryStructure3D;


/** BlockLatticeStructure3D is a interface class for defining dynamics on a BlockStructure3D.
 *  The pure virtual methods are wrapped by generic methods.
 */
template<typename T, template<typename U> class Lattice>
class BlockLatticeStructure3D : public BlockStructure3D, public SpatiallyExtendedObject3D {
public:
  BlockLatticeStructure3D(int nx, int ny, int nz) : BlockStructure3D(nx,ny,nz) {};
  ~BlockLatticeStructure3D() override { }
public:
  virtual void defineRho(BlockGeometryStructure3D<T>& blockGeometry, int material,
                         AnalyticalF3D<T,T>& rho);
  virtual void defineU(BlockGeometryStructure3D<T>& blockGeometry, int material,
                       AnalyticalF3D<T,T>& u);
  virtual void defineRhoU(BlockGeometryStructure3D<T>& blockGeometry,
                          int material, AnalyticalF3D<T,T>& rho, AnalyticalF3D<T,T>& u);
  virtual void definePopulations(BlockGeometryStructure3D<T>& blockGeometry, int material,
                                 AnalyticalF3D<T,T>& Pop);
  virtual void defineExternalField(BlockGeometryStructure3D<T>& blockGeometry,
                                   int material, int fieldBeginsAt, int sizeOfField,
                                   AnalyticalF3D<T,T>& field);
  virtual void defineExternalField(BlockIndicatorF3D<T>& indicator, int overlap,
                                   int fieldBeginsAt, int sizeOfField,
                                   AnalyticalF3D<T,T>& field);
  virtual void defineExternalField(BlockGeometryStructure3D<T>& blockGeometry,
                                   IndicatorF3D<T>& indicator, int fieldBeginsAt, int sizeOfField,
                                   AnalyticalF3D<T,T>& field);
  virtual void setExternalParticleField(BlockGeometryStructure3D<T>& blockGeometry, AnalyticalF3D<T,T>& velocity, ParticleIndicatorF3D<T,T>& sIndicator);

  virtual void iniEquilibrium(BlockGeometryStructure3D<T>& blockGeometry, int material,
                              AnalyticalF3D<T,T>& rho , AnalyticalF3D<T,T>& u);
  // pure virtual member functions
  virtual Cell<T,Lattice>& get(int iX, int iY, int iZ) =0;
  virtual Cell<T,Lattice> const& get(int iX, int iY, int iZ) const =0;
  virtual void initialize() =0;
  virtual void defineDynamics(int x0_, int x1_, int y0_, int y1_, int z0_,
                              int z1_, Dynamics<T,Lattice>* dynamics) =0;
  virtual void defineDynamics(int iX, int iY, int iZ, Dynamics<T,Lattice>* dynamics ) =0;
  virtual void defineDynamics(BlockGeometryStructure3D<T>& blockGeometry,
                              int material, Dynamics<T,Lattice>* dynamics)=0;
  virtual Dynamics<T,Lattice>* getDynamics(int iX, int iY, int iZ) = 0;
  virtual void collide(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_) =0;
  virtual void collide() =0;
  virtual void staticCollide(int x0_, int x1_, int y0_, int y1_, int z0_,
                             int z1_, BlockData3D<T,T> const& u) =0;
  virtual void staticCollide( BlockData3D<T,T> const& u) =0;
  virtual void stream(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_) =0;
  virtual void stream(bool periodic=false) =0;
  virtual void collideAndStream(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_) =0;
  virtual void collideAndStream(bool periodic=false) =0;
  virtual T computeAverageDensity(int x0_, int x1_, int y0_, int y1_, int z0_,
                                  int z1_) const =0;
  virtual T computeAverageDensity() const =0;
  virtual void computeStress(int iX, int iY, int iZ, T pi[util::TensorVal<Lattice<T> >::n]) = 0;
  virtual void stripeOffDensityOffset(int x0_, int x1_, int y0_, int y1_, int z0_,
                                      int z1_, T offset) =0;
  virtual void stripeOffDensityOffset(T offset) =0;
  virtual void forAll(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
                      WriteCellFunctional<T,Lattice> const& application) =0;
  virtual void forAll(WriteCellFunctional<T,Lattice> const& application) =0;
  virtual void addPostProcessor(PostProcessorGenerator3D<T,Lattice> const& ppGen) =0;
  virtual void resetPostProcessors() =0;
  virtual void postProcess(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_) =0;
  virtual void postProcess() =0;
  virtual void addLatticeCoupling(LatticeCouplingGenerator3D<T,Lattice> const& lcGen,
                                  std::vector<SpatiallyExtendedObject3D*> partners ) =0;
  virtual void executeCoupling(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_) =0;
  virtual void executeCoupling() =0;
  virtual void subscribeReductions(Reductor<T>& reductor) =0;
  virtual LatticeStatistics<T>& getStatistics() =0;
  virtual LatticeStatistics<T> const& getStatistics() const =0;
};

}  // namespace olb

#endif
