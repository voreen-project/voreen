/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Jonas Latt
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
 * Dynamics for a generic 2D block lattice view -- header file.
 */
#ifndef BLOCK_LATTICE_VIEW_2D_H
#define BLOCK_LATTICE_VIEW_2D_H

#include <vector>
#include "blockLatticeStructure2D.h"
#include "geometry/blockGeometry2D.h"
#include "geometry/blockGeometryStatistics2D.h"

namespace olb {

/// A rectangular extract from a given BlockLatticeStructure.
/** Attention:
 * - This class can only be applied to already existing BlockLattices.
 * - The postProcessors of the original BlockLattice are not called any
 *   more. New, appropriate PostProcessors are generated instead.
 */
template<typename T, template<typename U> class Lattice>
class BlockLatticeView2D : public BlockLatticeStructure2D<T,Lattice> {
public:
  BlockLatticeView2D(BlockLatticeStructure2D<T,Lattice>& originalLattice_);
  BlockLatticeView2D(BlockLatticeStructure2D<T,Lattice>& originalLattice_,
                     int x0_, int x1_, int y0_, int y1_);
  ~BlockLatticeView2D() override;
  BlockLatticeView2D(BlockLatticeView2D const& rhs);
  BlockLatticeView2D<T,Lattice>& operator=
  (BlockLatticeView2D<T,Lattice> const& rhs);
  void swap(BlockLatticeView2D<T,Lattice>& rhs);

  Cell<T,Lattice>& get(int iX, int iY) override;
  Cell<T,Lattice> const& get(int iX, int iY) const override;
  void initialize() override;
  void defineDynamics (
    int x0_, int x1_, int y0_, int y1_,
    Dynamics<T,Lattice>* dynamics ) override;
  void defineDynamics(int iX, int iY, Dynamics<T,Lattice>* dynamics ) override;
  Dynamics<T,Lattice>* getDynamics(int iX, int iY) override;
  void collide(int x0_, int x1_, int y0_, int y1_) override;
  void collide() override;
  /*virtual void staticCollide (int x0, int x1, int y0, int y1,
                              TensorFieldBase2D<T,2> const& u);
  virtual void staticCollide (TensorFieldBase2D<T,2> const& u);*/
  void stream(int x0_, int x1_, int y0_, int y1_) override;
  void collideAndStream(int x0_, int x1_, int y0_, int y1_) override;
  void stream(bool periodic=false) override;
  void collideAndStream(bool periodic=false) override;
  T computeAverageDensity(int x0_, int x1_, int y0_, int y1_) const override;
  T computeAverageDensity() const override;
  void computeStress(int iX, int iY, T pi[util::TensorVal<Lattice<T> >::n]) override;
  void stripeOffDensityOffset (
    int x0_, int x1_, int y0_, int y1_, T offset ) override;
  void stripeOffDensityOffset(T offset) override;
  void forAll(int x0_, int x1_, int y0_, int y1_,
                      WriteCellFunctional<T,Lattice> const& application) override;
  void forAll(WriteCellFunctional<T,Lattice> const& application) override;
  void addPostProcessor (
    PostProcessorGenerator2D<T,Lattice> const& ppGen) override;
  void resetPostProcessors() override;
  void postProcess(int x0_, int x1_, int y0_, int y1_) override;
  void postProcess() override;
  void addLatticeCoupling (
    LatticeCouplingGenerator2D<T,Lattice> const& lcGen,
    std::vector<SpatiallyExtendedObject2D*> partners ) override;
  void executeCoupling(int x0_, int x1_, int y0_, int y1_) override;
  void executeCoupling() override;
  void subscribeReductions(Reductor<T>& reductor) override;
  LatticeStatistics<T>& getStatistics() override;
  LatticeStatistics<T> const& getStatistics() const override;
private:
  BlockLatticeStructure2D<T,Lattice>  *originalLattice;
  int                          x0, y0;
};

}  // namespace olb

#endif
