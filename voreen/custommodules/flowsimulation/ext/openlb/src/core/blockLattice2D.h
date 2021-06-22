/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006-2008 Jonas Latt
 *                2008-2020 Mathias Krause
 *                2020 Adrian Kummerlaender
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
 * The dynamics of a 2D block lattice -- header file.
 */
#ifndef BLOCK_LATTICE_2D_H
#define BLOCK_LATTICE_2D_H

#include <vector>
#include <forward_list>

#include "olbDebug.h"
#include "postProcessing.h"
#include "blockLatticeStructure2D.h"
#include "core/cell.h"
#include "latticeStatistics.h"
#include "serializer.h"
#include "blockStaticPopulationD2D.h"
#include "blockStaticFieldsD2D.h"
#include "blockDynamicFieldsD2D.h"
#include "blockDynamicsMap.h"


namespace olb {

template<typename T> class BlockGeometryStructure2D;
template<typename T, typename DESCRIPTOR> struct Dynamics;
template<typename T> class BlockIndicatorF2D;

/// Regular lattice for highly efficient 2D LB dynamics.
template<typename T, typename DESCRIPTOR>
class BlockLattice2D final : public BlockLatticeStructure2D<T,DESCRIPTOR>, public Serializable {
private:
  /// Population data of a lattice with DESCRIPTOR structure
  BlockStaticPopulationD2D<T,DESCRIPTOR> _staticPopulationD;
  /// Field data of DESCRIPTOR provided fields
  BlockStaticFieldsD2D<T,DESCRIPTOR>     _staticFieldsD;
  /// Field data for runtime-declared fields not in DESCRIPTOR
  BlockDynamicFieldsD2D<T,DESCRIPTOR>    _dynamicFieldsD;
  /// Assignments of dynamics instances to spatial cell locations
  BlockDynamicsMap<T,DESCRIPTOR> _dynamicsMap;

  /// List of all post processors to be applied after streaming
  std::vector<PostProcessor2D<T,DESCRIPTOR>*> _postProcessors;
  /// List of all coupling post processors
  std::vector<PostProcessor2D<T,DESCRIPTOR>*> _latticeCouplings;

#ifdef PARALLEL_MODE_OMP
  LatticeStatistics<T> **_statistics;
#else
  LatticeStatistics<T> *_statistics;
#endif

public:
  /// Construction of an nx_ by ny_ lattice
  BlockLattice2D(int nx, int ny);
  /// Destruction of the lattice
  ~BlockLattice2D() override;
  /// Copy construction
  BlockLattice2D(BlockLattice2D<T,DESCRIPTOR> const& rhs) = delete;
  // Move constructor
  BlockLattice2D(BlockLattice2D&&) = default;
  /// Copy assignment
  BlockLattice2D& operator=(BlockLattice2D<T,DESCRIPTOR> const& rhs) = delete;

  BlockStaticPopulationD2D<T,DESCRIPTOR>& getStaticPopulationD() {
    return _staticPopulationD;
  }
  BlockStaticFieldsD2D<T,DESCRIPTOR>& getStaticFieldsD() {
    return _staticFieldsD;
  }
  BlockDynamicFieldsD2D<T,DESCRIPTOR>& getDynamicFieldsD() {
    return _dynamicFieldsD;
  }

  template<typename FIELD>
  FieldArrayD<T,DESCRIPTOR,FIELD>& getDynamicFieldArray() {
    return _dynamicFieldsD.template getFieldArray<FIELD>();
  }

  Cell<T,DESCRIPTOR> get(std::size_t iCell)
  {
    return Cell<T,DESCRIPTOR>(
      _staticPopulationD, _staticFieldsD, _dynamicFieldsD, _dynamicsMap, iCell);
  }

  /// Read/write access to lattice cells
  Cell<T,DESCRIPTOR> get(int iX, int iY) override
  {
    return get(this->getCellId(iX,iY));
  }
  /// Read/write access to lattice cells
  Cell<T,DESCRIPTOR> get(int latticeR[]) override
  {
    return get(latticeR[0], latticeR[1]);
  }
  /// Read only access to lattice cells
  ConstCell<T,DESCRIPTOR> get(int iX, int iY) const override
  {
    return ConstCell<T,DESCRIPTOR>(
      _staticPopulationD, _staticFieldsD, _dynamicFieldsD, _dynamicsMap, this->getCellId(iX,iY));
  }

  T& getPop(std::size_t iCell, unsigned iPop) override {
    return *_staticPopulationD.getPopulationPointer(iPop, iCell);
  }

  T& getPop(int iX, int iY, unsigned iPop) override {
    return *_staticPopulationD.getPopulationPointer(iPop, this->getCellId(iX, iY));
  }

  /// Initialize the lattice cells to become ready for simulation
  void initialize() override;
  /// Get the dynamics for the specific point
  Dynamics<T,DESCRIPTOR>* getDynamics(int iX, int iY) override;
  /// Define the dynamics on a rectangular domain
  void defineDynamics (int x0, int x1, int y0, int y1,
                       Dynamics<T,DESCRIPTOR>* dynamics ) override;
  /// Define the dynamics on a lattice site
  void defineDynamics(int iX, int iY, Dynamics<T,DESCRIPTOR>* dynamics) override;
  /// Define the dynamics on a domain described by an indicator
  void defineDynamics(BlockIndicatorF2D<T>& indicator,
                      Dynamics<T,DESCRIPTOR>* dynamics);
  /// Define the dynamics on a domain described by a material number
  void defineDynamics(BlockGeometryStructure2D<T>& blockGeometry, int material,
                      Dynamics<T,DESCRIPTOR>* dynamics);

  /// Apply collision step to a rectangular domain
  void collide(int x0, int x1, int y0, int y1) override;
  /// Perform streaming on the whole domain
  void stream();
  /// Apply collide and stream on a rectangular domain
  void collideAndStream(int x0, int x1, int y0, int y1) override;
  /// Apply collision step to the whole domain
  void collide() override;
  /// Apply collide and stream to the whole domain
  void collideAndStream() override;

  /// Compute the average density within a rectangular domain
  T computeAverageDensity(int x0, int x1, int y0, int y1) const override;
  /// Compute the average density within the whole domain
  T computeAverageDensity() const override;
  /// Compute components of the stress tensor on the cell.
  void computeStress(int iX, int iY, T pi[util::TensorVal<DESCRIPTOR >::n]) override;
  /// Subtract a constant offset from the density within the whole domain
  void stripeOffDensityOffset (
    int x0, int x1, int y0, int y1, T offset ) override;
  /// Subtract a constant offset from the density within a rect. domain
  void stripeOffDensityOffset(T offset) override;
  /// Apply an operation to all cells of a sub-domain
  void forAll(int x0_, int x1_, int y0_, int y1_,
              WriteCellFunctional<T,DESCRIPTOR> const& application) override;
  /// Apply an operation to all cells
  void forAll(WriteCellFunctional<T,DESCRIPTOR> const& application) override;
  /// Add a non-local post-processing step
  void addPostProcessor (    PostProcessorGenerator2D<T,DESCRIPTOR> const& ppGen ) override;
  /// Clean up all non-local post-processing steps
  void resetPostProcessors() override;
  /// Execute post-processing on a sub-lattice
  void postProcess(int x0_, int x1_, int y0_, int y1_) override;
  /// Execute post-processing steps
  void postProcess() override;
  /// Add a non-local post-processing step which couples together lattices
  void addLatticeCoupling( LatticeCouplingGenerator2D<T,DESCRIPTOR> const& lcGen,
                           std::vector<SpatiallyExtendedObject2D*> partners ) override;
  /// Execute couplings on a sub-lattice
  void executeCoupling(int x0_, int x1_, int y0_, int y1_) override;
  /// Execute couplings
  void executeCoupling() override;
  /// Return a handle to the LatticeStatistics object
  LatticeStatistics<T>& getStatistics() override;
  /// Return a constant handle to the LatticeStatistics object
  LatticeStatistics<T> const& getStatistics() const override;

public:
  /// Number of data blocks for the serializable interface
  std::size_t getNblock() const override;
  /// Binary size for the serializer
  std::size_t getSerializableSize() const override;
  /// Return a pointer to the memory of the current block and its size for the serializable interface
  bool* getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode) override;

private:
  /// Release memory for post processors
  void clearPostProcessors();
  /// Release memory for lattice couplings
  void clearLatticeCouplings();
  void periodicEdge(int x0, int x1, int y0, int y1);
  void makePeriodic();
};

}  // namespace olb

#endif
