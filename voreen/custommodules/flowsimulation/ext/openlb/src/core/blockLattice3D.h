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
 * The dynamics of a 3D block lattice -- header file.
 */
#ifndef BLOCK_LATTICE_3D_H
#define BLOCK_LATTICE_3D_H

#include <vector>
#include "olbDebug.h"
#include "core/cell.h"
#include "postProcessing.h"
#include "blockLatticeStructure3D.h"
#include "geometry/blockGeometry3D.h"
#include "blockStaticPopulationD3D.h"
#include "blockStaticFieldsD3D.h"
#include "blockDynamicFieldsD3D.h"
#include "blockDynamicsMap.h"
#include "latticeStatistics.h"
#include "serializer.h"


namespace olb {

template<typename T> class BlockGeometryStructure3D;
template<typename T> class BlockIndicatorF3D;
template<typename T, typename DESCRIPTOR> struct Dynamics;


/// Regular lattice for highly efficient 3D LB dynamics.
template<typename T, typename DESCRIPTOR>
class BlockLattice3D final : public BlockLatticeStructure3D<T,DESCRIPTOR>, public Serializable {
private:
  /// Population data of a lattice with DESCRIPTOR structure
  BlockStaticPopulationD3D<T,DESCRIPTOR> _staticPopulationD;
  /// Field data of DESCRIPTOR provided fields
  BlockStaticFieldsD3D<T,DESCRIPTOR>     _staticFieldsD;
  /// Field data for runtime-declared fields not in DESCRIPTOR
  BlockDynamicFieldsD3D<T,DESCRIPTOR>    _dynamicFieldsD;
  /// Assignments of dynamics instances to spatial cell locations
  BlockDynamicsMap<T,DESCRIPTOR> _dynamicsMap;

  /// List of all post processors to be applied after streaming
  std::vector<PostProcessor3D<T,DESCRIPTOR>*> _postProcessors;
  /// List of all coupling post processors
  std::vector<PostProcessor3D<T,DESCRIPTOR>*> _latticeCouplings;

#ifdef PARALLEL_MODE_OMP
  LatticeStatistics<T> **statistics;
#else
  LatticeStatistics<T> *statistics;
#endif

public:
  /// Construction of an nx_ by ny_ by nz_ lattice
  BlockLattice3D(int nx_, int ny_, int nz_);
  /// Destruction of the lattice
  ~BlockLattice3D() override;
  /// Copy construction
  BlockLattice3D(BlockLattice3D<T,DESCRIPTOR> const& rhs) = delete;
  /// Move constructor
  BlockLattice3D(BlockLattice3D<T,DESCRIPTOR>&&) = default;
  /// Copy assignment
  BlockLattice3D& operator=(BlockLattice3D<T,DESCRIPTOR> const& rhs) = delete;

  BlockStaticPopulationD3D<T,DESCRIPTOR>& getStaticPopulationD() {
    return _staticPopulationD;
  }
  BlockStaticFieldsD3D<T,DESCRIPTOR>& getStaticFieldsD() {
    return _staticFieldsD;
  }
  BlockDynamicFieldsD3D<T,DESCRIPTOR>& getDynamicFieldsD() {
    return _dynamicFieldsD;
  }

  template<typename FIELD>
  FieldArrayD<T, DESCRIPTOR, FIELD>& getDynamicFieldArray() {
    return _dynamicFieldsD.template getFieldArray<FIELD>();
  }

  Cell<T,DESCRIPTOR> get(std::size_t iCell)
  {
    return Cell<T,DESCRIPTOR>(
      _staticPopulationD, _staticFieldsD, _dynamicFieldsD, _dynamicsMap, iCell);
  }
  /// Read/write access to lattice cells
  Cell<T,DESCRIPTOR> get(int iX, int iY, int iZ) override
  {
    return get(this->getCellId(iX,iY,iZ));
  }
  /// Read/write access to lattice cells
  Cell<T,DESCRIPTOR> get(const int latticeR[]) override
  {
    return get(latticeR[0], latticeR[1], latticeR[2]);
  }
  /// Read only access to lattice cells
  ConstCell<T,DESCRIPTOR> get(int iX, int iY, int iZ) const override
  {
  	return ConstCell<T,DESCRIPTOR>(
  	      _staticPopulationD, _staticFieldsD, _dynamicFieldsD, _dynamicsMap, this->getCellId(iX,iY,iZ));
  }

  T& getPop(std::size_t iCell, unsigned iPop) override {
    return *_staticPopulationD.getPopulationPointer(iPop, iCell);
  }

  T& getPop(int iX, int iY, int iZ, unsigned iPop) override {
    return *_staticPopulationD.getPopulationPointer(iPop, this->getCellId(iX, iY, iZ));
  }

  /// Initialize the lattice cells to become ready for simulation
  void initialize() override;

  /// Get the dynamics on a lattice site
  Dynamics<T,DESCRIPTOR>* getDynamics(int iX, int iY, int iZ) override;

  /// Define the dynamics on a lattice site
  void defineDynamics(int iX, int iY, int iZ, Dynamics<T,DESCRIPTOR>* dynamics) override;
  /// Define the dynamics on a 3D sub-box
  void defineDynamics(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
                      Dynamics<T,DESCRIPTOR>* dynamics ) override;
  /// Define the dynamics on a domain described by an indicator
  void defineDynamics(BlockIndicatorF3D<T>& indicator,
                      Dynamics<T,DESCRIPTOR>* dynamics) override;
  /// Define the dynamics on a domain with a particular material number
  void defineDynamics(BlockGeometryStructure3D<T>& blockGeometry, int material,
                      Dynamics<T,DESCRIPTOR>* dynamics) override;

  /// Apply collision step to a 3D sub-box
  void collide(int x0_, int x1_, int y0_, int y1_,
               int z0_, int z1_) override;
  /// Perform streaming on the whole domain
  void stream();
  /// Apply collide and stream step to a 3D sub-box
  void collideAndStream(int x0_, int x1_, int y0_, int y1_,
                        int z0_, int z1_) override;
  /// Apply collision step to the whole domain
  void collide() override;
  /// Apply collide and stream to the whole domain
  void collideAndStream() override;

  /// Compute the average density within a rectangular domain
  T computeAverageDensity(int x0_, int x1_, int y0_, int y1_,
                          int z0_, int z1_ ) const override;
  /// Compute the average density within the whole domain
  T computeAverageDensity() const override;
  /// Compute components of the stress tensor on the cell.
  void computeStress(int iX, int iY, int iZ, T pi[util::TensorVal<DESCRIPTOR >::n]) override;
  /// Subtract a constant offset from the density within the whole domain
  void stripeOffDensityOffset (
    int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, T offset ) override;
  /// Subtract a constant offset from the density within a rect. domain
  void stripeOffDensityOffset(T offset) override;
  /// Apply an operation to all cells of a sub-domain
  void forAll(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
              WriteCellFunctional<T,DESCRIPTOR> const& application) override;
  /// Apply an operation to all cells
  void forAll(WriteCellFunctional<T,DESCRIPTOR> const& application) override;
  /// Add a non-local post-processing step
  void addPostProcessor (
    PostProcessorGenerator3D<T,DESCRIPTOR> const& ppGen ) override;
  /// Clean up all non-local post-processing steps
  void resetPostProcessors() override;
  /// Execute post-processing on a sub-lattice
  void postProcess(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_) override;
  /// Execute post-processing steps
  void postProcess() override;
  /// Add a non-local post-processing step
  void addLatticeCoupling (
    LatticeCouplingGenerator3D<T,DESCRIPTOR> const& lcGen,
    std::vector<SpatiallyExtendedObject3D*> partners ) override;
  /// Execute couplings on a sub-lattice
  void executeCoupling(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_) override;
  /// Execute couplings steps
  void executeCoupling() override;
  /// Return a handle to the LatticeStatistics object
  LatticeStatistics<T>& getStatistics() override;
  /// Return a constant handle to the LatticeStatistics object
  LatticeStatistics<T> const& getStatistics() const override;

  /// Apply streaming step to bulk (non-boundary) cells
  void bulkStream(int x0, int x1, int y0, int y1, int z0, int z1);
  /// Apply streaming step to boundary cells
  void boundaryStream (
    int lim_x0, int lim_x1, int lim_y0, int lim_y1,
    int lim_z0, int lim_z1,
    int x0, int x1, int y0, int y1, int z0, int z1 );

  /// Number of data blocks for the serializable interface
  std::size_t getNblock() const override;
  /// Binary size for the serializer
  std::size_t getSerializableSize() const override;
  /// Return a pointer to the memory of the current block and its size for the serializable interface
  bool* getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode) override;

private:

  /// Release memory for post processors
  void clearPostProcessors();
  /// Release memory for post processors
  void clearLatticeCouplings();
  /// Make the lattice periodic in all directions
  void makePeriodic();

  void periodicSurface(int x0, int x1, int y0, int y1, int z0, int z1);

};

}  // namespace olb

#endif
