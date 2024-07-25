/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Claudius Holeksa
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

#pragma once

#include "dynamics/descriptorField.h"
#include "dynamics/freeSurfaceHelpers.h"
#include "core/postProcessing.h"
#include "core/blockLattice.h"

namespace olb {
class FreeSurface3D {
public:

  template<typename T, typename DESCRIPTOR>
  static bool isCellType(BlockLattice<T, DESCRIPTOR>& blockLattice, int X, int Y, int Z, const FreeSurface::Type& type);

  template<typename T, typename DESCRIPTOR>
  static bool hasCellFlags(BlockLattice<T, DESCRIPTOR>& blockLattice, int X, int Y, int Z, const FreeSurface::Flags& flags);

  template<typename T, typename DESCRIPTOR>
  static bool hasNeighbour (BlockLattice<T,DESCRIPTOR>& blockLattice, int X, int Y, int Z, const FreeSurface::Type& type);
  
  template<typename T, typename DESCRIPTOR>
  static bool hasNeighbourFlags (BlockLattice<T,DESCRIPTOR>& blockLattice, int X, int Y, int Z, const FreeSurface::Flags& flags);

  template<typename T, typename DESCRIPTOR>
  static std::array<T,DESCRIPTOR::d> computeInterfaceNormal(BlockLattice<T,DESCRIPTOR>& blockLattice, int X, int Y, int Z);


  template<typename T, typename DESCRIPTOR>
  static std::array<T,DESCRIPTOR::d> computeParkerYoungInterfaceNormal(BlockLattice<T,DESCRIPTOR>& blockLattice, int X, int Y, int Z);

  template<typename T, size_t S>
  static std::array<T,S> solvePivotedLU(std::array<std::array<T,S>,S>& matrix, const std::array<T,S>& b, size_t N = S);

  template<typename T, typename DESCRIPTOR>
  static T calculateCubeOffset(T volume, const std::array<T,DESCRIPTOR::d>& normal);
  template<typename T, typename DESCRIPTOR>
  static T calculateCubeOffsetOpt(T volume, const std::array<T,DESCRIPTOR::d>& normal);

  template<typename T, typename DESCRIPTOR>
  static T getClampedEpsilon(BlockLattice<T,DESCRIPTOR>& blockLattice, int X, int Y, int Z);

  template<typename T, typename DESCRIPTOR>
  static T getClampedEpsilonCorrected(BlockLattice<T,DESCRIPTOR>& blockLattice, int X, int Y, int Z);

  template<typename T, typename DESCRIPTOR>
  static T calculateSurfaceTensionCurvature(BlockLattice<T,DESCRIPTOR>& blockLattice, int X, int Y, int Z);

  template<typename T, typename DESCRIPTOR>
  static T plicInverse(T d, const std::array<T,DESCRIPTOR::d>& normal);

  template<typename T, typename DESCRIPTOR>
  static bool isHealthyInterface(BlockLattice<T,DESCRIPTOR>& blockLattice, int X, int Y, int Z);


  template<typename T, typename DESCRIPTOR>
  static void setCellType(BlockLattice<T,DESCRIPTOR>& blockLattice, int X, int Y, int Z, const FreeSurface::Type& type);

  template<typename T, typename DESCRIPTOR>
  static void setCellFlags(BlockLattice<T,DESCRIPTOR>& blockLattice, int X, int Y, int Z, const FreeSurface::Flags& flags);

  struct NeighbourInfo {
    bool has_fluid_neighbours = false;
    bool has_gas_neighbours = false;
    size_t interface_neighbours = 0;
  };

  template<typename T, typename DESCRIPTOR>
  static NeighbourInfo getNeighbourInfo(BlockLattice<T,DESCRIPTOR>& blockLattice, int X, int Y, int Z);

  /**
   * Helper class meant to make the necessary variables more accesible and interchangeable
   */
  template<typename T, typename DESCRIPTOR>
  struct Variables {
    static_assert(DESCRIPTOR::d == 3, "Descriptor doesn't fit the helper class it is in");

    bool drop_isolated_cells;
    T transition;
    T lonely_threshold;

    bool has_surface_tension = true;
    T surface_tension_parameter = 0.1;
    T force_conversion_factor = 1.0;

    T lattice_size = 1.0;
  };
};

/**
 * Free Surface Processor 1-3
 * Mass Flow 
 * Cleans up leftover flags from the previous simulation step.
 * This post processor is responsible for the calculation of exchange mass with the help of the distribution functions.
 * Replaces incoming DFs by calculating equilibrium functions and using the laplace pressure to include surface tension.
 * Marks cells which may be changed at the last step.
 * This whole step should be included in the collideAndStream step, though heavy modification of openlb would be necessary
 */
template<typename T, typename DESCRIPTOR>
class FreeSurfaceMassFlowPostProcessor3D final : public LocalPostProcessor3D<T, DESCRIPTOR> {
public:
  FreeSurfaceMassFlowPostProcessor3D(const FreeSurface3D::Variables<T,DESCRIPTOR>& vars_);
  FreeSurfaceMassFlowPostProcessor3D(int x0, int x1, int y0, int y1, int z0, int z1, const FreeSurface3D::Variables<T,DESCRIPTOR>& vars_);

  int extent() const override
  {
    return 1;
  }

  int extent(int whichDirection) const override
  {
    (void)whichDirection;
    return 1;
  }

  void process(BlockLattice<T,DESCRIPTOR>& blockLattice) override;
  void processSubDomain(BlockLattice<T,DESCRIPTOR>& blockLattice,
                        int x0_, int x1_, int y0_, int y1_, int z0, int z1) override;
private:
  int x0, x1, y0, y1, z0, z1;
  const FreeSurface3D::Variables<T,DESCRIPTOR>& vars;
};

/*
 * Free Surface Processor 4
 * ToFluid  
 * Converts cells to interface from gas if a neighbouring cell was converted to a fluid cell
 */
template<typename T, typename DESCRIPTOR>
class FreeSurfaceToFluidCellConversionPostProcessor3D final : public LocalPostProcessor3D<T, DESCRIPTOR> {
public:
  FreeSurfaceToFluidCellConversionPostProcessor3D(const FreeSurface3D::Variables<T,DESCRIPTOR>& vars_);
  FreeSurfaceToFluidCellConversionPostProcessor3D(int x0, int x1, int y0, int y1, int z0, int z1, const FreeSurface3D::Variables<T,DESCRIPTOR>& vars_);

  int extent() const override
  {
    return 1;
  }
  int extent(int whichDirection) const override
  {
    (void)whichDirection;
    return 1;
  }
  void process(BlockLattice<T,DESCRIPTOR>& blockLattice) override;
  void processSubDomain(BlockLattice<T,DESCRIPTOR>& blockLattice,
                        int x0_, int x1_, int y0_, int y1_, int z0, int z1) override;
private:
  int x0, x1, y0, y1, z0, z1;
  const FreeSurface3D::Variables<T,DESCRIPTOR>& vars;
};

/**
 * Free Surface Processor 5
 * ToGas
 * Converts cells to interface from fluid if a neighbouring cell was converted to a gas cell
 */
template<typename T, typename DESCRIPTOR>
class FreeSurfaceToGasCellConversionPostProcessor3D final : public LocalPostProcessor3D<T, DESCRIPTOR> {
public:
  FreeSurfaceToGasCellConversionPostProcessor3D(const FreeSurface3D::Variables<T,DESCRIPTOR>& vars_);
  FreeSurfaceToGasCellConversionPostProcessor3D(int x0, int x1, int y0, int y1, int z0, int z1, const FreeSurface3D::Variables<T,DESCRIPTOR>& vars_);

  int extent() const override
  {
    return 1;
  }
  int extent(int whichDirection) const override
  {
    (void)whichDirection;
    return 1;
  }
  void process(BlockLattice<T,DESCRIPTOR>& blockLattice) override;
  void processSubDomain(BlockLattice<T,DESCRIPTOR>& blockLattice,
                        int x0_, int x1_, int y0_, int y1_, int z0, int z1) override;
private:
  int x0, x1, y0, y1, z0, z1;
  const FreeSurface3D::Variables<T,DESCRIPTOR>& vars;
};

/**
 * Free Surface Processor 6
 * Calculates mass excess from the cell type conversions and distributes them to neighbouring interface cells
 * Keeps mass local if no neighbour exists until an interface reappears at this position
 */
template<typename T, typename DESCRIPTOR>
class FreeSurfaceMassExcessPostProcessor3D final : public LocalPostProcessor3D<T, DESCRIPTOR> {
public:
  FreeSurfaceMassExcessPostProcessor3D(const FreeSurface3D::Variables<T,DESCRIPTOR>& vars_);
  FreeSurfaceMassExcessPostProcessor3D(int x0, int x1, int y0, int y1, int z0, int z1, const FreeSurface3D::Variables<T,DESCRIPTOR>& vars_);

  int extent() const override
  {
    return 1;
  }
  int extent(int whichDirection) const override
  {
    (void)whichDirection;
    return 1;
  }
  void process(BlockLattice<T,DESCRIPTOR>& blockLattice) override;
  void processSubDomain(BlockLattice<T,DESCRIPTOR>& blockLattice,
                        int x0_, int x1_, int y0_, int y1_, int z0_, int z1_) override;
private:
  int x0, x1, y0, y1, z0, z1;
  const FreeSurface3D::Variables<T,DESCRIPTOR>& vars;
};

/**
 * Free Surface Processor 7
 * Finishes up left over cell conversions and prepares the state for the next simulation step
 */
template<typename T, typename DESCRIPTOR>
class FreeSurfaceFinalizeConversionPostProcessor3D final : public LocalPostProcessor3D<T, DESCRIPTOR> {
public:
  FreeSurfaceFinalizeConversionPostProcessor3D(const FreeSurface3D::Variables<T,DESCRIPTOR>& vars_);
  FreeSurfaceFinalizeConversionPostProcessor3D(int x0, int x1, int y0, int y1, int z0, int z1, const FreeSurface3D::Variables<T,DESCRIPTOR>& vars_);

  int extent() const override
  {
    return 1;
  }
  int extent(int whichDirection) const override
  {
    (void)whichDirection;
    return 1;
  }
  void process(BlockLattice<T,DESCRIPTOR>& blockLattice) override;
  void processSubDomain(BlockLattice<T,DESCRIPTOR>& blockLattice,
                        int x0_, int x1_, int y0_, int y1_, int z0, int z1) override;
private:
  int x0, x1, y0, y1, z0, z1;
  const FreeSurface3D::Variables<T,DESCRIPTOR>& vars;
};


/// Generator class for the PostProcessors tracking the interface.

/*
* Generator 1, 2, 3
*/
template<typename T, typename DESCRIPTOR>
class FreeSurfaceMassFlowGenerator3D final : public PostProcessorGenerator3D<T,DESCRIPTOR> {
public:
  /// \param[in] alpha_ - dummy parameter. [lattice units]
  FreeSurfaceMassFlowGenerator3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, const FreeSurface3D::Variables<T,DESCRIPTOR>& vars_);
  /// \param[in] alpha_ - dummy parameter. [lattice units]
  FreeSurfaceMassFlowGenerator3D(const FreeSurface3D::Variables<T,DESCRIPTOR>& vars_);
  PostProcessor3D<T,DESCRIPTOR>* generate() const override;
  PostProcessorGenerator3D<T,DESCRIPTOR>* clone() const override;
private:
  const FreeSurface3D::Variables<T,DESCRIPTOR>& vars;
};
/*
* Generator 4
*/
template<typename T, typename DESCRIPTOR>
class FreeSurfaceToFluidCellConversionGenerator3D final : public PostProcessorGenerator3D<T,DESCRIPTOR> {
public:
  /// \param[in] alpha_ - dummy parameter. [lattice units]
  FreeSurfaceToFluidCellConversionGenerator3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, const FreeSurface3D::Variables<T,DESCRIPTOR>& vars_);
  /// \param[in] alpha_ - dummy parameter. [lattice units]
  FreeSurfaceToFluidCellConversionGenerator3D(const FreeSurface3D::Variables<T,DESCRIPTOR>& vars_);
  PostProcessor3D<T,DESCRIPTOR>* generate() const override;
  PostProcessorGenerator3D<T,DESCRIPTOR>* clone() const override;
private:
  const FreeSurface3D::Variables<T,DESCRIPTOR>& vars;
};

/*
* Generator 5
*/
template<typename T, typename DESCRIPTOR>
class FreeSurfaceToGasCellConversionGenerator3D final : public PostProcessorGenerator3D<T,DESCRIPTOR> {
public:
  /// \param[in] alpha_ - dummy parameter. [lattice units]
  FreeSurfaceToGasCellConversionGenerator3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, const FreeSurface3D::Variables<T,DESCRIPTOR>& vars_);
  /// \param[in] alpha_ - dummy parameter. [lattice units]
  FreeSurfaceToGasCellConversionGenerator3D(const FreeSurface3D::Variables<T,DESCRIPTOR>& vars_);
  PostProcessor3D<T,DESCRIPTOR>* generate() const override;
  PostProcessorGenerator3D<T,DESCRIPTOR>* clone() const override;
private:
  const FreeSurface3D::Variables<T,DESCRIPTOR>& vars;
};

/*
* Generator 6
*/
template<typename T, typename DESCRIPTOR>
class FreeSurfaceMassExcessGenerator3D final : public PostProcessorGenerator3D<T,DESCRIPTOR> {
public:
  /// \param[in] alpha_ - dummy parameter. [lattice units]
  FreeSurfaceMassExcessGenerator3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, const FreeSurface3D::Variables<T,DESCRIPTOR>& vars_);
  /// \param[in] alpha_ - dummy parameter. [lattice units]
  FreeSurfaceMassExcessGenerator3D(const FreeSurface3D::Variables<T,DESCRIPTOR>& vars_);
  PostProcessor3D<T,DESCRIPTOR>* generate() const override;
  PostProcessorGenerator3D<T,DESCRIPTOR>* clone() const override;
private:
  const FreeSurface3D::Variables<T,DESCRIPTOR>& vars;
};

/*
* Generator 7
*/
template<typename T, typename DESCRIPTOR>
class FreeSurfaceFinalizeConversionGenerator3D final : public PostProcessorGenerator3D<T,DESCRIPTOR> {
public:
  /// \param[in] alpha_ - dummy parameter. [lattice units]
  FreeSurfaceFinalizeConversionGenerator3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, const FreeSurface3D::Variables<T,DESCRIPTOR>& vars_);
  /// \param[in] alpha_ - dummy parameter. [lattice units]
  FreeSurfaceFinalizeConversionGenerator3D(const FreeSurface3D::Variables<T,DESCRIPTOR>& vars_);
  PostProcessor3D<T,DESCRIPTOR>* generate() const override;
  PostProcessorGenerator3D<T,DESCRIPTOR>* clone() const override;
private:
  const FreeSurface3D::Variables<T,DESCRIPTOR>& vars;
};

/*
* Setup helper
*/

template<typename T, typename DESCRIPTOR>
class FreeSurface3DSetup {
public:
private:
  const FreeSurface3D::Variables<T,DESCRIPTOR> vars;
  SuperLattice<T, DESCRIPTOR>& sLattice;

  // BlockPostProcessors
  // Step 1
  std::unique_ptr<FreeSurfaceMassFlowGenerator3D<T,DESCRIPTOR>> mass_flow;
  // Step 2
  //std::unique_ptr<FreeSurfaceInterfaceGenerator3D<T,DESCRIPTOR>> interface;
  // Step 3
  //std::unique_ptr<FreeSurfaceCellConversionPreparationGenerator3D<T,DESCRIPTOR>> cell_conversion_prep;
  // Step 4
  std::unique_ptr<FreeSurfaceToFluidCellConversionGenerator3D<T,DESCRIPTOR>> to_fluid;
  // Step 5
  std::unique_ptr<FreeSurfaceToGasCellConversionGenerator3D<T,DESCRIPTOR>> to_gas;
  // Step 6
  std::unique_ptr<FreeSurfaceMassExcessGenerator3D<T,DESCRIPTOR>> mass_excess;
  // Step 7
  std::unique_ptr<FreeSurfaceFinalizeConversionGenerator3D<T,DESCRIPTOR>> finalize_conversion;

  // SuperPostProcessors
  // Corresponding to the local block processors
public:
  FreeSurface3DSetup(
    SuperLattice<T, DESCRIPTOR>& sLattice,
    const FreeSurface3D::Variables<T,DESCRIPTOR>& vars);

  void addPostProcessor();
};

}
