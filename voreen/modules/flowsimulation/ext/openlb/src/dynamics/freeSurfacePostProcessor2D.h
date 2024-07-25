/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2020 Claudius Holeksa, Robin Trunk
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

#ifndef FREE_SURFACE_POST_PROCESSOR_2D_H
#define FREE_SURFACE_POST_PROCESSOR_2D_H

#include "dynamics/freeSurfaceHelpers.h"
#include "core/postProcessing.h"
#include "core/blockLattice.h"
#include "core/superLattice.h"

#include <array>
#include <memory>

/* \file
 * PostProcessor classes organising the interface tracking and mass distibution for a
 * free surface model.
 *
 * Description how PostProcessors are applied.
 */

namespace olb {
/**
 * Free surface 2D helper. Since the data structure and the necessary functions are quite complex, we provide
 * a class which makes the
 */
class FreeSurface2D {
public:

  template <typename CELL>
  static bool isCellType(CELL& cell, const FreeSurface::Type& type) any_platform;

  template <typename CELL>
  static bool hasCellFlags(CELL& cell, const FreeSurface::Flags& type) any_platform;

  template <typename CELL>
  static bool hasNeighbour(CELL& cell, const FreeSurface::Type& type) any_platform;

  template <typename CELL>
  static bool hasNeighbourFlags(CELL& cell, const FreeSurface::Flags& flags) any_platform;

  template <typename CELL, typename V=typename CELL::value_t>
  static V getClampedEpsilon(CELL& cell) any_platform;

  template <typename CELL, typename V=typename CELL::value_t>
  static Vector<V,CELL::descriptor_t::d> computeInterfaceNormal(CELL& cell) any_platform;

  template <typename CELL, typename V=typename CELL::value_t>
  static Vector<V,CELL::descriptor_t::d> computeParkerYoungInterfaceNormal(CELL& cell) any_platform;

  template<typename T, typename DESCRIPTOR>
  static T calculateCubeOffset(T volume, const Vector<T,DESCRIPTOR::d>& normal) any_platform;

  template<typename T, typename DESCRIPTOR>
  static T plicInverse(T d, const Vector<T,DESCRIPTOR::d>& normal) any_platform;

  template <typename CELL, typename V=typename CELL::value_t>
  static V calculateSurfaceTensionCurvature(CELL& cell) any_platform;

  template <typename CELL, typename V=typename CELL::value_t>
  static bool isHealthyInterface(CELL& cell) any_platform;

  template <typename CELL, typename V=typename CELL::value_t>
  static void setCellType(CELL& cell, const FreeSurface::Type& type) any_platform;

  template <typename CELL, typename V=typename CELL::value_t>
  static void setCellFlags(CELL& cell, const FreeSurface::Flags& flags) any_platform;

  struct NeighbourInfo {
    bool has_fluid_neighbours = false;
    bool has_gas_neighbours = false;
    size_t interface_neighbours = 0;
  };

  template <typename CELL>
  static NeighbourInfo getNeighbourInfo(CELL& cell) any_platform;

  /**
   * Helper class meant to make the necessary variables more accesible and interchangeable
   */
  template <typename T>
  struct Variables {
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
class FreeSurfaceMassFlowPostProcessor2D {
public:
  using ParametersD = FreeSurface2D::Variables<T>;

  static constexpr OperatorScope scope = OperatorScope::PerCellWithParameters;

  int getPriority() const {
    return 1;
  }

  template <typename CELL, typename PARAMETERS>
  void apply(CELL& cell, PARAMETERS& parameters) any_platform;

};

/*
 * Free Surface Processor 4
 * ToFluid
 * Converts cells to interface from gas if a neighbouring cell was converted to a fluid cell
 */
template<typename T, typename DESCRIPTOR>
class FreeSurfaceToFluidCellConversionPostProcessor2D  {
public:
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  int getPriority() const {
    return 4;
  }

  template <typename CELL>
  void apply(CELL& cell) any_platform;

};

/**
 * Free Surface Processor 5
 * ToGas
 * Converts cells to interface from fluid if a neighbouring cell was converted to a gas cell
 */
template<typename T, typename DESCRIPTOR>
class FreeSurfaceToGasCellConversionPostProcessor2D  {
public:
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  int getPriority() const {
    return 5;
  }

  template <typename CELL>
  void apply(CELL& cell) any_platform;

};

/**
 * Free Surface Processor 6
 * Calculates mass excess from the cell type conversions and distributes them to neighbouring interface cells
 * Keeps mass local if no neighbour exists until an interface reappears at this position
 */
template<typename T, typename DESCRIPTOR>
class FreeSurfaceMassExcessPostProcessor2D  {
public:
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  int getPriority() const {
    return 6;
  }

  template <typename CELL>
  void apply(CELL& cell) any_platform;

};

/**
 * Free Surface Processor 7
 * Finishes up left over cell conversions and prepares the state for the next simulation step
 */
template<typename T, typename DESCRIPTOR>
class FreeSurfaceFinalizeConversionPostProcessor2D {
public:
  static constexpr OperatorScope scope = OperatorScope::PerCell;

  int getPriority() const {
    return 7;
  }

  template <typename CELL>
  void apply(CELL& cell) any_platform;

};


/// Generator class for the PostProcessors tracking the interface.

/*
* Setup helper
*/

template<typename T, typename DESCRIPTOR>
class FreeSurface2DSetup {
public:
private:
  const FreeSurface2D::Variables<T> vars;
  SuperLattice<T, DESCRIPTOR>& sLattice;

  // SuperPostProcessors
  // Corresponding to the local block processors
public:
  FreeSurface2DSetup(
    SuperLattice<T, DESCRIPTOR>& sLattice,
    const FreeSurface2D::Variables<T>& vars);

  void addPostProcessor();
};
}

#endif
