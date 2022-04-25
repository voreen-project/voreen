/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2020,2021 Claudius Holeksa, Robin Trunk
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

#ifndef FREE_SURFACE_POST_PROCESSOR_2D_HH
#define FREE_SURFACE_POST_PROCESSOR_2D_HH

#include "freeSurfacePostProcessor2D.h"
#include "core/blockLattice.h"

#include <cfenv>

namespace olb {

template <typename CELL>
bool FreeSurface2D::isCellType(CELL& cell, const FreeSurface::Type& type) {
  return cell.template getField<FreeSurface::CELL_TYPE>() == type;
}

template <typename CELL>
bool FreeSurface2D::hasCellFlags(CELL& cell, const FreeSurface::Flags& flags) {
  return static_cast<bool>(cell.template getField<FreeSurface::CELL_FLAGS>() & flags);
}

template <typename CELL>
bool FreeSurface2D::hasNeighbour(CELL& cell, const FreeSurface::Type& type) {
  for(int i = -1; i <= 1; ++i ){
    for(int j = -1; j <= 1; ++j){
      if( i == 0 && j == 0){
        continue;
      }
      auto cellC = cell.neighbor({i,j});
      if(isCellType(cellC, type)){
        return true;
      }
    }
  }
  return false;
}

template <typename CELL>
bool FreeSurface2D::hasNeighbourFlags(CELL& cell, const FreeSurface::Flags& flags) {
  for(int i = -1; i <= 1; ++i ){
    for(int j = -1; j <= 1; ++j){
      if( i == 0 && j == 0){
        continue;
      }
      auto cellC = cell.neighbor({i,j});
      if (hasCellFlags(cellC, flags)) {
        return true;
      }
    }
  }
  return false;
}

template <typename CELL, typename V>
V FreeSurface2D::getClampedEpsilon(CELL& cell) {
  V epsilon = cell.template getField<FreeSurface::EPSILON>();
  return util::max(0., util::min(1., epsilon));
}

template <typename CELL, typename V>
Vector<V,CELL::descriptor_t::d> FreeSurface2D::computeInterfaceNormal(CELL& cell) {
  Vector<V,CELL::descriptor_t::d> normal{};
  normal[0] = 0.5 * (  getClampedEpsilon(cell.neighbor({-1, 0}))
                     - getClampedEpsilon(cell.neighbor({ 1, 0})));
  normal[1] = 0.5 * (  getClampedEpsilon(cell.neighbor({ 0,-1}))
                     - getClampedEpsilon(cell.neighbor({ 0, 1})));
  return normal;
}

template <typename CELL, typename V>
Vector<V,CELL::descriptor_t::d> FreeSurface2D::computeParkerYoungInterfaceNormal(CELL& cell) {
  Vector<V,CELL::descriptor_t::d> normal{};
  for(int i = -1; i <= 1; ++i){
    for(int j = -1; j <= 1; ++j){
      if(i == 0 && j == 0){
        continue;
      }

      int omega_weight = 1;
      if(i != 0){
        omega_weight *= 2;
      }
      if(j != 0){
        omega_weight *= 2;
      }
      omega_weight /= 2;

      auto cellC = cell.neighbor({i,j});
      V epsilon = getClampedEpsilon(cellC);

      normal[0] -= omega_weight * (i * epsilon);
      normal[1] -= omega_weight * (j * epsilon);
    }
  }
  return normal;
}

template<typename T, typename DESCRIPTOR>
T FreeSurface2D::plicInverse(T d_o, const Vector<T,DESCRIPTOR::d>& normal) {
  const T n1 = util::min(util::abs(normal[0]), util::abs(normal[1]));
  const T n2 = util::max(util::abs(normal[0]), util::abs(normal[1]));

  const T d = 0.5 * (n1+n2) - util::abs(d_o);

  T vol;
  if(d < n1){
    vol = d * d / (2. * n1 * n2);
  }
  if(d >= n1){
    vol = d / n2 - n1 / (2. * n2);
  }

  return std::copysign(0.5 - vol, d_o) + 0.5;
}

namespace {

template<typename T, typename DESCRIPTOR>
std::enable_if_t<DESCRIPTOR::d == 2,T>
offsetHelper(T volume, const Vector<T,DESCRIPTOR::d>& sorted_normal) any_platform {
  T d2 = volume * sorted_normal[1] + 0.5 * sorted_normal[0];
  if(d2 >= sorted_normal[0]){
    return d2;
  }

  T d1 = util::sqrt(2. * sorted_normal[0] * sorted_normal[1] * volume);
  // if ( d1 < sorted_normal[0] )
  return d1;
}

}

template<typename T, typename DESCRIPTOR>
T FreeSurface2D::calculateCubeOffset(T volume, const Vector<T,DESCRIPTOR::d>& normal) {
  Vector<T,DESCRIPTOR::d> abs_normal{
    util::abs(normal[0]),
    util::abs(normal[1])
  };

  T volume_symmetry = 0.5 - std::abs(volume - 0.5);

  Vector<T,DESCRIPTOR::d> sorted_normal{
    util::max(util::min(abs_normal[0], abs_normal[1]), 1e-5),
    util::max(abs_normal[0], abs_normal[1])
  };

  T d = offsetHelper<T,DESCRIPTOR>(volume_symmetry, sorted_normal);

  return std::copysign(d - 0.5 * (sorted_normal[0] + sorted_normal[1]), volume - 0.5);
}

template <typename CELL>
FreeSurface2D::NeighbourInfo FreeSurface2D::getNeighbourInfo(CELL& cell) {
  NeighbourInfo info{};
  using DESCRIPTOR = typename CELL::descriptor_t;

  for(int iPop = 1; iPop < DESCRIPTOR::q; ++iPop){
    auto cellC = cell.neighbor({descriptors::c<DESCRIPTOR>(iPop, 0),
                                descriptors::c<DESCRIPTOR>(iPop, 1)});

    if(isCellType(cellC, FreeSurface::Type::Gas)){
      info.has_gas_neighbours = true;
    }
    else if(isCellType(cellC, FreeSurface::Type::Fluid)){
      info.has_fluid_neighbours = true;
    }
    else if(isCellType(cellC, FreeSurface::Type::Interface)){
      ++info.interface_neighbours;
    }
  }
  return info;
}

template <typename CELL, typename V>
V FreeSurface2D::calculateSurfaceTensionCurvature(CELL& cell) {
  auto normal = FreeSurface2D::computeParkerYoungInterfaceNormal(cell);

  using DESCRIPTOR = typename CELL::descriptor_t;
  {
    V norm = 0.;
    for(size_t i = 0; i < DESCRIPTOR::d; ++i){
      norm += normal[i] * normal[i];
    }

    norm = util::sqrt(norm);

    if(norm < 1e-6){
      return 0.;
    }

    for(size_t i = 0; i <DESCRIPTOR::d; ++i){
      normal[i] /= norm;
    }
  }

  // Rotation matrix is
  // ( n1 | -n0 )
  // ( n0 |  n1 )

  // It is 2 because of the amount of fitting parameters. Not because of the dimension
  constexpr size_t S = 2;
  std::array<std::array<V,S>, S> lq_matrix;
  std::array<V,S> b_rhs;
  for(size_t i = 0; i < S; ++i){
    for(size_t j = 0; j < S; ++j){
      lq_matrix[i][j] = 0.;
    }
    b_rhs[i] = 0.;
  }

  // Offset for the plic correction
  V origin_offset = 0.;
  {
    V fill_level = getClampedEpsilon(cell);
    origin_offset = FreeSurface2D::calculateCubeOffset<V,DESCRIPTOR>(fill_level, normal);
  }

  // The amount of neighbouring interfaces. if less are available we will solve a reduced curve by setting the less important parameters to zero
  std::size_t healthy_interfaces = 0;

  for(int i = -1; i <= 1; ++i){
    for(int j = -1 ; j <= 1; ++j){
      if(i == 0 && j == 0){
        continue;
      }

      auto cellC = cell.neighbor({i,j});

      if(   !FreeSurface2D::isCellType(cellC, FreeSurface::Type::Interface)
         || !FreeSurface2D::hasNeighbour(cellC, FreeSurface::Type::Gas)) {
        continue;
      }

      ++healthy_interfaces;

      V fill_level = getClampedEpsilon(cellC);

      V cube_offset = FreeSurface2D::calculateCubeOffset<V,DESCRIPTOR>(fill_level, normal);

      V x_pos = i;
      V y_pos = j;

      // Rotation
      V rot_x_pos = x_pos * normal[1] - y_pos * normal[0];
      V rot_y_pos = x_pos * normal[0] + y_pos * normal[1] + (cube_offset - origin_offset);

      V rot_x_pos_2 = rot_x_pos * rot_x_pos;
      V rot_x_pos_3 = rot_x_pos_2 * rot_x_pos;
      V rot_x_pos_4 = rot_x_pos_3 * rot_x_pos;

      lq_matrix[1][1] += rot_x_pos_2;
      lq_matrix[1][0] += rot_x_pos_3;
      lq_matrix[0][0] += rot_x_pos_4;

      b_rhs[0] += rot_x_pos_2*(rot_y_pos);
      b_rhs[1] += rot_x_pos*(rot_y_pos);
    }
  }

  lq_matrix[0][1] = lq_matrix[1][0];

  // Thikonov regularization parameter
  V alpha = 0.0;
  for(size_t i = 0; i < DESCRIPTOR::d; ++i){
    lq_matrix[i][i] += alpha;
  }

  // It is 2 because of the fitting parameters. Not dependent on the dimension
  std::array<V,S> solved_fit = FreeSurface::solvePivotedLU<V,S>(lq_matrix, b_rhs, healthy_interfaces);

  // signed curvature -> kappa = y'' / ( (1 + y'²)^(3/2) )
  V denom = std::sqrt(1. + solved_fit[1]*solved_fit[1]);
  denom = denom * denom * denom;
  V curvature = 2.*solved_fit[0] / denom;
  return util::max(-1., util::min(1., curvature));
}

template <typename CELL, typename V>
void FreeSurface2D::setCellType(CELL& cell, const FreeSurface::Type& type) {
  cell.template setField<FreeSurface::CELL_TYPE>(type);
}

template <typename CELL, typename V>
void FreeSurface2D::setCellFlags(CELL& cell, const FreeSurface::Flags& flags){
  cell.template setField<FreeSurface::CELL_FLAGS>(flags);
}

template <typename CELL, typename V>
bool FreeSurface2D::isHealthyInterface(CELL& cell) {
  bool has_fluid_neighbours = false;
  bool has_gas_neighbours = false;

  if(!isCellType(cell, FreeSurface::Type::Interface)){
    return false;
  }

  using DESCRIPTOR = typename CELL::descriptor_t;
  for(int iPop = 1; iPop < DESCRIPTOR::q; ++iPop){
    auto cellC = cell.neighbor({descriptors::c<DESCRIPTOR>(iPop, 0),
                                descriptors::c<DESCRIPTOR>(iPop, 1)});
    if(isCellType(cellC, FreeSurface::Type::Gas)){
      has_gas_neighbours = true;
      if(has_fluid_neighbours){
        return true;
      }
    }
    else if(isCellType(cellC, FreeSurface::Type::Fluid)){
      has_fluid_neighbours = true;
      if(has_gas_neighbours){
        return true;
      }
    }
  }
  return false;
}

// LocalProcessor 1

// Read

// range 0
// Mass
// CellType

// range 1
// DFs
// Epsilon
// CellType

// Write (always range 0)
// Mass
// CellFlags
// DFs (only replacing from incoming gas stream)

template <typename T, typename DESCRIPTOR>
template <typename CELL, typename PARAMETERS>
void FreeSurfaceMassFlowPostProcessor2D<T, DESCRIPTOR>::apply(CELL& cell, PARAMETERS& vars) {
  /*
  * Minor "hack". Remove all cell flags here, because it is needed in the last processor due to pulling steps in processor 6 and 7
  */
  FreeSurface2D::setCellFlags(cell, static_cast<FreeSurface::Flags>(0));

  /*
  * This processor only works on interface types
  */
  /*if(FreeSurface2D::isCellType(cell, FreeSurface::Type::Fluid )){

    T mass_tmp = blockLattice.get(iX, iY).template getField<FreeSurface::MASS>();
    for(int iPop = 1; iPop < DESCRIPTOR::q; ++iPop){
      int iXc = iX + descriptors::c<DESCRIPTOR>(iPop, 0);
      int iYc = iY + descriptors::c<DESCRIPTOR>(iPop, 1);
      int iPop_op = descriptors::opposite<DESCRIPTOR>(iPop);
      mass_tmp += blockLattice.get(iX, iY)[iPop_op] - blockLattice.get(iXc, iYc)[iPop];
    }
    blockLattice.get(iX, iY).template setField<FreeSurface::MASS>(mass_tmp);
  }
  else */
  if (FreeSurface2D::isCellType(cell, FreeSurface::Type::Interface )) {
    T mass_tmp = cell.template getField<FreeSurface::MASS>();

    FreeSurface2D::NeighbourInfo neighbour_info = FreeSurface2D::getNeighbourInfo(cell);

    for (int iPop = 1; iPop < DESCRIPTOR::q; ++iPop){
      auto cellC = cell.neighbor({descriptors::c<DESCRIPTOR>(iPop, 0),
                                  descriptors::c<DESCRIPTOR>(iPop, 1)});
      int iPop_op = descriptors::opposite<DESCRIPTOR>(iPop);

      /*
      * Iterate over neighbours and perform a mass exchange. Interface to fluid are simple cases.
      * Interface to interface has to be symmetric and multiple cases are for artifact reduction
      * from Thuerey
      * Added a distinction between the amount of interface nodes. Weight consideration seems to cause artifacts.
      */
      if( FreeSurface2D::isCellType(cellC, FreeSurface::Type::Fluid)) {
        mass_tmp += cell[iPop_op] - cellC[iPop];
      } else if ( FreeSurface2D::isCellType(cellC, FreeSurface::Type::Interface)) {
        FreeSurface2D::NeighbourInfo neighbour_neighbour_info = FreeSurface2D::getNeighbourInfo(cellC);

        T mass_flow = 0.;

        if( !neighbour_info.has_fluid_neighbours){
          if(!neighbour_neighbour_info.has_fluid_neighbours){
            if(neighbour_info.interface_neighbours < neighbour_neighbour_info.interface_neighbours){
              mass_flow = -cellC[iPop];// - descriptors::t<T,DESCRIPTOR>(iPop);
            }else if(neighbour_info.interface_neighbours > neighbour_neighbour_info.interface_neighbours){
              mass_flow = cell[iPop_op];// + descriptors::t<T,DESCRIPTOR>(iPop_op);
            }else{
              mass_flow = cell[iPop_op] - cellC[iPop];
            }
          }else {
            mass_flow = -cellC[iPop];// - descriptors::t<T,DESCRIPTOR>(iPop);
          }
        }else if(!neighbour_info.has_gas_neighbours){
          if(!neighbour_neighbour_info.has_gas_neighbours){
            if(neighbour_info.interface_neighbours < neighbour_neighbour_info.interface_neighbours){
              mass_flow = cell[iPop_op];// + descriptors::t<T,DESCRIPTOR>(iPop_op);
            }else if(neighbour_info.interface_neighbours > neighbour_neighbour_info.interface_neighbours){
              mass_flow = -cellC[iPop];// - descriptors::t<T,DESCRIPTOR>(iPop);
            }else{
              mass_flow = cell[iPop_op] - cellC[iPop];
            }
          }else {
            mass_flow = cell[iPop_op];// + descriptors::t<T,DESCRIPTOR>(iPop_op);
          }
        }else {
          if(!neighbour_neighbour_info.has_fluid_neighbours){
            mass_flow = cell[iPop_op];// + descriptors::t<T,DESCRIPTOR>(iPop_op);
          }else if(!neighbour_neighbour_info.has_gas_neighbours){
            mass_flow = -cellC[iPop];// - descriptors::t<T,DESCRIPTOR>(iPop);
          }else {
            mass_flow = cell[iPop_op] - cellC[iPop];
          }
        }

        /*
        * Exchange depends on how filled the interfaces are
        */
        mass_tmp += mass_flow * 0.5 * (FreeSurface2D::getClampedEpsilon(cell) + FreeSurface2D::getClampedEpsilon(cellC));
      }
    }

    cell.template setField<FreeSurface::MASS>(mass_tmp);

    // Former 2 Step

    // Because I need the distribution functions of the last step I will write results in a temporary
    // array, before copying it back into the DFs

    Vector<T,DESCRIPTOR::q> dfs;

    T curvature = 0.;

    if(vars.has_surface_tension){
      FreeSurface2D::NeighbourInfo info = FreeSurface2D::getNeighbourInfo(cell);
      if(info.has_gas_neighbours){
        curvature = FreeSurface2D::calculateSurfaceTensionCurvature(cell);
      }

    }

    // Gas pressure adjusting
    T gas_pressure = 1. - 6. * vars.surface_tension_parameter * curvature;

    for(int iPop=1; iPop < DESCRIPTOR::q; iPop++) {
      auto cellC = cell.neighbor({descriptors::c<DESCRIPTOR>(iPop, 0),
                                  descriptors::c<DESCRIPTOR>(iPop, 1)});
      int iPop_op = descriptors::opposite<DESCRIPTOR>(iPop);

      /*
      * Gas replacement
      */
      if ( FreeSurface2D::isCellType(cellC, FreeSurface::Type::Gas )) {
        Vector<T, DESCRIPTOR::d> u_vel = cell.template getField<FreeSurface::PREVIOUS_VELOCITY>();
        T u[DESCRIPTOR::d];
        for(size_t u_i = 0; u_i < DESCRIPTOR::d; ++u_i){
          u[u_i] = u_vel[u_i];
        }
        T uSqr = util::normSqr<T,DESCRIPTOR::d>(u);
        dfs[iPop_op] = equilibrium<DESCRIPTOR>::secondOrder(iPop, gas_pressure, u)
          + equilibrium<DESCRIPTOR>::secondOrder(iPop_op, gas_pressure, u)
          - cellC[iPop];
      }else {
        dfs[iPop_op] = cell[iPop_op];
      }
    }

    for(int iPop=1; iPop<DESCRIPTOR::q; iPop++) {
      cell[iPop] = dfs[iPop];
    }

    // Former 3 Step
    /*
    * Based on the mass calculation, flag this interface cell as toFluid or toGas if set boundaries are met
    */
    T rho = cell.computeRho();

    // Check if transition needs to happen.
    if ( mass_tmp < -vars.transition * rho || (mass_tmp < vars.lonely_threshold * rho && !neighbour_info.has_fluid_neighbours) ){
      FreeSurface2D::setCellFlags(cell, FreeSurface::Flags::ToGas);
    }
    else if ( mass_tmp > (1. + vars.transition)*rho  || ( mass_tmp > (1-vars.lonely_threshold) * rho && !neighbour_info.has_gas_neighbours) ){
      FreeSurface2D::setCellFlags(cell, FreeSurface::Flags::ToFluid);
    }else if(vars.drop_isolated_cells && (neighbour_info.interface_neighbours == 0)){
      if(!neighbour_info.has_gas_neighbours){
        FreeSurface2D::setCellFlags(cell, FreeSurface::Flags::ToFluid);
      }else if(!neighbour_info.has_fluid_neighbours){
        //FreeSurface2D::setCellFlags(cell, FreeSurface::Flags::ToGas);
      }
    }
  }
}


// Processor 4

// Read
// range 0
// CellType
// CellFlags

// range 1
// CellType
// CellFlags
// DFs

// Write (always range 0)
// CellFlags
// DFs
template <typename T, typename DESCRIPTOR>
template <typename CELL>
void FreeSurfaceToFluidCellConversionPostProcessor2D<T, DESCRIPTOR>::apply(CELL& cell) {
  /*
  * Initializing new interface cells with DFs from surrounding fluid and interface cells
  * Just takes the arithmetic average.
  */
  if(FreeSurface2D::isCellType(cell, FreeSurface::Type::Gas)){
    if(FreeSurface2D::hasNeighbourFlags(cell, FreeSurface::Flags::ToFluid)){
      FreeSurface2D::setCellFlags(cell, FreeSurface::Flags::NewInterface);
      T rho_avg = 0.;
      T u_avg[DESCRIPTOR::d] = {0., 0.};

      std::size_t ctr = 0;

      for(int iPop=1; iPop<DESCRIPTOR::q; iPop++) {
        auto cellC = cell.neighbor({descriptors::c<DESCRIPTOR>(iPop,0),
                                    descriptors::c<DESCRIPTOR>(iPop,1)});

        if (FreeSurface2D::isCellType(cellC, FreeSurface::Type::Fluid) || FreeSurface2D::isCellType(cellC, FreeSurface::Type::Interface)){
          T rho_tmp = 0.;
          T u_tmp[DESCRIPTOR::d] = {0., 0.};
          ++ctr;
          cellC.computeRhoU(rho_tmp, u_tmp);
          rho_avg += rho_tmp;
          for(size_t i = 0; i < DESCRIPTOR::d; ++i){
            u_avg[i] += u_tmp[i];
          }
        }
      }

      if(ctr > 0){
        rho_avg /= static_cast<T>(ctr);
        for(size_t i = 0; i < DESCRIPTOR::d; ++i){
          u_avg[i] /= static_cast<T>(ctr);
        }
      }

      cell.iniEquilibrium(rho_avg, u_avg);
    }
  }else if(FreeSurface2D::hasCellFlags(cell, FreeSurface::Flags::ToGas)){
    /*
    * If a toGas cell has a neighbouring toFluid cell, unset the toGas flag
    */
    if(FreeSurface2D::hasNeighbourFlags(cell, FreeSurface::Flags::ToFluid)){
      FreeSurface2D::setCellFlags(cell, static_cast<FreeSurface::Flags>(0));
    }
  }
}

// LocalProcessor 5

// Read
// range 0
// CellType
// DFs

// range 1
// CellFlags

// Write (always range 0)
// CellFlags

template <typename T, typename DESCRIPTOR>
template <typename CELL>
void FreeSurfaceToGasCellConversionPostProcessor2D<T, DESCRIPTOR>::apply(CELL& cell) {
  /*
  * For the to be converted toGas cells, set the neighbours to interface cells
  */
  if (   FreeSurface2D::isCellType(cell, FreeSurface::Type::Fluid)
      && FreeSurface2D::hasNeighbourFlags(cell, FreeSurface::Flags::ToGas)){
    FreeSurface2D::setCellFlags(cell, FreeSurface::Flags::NewInterface);
    T rho = cell.computeRho();
    cell.template setField<FreeSurface::MASS>(rho);
  }
}

// LocalProcessor 6

// Read
// range 0
// CellType
// CellFlags
// DFs
// Mass

// range 1
// Epsilon
// CellType
// CellFlags

// Write (always range 0)
// Mass
// TempMassExchange
template <typename T, typename DESCRIPTOR>
template <typename CELL>
void FreeSurfaceMassExcessPostProcessor2D<T, DESCRIPTOR>::apply(CELL& cell) {
  if( !FreeSurface2D::isCellType(cell, FreeSurface::Type::Interface) ){
    return;
  }

  T rho = cell.computeRho();
  T mass = cell.template getField<FreeSurface::MASS>( );
  T mass_excess = 0.;

  auto normal = FreeSurface2D::computeParkerYoungInterfaceNormal(cell);
  // redistribute excess mass

  /// @hint EPSILON of neighbours used here
  /// @hint Mass can be set in this processor, but not epsilon since it is needed for the normal computation. epsilon is set in the next processor
  /// Became untrue due to code section removal, but epsilon is still set in the next part because of pushed mass excess
  if(FreeSurface2D::hasCellFlags(cell, FreeSurface::Flags::ToGas)){
    mass_excess = mass;
    cell.template setField<FreeSurface::MASS>( 0. );
    normal = {-normal[0], -normal[1]};
  } else if (FreeSurface2D::hasCellFlags(cell, FreeSurface::Flags::ToFluid)) {
    mass_excess = mass - rho;
    cell.template setField<FreeSurface::MASS>( rho );
  } else {
    return;
  }

  std::array<T,DESCRIPTOR::q> products;
  products[0] = 0.;
  T product_sum = 0.;
  std::size_t product_total = 0;

  for(int iPop=1; iPop<DESCRIPTOR::q; iPop++) {
    auto cellC = cell.neighbor({descriptors::c<DESCRIPTOR>(iPop,0),
                                descriptors::c<DESCRIPTOR>(iPop,1)});
    products[iPop] = 0.;

    // Thuerey Paper says we can't use new interface cells
    // or flagged cells
    // But surface tension showed us that it has anisotropic effects
    if( (FreeSurface2D::isCellType(cellC, FreeSurface::Type::Interface) && (!FreeSurface2D::hasCellFlags(cellC,
          static_cast<FreeSurface::Flags>(255)) /*|| FreeSurface2D::hasCellFlags(cellC, FreeSurface::Flags::ToFluid)*/ ))
        /*|| FreeSurface2D::isCellType(cellC, FreeSurface::Type::Fluid)*/
        ){
      products[iPop] = (normal[0] * descriptors::c<DESCRIPTOR>(iPop, 0) + normal[1] * descriptors::c<DESCRIPTOR>(iPop,1));
      if(products[iPop] <= 0){
        products[iPop] = 0.;
      }
      ++product_total;
      product_sum += products[iPop];
    }
  }

  /* Prepare Mass excess push */
  Vector<T,DESCRIPTOR::q> mass_excess_vector{};
  mass_excess_vector[0] = 0.;
  /*
  if(product_sum > 0){
    T fraction = 1./ product_sum;

    for(int iPop = 1; iPop < DESCRIPTOR::q; ++iPop){
      mass_excess_vector[iPop] = mass_excess * products[iPop] * fraction;
    }
    cell.template setField<FreeSurface::TEMP_MASS_EXCHANGE>( mass_excess_vector );
  }
  else*/
  if (product_total > 0) {
    T product_fraction = 1. / product_total;
    for(int iPop=1; iPop < DESCRIPTOR::q; iPop++) {
      auto cellC = cell.neighbor({descriptors::c<DESCRIPTOR>(iPop,0),
                                  descriptors::c<DESCRIPTOR>(iPop,1)});
      if( (FreeSurface2D::isCellType(cellC, FreeSurface::Type::Interface) && (!FreeSurface2D::hasCellFlags(cellC,
          static_cast<FreeSurface::Flags>(255)) /*|| FreeSurface2D::hasCellFlags(cellC, FreeSurface::Flags::ToFluid)*/ ))
          /*|| FreeSurface2D::isCellType(cellC, FreeSurface::Type::Fluid)*/
          ){
        mass_excess_vector[iPop] = mass_excess * product_fraction;
      } else {
        mass_excess_vector[iPop] = 0.;
      }
    }
    cell.template setField<FreeSurface::TEMP_MASS_EXCHANGE>( mass_excess_vector );
  } else {
    mass_excess_vector[0] = mass_excess;
    for(int iPop=1; iPop < DESCRIPTOR::q; iPop++) {
      mass_excess_vector[iPop] = 0.;
    }
    cell.template setField<FreeSurface::TEMP_MASS_EXCHANGE>( mass_excess_vector );
  }
}


// LocalProcessor 7

// Read
// range 0
// CellFlags
// CellType
// DFs
// Mass
// Epsilon

// range 1
// CellFlags
// TempMassExchange

// Write (always range 0)
// Epsilon
// CellType
// Mass

template <typename T, typename DESCRIPTOR>
template <typename CELL>
void FreeSurfaceFinalizeConversionPostProcessor2D<T, DESCRIPTOR>::apply(CELL& cell) {
  /* Convert flagged cells to appropriate cell types */
  FreeSurface::Flags flags = static_cast<FreeSurface::Flags>(cell.template getField<FreeSurface::CELL_FLAGS>());

  switch(flags){
    case FreeSurface::Flags::ToFluid:
    {
      /// @hint moved flag removal to processor 1 without any negative effects
      FreeSurface2D::setCellType(cell, FreeSurface::Type::Fluid);
      cell.template setField<FreeSurface::EPSILON>( 1. );
      T mass_tmp = cell.template getField<FreeSurface::MASS>();
      mass_tmp += cell.template getField<FreeSurface::TEMP_MASS_EXCHANGE>()[0];
      cell.template setField<FreeSurface::MASS>(mass_tmp);

    }
    break;
    case FreeSurface::Flags::ToGas:
    {
      /// @hint moved flag removal to processor 1 without any negative effects
      FreeSurface2D::setCellType(cell, FreeSurface::Type::Gas);
      cell.template setField<FreeSurface::EPSILON>( 0. );
      T mass_tmp = cell.template getField<FreeSurface::MASS>();
      mass_tmp += cell.template getField<FreeSurface::TEMP_MASS_EXCHANGE>()[0];
      cell.template setField<FreeSurface::MASS>(mass_tmp);

    }
    break;
    case FreeSurface::Flags::NewInterface:
    {
      FreeSurface2D::setCellType(cell, FreeSurface::Type::Interface);
    }
    break;
    default:
    break;
  }

  FreeSurface::Type type = static_cast<FreeSurface::Type>(cell.template getField<FreeSurface::CELL_TYPE>());

  /* Collection of mass excess in a pulling step */
  switch(type){
    case FreeSurface::Type::Interface:
    {
      T collected_excess = 0.;
      for(int iPop = 1; iPop < DESCRIPTOR::q; ++iPop){
        auto cellC = cell.neighbor({descriptors::c<DESCRIPTOR>(iPop,0),
                                    descriptors::c<DESCRIPTOR>(iPop,1)});
        auto tempMassExchange = cellC.template getFieldPointer<FreeSurface::TEMP_MASS_EXCHANGE>();

        if (FreeSurface2D::hasCellFlags(cellC,
                FreeSurface::Flags::ToFluid
              | FreeSurface::Flags::ToGas)){
            int iPop_op = descriptors::opposite<DESCRIPTOR>(iPop);
            collected_excess += tempMassExchange[iPop_op];
          }
      }

      T mass_tmp = cell.template getField<FreeSurface::MASS>();

      mass_tmp += collected_excess;
      T rho;
      T u_tmp[DESCRIPTOR::d];
      cell.computeRhoU(rho, u_tmp);

      Vector<T,DESCRIPTOR::d> u_vel{u_tmp[0], u_tmp[1]};

      cell.template setField<FreeSurface::EPSILON>( mass_tmp / rho );
      cell.template setField<FreeSurface::MASS>(mass_tmp);
      cell.template setField<FreeSurface::PREVIOUS_VELOCITY>(u_vel);
    }
    break;
    case FreeSurface::Type::Fluid:
    {
      T collected_excess = 0.;

      for(int iPop = 1; iPop < DESCRIPTOR::q; ++iPop){
        auto cellC = cell.neighbor({descriptors::c<DESCRIPTOR>(iPop,0),
                                    descriptors::c<DESCRIPTOR>(iPop,1)});
        auto tempMassExchange = cellC.template getFieldPointer<FreeSurface::TEMP_MASS_EXCHANGE>();
        if (FreeSurface2D::hasCellFlags(cellC,
              FreeSurface::Flags::ToFluid
            | FreeSurface::Flags::ToGas)) {
          int iPop_op = descriptors::opposite<DESCRIPTOR>(iPop);
          collected_excess += tempMassExchange[iPop_op];
        }
      }

      T mass_tmp = cell.template getField<FreeSurface::MASS>();
      mass_tmp += collected_excess;
      cell.template setField<FreeSurface::MASS>(mass_tmp);
    }
    break;
    default:
    break;
  }
}

// Setup

template<typename T, typename DESCRIPTOR>
FreeSurface2DSetup<T,DESCRIPTOR>::FreeSurface2DSetup(SuperLattice<T, DESCRIPTOR>& sLattice,
    const FreeSurface2D::Variables<T>& vars_)
:
  vars{vars_},
  sLattice{sLattice}
{}

template<typename T, typename DESCRIPTOR>
void FreeSurface2DSetup<T,DESCRIPTOR>::addPostProcessor(){
  sLattice.template addPostProcessor<FreeSurface::Stage0>(
    meta::id<FreeSurfaceMassFlowPostProcessor2D<T,DESCRIPTOR>>{});
  sLattice.template addPostProcessor<FreeSurface::Stage1>(
    meta::id<FreeSurfaceToFluidCellConversionPostProcessor2D<T,DESCRIPTOR>>{});
  sLattice.template addPostProcessor<FreeSurface::Stage2>(
    meta::id<FreeSurfaceToGasCellConversionPostProcessor2D<T,DESCRIPTOR>>{});
  sLattice.template addPostProcessor<FreeSurface::Stage3>(
    meta::id<FreeSurfaceMassExcessPostProcessor2D<T,DESCRIPTOR>>{});
  sLattice.template addPostProcessor<FreeSurface::Stage4>(
    meta::id<FreeSurfaceFinalizeConversionPostProcessor2D<T,DESCRIPTOR>>{});

  for (int iC=0; iC < sLattice.getLoadBalancer().size(); ++iC) {
    sLattice.getBlock(iC)
            .template getData<TrivialParameters<FreeSurfaceMassFlowPostProcessor2D<T,DESCRIPTOR>>>() = vars;
  }

  {
    // Communicate DFs, Epsilon and Cell Types
    auto& communicator = sLattice.getCommunicator(FreeSurface::Stage0());
    communicator.requestOverlap(2);
    communicator.template requestField<FreeSurface::EPSILON>();
    communicator.template requestField<FreeSurface::CELL_TYPE>();
    communicator.template requestField<descriptors::POPULATION>();
    communicator.exchangeRequests();
  }

  {
    // Communicate DFs, Cell Flags
    auto& communicator = sLattice.getCommunicator(FreeSurface::Stage1());
    communicator.requestOverlap(2);
    communicator.template requestField<FreeSurface::CELL_FLAGS>();
    communicator.template requestField<descriptors::POPULATION>();
    communicator.exchangeRequests();
  }

  {
    // Communicate Cell Flags
    auto& communicator = sLattice.getCommunicator(FreeSurface::Stage2());
    communicator.requestOverlap(2);
    communicator.template requestField<FreeSurface::CELL_FLAGS>();
    communicator.exchangeRequests();
  }

  {
    // Communicate Cell Flags
    auto& communicator = sLattice.getCommunicator(FreeSurface::Stage3());
    communicator.requestOverlap(2);
    communicator.template requestField<FreeSurface::CELL_FLAGS>();
    communicator.exchangeRequests();
  }

  {
    // Communicate TempMassExchange
    auto& communicator = sLattice.getCommunicator(FreeSurface::Stage4());
    communicator.requestOverlap(2);
    communicator.template requestField<FreeSurface::TEMP_MASS_EXCHANGE>();
    communicator.exchangeRequests();
  }

  sLattice.scheduleCustomPostProcessing([](SuperLattice<T,DESCRIPTOR>& lattice) {
    lattice.executePostProcessors(FreeSurface::Stage0());
    lattice.executePostProcessors(FreeSurface::Stage1());
    lattice.executePostProcessors(FreeSurface::Stage2());
    lattice.executePostProcessors(FreeSurface::Stage3());
    lattice.executePostProcessors(FreeSurface::Stage4());
  });
}

}
#endif
