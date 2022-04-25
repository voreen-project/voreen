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

#ifndef FREE_SURFACE_POST_PROCESSOR_3D_HH
#define FREE_SURFACE_POST_PROCESSOR_3D_HH

#include "freeSurfacePostProcessor3D.h"
#include "core/blockLattice.h"

#include <cmath>
#include <algorithm>

namespace olb {
template<typename T, typename DESCRIPTOR>
bool FreeSurface3D::isCellType(BlockLattice<T, DESCRIPTOR>& blockLattice, int X, int Y, int Z, const FreeSurface::Type& type){
  return blockLattice.get(X,Y,Z).template getField<FreeSurface::CELL_TYPE>() == type;
}

template<typename T, typename DESCRIPTOR>
bool FreeSurface3D::hasCellFlags(BlockLattice<T, DESCRIPTOR>& blockLattice, int X, int Y, int Z, const FreeSurface::Flags& flags){
  return static_cast<bool>(blockLattice.get(X,Y,Z).template getField<FreeSurface::CELL_FLAGS>() & flags);
}

template<typename T, typename DESCRIPTOR>
bool FreeSurface3D::hasNeighbour(BlockLattice<T,DESCRIPTOR>& blockLattice, int X, int Y, int Z, const FreeSurface::Type& type){
  for(int i = -1; i <= 1; ++i ){
    for(int j = -1; j <= 1; ++j){
      for(int k = -1; k <= 1; ++k){
        if( i == 0 && j == 0 && k == 0){
          continue;
        }

        if(isCellType<T,DESCRIPTOR>(blockLattice, X+i, Y+j, Z+k, type)){
          return true;
        }
      }
    }
  }
  return false;
}

template<typename T, typename DESCRIPTOR>
T FreeSurface3D::getClampedEpsilon(BlockLattice<T,DESCRIPTOR>& blockLattice, int X, int Y, int Z){
  T epsilon = blockLattice.get( X, Y, Z ).template getField<FreeSurface::EPSILON>();

  return std::max(0., std::min(1., epsilon));
}

template<typename T, typename DESCRIPTOR>
T FreeSurface3D::getClampedEpsilonCorrected(BlockLattice<T,DESCRIPTOR>& blockLattice, int X, int Y, int Z){
  if(isCellType(blockLattice, X, Y, Z, FreeSurface::Type::Interface) && !hasNeighbour(blockLattice, X, Y, Z, FreeSurface::Type::Gas)){
    return 1.0;
  }
  T epsilon = blockLattice.get( X, Y, Z ).template getField<FreeSurface::EPSILON>();

  T clamped = std::max(0., std::min(1., epsilon));

  return clamped;
}

template<typename T, typename DESCRIPTOR>
std::array<T,DESCRIPTOR::d> FreeSurface3D::computeParkerYoungInterfaceNormal(BlockLattice<T,DESCRIPTOR>& blockLattice, int X, int Y, int Z){
  std::array<T,DESCRIPTOR::d> normal;
  for(size_t dim = 0; dim < DESCRIPTOR::d; ++dim){
    normal[dim] = 0;
  }

  for(int i = -1; i <= 1; ++i){
    for(int j = -1; j <= 1; ++j){
      for(int k = -1; k <= 1; ++k){
        if(i == 0 && j == 0 && k == 0){
          continue;
        }
        int iXc = X + i;
        int iYc = Y + j;
        int iZc = Z + k;

        int omega_weight = 1;
        if(i != 0){
          omega_weight *= 2;
        }
        if(j != 0){
          omega_weight *= 2;
        }
        if(k != 0){
          omega_weight *= 2;
        }
        omega_weight /= 2;

        T epsilon = getClampedEpsilon(blockLattice, iXc, iYc, iZc);
        
        normal[0] -= omega_weight * (i * epsilon);
        normal[1] -= omega_weight * (j * epsilon);
        normal[2] -= omega_weight * (k * epsilon);
      }
    }
  }

  return normal;
}


template<typename T, typename DESCRIPTOR>
bool FreeSurface3D::hasNeighbourFlags(BlockLattice<T,DESCRIPTOR>& blockLattice, int X, int Y, int Z, const FreeSurface::Flags& flags){
  for(int i = -1; i <= 1; ++i ){
    for(int j = -1; j <= 1; ++j){
      for(int k = -1; k <= 1; ++k){
        if( i == 0 && j == 0 && k == 0){
          continue;
        }

        if(hasCellFlags<T,DESCRIPTOR>(blockLattice, X+i, Y+j, Z+k, flags)){
          return true;
        }
      }
    }
  }
  return false;
}

template<typename T, typename DESCRIPTOR>
std::array<T,DESCRIPTOR::d> FreeSurface3D::computeInterfaceNormal(BlockLattice<T,DESCRIPTOR>& blockLattice, int X, int Y, int Z){
  std::array<T,DESCRIPTOR::d> normal;

  normal[0] = 0.5 * (blockLattice.get( X-1, Y, Z ).template getField<FreeSurface::EPSILON>()
              - blockLattice.get( X+1, Y, Z ).template getField<FreeSurface::EPSILON>());
  normal[1] = 0.5 * (blockLattice.get( X, Y-1, Z ).template getField<FreeSurface::EPSILON>()
              - blockLattice.get( X, Y+1, Z ).template getField<FreeSurface::EPSILON>());
  normal[2] = 0.5 * (blockLattice.get( X, Y, Z-1 ).template getField<FreeSurface::EPSILON>()
              - blockLattice.get( X, Y, Z+1 ).template getField<FreeSurface::EPSILON>());
  return normal;
}

namespace {

// A lot of magic numbers are happening here. Optimized algorithm taken from Moritz Lehmann
template<typename T, typename DESCRIPTOR>
std::enable_if_t<DESCRIPTOR::d == 3,T>
offsetHelper(T volume, const std::array<T,DESCRIPTOR::d>& sorted_normal){
  T sn0_plus_sn1 = sorted_normal[0] + sorted_normal[1];
  T sn0_times_sn1 = sorted_normal[0] * sorted_normal[1];
  T sn2_volume = sorted_normal[2] * volume;

  T min_sn0_plus_sn1_and_sn2 = std::min(sn0_plus_sn1, sorted_normal[2]);

  T d5 = sn2_volume + 0.5 * sn0_plus_sn1;
  if(d5 > min_sn0_plus_sn1_and_sn2 && d5 <= sorted_normal[2]){
    return d5;
  }

  T d2 = 0.5 * sorted_normal[0] + 0.28867513 * std::sqrt( std::max(0., 24. * sorted_normal[1] * sn2_volume - sorted_normal[0]*sorted_normal[0]) );

  if(d2 > sorted_normal[0] && d2 <= sorted_normal[1]){
    return d2;
  }

  T d1  = std::cbrt(6.0 * sn0_times_sn1 * sn2_volume);
  if(d1 <= sorted_normal[0]){
    return d1;
  }

  T x3 = 81.0  * sn0_times_sn1 * (sn0_plus_sn1 - 2. * sn2_volume);
  T y3 = std::sqrt(std::max(0., 23328. * sn0_times_sn1*sn0_times_sn1*sn0_times_sn1 - x3*x3 ));
  T u3 = std::cbrt(x3*x3 + y3*y3);
  T d3 = sn0_plus_sn1 - (7.5595264 * sn0_times_sn1 + 0.26456684 * u3) * (1./std::sqrt(u3)) * std::sin(0.5235988 - 0.3333334 * std::atan(y3 / x3));
  if(d3 > sorted_normal[1] && d3 <= min_sn0_plus_sn1_and_sn2){
    return d3;
  }

  T t4 = 9. * std::pow(sn0_plus_sn1 + sorted_normal[2], 2) - 18.;
  T x4 = std::max(sn0_times_sn1 * sorted_normal[2] * (324. - 648. * volume), 1.1754944e-38);
  T y4 = std::sqrt(std::max(4. * t4*t4*t4 - x4*x4, 0.));
  T u4 = std::cbrt(x4*x4 + y4*y4);
  T d4  = 0.5 * (sn0_plus_sn1 + sorted_normal[2]) - (0.20998684 * t4 + 0.13228342 * u4) * (1./std::sqrt(u4)) * std::sin(0.5235988- 0.3333334 * std::atan(y4/x4));

  return d4;
}

// A lot of magic numbers are happening here. Optimized algorithm taken from Moritz Lehmann
template<typename T, typename DESCRIPTOR>
T offsetHelperOpt(T vol, const std::array<T,DESCRIPTOR::d>& sn){
  const T sn0_p_sn1 = sn[0] + sn[1];
  const T sn2_t_V = sn[2] * vol;

  if(sn0_p_sn1 <= 2. * sn2_t_V){
    return sn2_t_V + 0.5 * sn0_p_sn1;
  }

  const T sq_sn0 = std::pow(sn[0],2), sn1_6 = 6. * sn[1], v1 = sq_sn0 / sn1_6;

  if(v1 <= sn2_t_V && sn2_t_V < v1 + 0.5 * (sn[1]-sn[0])){
    return 0.5 *(sn[0] + std::sqrt(sq_sn0 + 8.0 * sn[1] * (sn2_t_V - v1)));
  }

  const T v6 = sn[0] * sn1_6 * sn2_t_V;
  if(sn2_t_V < v1){
    return std::cbrt(v6);
  }

  const T v3 = sn[2] < sn0_p_sn1 ? (std::pow(sn[2],2) * (3. * sn0_p_sn1 - sn[2]) + sq_sn0 *(sn[0] - 3.0 * sn[2]) + std::pow(sn[1],2)*(sn[1]-3.0 * sn[2])) / (sn[0] * sn1_6) : 0.5 * sn0_p_sn1;

  const T sq_sn0_sq_sn1 = sq_sn0 + std::pow(sn[1],2), v6_cb_sn0_sn1 = v6 - std::pow(sn[0],3) - std::pow(sn[1],3);

  const bool case34 = sn2_t_V < v3;
  const T a = case34 ? v6_cb_sn0_sn1 : 0.5 * (v6_cb_sn0_sn1 - std::pow(sn[2], 3));
  const T b = case34 ? sq_sn0_sq_sn1 : 0.5 * (sq_sn0_sq_sn1 + std::pow(sn[2], 2));
  const T c = case34 ? sn0_p_sn1 : 0.5;
  const T t = std::sqrt(std::pow(c,2) - b);
  return c - 2.0 * t * std::sin(0.33333334 * std::asin((std::pow(c,3) - 0.5 * a - 1.5 * b * c) / std::pow(t,3)));
}

}

template<typename T, typename DESCRIPTOR>
T FreeSurface3D::calculateCubeOffset(T volume, const std::array<T,DESCRIPTOR::d>& normal){

  std::array<T, DESCRIPTOR::d> abs_normal{{
    std::abs(normal[0]),
    std::abs(normal[1]),
    std::abs(normal[2])
  }};

  T volume_symmetry = 0.5 - std::abs(volume - 0.5);

  std::array<T, DESCRIPTOR::d> sorted_normal{{
    std::max(std::min(std::min(abs_normal[0], abs_normal[1]), abs_normal[2]),1e-12),
    0.,
    std::max(std::max(abs_normal[0], abs_normal[1]), abs_normal[2])
  }};
  sorted_normal[1] = std::max((abs_normal[0] + abs_normal[1] + abs_normal[2]) - sorted_normal[0] - sorted_normal[2], 1e-12);

  T d = offsetHelper<T,DESCRIPTOR>(volume_symmetry, sorted_normal);

  return std::copysign(d-0.5*(sorted_normal[0] + sorted_normal[1] + sorted_normal[2]), volume - 0.5);
}

template<typename T, typename DESCRIPTOR>
T FreeSurface3D::calculateCubeOffsetOpt(T volume, const std::array<T,DESCRIPTOR::d>& normal){

  std::array<T, DESCRIPTOR::d> abs_normal{{
    std::abs(normal[0]),
    std::abs(normal[1]),
    std::abs(normal[2])
  }};
  T a_l1 = abs_normal[0] + abs_normal[1] + abs_normal[2];

  T volume_symmetry = 0.5 - std::abs(volume - 0.5);

  std::array<T, DESCRIPTOR::d> sorted_normal{{
    std::min(std::min(abs_normal[0], abs_normal[1]), abs_normal[2]) / a_l1,
    0.,
    std::max(std::max(abs_normal[0], abs_normal[1]), abs_normal[2]) / a_l1
  }};
  sorted_normal[1] = std::max(1. - sorted_normal[0] - sorted_normal[2], 0.);

  T d = offsetHelperOpt<T,DESCRIPTOR>(volume_symmetry, sorted_normal);

  return a_l1 * std::copysign(0.5 - d, volume - 0.5);
}

template<typename T, typename DESCRIPTOR>
T FreeSurface3D::calculateSurfaceTensionCurvature(BlockLattice<T,DESCRIPTOR>& blockLattice, int X, int Y, int Z){
  // This is b_z
  std::array<T,DESCRIPTOR::d> normal = FreeSurface3D::computeParkerYoungInterfaceNormal<T,DESCRIPTOR>(blockLattice, X, Y, Z);

  {
    T norm = 0.;
    for(size_t i = 0; i < DESCRIPTOR::d; ++i){
      norm += normal[i] * normal[i];
    }

    norm = std::sqrt(norm);

    if(norm < 1e-12){
      return 0.;
    }

    for(size_t i = 0; i <DESCRIPTOR::d; ++i){
      normal[i] /= norm;
    }
  }

  std::array<T,DESCRIPTOR::d> r_vec{
    0.56270900, 0.32704452, 0.75921047
  };
  /*
  std::array<T,DESCRIPTOR::d> r_vec{
    0.,0.,1.
  };
  */
  std::array<std::array<T,DESCRIPTOR::d>,DESCRIPTOR::d> rotation{{
    {{0., 0., 0.}},
    //{{normal[1], -normal[0], 0.}},
    {{normal[1] * r_vec[2] - normal[2] * r_vec[1], normal[2] * r_vec[0] - normal[0] * r_vec[2], normal[0] * r_vec[1] - normal[1] * r_vec[0]}},
    {{normal[0], normal[1], normal[2]}}
  }};

  // Cross product with (0,0,1) x normal
  // This is b_y

  // (normal[0], normal[1], normal[2])
  
  T cross_norm = 0.;
  for(size_t i = 0; i < DESCRIPTOR::d; ++i){
    cross_norm += rotation[1][i] * rotation[1][i];
  }

  // If too close too each other use the identity matrix
  if(cross_norm > 1e-6){

    cross_norm = std::sqrt(cross_norm);

    for(size_t i = 0; i <DESCRIPTOR::d; ++i){
      rotation[1][i] /= cross_norm;
    }
  }else {

    rotation[1] = {{
      -normal[2],
      0.,
      normal[0]
    }};

    cross_norm = 0.;
    for(size_t i = 0; i < DESCRIPTOR::d; ++i){
      cross_norm += rotation[1][i] * rotation[1][i];
    }

    cross_norm = std::sqrt(cross_norm);

    for(size_t i = 0; i <DESCRIPTOR::d; ++i){
      rotation[1][i] /= cross_norm;
    }
  }


  // Cross product of ((0,0,1) x normal / | (0,0,1) x normal |) x normal
  // This is b_x
  rotation[0] = {{
    rotation[1][1] * normal[2] - rotation[1][2] * normal[1],
    rotation[1][2] * normal[0] - rotation[1][0] * normal[2],
    rotation[1][0] * normal[1] - rotation[1][1] * normal[0]
  }};

  // These three form a matrix and are entered into each row
  // ( b_x )
  // ( b_y )
  // ( b_z )

  constexpr size_t S = 5;
  std::array<std::array<T,S>, S> lq_matrix;
  std::array<T,S> b_rhs;
  for(size_t i = 0; i < S; ++i){
    for(size_t j = 0; j < S; ++j){
      lq_matrix[i][j] = 0.;
    }
    b_rhs[i] = 0.;
  }
  T origin_offset = 0.;
  {
    T fill_level = getClampedEpsilon(blockLattice, X, Y, Z);
    origin_offset = FreeSurface3D::calculateCubeOffsetOpt<T,DESCRIPTOR>(fill_level, normal);
  }

  size_t healthy_interfaces = 0;
  for(int i = -1; i <= 1; ++i){
    for(int j = -1 ; j <= 1; ++j){
      for(int k = -1 ; k <= 1; ++k){
        if( i == 0 && j == 0 && k == 0){
          continue;
        }

        int iXc = X + i;
        int iYc = Y + j;
        int iZc = Z + k;

        if(!FreeSurface3D::isCellType(blockLattice, iXc, iYc, iZc, FreeSurface::Type::Interface) || !FreeSurface3D::hasNeighbour(blockLattice, iXc, iYc, iZc, FreeSurface::Type::Gas)){
          continue;
        }

        ++healthy_interfaces;

        T fill_level = getClampedEpsilon(blockLattice, iXc, iYc, iZc);

        T cube_offset = FreeSurface3D::calculateCubeOffsetOpt<T,DESCRIPTOR>(fill_level, normal);

        std::array<T,DESCRIPTOR::d> pos{static_cast<T>(i),static_cast<T>(j),static_cast<T>(k)};
        std::array<T,DESCRIPTOR::d> r_pos{0.,0.,cube_offset - origin_offset};

        for(size_t a = 0; a < DESCRIPTOR::d; ++a){
          for(size_t b = 0; b < DESCRIPTOR::d; ++b){
            r_pos[a] += rotation[a][b] * pos[b];
          }
        }

        T r_x_2 = r_pos[0] * r_pos[0];
        T r_x_3 = r_x_2 * r_pos[0];
        T r_x_4 = r_x_3 * r_pos[0];

        T r_y_2 = r_pos[1] * r_pos[1];
        T r_y_3 = r_y_2 * r_pos[1];
        T r_y_4 = r_y_3 * r_pos[1];

        T r_x_2_y_2 = r_x_2 * r_y_2;
        T r_x_3_y = r_x_3 * r_pos[1];
        T r_x_2_y = r_x_2 * r_pos[1];

        T r_x_y_3 = r_pos[0] * r_y_3;
        T r_x_y_2 = r_pos[0] * r_y_2;

        T r_x_y = r_pos[0] * r_pos[1];

        lq_matrix[0][0] += r_x_4;
        lq_matrix[1][1] += r_y_4;
        lq_matrix[2][2] += r_x_2_y_2;
        lq_matrix[3][3] += r_x_2;
        lq_matrix[4][4] += r_y_2;

        // skip [1][0] copy later from [2][2]
        lq_matrix[2][0] += r_x_3_y;
        lq_matrix[3][0] += r_x_3;
        lq_matrix[4][0] += r_x_2_y;

        lq_matrix[2][1] += r_x_y_3;
        lq_matrix[3][1] += r_x_y_2;
        lq_matrix[4][1] += r_y_3;

        // skip [3][2] copy from [4][0]
        // skip [4][2] copy from [3][1]

        lq_matrix[4][3] += r_x_y;

        b_rhs[0] +=  r_x_2 * r_pos[2];
        b_rhs[1] +=  r_y_2 * r_pos[2];
        b_rhs[2] +=  r_x_y * r_pos[2];
        b_rhs[3] +=  r_pos[0] * r_pos[2];
        b_rhs[4] +=  r_pos[1] * r_pos[2];
      }
    }
  }

  lq_matrix[1][0] = lq_matrix[2][2];
  lq_matrix[3][2] = lq_matrix[4][0];
  lq_matrix[4][2] = lq_matrix[3][1];

  for(size_t i = 0; i < S; ++i){
    for(size_t j = i + 1; j < S; ++j){
      lq_matrix[i][j] = lq_matrix[j][i];
    }
  }

  // Consider using Thikonov regularization?
  //T alpha = 1e-8;
  T alpha = 0.0;
  for(size_t i = 0; i < S; ++i){
    lq_matrix[i][i] += alpha;
  }

  std::array<T,S> solved_fit = FreeSurface::solvePivotedLU<T,S>(lq_matrix, b_rhs, healthy_interfaces);

  T denom = std::sqrt(1. + solved_fit[3]*solved_fit[3] + solved_fit[4]*solved_fit[4]);
  denom = denom * denom * denom;
  T curvature = ( (1.+solved_fit[4]*solved_fit[4]) * solved_fit[0] + (1. + solved_fit[3]*solved_fit[3] ) * solved_fit[1] - solved_fit[3] * solved_fit[4] * solved_fit[2] ) / denom;

  return std::max(-1., std::min(1., curvature));
}

template<typename T, typename DESCRIPTOR>
T FreeSurface3D::plicInverse(T d_o, const std::array<T,DESCRIPTOR::d>& normal){
  const T n1 = std::min(std::min(std::abs(normal[0]),std::abs(normal[1])),std::abs(normal[2]));
  const T n3 = std::max(std::max(std::abs(normal[0]), std::abs(normal[1])), std::abs(normal[2]));
  const T n2 = std::abs(normal[0])+std::abs(normal[1])+std::abs(normal[2])-n1-n3;

  const T d = 0.5 * (n1+n2+n3) - std::abs(d_o);

  T vol;
  if(std::min(n1+n2,n3) <= d && d <= n3){
    vol = (d-0.5 *(n1+n2))/n3;
  }else if(d < n1){
    vol = std::pow(d,3) / (6. * n1 * n2 * n3);
  }else if(d <= n2){
    vol = (3.0 * d * (d-n1) + std::pow(n1,2))/(6. * n2 * n3);
  }else {
    vol = (std::pow(d,3) - std::pow(d-n1,3) - std::pow(d-n2,3) - std::pow(std::max(0., d-n3),3)) / (6. * n1* n2 * n3);
  }

  return std::copysign(0.5 - vol, d_o) + 0.5;
}

template<typename T, typename DESCRIPTOR>
FreeSurface3D::NeighbourInfo FreeSurface3D::getNeighbourInfo(BlockLattice<T,DESCRIPTOR>& blockLattice, int X, int Y, int Z){
  NeighbourInfo info;

  for(int iPop = 1; iPop < DESCRIPTOR::q; ++iPop){
    int iXc = X + descriptors::c<DESCRIPTOR>(iPop, 0);
    int iYc = Y + descriptors::c<DESCRIPTOR>(iPop, 1);
    int iZc = Z + descriptors::c<DESCRIPTOR>(iPop, 2);

    if(isCellType(blockLattice, iXc, iYc, iZc, FreeSurface::Type::Gas)){
      info.has_gas_neighbours = true;
    }
    else if(isCellType(blockLattice, iXc, iYc, iZc, FreeSurface::Type::Fluid)){
      info.has_fluid_neighbours = true;
    }
    else if(isCellType(blockLattice, iXc, iYc, iZc, FreeSurface::Type::Interface)){
      ++info.interface_neighbours;
    }
  }
  return info;
}

template<typename T, typename DESCRIPTOR>
bool FreeSurface3D::isHealthyInterface(BlockLattice<T,DESCRIPTOR>& blockLattice, int X, int Y, int Z){
  bool has_fluid_neighbours = false;
  bool has_gas_neighbours = false;

  if(!isCellType(blockLattice, X, Y, Z, FreeSurface::Type::Interface)){
    return false;
  }

  for(int iPop = 1; iPop < DESCRIPTOR::q; ++iPop){
    int iXc = X + descriptors::c<DESCRIPTOR>(iPop, 0);
    int iYc = Y + descriptors::c<DESCRIPTOR>(iPop, 1);
    int iZc = Z + descriptors::c<DESCRIPTOR>(iPop, 2);

    if(isCellType(blockLattice, iXc, iYc, iZc, FreeSurface::Type::Gas)){
      has_gas_neighbours = true;

      if(has_fluid_neighbours){
        return true;
      }
    }
    else if(isCellType(blockLattice, iXc, iYc, iZc, FreeSurface::Type::Fluid)){
      has_fluid_neighbours = true;

      if(has_gas_neighbours){
        return true;
      }
    }
  }
  return false;
}

template<typename T, typename DESCRIPTOR>
void FreeSurface3D::setCellType(BlockLattice<T,DESCRIPTOR>& blockLattice, int X, int Y, int Z, const FreeSurface::Type& type){
  blockLattice.get(X,Y,Z).template setField<FreeSurface::CELL_TYPE>(type);
}

template<typename T, typename DESCRIPTOR>
void FreeSurface3D::setCellFlags(BlockLattice<T,DESCRIPTOR>& blockLattice, int X, int Y, int Z, const FreeSurface::Flags& flags){
  blockLattice.get(X,Y,Z).template setField<FreeSurface::CELL_FLAGS>(flags);

}

// LocalProcessor 1

template<typename T, typename DESCRIPTOR>
FreeSurfaceMassFlowPostProcessor3D<T, DESCRIPTOR>::FreeSurfaceMassFlowPostProcessor3D(
  int _x0, int _x1, int _y0, int _y1, int _z0, int _z1,
  const FreeSurface3D::Variables<T,DESCRIPTOR>& vars_):
  x0{_x0}, x1{_x1}, y0{_y0}, y1{_y1}, z0{_z0}, z1{_z1},
  vars{vars_}
{
  this->getName() = "FreeSurfaceMassFlowPostProcessor3D";
  this->_priority = 1;
}

template<typename T, typename DESCRIPTOR>
FreeSurfaceMassFlowPostProcessor3D<T, DESCRIPTOR>::FreeSurfaceMassFlowPostProcessor3D(const FreeSurface3D::Variables<T,DESCRIPTOR>& vars_):
  FreeSurfaceMassFlowPostProcessor3D{0,0,0,0,0,0, vars_}{}

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
// DFs (only replacing from gas cells)

template <typename T, typename DESCRIPTOR>
void FreeSurfaceMassFlowPostProcessor3D<T, DESCRIPTOR>::processSubDomain(
  BlockLattice<T, DESCRIPTOR>& blockLattice, int x0_, int x1_, int y0_, int y1_, int z0_, int z1_){
  int newX0, newX1, newY0, newY1, newZ0, newZ1;
  if (util::intersect(x0, x1, y0, y1, z0, z1, x0_, x1_, y0_, y1_, z0_, z1_, newX0, newX1, newY0, newY1, newZ0, newZ1)) {
    for (int iX = newX0; iX <= newX1; ++iX) {
      for (int iY = newY0; iY <= newY1; ++iY) {
        for (int iZ = newZ0; iZ <= newZ1; ++iZ) {
          /*
          * Minor "hack". Remove all cell flags here, because it is needed in the last processor due to pulling steps in processor 6 and 7
          */
          FreeSurface3D::setCellFlags(blockLattice, iX, iY, iZ, static_cast<FreeSurface::Flags>(0));

          // I'm leaving this in because the algorithm is very unstable and should someone decide to turn on fluid mass exchange they are free to do so.
          // Though it makes the algorithm more unstable in my test cases since more cells are involved.
          /*if(FreeSurface3D::isCellType(blockLattice, iX, iY, iZ, FreeSurface::Type::Fluid )){

            T mass_tmp = blockLattice.get(iX, iY, iZ).template getField<FreeSurface::MASS>();
            for(int iPop = 1; iPop < DESCRIPTOR::q; ++iPop){
              int iXc = iX + descriptors::c<DESCRIPTOR>(iPop, 0);
              int iYc = iY + descriptors::c<DESCRIPTOR>(iPop, 1);
              int iZc = iZ + descriptors::c<DESCRIPTOR>(iPop, 2);
              int iPop_op = descriptors::opposite<DESCRIPTOR>(iPop);
              
              mass_tmp += blockLattice.get(iX, iY, iZ)[iPop_op] - blockLattice.get(iXc, iYc, iZc)[iPop];
            }
            blockLattice.get(iX, iY, iZ).template setField<FreeSurface::MASS>(mass_tmp);
          } else */
          /*
          * This processor only works on interface types
          */
          if( FreeSurface3D::isCellType(blockLattice, iX, iY, iZ, FreeSurface::Type::Interface )){

            T mass_tmp = blockLattice.get(iX, iY, iZ).template getField<FreeSurface::MASS>();

            FreeSurface3D::NeighbourInfo neighbour_info = FreeSurface3D::getNeighbourInfo(blockLattice, iX, iY, iZ);

            for(int iPop = 1; iPop < DESCRIPTOR::q; ++iPop){
              int iXc = iX + descriptors::c<DESCRIPTOR>(iPop, 0);
              int iYc = iY + descriptors::c<DESCRIPTOR>(iPop, 1);
              int iZc = iZ + descriptors::c<DESCRIPTOR>(iPop, 2);
              int iPop_op = descriptors::opposite<DESCRIPTOR>(iPop);

              /*
              * Iterate over neighbours and perform a mass exchange. Interface to fluid are simple cases.
              * Interface to interface has to be symmetric and multiple cases are for artifact reduction
              * from Thuerey
              */

              if( FreeSurface3D::isCellType(blockLattice, iXc, iYc, iZc, FreeSurface::Type::Fluid)) {
                mass_tmp += blockLattice.get(iX, iY, iZ)[iPop_op] - blockLattice.get(iXc, iYc, iZc)[iPop];
              }else if ( FreeSurface3D::isCellType(blockLattice, iXc, iYc, iZc, FreeSurface::Type::Interface)){
                FreeSurface3D::NeighbourInfo neighbour_neighbour_info = FreeSurface3D::getNeighbourInfo(blockLattice, iXc, iYc, iZc);
                T mass_flow = 0.;

                if( !neighbour_info.has_fluid_neighbours){
                  if(!neighbour_neighbour_info.has_fluid_neighbours){
                    if(neighbour_info.interface_neighbours < neighbour_neighbour_info.interface_neighbours){
                      mass_flow = -blockLattice.get(iXc, iYc, iZc)[iPop];// - descriptors::t<T,DESCRIPTOR>(iPop);
                    }else if(neighbour_info.interface_neighbours > neighbour_neighbour_info.interface_neighbours){
                      mass_flow = blockLattice.get(iX, iY, iZ)[iPop_op];// + descriptors::t<T,DESCRIPTOR>(iPop_op);
                    }else{
                      mass_flow = blockLattice.get(iX, iY, iZ)[iPop_op] - blockLattice.get(iXc, iYc, iZc)[iPop];
                    }
                  }else {
                    mass_flow = -blockLattice.get(iXc, iYc, iZc)[iPop];// - descriptors::t<T,DESCRIPTOR>(iPop);
                  }
                }else if(!neighbour_info.has_gas_neighbours){
                  if(!neighbour_neighbour_info.has_gas_neighbours){
                    if(neighbour_info.interface_neighbours < neighbour_neighbour_info.interface_neighbours){
                      mass_flow = blockLattice.get(iX, iY, iZ)[iPop_op];// + descriptors::t<T,DESCRIPTOR>(iPop_op);
                    }else if(neighbour_info.interface_neighbours > neighbour_neighbour_info.interface_neighbours){
                      mass_flow = -blockLattice.get(iXc, iYc, iZc)[iPop];// - descriptors::t<T,DESCRIPTOR>(iPop);
                    }else{
                      mass_flow = blockLattice.get(iX, iY, iZ)[iPop_op] - blockLattice.get(iXc, iYc, iZc)[iPop];
                    }
                  }else {
                    mass_flow = blockLattice.get(iX, iY, iZ)[iPop_op];// + descriptors::t<T,DESCRIPTOR>(iPop_op);
                  }
                }else {
                  if(!neighbour_neighbour_info.has_fluid_neighbours){
                    mass_flow = blockLattice.get(iX, iY, iZ)[iPop_op];// + descriptors::t<T,DESCRIPTOR>(iPop_op);
                  }else if(!neighbour_neighbour_info.has_gas_neighbours){
                    mass_flow = -blockLattice.get(iXc, iYc, iZc)[iPop];// - descriptors::t<T,DESCRIPTOR>(iPop);
                  }else {
                    mass_flow = blockLattice.get(iX, iY, iZ)[iPop_op] - blockLattice.get(iXc, iYc, iZc)[iPop];
                  }
                }

                /*
                * Exchange depends on how filled the interfaces are
                */
                mass_tmp += mass_flow * 0.5 * (FreeSurface3D::getClampedEpsilon(blockLattice, iX, iY, iZ) + FreeSurface3D::getClampedEpsilon(blockLattice, iXc, iYc, iZc));

              }
            }

            blockLattice.get(iX, iY, iZ).template setField<FreeSurface::MASS>(mass_tmp);


            // Former 2 Step

            // Because I need the distribution functions of the last step I will write results in a temporary
            // array, before copying it back into the DFs

            std::array<T, DESCRIPTOR::q> dfs;

            T curvature = 0.;

            if(vars.has_surface_tension){
              /*
              *
              */
              if(neighbour_info.has_gas_neighbours){
                curvature = FreeSurface3D::calculateSurfaceTensionCurvature(blockLattice, iX, iY, iZ);
              }
            }

            T gas_pressure = 1. - 6. * vars.surface_tension_parameter * curvature;

            // std::array<T,DESCRIPTOR::d> normal = FreeSurface3D::computeInterfaceNormal(blockLattice, iX, iY, iZ);

            for(int iPop=1; iPop < DESCRIPTOR::q; iPop++) {

              int iXc = iX + descriptors::c<DESCRIPTOR>(iPop,0);
              int iYc = iY + descriptors::c<DESCRIPTOR>(iPop,1);
              int iZc = iZ + descriptors::c<DESCRIPTOR>(iPop,2);
              int iPop_op = descriptors::opposite<DESCRIPTOR>(iPop);

              /*
              * Replace incoming streaming distributing functions from gas cells
              */

              if ( FreeSurface3D::isCellType(blockLattice, iXc, iYc, iZc, FreeSurface::Type::Gas )) {
                // assume atmosphere pressure as 1
                Vector<T, DESCRIPTOR::d> u_vel = blockLattice.get(iX,iY,iZ).template getField<FreeSurface::PREVIOUS_VELOCITY>();
                T u[DESCRIPTOR::d];
                for(size_t u_i = 0; u_i < DESCRIPTOR::d; ++u_i){
                  u[u_i] = u_vel[u_i];
                }
                util::normSqr<T,DESCRIPTOR::d>(u);
                dfs[iPop_op] = equilibrium<DESCRIPTOR>::secondOrder(iPop, gas_pressure, u)
                                                   + equilibrium<DESCRIPTOR>::secondOrder(iPop_op, gas_pressure, u)
                                                   - blockLattice.get(iXc,iYc,iZc)[iPop];
              }else {
                dfs[iPop_op] = blockLattice.get(iX, iY, iZ)[iPop_op];
              }
            }

            for(int iPop=1; iPop<DESCRIPTOR::q; iPop++) {

              blockLattice.get(iX,iY,iZ)[iPop] = dfs[iPop];

            }

            // Former 3 Step
            /*
            * Based on the mass calculation, flag this interface cell as toFluid or toGas if set boundaries are met
            */
            T rho = blockLattice.get(iX,iY,iZ).computeRho();

            // Check if transition needs to happen.
            if ( mass_tmp < -vars.transition * rho || (mass_tmp < vars.lonely_threshold * rho && !neighbour_info.has_fluid_neighbours) ){
              FreeSurface3D::setCellFlags(blockLattice, iX, iY, iZ, FreeSurface::Flags::ToGas);
            }
            else if ( mass_tmp > (1. + vars.transition)*rho  || ( mass_tmp > (1-vars.lonely_threshold) * rho && !neighbour_info.has_gas_neighbours) ){
              FreeSurface3D::setCellFlags(blockLattice, iX, iY, iZ, FreeSurface::Flags::ToFluid);
            }else if(vars.drop_isolated_cells && neighbour_info.interface_neighbours == 0){
              if(!neighbour_info.has_gas_neighbours){
                FreeSurface3D::setCellFlags(blockLattice, iX, iY, iZ, FreeSurface::Flags::ToFluid);
              }else if(!neighbour_info.has_fluid_neighbours){
                //FreeSurface3D::setCellFlags(blockLattice, iX, iY, FreeSurface::Flags::ToGas);
              }
            }
          }
        }
      }
    }
  }
}

template <typename T, typename DESCRIPTOR>
void FreeSurfaceMassFlowPostProcessor3D<T, DESCRIPTOR>::process(
  BlockLattice<T, DESCRIPTOR>& blockLattice
){
  processSubDomain(blockLattice, x0, x1, y0, y1, z0, z1);
}

// Processor 4
template<typename T, typename DESCRIPTOR>
FreeSurfaceToFluidCellConversionPostProcessor3D<T,DESCRIPTOR>::FreeSurfaceToFluidCellConversionPostProcessor3D(int _x0, int _x1, int _y0, int _y1, int _z0, int _z1, const FreeSurface3D::Variables<T,DESCRIPTOR>& vars_):
  x0{_x0},
  x1{_x1},
  y0{_y0},
  y1{_y1},
  z0{_z0},
  z1{_z1},
  vars{vars_}
{
  this->getName() = "FreeSurfaceToFluidCellConversionPostProcessor3D";
  this->_priority = 4;
}
template<typename T, typename DESCRIPTOR>
FreeSurfaceToFluidCellConversionPostProcessor3D<T,DESCRIPTOR>::FreeSurfaceToFluidCellConversionPostProcessor3D(const FreeSurface3D::Variables<T,DESCRIPTOR>& vars_):
  FreeSurfaceToFluidCellConversionPostProcessor3D{0,0,0,0,0,0, vars_}
{}

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
// DFs ( only toFluid )
template <typename T, typename DESCRIPTOR>
void FreeSurfaceToFluidCellConversionPostProcessor3D<T, DESCRIPTOR>::processSubDomain(
  BlockLattice<T, DESCRIPTOR>& blockLattice, int x0_, int x1_, int y0_, int y1_, int z0_, int z1_){
  int newX0, newX1, newY0, newY1, newZ0, newZ1;

  if(util::intersect(x0, x1, y0, y1, z0, z1, x0_, x1_, y0_, y1_, z0_, z1_, newX0, newX1, newY0, newY1, newZ0, newZ1)){

    for(int iX = newX0; iX <= newX1; ++iX){
      for(int iY = newY0; iY <= newY1; ++iY){
        for(int iZ = newZ0; iZ <= newZ1; ++iZ){

          /*
          * Initializing new interface cells with DFs from surrounding fluid and interface cells
          * Just takes the arithmetic average.
          */

          /*
          * If Gas has a ToFluid cell next to it then set a new interface flag
          */
          if(FreeSurface3D::isCellType(blockLattice, iX, iY, iZ, FreeSurface::Type::Gas)){
            if(FreeSurface3D::hasNeighbourFlags(blockLattice, iX, iY, iZ, FreeSurface::Flags::ToFluid)){
              FreeSurface3D::setCellFlags(blockLattice, iX, iY, iZ, FreeSurface::Flags::NewInterface);
              T rho_avg = 0.;
              T u_avg[DESCRIPTOR::d] = {0., 0.};

              size_t ctr = 0;

              for(int iPop=1; iPop<DESCRIPTOR::q; iPop++) {
                int iXc = iX+descriptors::c<DESCRIPTOR>(iPop,0);
                int iYc = iY+descriptors::c<DESCRIPTOR>(iPop,1);
                int iZc = iZ+descriptors::c<DESCRIPTOR>(iPop,2);

                if(FreeSurface3D::isCellType(blockLattice, iXc, iYc, iZc, FreeSurface::Type::Fluid) || FreeSurface3D::isCellType(blockLattice, iXc, iYc, iZc, FreeSurface::Type::Interface)){
                 
                  T rho_tmp = 0.;
                  T u_tmp[DESCRIPTOR::d] = {0., 0.};
                  ++ctr;
                  blockLattice.get(iXc,iYc,iZc).computeRhoU(rho_tmp, u_tmp);
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
              
              blockLattice.get(iX,iY,iZ).iniEquilibrium(rho_avg, u_avg);
            }
          }else if(FreeSurface3D::hasCellFlags(blockLattice, iX, iY, iZ, FreeSurface::Flags::ToGas)){
            /*
            * If a toGas cell has a neighbouring toFluid cell, unset the toGas flag
            */
            if(FreeSurface3D::hasNeighbourFlags(blockLattice, iX, iY, iZ, FreeSurface::Flags::ToFluid)){
              FreeSurface3D::setCellFlags(blockLattice, iX, iY, iZ, static_cast<FreeSurface::Flags>(0));
            }
          }
        }
      }
    }
  }
}

template <typename T, typename DESCRIPTOR>
void FreeSurfaceToFluidCellConversionPostProcessor3D<T, DESCRIPTOR>::process(
  BlockLattice<T, DESCRIPTOR>& blockLattice
){
  processSubDomain(blockLattice, x0, x1, y0, y1, z0, z1);
}

// LocalProcessor 5
template<typename T, typename DESCRIPTOR>
FreeSurfaceToGasCellConversionPostProcessor3D<T,DESCRIPTOR>::FreeSurfaceToGasCellConversionPostProcessor3D(
  int _x0, int _x1, int _y0, int _y1, int _z0, int _z1,
  const FreeSurface3D::Variables<T,DESCRIPTOR>& vars_):
  x0{_x0},
  x1{_x1},
  y0{_y0},
  y1{_y1},
  z0{_z0},
  z1{_z1},
  vars{vars_}
{
  this->getName() = "FreeSurfaceToGasCellConversionPostProcessor3D";
  this->_priority = 5;
}
template<typename T, typename DESCRIPTOR>
FreeSurfaceToGasCellConversionPostProcessor3D<T,DESCRIPTOR>::FreeSurfaceToGasCellConversionPostProcessor3D(const FreeSurface3D::Variables<T,DESCRIPTOR>& vars_):
  FreeSurfaceToGasCellConversionPostProcessor3D{0,0,0,0,0,0, vars_}
{}

// Read
// range 0
// CellType
// DFs

// range 1
// CellFlags (only to gas flags on interface cells)

// Write (always range 0)
// CellFlags (only on fluid cells)
// Mass

template <typename T, typename DESCRIPTOR>
void FreeSurfaceToGasCellConversionPostProcessor3D<T, DESCRIPTOR>::processSubDomain(
  BlockLattice<T, DESCRIPTOR>& blockLattice, int x0_, int x1_, int y0_, int y1_, int z0_, int z1_){
  int newX0, newX1, newY0, newY1, newZ0, newZ1;

  if(util::intersect(x0, x1, y0, y1, z0, z1, x0_, x1_, y0_, y1_, z0_, z1_, newX0, newX1, newY0, newY1, newZ0, newZ1)){

    for(int iX = newX0; iX <= newX1; ++iX){
      for(int iY = newY0; iY <= newY1; ++iY){
        for(int iZ = newZ0; iZ <= newZ1; ++iZ){
        
          // Parallel version. Write locally.
          /*
          * For the to be converted toGas cells, set the neighbours to interface cells
          */
          if(FreeSurface3D::isCellType(blockLattice, iX, iY, iZ, FreeSurface::Type::Fluid)){
            if(FreeSurface3D::hasNeighbourFlags(blockLattice, iX, iY, iZ, FreeSurface::Flags::ToGas)){
              FreeSurface3D::setCellFlags(blockLattice, iX, iY, iZ, FreeSurface::Flags::NewInterface);

              T rho = blockLattice.get(iX,iY,iZ).computeRho();
              blockLattice.get(iX, iY, iZ).template setField<FreeSurface::MASS>(rho);
            }
          }
        }
      }
    }
  }

}

template <typename T, typename DESCRIPTOR>
void FreeSurfaceToGasCellConversionPostProcessor3D<T, DESCRIPTOR>::process(
  BlockLattice<T, DESCRIPTOR>& blockLattice
){
  processSubDomain(blockLattice, x0, x1, y0, y1, z0, z1);
}
// LocalProcessor 6
template<typename T, typename DESCRIPTOR>
FreeSurfaceMassExcessPostProcessor3D<T,DESCRIPTOR>::FreeSurfaceMassExcessPostProcessor3D(
  int _x0, int _x1, int _y0, int _y1, int _z0, int _z1, const FreeSurface3D::Variables<T,DESCRIPTOR>& vars_):
  x0{_x0},
  x1{_x1},
  y0{_y0},
  y1{_y1},
  z0{_z0},
  z1{_z1},
  vars{vars_}
{
  this->getName() = "FreeSurfaceMassExcessPostProcessor3D";
  this->_priority = 6;
}
template<typename T, typename DESCRIPTOR>
FreeSurfaceMassExcessPostProcessor3D<T,DESCRIPTOR>::FreeSurfaceMassExcessPostProcessor3D(const FreeSurface3D::Variables<T,DESCRIPTOR>& vars_):
  FreeSurfaceMassExcessPostProcessor3D{0,0,0,0,0,0, vars_}
{}

// Read
// range 0
// CellType
// CellFlags
// DFs
// Mass

// range 1
// CellType
// CellFlags (only to fluid flags on interface cells)

// Write (always range 0)
// Mass
// TempMassExchange
template <typename T, typename DESCRIPTOR>
void FreeSurfaceMassExcessPostProcessor3D<T, DESCRIPTOR>::processSubDomain(
  BlockLattice<T, DESCRIPTOR>& blockLattice, int x0_, int x1_, int y0_, int y1_, int z0_, int z1_){
  int newX0, newX1, newY0, newY1, newZ0, newZ1;

  if(util::intersect(x0, x1, y0, y1, z0, z1, x0_, x1_, y0_, y1_, z0_, z1_, newX0, newX1, newY0, newY1, newZ0, newZ1)){

    for(int iX = newX0; iX <= newX1; ++iX){
      for(int iY = newY0; iY <= newY1; ++iY){
        for(int iZ = newZ0; iZ <= newZ1; ++iZ){
          if( !FreeSurface3D::isCellType(blockLattice, iX, iY, iZ, FreeSurface::Type::Interface) ){
            continue;
          }

          T rho = blockLattice.get(iX,iY,iZ).computeRho();
          T mass = blockLattice.get(iX,iY,iZ).template getField<FreeSurface::MASS>( );
          T mass_excess = 0.;

          // redistribute excess mass

          /// @hint EPSILON of neighbours used here

          /// @hint Mass can be set in this processor, but not epsilon since it is needed for the normal computation. epsilon is set in the next processor
          if(FreeSurface3D::hasCellFlags(blockLattice, iX, iY, iZ, FreeSurface::Flags::ToGas)){

            mass_excess = mass;
            blockLattice.get(iX, iY, iZ).template setField<FreeSurface::MASS>( 0. );

          }else if (FreeSurface3D::hasCellFlags(blockLattice, iX, iY, iZ, FreeSurface::Flags::ToFluid)){

            mass_excess = mass - rho;
            blockLattice.get(iX, iY, iZ).template setField<FreeSurface::MASS>( rho );


          }else {
            continue;
          }

          size_t product_total = 0;

          for(int iPop=1; iPop<DESCRIPTOR::q; iPop++) {
            int iXc = iX+descriptors::c<DESCRIPTOR>(iPop,0);
            int iYc = iY+descriptors::c<DESCRIPTOR>(iPop,1);
            int iZc = iZ+descriptors::c<DESCRIPTOR>(iPop,2);

            // Thuerey Paper says we can't use new interface cells
            // or flagged 
            // And it's kind of logical. It has a state which is new and is not be deterministicly defined at this timestep
            // 

            if( (FreeSurface3D::isCellType(blockLattice, iXc, iYc, iZc, FreeSurface::Type::Interface) && (!FreeSurface3D::hasCellFlags(blockLattice, iXc, iYc, iZc,
              static_cast<FreeSurface::Flags>(255)) /*|| FreeSurface3D::hasCellFlags(blockLattice, iXc, iYc, iZc, FreeSurface::Flags::ToFluid)*/))
              /*|| FreeSurface3D::isCellType(blockLattice, iXc, iYc, iZc, FreeSurface::Type::Fluid)*/
            ){
              ++product_total;
            }
          }
          
          Vector<T,DESCRIPTOR::q> mass_excess_vector;
          mass_excess_vector[0] = 0.;
          
          if(product_total > 0){
            T product_fraction = 1. / product_total;

            for(int iPop=1; iPop < DESCRIPTOR::q; iPop++) {
              int iXc = iX+descriptors::c<DESCRIPTOR>(iPop,0);
              int iYc = iY+descriptors::c<DESCRIPTOR>(iPop,1);
              int iZc = iZ+descriptors::c<DESCRIPTOR>(iPop,2);
              if( (FreeSurface3D::isCellType(blockLattice, iXc, iYc, iZc, FreeSurface::Type::Interface) && (!FreeSurface3D::hasCellFlags(blockLattice, iXc, iYc, iZc,
              static_cast<FreeSurface::Flags>(255)) /*|| FreeSurface3D::hasCellFlags(blockLattice, iXc, iYc, iZc, FreeSurface::Flags::ToFluid)*/))
              /*|| FreeSurface3D::isCellType(blockLattice, iXc, iYc, iZc, FreeSurface::Type::Fluid)*/
                ){
                  mass_excess_vector[iPop] = mass_excess * product_fraction;
                }else{
                  mass_excess_vector[iPop] = 0.;
                }
            }
            blockLattice.get(iX, iY, iZ).template setField<FreeSurface::TEMP_MASS_EXCHANGE>( mass_excess_vector );
          }else {
            mass_excess_vector[0] = mass_excess;
            for(int iPop=1; iPop < DESCRIPTOR::q; iPop++) {
              mass_excess_vector[iPop] = 0.;
            }
            blockLattice.get(iX, iY, iZ).template setField<FreeSurface::TEMP_MASS_EXCHANGE>( mass_excess_vector );
          }
        }
      }
    }
  }

}

template <typename T, typename DESCRIPTOR>
void FreeSurfaceMassExcessPostProcessor3D<T, DESCRIPTOR>::process(
  BlockLattice<T, DESCRIPTOR>& blockLattice
){
  processSubDomain(blockLattice, x0, x1, y0, y1, z0, z1);
}

// LocalProcessor 7
template<typename T, typename DESCRIPTOR>
FreeSurfaceFinalizeConversionPostProcessor3D<T,DESCRIPTOR>::FreeSurfaceFinalizeConversionPostProcessor3D(
  int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, const FreeSurface3D::Variables<T,DESCRIPTOR>& vars_):
  x0{x0_},
  x1{x1_},
  y0{y0_},
  y1{y1_},
  z0{z0_},
  z1{z1_},
  vars{vars_}
{
  this->getName() = "FreeSurfaceFinalizeConversionPostProcessor3D";
  this->_priority = 7;
}

template<typename T, typename DESCRIPTOR>
FreeSurfaceFinalizeConversionPostProcessor3D<T,DESCRIPTOR>::FreeSurfaceFinalizeConversionPostProcessor3D(const FreeSurface3D::Variables<T,DESCRIPTOR>& vars_):
  FreeSurfaceFinalizeConversionPostProcessor3D{0,0,0,0,0,0, vars_}
{}

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
// DFs

template <typename T, typename DESCRIPTOR>
void FreeSurfaceFinalizeConversionPostProcessor3D<T, DESCRIPTOR>::processSubDomain(
  BlockLattice<T, DESCRIPTOR>& blockLattice, int x0_, int x1_, int y0_, int y1_, int z0_, int z1_){
  int newX0, newX1, newY0, newY1, newZ0, newZ1;

  if(util::intersect(x0, x1, y0, y1, z0, z1, x0_, x1_, y0_, y1_, z0_, z1_, newX0, newX1, newY0, newY1, newZ0, newZ1)){

    for(int iX = newX0; iX <= newX1; ++iX){
      for(int iY = newY0; iY <= newY1; ++iY){
        for(int iZ = newZ0; iZ <= newZ1; ++iZ){

          FreeSurface::Flags flags = static_cast<FreeSurface::Flags>(blockLattice.get(iX,iY, iZ).template getField<FreeSurface::CELL_FLAGS>());

          switch(flags){
            case FreeSurface::Flags::ToFluid:
            {
              /// @hint moved flag removal to processor 1 without any negative effects
              FreeSurface3D::setCellType(blockLattice, iX, iY, iZ, FreeSurface::Type::Fluid);
              blockLattice.get(iX, iY, iZ).template setField<FreeSurface::EPSILON>( 1. );
              T mass_tmp = blockLattice.get(iX, iY, iZ).template getField<FreeSurface::MASS>();
              mass_tmp += blockLattice.get(iX, iY, iZ).template getField<FreeSurface::TEMP_MASS_EXCHANGE>()[0];
              blockLattice.get(iX, iY, iZ).template setField<FreeSurface::MASS>(mass_tmp);
            }
            break;
            case FreeSurface::Flags::ToGas:
            {
              /// @hint moved flag removal to processor 1 without any negative effects
              FreeSurface3D::setCellType(blockLattice, iX, iY, iZ, FreeSurface::Type::Gas);
              blockLattice.get(iX, iY, iZ).template setField<FreeSurface::EPSILON>( 0. );
              T mass_tmp = blockLattice.get(iX, iY, iZ).template getField<FreeSurface::MASS>();
              mass_tmp += blockLattice.get(iX, iY, iZ).template getField<FreeSurface::TEMP_MASS_EXCHANGE>()[0];
              blockLattice.get(iX, iY, iZ).template setField<FreeSurface::MASS>(mass_tmp);
            }
            break;
            case FreeSurface::Flags::NewInterface:
            {
              FreeSurface3D::setCellType(blockLattice, iX, iY, iZ, FreeSurface::Type::Interface);
            }
            break;
          }

          FreeSurface::Type type = static_cast<FreeSurface::Type>(blockLattice.get(iX,iY,iZ).template getField<FreeSurface::CELL_TYPE>());

          switch(type){
            case FreeSurface::Type::Interface:
            {
              T collected_excess = 0.; // Doesnt happen // blockLattice.get(iX, iY, iZ).template getField<FreeSurface::TEMP_MASS_EXCHANGE> [0];
              if(type == FreeSurface::Type::Interface){

                for(int iPop = 1; iPop < DESCRIPTOR::q; ++iPop){
                    int iXc = iX+descriptors::c<DESCRIPTOR>(iPop,0);
                    int iYc = iY+descriptors::c<DESCRIPTOR>(iPop,1);
                    int iZc = iZ+descriptors::c<DESCRIPTOR>(iPop,2);

                    if(FreeSurface3D::hasCellFlags(blockLattice, iXc, iYc, iZc,
                          FreeSurface::Flags::ToFluid
                        | FreeSurface::Flags::ToGas)
                    ){

                      int iPop_op = descriptors::opposite<DESCRIPTOR>(iPop);

                      collected_excess += blockLattice.get(iXc, iYc, iZc).template getField<FreeSurface::TEMP_MASS_EXCHANGE>()[iPop_op];
                    }


                }
              }
            
              T mass_tmp = blockLattice.get(iX, iY, iZ).template getField<FreeSurface::MASS>();

              mass_tmp += collected_excess;

              T rho;
              T u_tmp[DESCRIPTOR::d];
              blockLattice.get(iX, iY, iZ).computeRhoU(rho, u_tmp);
              Vector<T,DESCRIPTOR::d> u_vel{u_tmp[0], u_tmp[1], u_tmp[2]};

              blockLattice.get(iX, iY, iZ).template setField<FreeSurface::EPSILON>( mass_tmp / rho );
              blockLattice.get(iX, iY, iZ).template setField<FreeSurface::MASS>(mass_tmp);
           
              blockLattice.get(iX, iY, iZ).template setField<FreeSurface::PREVIOUS_VELOCITY>(u_vel);
              
            }
            break;
            case FreeSurface::Type::Fluid:
            {
              T collected_excess = 0.;

              for(int iPop = 1; iPop < DESCRIPTOR::q; ++iPop){
                  int iXc = iX+descriptors::c<DESCRIPTOR>(iPop,0);
                  int iYc = iY+descriptors::c<DESCRIPTOR>(iPop,1);
                  int iZc = iZ+descriptors::c<DESCRIPTOR>(iPop,2);

                  if(FreeSurface3D::hasCellFlags(blockLattice, iXc, iYc, iZc,
                        FreeSurface::Flags::ToFluid
                      | FreeSurface::Flags::ToGas)
                  ){
                    int iPop_op = descriptors::opposite<DESCRIPTOR>(iPop);
                    collected_excess += blockLattice.get(iXc, iYc, iZc).template getField<FreeSurface::TEMP_MASS_EXCHANGE>()[iPop_op];
                  }
              }
            
              T mass_tmp = blockLattice.get(iX, iY, iZ).template getField<FreeSurface::MASS>();
              mass_tmp += collected_excess;
              blockLattice.get(iX, iY, iZ).template setField<FreeSurface::MASS>(mass_tmp);
            }
            break;
            case FreeSurface::Type::Gas:
            case FreeSurface::Type::Solid:
            break;
            default:
            break;
          }
        }
      }
    }
  }
}

template <typename T, typename DESCRIPTOR>
void FreeSurfaceFinalizeConversionPostProcessor3D<T, DESCRIPTOR>::process(
  BlockLattice<T, DESCRIPTOR>& blockLattice
){
  processSubDomain(blockLattice, x0, x1, y0, y1, z0, z1);
}

/*
*
*/

// Generator 1
template<typename T, typename DESCRIPTOR>
FreeSurfaceMassFlowGenerator3D<T, DESCRIPTOR>::FreeSurfaceMassFlowGenerator3D(
  int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, const FreeSurface3D::Variables<T,DESCRIPTOR>& vars_
)
  : PostProcessorGenerator3D<T, DESCRIPTOR>(x0_, x1_, y0_, y1_, z0_, z1_),
  vars{vars_}
{}


template<typename T, typename DESCRIPTOR>
FreeSurfaceMassFlowGenerator3D<T, DESCRIPTOR>::FreeSurfaceMassFlowGenerator3D(const FreeSurface3D::Variables<T,DESCRIPTOR>& vars_):
  FreeSurfaceMassFlowGenerator3D<T,DESCRIPTOR>(0,0,0,0,0,0, vars_){}


template<typename T, typename DESCRIPTOR>
PostProcessor3D<T,DESCRIPTOR>* FreeSurfaceMassFlowGenerator3D<T, DESCRIPTOR>::generate() const {
  return new FreeSurfaceMassFlowPostProcessor3D<T,DESCRIPTOR>(this->x0, this->x1, this->y0, this->y1, this->z0, this->z1, this->vars);
}

template<typename T, typename DESCRIPTOR>
PostProcessorGenerator3D<T,DESCRIPTOR>* FreeSurfaceMassFlowGenerator3D<T, DESCRIPTOR>::clone() const {
  return new FreeSurfaceMassFlowGenerator3D<T,DESCRIPTOR>(*this);
}

// Generator 4
template<typename T, typename DESCRIPTOR>
FreeSurfaceToFluidCellConversionGenerator3D<T, DESCRIPTOR>::FreeSurfaceToFluidCellConversionGenerator3D(
  int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, const FreeSurface3D::Variables<T,DESCRIPTOR>& vars_
)
  : PostProcessorGenerator3D<T, DESCRIPTOR>(x0_, x1_, y0_, y1_, z0_, z1_),
  vars{vars_}
{}


template<typename T, typename DESCRIPTOR>
FreeSurfaceToFluidCellConversionGenerator3D<T, DESCRIPTOR>::FreeSurfaceToFluidCellConversionGenerator3D(const FreeSurface3D::Variables<T,DESCRIPTOR>& vars_):
  FreeSurfaceToFluidCellConversionGenerator3D<T,DESCRIPTOR>(0,0,0,0,0,0, vars_){}


template<typename T, typename DESCRIPTOR>
PostProcessor3D<T,DESCRIPTOR>* FreeSurfaceToFluidCellConversionGenerator3D<T, DESCRIPTOR>::generate() const {
  return new FreeSurfaceToFluidCellConversionPostProcessor3D<T,DESCRIPTOR>(this->x0, this->x1, this->y0, this->y1, this->z0, this->z1, this->vars);
}

template<typename T, typename DESCRIPTOR>
PostProcessorGenerator3D<T,DESCRIPTOR>* FreeSurfaceToFluidCellConversionGenerator3D<T, DESCRIPTOR>::clone() const {
  return new FreeSurfaceToFluidCellConversionGenerator3D<T,DESCRIPTOR>(*this);
}

// Generator 5
template<typename T, typename DESCRIPTOR>
FreeSurfaceToGasCellConversionGenerator3D<T, DESCRIPTOR>::FreeSurfaceToGasCellConversionGenerator3D(
  int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, const FreeSurface3D::Variables<T,DESCRIPTOR>& vars_
)
  : PostProcessorGenerator3D<T, DESCRIPTOR>(x0_, x1_, y0_, y1_, z0_, z1_),
  vars{vars_}
{}


template<typename T, typename DESCRIPTOR>
FreeSurfaceToGasCellConversionGenerator3D<T, DESCRIPTOR>::FreeSurfaceToGasCellConversionGenerator3D(const FreeSurface3D::Variables<T,DESCRIPTOR>& vars_):
  FreeSurfaceToGasCellConversionGenerator3D<T,DESCRIPTOR>(0,0,0,0,0,0, vars_){}


template<typename T, typename DESCRIPTOR>
PostProcessor3D<T,DESCRIPTOR>* FreeSurfaceToGasCellConversionGenerator3D<T, DESCRIPTOR>::generate() const {
  return new FreeSurfaceToGasCellConversionPostProcessor3D<T,DESCRIPTOR>(this->x0, this->x1, this->y0, this->y1, this->z0, this->z1, this->vars);
}

template<typename T, typename DESCRIPTOR>
PostProcessorGenerator3D<T,DESCRIPTOR>* FreeSurfaceToGasCellConversionGenerator3D<T, DESCRIPTOR>::clone() const {
  return new FreeSurfaceToGasCellConversionGenerator3D<T,DESCRIPTOR>(*this);
}

// Generator 6
template<typename T, typename DESCRIPTOR>
FreeSurfaceMassExcessGenerator3D<T, DESCRIPTOR>::FreeSurfaceMassExcessGenerator3D(
  int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, const FreeSurface3D::Variables<T,DESCRIPTOR>& vars_
)
  : PostProcessorGenerator3D<T, DESCRIPTOR>(x0_, x1_, y0_, y1_, z0_, z1_),
  vars{vars_}
{}


template<typename T, typename DESCRIPTOR>
FreeSurfaceMassExcessGenerator3D<T, DESCRIPTOR>::FreeSurfaceMassExcessGenerator3D(const FreeSurface3D::Variables<T,DESCRIPTOR>& vars_):
  FreeSurfaceMassExcessGenerator3D<T,DESCRIPTOR>(0,0,0,0,0,0,vars_){}


template<typename T, typename DESCRIPTOR>
PostProcessor3D<T,DESCRIPTOR>* FreeSurfaceMassExcessGenerator3D<T, DESCRIPTOR>::generate() const {
  return new FreeSurfaceMassExcessPostProcessor3D<T,DESCRIPTOR>(this->x0, this->x1, this->y0, this->y1, this->z0, this->z1, this->vars);
}

template<typename T, typename DESCRIPTOR>
PostProcessorGenerator3D<T,DESCRIPTOR>* FreeSurfaceMassExcessGenerator3D<T, DESCRIPTOR>::clone() const {
  return new FreeSurfaceMassExcessGenerator3D<T,DESCRIPTOR>(*this);
}

// Generator 7
template<typename T, typename DESCRIPTOR>
FreeSurfaceFinalizeConversionGenerator3D<T, DESCRIPTOR>::FreeSurfaceFinalizeConversionGenerator3D(
  int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, const FreeSurface3D::Variables<T,DESCRIPTOR>& vars_
)
  : PostProcessorGenerator3D<T, DESCRIPTOR>(x0_, x1_, y0_, y1_, z0_, z1_),
  vars{vars_}
{}

template<typename T, typename DESCRIPTOR>
FreeSurfaceFinalizeConversionGenerator3D<T, DESCRIPTOR>::FreeSurfaceFinalizeConversionGenerator3D(const FreeSurface3D::Variables<T,DESCRIPTOR>& vars_):
  FreeSurfaceFinalizeConversionGenerator3D<T,DESCRIPTOR>(0,0,0,0,0,0, vars_){}


template<typename T, typename DESCRIPTOR>
PostProcessor3D<T,DESCRIPTOR>* FreeSurfaceFinalizeConversionGenerator3D<T, DESCRIPTOR>::generate() const {
  return new FreeSurfaceFinalizeConversionPostProcessor3D<T,DESCRIPTOR>(this->x0, this->x1, this->y0, this->y1, this->z0, this->z1, this->vars);
}

template<typename T, typename DESCRIPTOR>
PostProcessorGenerator3D<T,DESCRIPTOR>* FreeSurfaceFinalizeConversionGenerator3D<T, DESCRIPTOR>::clone() const {
  return new FreeSurfaceFinalizeConversionGenerator3D<T,DESCRIPTOR>(*this);
}

template<typename T, typename DESCRIPTOR>
FreeSurface3DSetup<T,DESCRIPTOR>::FreeSurface3DSetup(SuperLattice<T, DESCRIPTOR>& sLattice,
    const FreeSurface3D::Variables<T,DESCRIPTOR>& vars_)
:
  vars{vars_},
  sLattice{sLattice},
  mass_flow{std::make_unique<FreeSurfaceMassFlowGenerator3D<T,DESCRIPTOR>>(vars)},
  to_fluid{std::make_unique<FreeSurfaceToFluidCellConversionGenerator3D<T,DESCRIPTOR>>(vars)},
  to_gas{std::make_unique<FreeSurfaceToGasCellConversionGenerator3D<T,DESCRIPTOR>>(vars)},
  mass_excess{std::make_unique<FreeSurfaceMassExcessGenerator3D<T,DESCRIPTOR>>(vars)},
  finalize_conversion{std::make_unique<FreeSurfaceFinalizeConversionGenerator3D<T,DESCRIPTOR>>(vars)}
{}

template<typename T, typename DESCRIPTOR>
void FreeSurface3DSetup<T,DESCRIPTOR>::addPostProcessor(){
  sLattice.template addPostProcessor<FreeSurface::Stage0>(*mass_flow);
  sLattice.template addPostProcessor<FreeSurface::Stage1>(*to_fluid);
  sLattice.template addPostProcessor<FreeSurface::Stage2>(*to_gas);
  sLattice.template addPostProcessor<FreeSurface::Stage3>(*mass_excess);
  sLattice.template addPostProcessor<FreeSurface::Stage4>(*finalize_conversion);

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
