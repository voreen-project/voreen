/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Jonas Latt, 2020 Adrian Kummerlaender
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
 * Set of functions commonly used in LB computations
 *  -- header file
 */
#ifndef UTIL_H
#define UTIL_H

#include <sstream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <type_traits>

#include "core/baseType.h"
#include "core/vector.h"
#include "utilities/vectorHelpers.h"
#include "meta.h"
#include "utilities/omath.h"
#include "utilities/oalgorithm.h"
#include "dynamics/descriptorFunction.h"

namespace olb {

namespace util {

// *INDENT-OFF*

template<typename T> T norm(const std::vector<T>& a);

template <typename T>
inline int sign(T val)
{
  return (0 < val) - (val < 0);
}

template <typename T>
inline bool aligned_to_x(const std::vector<T>& vec)
{
  return (vec[0]!=0 and vec[1]==0 and vec[2]==0);
}

template <typename T>
inline bool aligned_to_y(const std::vector<T>& vec)
{
  return (vec[0]==0 and vec[1]!=0 and vec[2]==0);
}

template <typename T>
inline bool aligned_to_z(const std::vector<T>& vec)
{
  return (vec[0]==0 and vec[1]==0 and vec[2]!=0);
}

template <typename T>
inline bool aligned_to_grid(const std::vector<T>& vec)
{
  return (aligned_to_x<T>(vec) or
          aligned_to_y<T>(vec) or
          aligned_to_z<T>(vec));
}





inline bool intersect (
  int x0, int x1, int y0, int y1,
  int x0_, int x1_, int y0_, int y1_,
  int& newX0, int& newX1, int& newY0, int& newY1 )
{
  newX0 = util::max(x0,x0_);
  newY0 = util::max(y0,y0_);

  newX1 = util::min(x1,x1_);
  newY1 = util::min(y1,y1_);

  return newX1>=newX0 && newY1>=newY0;
}

inline bool intersect (
  int x0, int x1, int y0, int y1, int z0, int z1,
  int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
  int& newX0, int& newX1, int& newY0, int& newY1, int& newZ0, int& newZ1 )
{
  newX0 = util::max(x0,x0_);
  newY0 = util::max(y0,y0_);
  newZ0 = util::max(z0,z0_);

  newX1 = util::min(x1,x1_);
  newY1 = util::min(y1,y1_);
  newZ1 = util::min(z1,z1_);

  return newX1>=newX0 && newY1>=newY0 && newZ1>=newZ0;
}

inline bool contained(int x, int y,
                      int x0, int x1, int y0, int y1)
{
  return x>=x0 && x<=x1 &&
         y>=y0 && y<=y1;
}

inline bool contained(int x, int y, int z,
                      int x0, int x1, int y0, int y1, int z0, int z1)
{
  return x>=x0 && x<=x1 &&
         y>=y0 && y<=y1 &&
         z>=z0 && z<=z1;
}


template<typename T>
T sqr(T arg)
{
  return arg*arg;
}

/// Compute norm square of a d-dimensional vector
template<typename ARRAY_LIKE, unsigned D>
auto normSqr(const ARRAY_LIKE& u) any_platform
{
  auto uSqr = decltype(u[0]){};
  for (unsigned iD=0; iD < D; ++iD) {
    uSqr += u[iD]*u[iD];
  }
  return uSqr;
}

template<unsigned D, typename ARRAY_LIKE>
auto norm(const ARRAY_LIKE& u)
{
  return sqrt(normSqr<ARRAY_LIKE,D>(u));
}

template<typename T, unsigned D>
T normSqr(const T* u) any_platform
{
  auto uSqr = T{};
  for (unsigned iD=0; iD < D; ++iD) {
    uSqr += u[iD]*u[iD];
  }
  return uSqr;
}

template<typename T>
T normSqr(std::initializer_list<T> data)
{
  return std::inner_product(data.begin(), data.end(), data.begin(), T(0));
}


/// Compute norm square of a d-dimensional vector
template<typename T, unsigned D, typename IMPL>
T normSqr(const ScalarVector<T,D,IMPL>& u) any_platform
{
  T uSqr = T();
  for (unsigned iD=0; iD < D; ++iD) {
    uSqr += u[iD]*u[iD];
  }
  return uSqr;
}

template<typename T, int d>
T scalarProduct(const T u1[d], const T u2[d])
{
  T prod = T();
  for (int iD=0; iD<d; ++iD) {
    prod += u1[iD]*u2[iD];
  }
  return prod;
}

template<typename T>
T scalarProduct(const std::vector<T>& u1, const std::vector<T>& u2)
{
  T prod = T();
  if (u1.size() == u2.size()) {
    for (int iD=0; iD<u1.size(); ++iD) {
      prod += u1[iD]*u2[iD];
    }
  }
  return prod;
}

/// Compute number of elements of a symmetric d-dimensional tensor
template <typename DESCRIPTORBASE>
struct TensorVal {
  static constexpr int n = (DESCRIPTORBASE::d*(DESCRIPTORBASE::d+1))/2; ///< result stored in n
};

/// Compute the opposite of a given direction
template <typename DESCRIPTORBASE> inline int opposite(int iPop)
{
  return descriptors::opposite<DESCRIPTORBASE>(iPop);
}

/// Return array of population indices where c[iVel] == value
template <typename DESCRIPTOR, unsigned iVel, int value>
constexpr auto populationsContributingToVelocity() any_platform {
  constexpr auto velocity_matches_value = [](unsigned iPop) -> bool {
    return descriptors::c<DESCRIPTOR>(iPop,iVel) == value;
  };
  return meta::array_from_index_sequence(
    descriptors::filter_population_indices<DESCRIPTOR>(velocity_matches_value));
};

/// Return array of population indices where c[iPop][iD] == NORMAL[iD]
template <typename DESCRIPTOR, int... NORMAL>
constexpr auto populationsContributingToDirection() any_platform {
  constexpr auto velocity_matches_direction = [](unsigned iPop) -> bool {
    return meta::indexed_pack_contains<NORMAL...>([iPop](unsigned iD, int x) -> bool {
      return x != 0 && x == descriptors::c<DESCRIPTOR>(iPop,iD);
    });
  };
  return meta::array_from_index_sequence(
    descriptors::filter_population_indices<DESCRIPTOR>(velocity_matches_direction));
};

template <typename DESCRIPTORBASE>
int findVelocity(const int v[DESCRIPTORBASE::d]) any_platform
{
  for (int iPop=0; iPop<DESCRIPTORBASE::q; ++iPop) {
    bool fit = true;
    for (int iD=0; iD<DESCRIPTORBASE::d; ++iD) {
      if (descriptors::c<DESCRIPTORBASE>(iPop,iD) != v[iD]) {
        fit = false;
        break;
      }
    }
    if (fit) {
      return iPop;
    }
  }
  return DESCRIPTORBASE::q;
}

/// Compute opposites of wall-incoming population indices
template <typename DESCRIPTOR, int direction, int orientation>
constexpr auto subIndexOutgoing() {
  // Should be expressed in terms of populationsContributingToVelocity
  // but std::array / std::index_sequence distinction makes this ugly
  // To be revisitited
  constexpr auto velocity_matches_value = [](unsigned iPop) {
    return descriptors::c<DESCRIPTOR>(iPop,direction) == orientation;
  };
  constexpr auto opposite_population = [](unsigned iPop) {
    return descriptors::opposite<DESCRIPTOR>(iPop);
  };
  return meta::array_from_index_sequence(
    meta::map_index_sequence(opposite_population,
                                        descriptors::filter_population_indices<DESCRIPTOR>(velocity_matches_value)));
}

///finds all the remaining indexes of a lattice given some other indexes
template <typename DESCRIPTORBASE>
std::vector<int> remainingIndexes(const std::vector<int> &indices)
{
  std::vector<int> remaining;
  for (int iPop = 0; iPop < DESCRIPTORBASE::q; ++iPop) {
    bool found = false;
    for (unsigned jPop = 0; jPop < indices.size(); ++jPop) {
      if (indices[jPop] == iPop) {
        found = true;
      }
    }
    if (!found) {
      remaining.push_back(iPop);
    }
  }
  return remaining;
}

template <typename DESCRIPTORBASE, int plane, int normal1, int normal2>
class SubIndexOutgoing3DonEdges {
private:
  SubIndexOutgoing3DonEdges()
  {
    int normalX,normalY,normalZ;
    typedef DESCRIPTORBASE L;

    switch (plane) {
    case 0: {
      normalX=0;
      if (normal1==1) {
        normalY= 1;
      }
      else {
        normalY=-1;
      }
      if (normal2==1) {
        normalZ= 1;
      }
      else {
        normalZ=-1;
      }
    }
    case 1: {
      normalY=0;
      if (normal1==1) {
        normalX= 1;
      }
      else {
        normalX=-1;
      }
      if (normal2==1) {
        normalZ= 1;
      }
      else {
        normalZ=-1;
      }
    }
    case 2: {
      normalZ=0;
      if (normal1==1) {
        normalX= 1;
      }
      else {
        normalX=-1;
      }
      if (normal2==1) {
        normalY= 1;
      }
      else {
        normalY=-1;
      }
    }
    }

    // add zero velocity
    //knownIndexes.push_back(0);
    // compute scalar product with boundary normal for all other velocities
    for (int iP=1; iP<L::q; ++iP) {
      if (descriptors::c<L>(iP,0)*normalX + descriptors::c<L>(iP,1)*normalY + descriptors::c<L>(iP,2)*normalZ<0) {
        indices.push_back(iP);
      }
    }
  }
  std::vector<int> indices;

  template <typename DESCRIPTORBASE_,  int plane_, int normal1_, int normal2_>
  friend std::vector<int> const& subIndexOutgoing3DonEdges();
};

template <typename DESCRIPTORBASE,  int plane, int normal1, int normal2>
std::vector<int> const& subIndexOutgoing3DonEdges()
{
  static SubIndexOutgoing3DonEdges<DESCRIPTORBASE,  plane, normal1, normal2> subIndexOutgoing3DonEdgesSingleton;
  return subIndexOutgoing3DonEdgesSingleton.indices;
}

// For 2D Corners
template <typename DESCRIPTORBASE, int normalX, int normalY>
class SubIndexOutgoing2DonCorners {
private:
  SubIndexOutgoing2DonCorners()
  {
    typedef DESCRIPTORBASE L;

    // add zero velocity
    //knownIndexes.push_back(0);
    // compute scalar product with boundary normal for all other velocities
    for (int iPop=1; iPop<L::q; ++iPop) {
      if (descriptors::c<L>(iPop,0)*normalX + descriptors::c<L>(iPop,1)*normalY<0) {
        indices.push_back(iPop);
      }
    }
  }

  std::vector<int> indices;

  template <typename DESCRIPTORBASE_,  int normalX_, int normalY_>
  friend std::vector<int> const& subIndexOutgoing2DonCorners();
};

template <typename DESCRIPTORBASE,  int normalX, int normalY>
std::vector<int> const& subIndexOutgoing2DonCorners()
{
  static SubIndexOutgoing2DonCorners<DESCRIPTORBASE, normalX, normalY> subIndexOutgoing2DonCornersSingleton;
  return subIndexOutgoing2DonCornersSingleton.indices;
}

// For 3D Corners
template <typename DESCRIPTORBASE, int normalX, int normalY, int normalZ>
class SubIndexOutgoing3DonCorners {
private:
  SubIndexOutgoing3DonCorners()
  {
    typedef DESCRIPTORBASE L;

    // add zero velocity
    //knownIndexes.push_back(0);
    // compute scalar product with boundary normal for all other velocities
    for (int iP=1; iP<L::q; ++iP) {
      if (descriptors::c<L>(iP,0)*normalX + descriptors::c<L>(iP,1)*normalY + descriptors::c<L>(iP,2)*normalZ<0) {
        indices.push_back(iP);
      }
    }
  }

  std::vector<int> indices;

  template <typename DESCRIPTORBASE_,  int normalX_, int normalY_, int normalZ_>
  friend std::vector<int> const& subIndexOutgoing3DonCorners();
};

template <typename DESCRIPTORBASE,  int normalX, int normalY, int normalZ>
std::vector<int> const& subIndexOutgoing3DonCorners()
{
  static SubIndexOutgoing3DonCorners<DESCRIPTORBASE, normalX, normalY, normalZ> subIndexOutgoing3DonCornersSingleton;
  return subIndexOutgoing3DonCornersSingleton.indices;
}

/// Util Function for Wall Model of Malaspinas
/// get link with smallest angle to a vector
template <typename T, typename DESCRIPTOR>
int get_nearest_link(const std::vector<T>& vec)
{
  T max=-1;
  int max_index = 0;
  for (int iQ=1; iQ<DESCRIPTOR::q; ++iQ) {
    std::vector<T> c_i(DESCRIPTOR::c[iQ], DESCRIPTOR::c[iQ]+3);
    T tmp = util::scalarProduct<T>(c_i, vec)/util::norm(c_i);
    if (tmp > max) {
      max = tmp;
      max_index = iQ;
    }
  }
  return max_index;
}

namespace tensorIndices2D {
enum { xx=0, xy=1, yy=2 };
}

namespace tensorIndices3D {
enum { xx=0, xy=1, xz=2, yy=3, yz=4, zz=5 };
}

/// compute lattice density from lattice pressure
template <typename T, typename DESCRIPTOR>
T densityFromPressure( T latticePressure )
{
  // rho = p / c_s^2 + 1
  return latticePressure * descriptors::invCs2<T,DESCRIPTOR>() + 1.0;
}

/// compute lattice pressure from lattice density
template <typename T, typename DESCRIPTOR>
T pressureFromDensity( T latticeDensity )
{
  // p = (rho - 1) * c_s^2
  return (latticeDensity - 1.0) / descriptors::invCs2<T,DESCRIPTOR>();
}

}  // namespace util

template <typename T>
std::enable_if_t<std::is_arithmetic<T>::type::value, T> abs(T x)
{
  return util::fabs(x);
}

// *INDENT-OFF*

}  // namespace olb

#endif
