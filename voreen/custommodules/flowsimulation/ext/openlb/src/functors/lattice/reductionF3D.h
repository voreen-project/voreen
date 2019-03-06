/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012-2017 Lukas Baron, Tim Dornieden, Mathias J. Krause,
 *  Albert Mink, Benjamin Förster, Adrian Kummerländer
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

#ifndef REDUCTION_F_3D_H
#define REDUCTION_F_3D_H

#include "functors/analytical/analyticalF.h"
#include "blockBaseF3D.h"
#include "superBaseF3D.h"
#include "geometry/cuboidGeometry3D.h"
#include "geometry/blockGeometry3D.h"
#include "geometry/superGeometry3D.h"

namespace olb {


/// Functor used to convert analytical functions to lattice functions
/**
 *  Input functions are interpreted as SI->SI units, the resulting lattice
 *  function will map lattice->lattice units
 *
 *  Maintains block level BlockLatticeFfromAnalyticalF3D functors.
 */
template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticeFfromAnalyticalF3D final : public SuperLatticeF3D<T,DESCRIPTOR> {
protected:
  AnalyticalF3D<T,T>& _f;
public:
  /**
   * \param f        Analytical functor to be converted into a lattice functor
   * \param sLattice Lattice reference required for conversion and block functor construction
   **/
  SuperLatticeFfromAnalyticalF3D(AnalyticalF3D<T,T>&           f,
                                 SuperLattice3D<T,DESCRIPTOR>& sLattice);
  bool operator() (T output[], const int input[]) override;
};


/// Block level functor for conversion of analytical to lattice functors.
/**
 * Instances are contained in SuperLatticeFfromAnalyticalF3D::_blockF.
 **/
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticeFfromAnalyticalF3D final : public BlockLatticeF3D<T,DESCRIPTOR> {
protected:
  AnalyticalF3D<T,T>& _f;
  Cuboid3D<T>         _cuboid;
  const int           _overlap;
public:
  /**
   * \param f       Analytical functor to be converted into a lattice functor
   * \param lattice Block lattice structure required for BlockLatticeF3D construction
   * \param cuboid  Cuboid reference required for input parameter conversion
   * \param overlap Block lattice overlap for input conversion
   **/
  BlockLatticeFfromAnalyticalF3D(AnalyticalF3D<T,T>&                    f,
                                 BlockLatticeStructure3D<T,DESCRIPTOR>& lattice,
                                 Cuboid3D<T>&                           cuboid,
                                 int                                    overlap);
  bool operator() (T output[], const int input[]) override;
};

//////////// not yet working // symbolically ///////////////////
////////////////////////////////////////////////
template <typename T, template <typename U> class DESCRIPTOR>
class SmoothBlockIndicator3D final : public BlockDataF3D<T,T> {
protected:
  IndicatorF3D<T>&  _f;
  T _h;
public:
  SmoothBlockIndicator3D(IndicatorF3D<T>& f, T h);
  //bool operator() (T output[], const int input[]);
};

// TODO: comment code
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticeInterpPhysVelocity3Degree3D final : public
  BlockLatticeF3D<T,DESCRIPTOR> {
protected:
  UnitConverter<T,DESCRIPTOR>& _conv;
  Cuboid3D<T>* _cuboid;
  int _overlap;
  int _range;
public:
  BlockLatticeInterpPhysVelocity3Degree3D(
    BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
    UnitConverter<T,DESCRIPTOR>& conv, Cuboid3D<T>* c, int overlap, int range);
  BlockLatticeInterpPhysVelocity3Degree3D(
    const BlockLatticeInterpPhysVelocity3Degree3D<T,DESCRIPTOR>& rhs);
  bool operator() (T output[], const int input[]) override
  {
    return false;
  }
  void operator() (T output[], const T input[]);
};

// TODO: comment code
template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticeInterpPhysVelocity3Degree3D final : public
  SuperLatticeF3D<T,DESCRIPTOR> {
private:
  std::vector<BlockLatticeInterpPhysVelocity3Degree3D<T,DESCRIPTOR>* > _bLattices;
public:
  SuperLatticeInterpPhysVelocity3Degree3D(
    SuperLattice3D<T,DESCRIPTOR>& sLattice, UnitConverter<T,DESCRIPTOR>& conv,
    int range=1);
  bool operator() (T output[], const int input[]) override
  {
    return 0;
  }
  void operator()(T output[], const T input[], const int iC);
};

// TODO: comment code
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticeInterpDensity3Degree3D final : public
  BlockLatticeF3D<T,DESCRIPTOR> {
protected:
  BlockGeometryStructure3D<T>& _blockGeometry;
  UnitConverter<T,DESCRIPTOR>& _conv;
  Cuboid3D<T>* _cuboid;
  int _overlap;
  int _range; // degree of interpolation can be changed (2,3,4,...)
public:
  BlockLatticeInterpDensity3Degree3D(
    BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
    BlockGeometryStructure3D<T>& blockGeometry,
    UnitConverter<T,DESCRIPTOR>& conv, Cuboid3D<T>* c, int overlap, int range);
  BlockLatticeInterpDensity3Degree3D(
    const BlockLatticeInterpDensity3Degree3D<T,DESCRIPTOR>& rhs);
  bool operator() (T output[], const int input[]) override
  {
    return false;
  }
  void operator() (T output[DESCRIPTOR<T>::q], const T input[3]);
};

// TODO: comment code
template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticeInterpDensity3Degree3D final : public
  SuperLatticeF3D<T,DESCRIPTOR> {
private:
  std::vector<BlockLatticeInterpDensity3Degree3D<T,DESCRIPTOR>* > _bLattices;
public:
  SuperLatticeInterpDensity3Degree3D(SuperLattice3D<T,DESCRIPTOR>& sLattice,
                                     SuperGeometry3D<T>& sGeometry,
                                     UnitConverter<T,DESCRIPTOR>& conv, int range=1);
  ~SuperLatticeInterpDensity3Degree3D() override;
  // range equals degree of interpolation and can be changed (2,3,4,...)
  bool operator() (T output[], const int input[]) override
  {
    return 0;
  }
  void operator()(T output[], const T input[], const int iC);
};

// TODO: comment code
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticeSmoothDiracDelta3D final : public
  BlockLatticeF3D<T,DESCRIPTOR> {
protected:
  UnitConverter<T,DESCRIPTOR>& _conv;
  Cuboid3D<T>* _cuboid;
public:
  BlockLatticeSmoothDiracDelta3D(BlockLattice3D<T,DESCRIPTOR>& blockLattice,
                                 UnitConverter<T,DESCRIPTOR>& conv, Cuboid3D<T>* c);
  BlockLatticeSmoothDiracDelta3D(
    const BlockLatticeSmoothDiracDelta3D<T,DESCRIPTOR>& rhs);
  bool operator() (T output[], const int input[]) override
  {
    return false;
  }
  void operator() (T delta[4][4][4], const T physPosP[3]);
};

// TODO: comment code
template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticeSmoothDiracDelta3D final : public
  SuperLatticeF3D<T,DESCRIPTOR> {
private:
  std::vector<BlockLatticeSmoothDiracDelta3D<T,DESCRIPTOR>* > _bLattices;
public:
  SuperLatticeSmoothDiracDelta3D(SuperLattice3D<T, DESCRIPTOR>& sLattice,
                                 UnitConverter<T,DESCRIPTOR>& conv,
                                 SuperGeometry3D<T>& superGeometry);
  ~SuperLatticeSmoothDiracDelta3D() override;
  bool operator()(T output[], const int input[]) override
  {
    return false;
  };
  void operator()(T delta[4][4][4], const T physPos[3], const int iC);
};


} // end namespace olb

#endif
