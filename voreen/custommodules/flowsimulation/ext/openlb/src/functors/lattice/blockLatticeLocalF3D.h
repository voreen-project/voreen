/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Lukas Baron, Mathias J. Krause, Albert Mink
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

#ifndef BLOCK_LATTICE_LOCAL_F_3D_H
#define BLOCK_LATTICE_LOCAL_F_3D_H

#include "blockBaseF3D.h"
#include "geometry/blockGeometry3D.h"
#include "indicator/indicatorBaseF3D.h"
#include "dynamics/smagorinskyBGKdynamics.h"

/** Note: Throughout the whole source code directory genericFunctions, the
 *  template parameters for i/o dimensions are:
 *           F: S^m -> T^n  (S=source, T=target)
 */

namespace olb {

template<typename T, template<typename U> class Lattice> class blockLatticeStructure3D;
template<typename T, typename S> class ParticleIndicatorF3D;

////////////////////////////////////////////////////////////////////////////////
//// if globIC is not on the local processor, the returned vector is empty//////
////////////////////////////////////////////////////////////////////////////////

/// functor returns pointwise f population on local lattices
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticeFpop3D final : public BlockLatticeF3D<T,DESCRIPTOR> {
public:
  BlockLatticeFpop3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice);
  bool operator() (T output[], const int input[]) override;
};

/// functor returns pointwise dissipation density on local lattices
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticeDissipation3D final : public BlockLatticeF3D<T,DESCRIPTOR> {
protected:
  const UnitConverter<T,DESCRIPTOR>& _converter;
public:
  BlockLatticeDissipation3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
                            const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};

/// functor returns pointwise dissipation density on local lattices
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticePhysDissipation3D final : public BlockLatticeF3D<T,DESCRIPTOR> {
protected:
  const int             _overlap;
  const UnitConverter<T, DESCRIPTOR>& _converter;
public:
  BlockLatticePhysDissipation3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
                                int overlap,
                                const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};

/// functor returns pointwise turbulent dissipation density on local lattices
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticeEffevtiveDissipation3D final : public BlockLatticeF3D<T,DESCRIPTOR> {
protected:
  const UnitConverter<T,DESCRIPTOR>& _converter;
  T _smagoConst;
  LESDynamics<T, DESCRIPTOR>& _LESdynamics;
public:
  BlockLatticeEffevtiveDissipation3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
                                     const UnitConverter<T,DESCRIPTOR>& converter, T smagoConst,
                                     LESDynamics<T, DESCRIPTOR>& LESdynamics);
  bool operator() (T output[], const int input[]) override;
};

/// functor returns pointwise turbulent dissipation density on local lattices
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticePhysEffevtiveDissipation3D final : public BlockLatticeF3D<T,DESCRIPTOR> {
protected:
  const UnitConverter<T,DESCRIPTOR>& _converter;
  T _smagoConst;
  LESDynamics<T, DESCRIPTOR>& _LESdynamics;
public:
  BlockLatticePhysEffevtiveDissipation3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
                                         const UnitConverter<T,DESCRIPTOR>& converter, T smagoConst,
                                         LESDynamics<T, DESCRIPTOR>& LESdynamics);
  bool operator() (T output[], const int input[]) override;
};


/// functor returns pointwise density rho on local lattices
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticeDensity3D final : public BlockLatticeF3D<T,DESCRIPTOR> {
public:
  BlockLatticeDensity3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice);
  bool operator() (T output[], const int input[]) override;
};


/// functor returns pointwise velocity on local lattice
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticeVelocity3D final : public BlockLatticeF3D<T,DESCRIPTOR> {
public:
  BlockLatticeVelocity3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice);
  bool operator() (T output[], const int input[]) override;
};

/// functor returns pointwise external velocity (external field) on local lattice
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticeExternalVelocity3D final : public BlockLatticeF3D<T,DESCRIPTOR> {
public:
  BlockLatticeExternalVelocity3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice);
  bool operator() (T output[], const int input[]);
};

/// functor returns pointwise the material no. presenting the geometry on local lattice
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticeGeometry3D final : public BlockLatticeF3D<T,DESCRIPTOR> {
  BlockGeometryStructure3D<T>& _blockGeometry;
  int _material;
public:
  BlockLatticeGeometry3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
                         BlockGeometryStructure3D<T>& blockGeometry, int material = -1);
  bool operator() (T output[], const int input[]) override;
};

/// functor to get pointwise the rank no. + 1 on local lattice
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticeRank3D final : public BlockLatticeF3D<T,DESCRIPTOR> {
public:
  BlockLatticeRank3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice);
  bool operator() (T output[], const int input[]) override;
};

/// functor to get pointwise the cuboid no. + 1 on local lattice
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticeCuboid3D final : public BlockLatticeF3D<T,DESCRIPTOR> {
private:
  // holds cuboid nmb of current block
  int _iC;
public:
  BlockLatticeCuboid3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, int iC);
  bool operator() (T output[], const int input[]) override;
};


/// functor returns pointwise phys pressure from rho on local lattices
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticePhysPressure3D final : public BlockLatticePhysF3D<T,DESCRIPTOR> {
private:
  const int _overlap;
public:
  BlockLatticePhysPressure3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
                             int overlap,
                             const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};


/// functor returns pointwise phys velocity on local lattice
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticePhysVelocity3D final : public BlockLatticePhysF3D<T,DESCRIPTOR> {
private:
  const int  _overlap;
  const bool _print;
public:
  BlockLatticePhysVelocity3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
                             int overlap,
                             const UnitConverter<T,DESCRIPTOR>& converter,
                             bool print=false);
  bool operator() (T output[], const int input[]) override;
};

template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticePhysExternalVelocity3D final : public BlockLatticePhysF3D<T,DESCRIPTOR> {
public:
  BlockLatticePhysExternalVelocity3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
                                     const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]);
};

template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticePhysExternalPorosity3D final : public BlockLatticePhysF3D<T,DESCRIPTOR> {
private:
  const int _overlap;
public:
  BlockLatticePhysExternalPorosity3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
                                     int overlap,
                                     const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]);
};

template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticePhysExternalParticleVelocity3D final : public BlockLatticePhysF3D<T,DESCRIPTOR> {
public:
  BlockLatticePhysExternalParticleVelocity3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
      const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]);
};

template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticeExternal3D final : public BlockLatticeF3D<T,DESCRIPTOR> {
public:
  BlockLatticeExternal3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, int begin, int size);
  bool operator() (T output[], const int input[]);
private:
  int _start, _size;
};

template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticePhysExternal3D final : public BlockLatticeF3D<T,DESCRIPTOR> {
public:
  BlockLatticePhysExternal3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
                             T convFactorToPhysUnits, int offset, int size);
  bool operator() (T output[], const int input[]) override;
private:
  T _convFactorToPhysUnits;
  int _offset, _size;
};

/// functor returns pointwise strain rate on local lattice, s_ij = 1/2*(du_idr_j + du_jdr_i)
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticeStrainRate3D final : public BlockLatticePhysF3D<T,DESCRIPTOR> {
public:
  BlockLatticeStrainRate3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
                           const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};

/// functor returns pointwise phys strain rate on local lattice, s_ij = 1/2*(du_idr_j + du_jdr_i)
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticePhysStrainRate3D final : public BlockLatticePhysF3D<T,DESCRIPTOR> {
private:
  const int _overlap;
public:
  BlockLatticePhysStrainRate3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
                               int overlap,
                               const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]);
};

/// functor returns pointwise phys force acting on a boundary with a given material on local lattice
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticePhysBoundaryForce3D final : public BlockLatticePhysF3D<T,DESCRIPTOR> {
private:
  BlockGeometryStructure3D<T>& _blockGeometry;
  int _material;
public:
  BlockLatticePhysBoundaryForce3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
                                  BlockGeometryStructure3D<T>& blockGeometry, int material,
                                  const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};

/// functor returns pointwise phys force acting on a boundary with a given material on local lattice
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticePhysWallShearStress3D final : public BlockLatticePhysF3D<T,DESCRIPTOR> {
private:
  BlockGeometryStructure3D<T>& _blockGeometry;
  BlockLatticeStructure3D<T,DESCRIPTOR>& _blockLattice;
  int _material;
  const UnitConverter<T,DESCRIPTOR>& _converter;
  std::vector<std::vector<std::vector<std::vector<int>>>> _discreteNormal;
  std::vector<std::vector<std::vector<std::vector<T>>>> _normal;
  T _physFactor;
public:
  BlockLatticePhysWallShearStress3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
                                    BlockGeometryStructure3D<T>& blockGeometry, int material,
                                    const UnitConverter<T,DESCRIPTOR>& converter,
                                    IndicatorF3D<T>& indicator);
  bool operator() (T output[], const int input[]) override;
};


/**
 *  functor returns pointwise phys force acting on a boundary with a given material on local lattice
 *  see: Caiazzo, Junk: Boundary Forces in lattice Boltzmann: Analysis of MEA
 */
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticePhysCorrBoundaryForce3D final : public BlockLatticePhysF3D<T,DESCRIPTOR> {
private:
  BlockGeometryStructure3D<T>& _blockGeometry;
  int _material;
public:
  BlockLatticePhysCorrBoundaryForce3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
                                      BlockGeometry3D<T>& blockGeometry, int material,
                                      const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};


/// functor to get pointwise, lattice-dependent external field
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticeExternalField3D final : public BlockLatticeF3D<T,DESCRIPTOR> {
private:
  int _beginsAt;
  int _sizeOf;
public:
  BlockLatticeExternalField3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, int beginsAt, int sizeOf);
  bool operator() (T output[], const int input[]) override;
};

/**
 *  functor returns pointwise, lattice-dependent porosity values in [0,1]
 *  in combination with (Extended)PorousBGKdynamics: 0->solid, 1->fluid
 */
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticePorosity3D final : public BlockLatticeF3D<T,DESCRIPTOR> {
public:
  BlockLatticePorosity3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice);
  bool operator() (T output[], const int input[]);
};


/**
 *  functor to get pointwise mesh-independent permeability values in (0,inf)
 *  in combination with (Extended)PorousBGKdynamics
 *  note: result is cropped to 999999
 */
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticePhysPermeability3D final : public BlockLatticePhysF3D<T,DESCRIPTOR> {
public:
  BlockLatticePhysPermeability3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]);
};


/**
 *  functor to get pointwise mesh-independent permeability values in (0,inf)
 *  in combination with (Extended)PorousBGKdynamics
 *  note: result is cropped to 1
 */
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticePhysCroppedPermeability3D final : public BlockLatticePhysF3D<T,DESCRIPTOR> {
public:
  BlockLatticePhysCroppedPermeability3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]);
};


/*
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticePhysPermeability3D final : public BlockLatticePhysF3D<T,DESCRIPTOR> {
private:
  BlockGeometry3D<T>& _blockGeometry;
  int _material;
public:
  BlockLatticePhysPermeability3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
                                 BlockGeometry3D<T>& blockGeometry,
                                 int material, const UnitConverter<T>& converter);
  bool operator() (T output[], const int input[]);
};*/


/// functor returns pointwise -nu/K*u on the lattice, can be used with BlockSum3D as objective
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticePhysDarcyForce3D final : public BlockLatticePhysF3D<T,DESCRIPTOR> {
private:
  BlockGeometry3D<T>& _blockGeometry;
  int _material;
public:
  BlockLatticePhysDarcyForce3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
                               BlockGeometry3D<T>& blockGeometry,
                               int material, const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]);
};


/**
 *  functor to get a pointwise local average of a passed functor with a given
 *  material and radius on local lattice
 *  the output data must be of the same size and dimension like f
 */
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticeAverage3D final : public BlockLatticeF3D<T,DESCRIPTOR> {
private:
  BlockLatticeF3D<T,DESCRIPTOR>& _f;
  BlockGeometry3D<T>& _blockGeometry;
  int _material;
  T _radius;
public:
  BlockLatticeAverage3D(BlockLatticeF3D<T,DESCRIPTOR>& f,
                        BlockGeometry3D<T>& blockGeometry, int material,
                        T radius);
  bool operator() (T output[], const int input[]) override;
};


/// functor returns pointwise the l2-norm, e.g. of a velocity
template <typename T, template <typename U> class DESCRIPTOR>
class BlockEuklidNorm3D final : public BlockF3D<T> {
protected:
  BlockF3D<T>& _f;
public:
  BlockEuklidNorm3D(BlockF3D<T>& f);
  bool operator() (T output[], const int input[]) override;
};


template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticeInterpPhysVelocity3D final : public BlockLatticePhysF3D<T,DESCRIPTOR> {
protected:
  Cuboid3D<T>* _cuboid;
  int _overlap;
public:
  BlockLatticeInterpPhysVelocity3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
                                   const UnitConverter<T,DESCRIPTOR>& conv, Cuboid3D<T>* c, int overlap);
  BlockLatticeInterpPhysVelocity3D(const BlockLatticeInterpPhysVelocity3D<T,DESCRIPTOR>& rhs);
  bool operator() (T output[3], const int input[3]) override
  {
    return false;
  }
  void operator() (T output[3], const T input[3]);
};

//template <typename T, template <typename U> class DESCRIPTOR>
//class BlockLatticeInterpPhysVelocity3Degree3D final : public BlockLatticeF3D<T,DESCRIPTOR> {
//protected:
//  UnitConverter<T>& _conv;
//  Cuboid3D<T>* _cuboid;
//  int _overlap;
//public:
//  BlockLatticeInterpPhysVelocity3Degree3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, UnitConverter<T>& conv, Cuboid3D<T>* c, int overlap);
//  BlockLatticeInterpPhysVelocity3Degree3D(const BlockLatticeInterpPhysVelocity3Degree3D<T,DESCRIPTOR>& rhs);
//  bool operator() (T output[3], const int input[3]) {
//    return false;
//  }
//  void operator() (T output[3], const T input[3]);
//};
//
//template <typename T, template <typename U> class DESCRIPTOR>
//class BlockLatticeInterpDensity3Degree3D final : public BlockLatticeF3D<T,DESCRIPTOR> {
//protected:
//  UnitConverter<T>& _conv;
//  Cuboid3D<T>* _cuboid;
//  int _overlap;
//public:
//  BlockLatticeInterpDensity3Degree3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice, UnitConverter<T>& conv, Cuboid3D<T>* c, int overlap);
//  BlockLatticeInterpDensity3Degree3D(const BlockLatticeInterpDensity3Degree3D<T,DESCRIPTOR>& rhs);
//  bool operator() (T output[3], const int input[3]) {
//    return false;
//  }
//  void operator() (T output[3], const T input[3]);
//};

/// functor that returns forces acting on a particle surface, returns for every particle 0...2: translational force, 3...5: torque, 7: number of voxels
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticeMomentumExchangeForce3D final : public BlockLatticePhysF3D<T,DESCRIPTOR> {
private:
  BlockGeometryStructure3D<T>& _blockGeometry;
  std::vector<ParticleIndicatorF3D<T,T>* >& _vectorOfIndicator;
  int _overlap;
public:
  BlockLatticeMomentumExchangeForce3D(BlockLatticeStructure3D<T,DESCRIPTOR>& blockLattice,
                                      BlockGeometryStructure3D<T>& blockGeometry,
                                      std::vector<ParticleIndicatorF3D<T,T>* >& indicator,
                                      const UnitConverter<T,DESCRIPTOR>& converter, int overlap);
  bool operator() (T output[], const int input[]) override;
};

template <typename T, template <typename U> class DESCRIPTOR, template <typename V> class ThermalDESCRIPTOR>
class BlockLatticePhysTemperature3D final : public BlockLatticeThermalPhysF3D<T,DESCRIPTOR,ThermalDESCRIPTOR> {
public:
  BlockLatticePhysTemperature3D(BlockLatticeStructure3D<T,ThermalDESCRIPTOR>& blockLattice,
                                ThermalUnitConverter<T,DESCRIPTOR,ThermalDESCRIPTOR> const& converter);
  bool operator() (T output[], const int input[]);
};
} // end namespace olb

#endif
