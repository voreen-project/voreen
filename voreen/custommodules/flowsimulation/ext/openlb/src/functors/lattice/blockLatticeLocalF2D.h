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

#ifndef BLOCK_LATTICE_LOCAL_F_2D_H
#define BLOCK_LATTICE_LOCAL_F_2D_H

#include<vector>

#include "blockBaseF2D.h"
#include "geometry/blockGeometry2D.h"
#include "core/blockLattice2D.h"
#include "core/blockLatticeStructure2D.h"

/** Note: Throughout the whole source code directory genericFunctions, the
 *  template parameters for i/o dimensions are:
 *           F: S^m -> T^n  (S=source, T=target)
 */

namespace olb {

////////////////////////////////////////////////////////////////////////////////
//// if globIC is not on the local processor, the returned vector is empty//////
////////////////////////////////////////////////////////////////////////////////


/// BlockLatticeDissipation2D returns pointwise dissipation density on local lattices.
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticeDissipation2D final : public BlockLatticeF2D<T,DESCRIPTOR> {
protected:
  const UnitConverter<T,DESCRIPTOR>& _converter;
public:
  BlockLatticeDissipation2D(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice,
                            const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};

/// BlockLatticePhysDissipation2D returns pointwise physical dissipation density on local lattices.
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticePhysDissipation2D final : public BlockLatticePhysF2D<T,DESCRIPTOR> {
public:
  BlockLatticePhysDissipation2D(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice,
                                const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};

/// BlockLatticeDensity2D returns pointwise density rho on local lattices.
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticeDensity2D final : public BlockLatticeF2D<T,DESCRIPTOR> {
public:
  BlockLatticeDensity2D(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice);
  bool operator() (T output[], const int input[]) override;
};



/// BlockLatticeVelocity2D returns pointwise velocity on local lattices.
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticeVelocity2D final : public BlockLatticeF2D<T,DESCRIPTOR> {
public:
  BlockLatticeVelocity2D(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice);
  bool operator() (T output[], const int input[]) override;
};


/// BlockLatticeGeometry2D returns pointwise the material no. presenting the geometry on local lattice.
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticeGeometry2D final : public BlockLatticeF2D<T,DESCRIPTOR> {
private:
  BlockGeometryStructure2D<T>& _blockGeometry;
  int _material;
public:
  BlockLatticeGeometry2D(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice,
                         BlockGeometryStructure2D<T>& blockGeometry, int material = -1);
  bool operator() (T output[], const int input[]) override;
};


/// BlockLatticeRank2D returns pointwise the rank no. + 1 on local lattice.
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticeRank2D final : public BlockLatticeF2D<T,DESCRIPTOR> {
public:
  BlockLatticeRank2D(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice);
  bool operator() (T output[], const int input[]) override;
};


/// BlockLatticeCuboid2D returns pointwise the cuboid no. + 1 on local lattice.
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticeCuboid2D final : public BlockLatticeF2D<T,DESCRIPTOR> {
public:
  BlockLatticeCuboid2D(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice, const int iC);
  bool operator() (T output[], const int input[]) override;
private:
  // holds cuboid nmb of current block
  const int _iC;
};


/// BlockLatticePhysPressure2D returns pointwise phys pressure from rho on local lattices.
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticePhysPressure2D final : public BlockLatticePhysF2D<T,DESCRIPTOR> {
public:
  BlockLatticePhysPressure2D(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice,
                             const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};


/// BlockLatticePhysVelocity2D returns pointwise phys velocity on local lattice.
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticePhysVelocity2D final : public BlockLatticePhysF2D<T,DESCRIPTOR> {
public:
  BlockLatticePhysVelocity2D(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice,
                             const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};


template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticePhysExternalPorosity2D final : public BlockLatticePhysF2D<T,DESCRIPTOR> {
public:
  BlockLatticePhysExternalPorosity2D(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice,
                                     const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]);
};

template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticePhysExternalVelocity2D final : public BlockLatticePhysF2D<T,DESCRIPTOR> {
public:
  BlockLatticePhysExternalVelocity2D(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice,
                                     const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]);
};

template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticePhysExternalParticleVelocity2D final : public BlockLatticePhysF2D<T,DESCRIPTOR> {
public:
  BlockLatticePhysExternalParticleVelocity2D(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice,
      const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]);
};

/// BlockLatticeStrainRate2D returns pointwise strain rate on local lattice.
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticeStrainRate2D final : public BlockLatticePhysF2D<T,DESCRIPTOR> {
public:
  BlockLatticeStrainRate2D(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice,
                           const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]);
};


/// BlockLatticePhysStrainRate2D returns pointwise phys strain rate on local lattice.
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticePhysStrainRate2D final : public BlockLatticePhysF2D<T,DESCRIPTOR> {
public:
  BlockLatticePhysStrainRate2D(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice,
                               const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};

/// BlockLatticePhysBoundaryForce2D returns pointwise phys force acting on a boundary
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticePhysBoundaryForce2D final : public BlockLatticePhysF2D<T,DESCRIPTOR> {
private:
  BlockGeometryStructure2D<T>& _blockGeometry;
  int _material;
public:
  BlockLatticePhysBoundaryForce2D(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice,
                                  BlockGeometryStructure2D<T>& blockGeometry,
                                  int material, const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};

/// BlockLatticePhysBoundaryForce2D returns pointwise wall shear stress
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticePhysWallShearStress2D final : public BlockLatticePhysF2D<T,DESCRIPTOR> {
private:
  BlockGeometryStructure2D<T>& _blockGeometry;
  BlockLatticeStructure2D<T,DESCRIPTOR>& _blockLattice;
  int _material;
  const UnitConverter<T,DESCRIPTOR>& _converter;
  std::vector<std::vector<std::vector<int>>> _discreteNormal;
  T _physFactor;
public:
  BlockLatticePhysWallShearStress2D(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice,
                                    BlockGeometryStructure2D<T>& blockGeometry,
                                    int material, const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};


/**
 *  BlockLatticePhysCorrBoundaryForce2D returns pointwise phys force acting on a
 *  boundary with a given material on local lattice.
 *  see: Caiazzo, Junk: Boundary Forces in lattice Boltzmann: Analysis of MEA
 */
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticePhysCorrBoundaryForce2D final : public BlockLatticePhysF2D<T,DESCRIPTOR> {
private:
  BlockGeometry2D<T>& _blockGeometry;
  int _material;
public:
  BlockLatticePhysCorrBoundaryForce2D(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice,
                                      BlockGeometry2D<T>& blockGeometry, int material,
                                      const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};


/**
 *  BlockLatticePorosity2D returns pointwise, lattice-dependent porosity values in [0,1]
 *  in combination with (Extended)PorousBGKdynamics: 0->solid, 1->fluid.
 */
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticePorosity2D final : public BlockLatticePhysF2D<T,DESCRIPTOR> {
private:
  BlockGeometryStructure2D<T>& _blockGeometry;
  int _material;
public:
  BlockLatticePorosity2D(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice,
                         BlockGeometryStructure2D<T>& blockGeometry, int material,
                         const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};


/**
 *  BlockLatticePhysPermeability2D returns pointwise mesh-independent permeability
 *  values in (0,inf) in combination with (Extended)PorousBGKdynamics
 *  note: result is cropped to 999999.
 */
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticePhysPermeability2D final : public BlockLatticePhysF2D<T,DESCRIPTOR> {
private:
  BlockGeometryStructure2D<T>& _blockGeometry;
  int _material;
public:
  BlockLatticePhysPermeability2D(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice,
                                 BlockGeometryStructure2D<T>& blockGeometry,
                                 int material, const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};


/// BlockLatticePhysDarcyForce2D computes pointwise -nu/K*u on the lattice. can be used with BlockSum2D as objective
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticePhysDarcyForce2D final : public BlockLatticePhysF2D<T,DESCRIPTOR> {
private:
  BlockGeometry2D<T>& _blockGeometry;
  int _material;
public:
  BlockLatticePhysDarcyForce2D(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice,
                               BlockGeometry2D<T>& blockGeometry, int material,
                               const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};


/**
 *  BlockLatticeAverage2D returns pointwise local average of a passed functor with
 *  a given material and radius on local lattice.
 *  the output data must be of the same size and dimension like f.
 */
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticeAverage2D final : public BlockLatticeF2D<T,DESCRIPTOR> {
private:
  BlockLatticeF2D<T,DESCRIPTOR>& _f;
  BlockGeometry2D<T>& _blockGeometry;
  int _material;
  T _radius;
public:
  BlockLatticeAverage2D(BlockLatticeF2D<T,DESCRIPTOR>& f,
                        BlockGeometry2D<T>& blockGeometry, int material, T radius);
  bool operator() (T output[], const int input[]) override;
};


///  BlockL2Norm2D returns pointwise the l2-norm, e.g. of a velocity.
template <typename T, template <typename U> class DESCRIPTOR>
class BlockEuklidNorm2D final : public BlockF2D<T> {
private:
  BlockF2D<T>& _f;
public:
  BlockEuklidNorm2D(BlockF2D<T>& f);
  bool operator() (T output[], const int input[]) override;
};

/// functor that returns forces acting on a particle surface, returns for every particle 0,1: translational force, 3: torque, 4: number of voxels
template <typename T, template <typename U> class DESCRIPTOR>
class BlockLatticeMomentumExchangeForce2D final : public BlockLatticePhysF2D<T,DESCRIPTOR> {
private:
  BlockGeometryStructure2D<T>& _blockGeometry;
  std::vector<ParticleIndicatorF2D<T,T>* >& _vectorOfIndicator;
  int _overlap;
public:
  BlockLatticeMomentumExchangeForce2D(BlockLatticeStructure2D<T,DESCRIPTOR>& blockLattice,
                                      BlockGeometryStructure2D<T>& blockGeometry,
                                      std::vector<ParticleIndicatorF2D<T,T>* >& indicator,
                                      const UnitConverter<T,DESCRIPTOR>& converter, int overlap);
  bool operator() (T output[], const int input[]) override;
};

/// BlockLatticePhysTemperature2D returns pointwise phys temperature from rho on local lattices.
template <typename T, template <typename U> class DESCRIPTOR, template <typename V> class ThermalDESCRIPTOR>
class BlockLatticePhysTemperature2D final : public BlockLatticeThermalPhysF2D<T,DESCRIPTOR,ThermalDESCRIPTOR> {
public:
  BlockLatticePhysTemperature2D(BlockLatticeStructure2D<T,ThermalDESCRIPTOR>& blockLattice,
                                ThermalUnitConverter<T,DESCRIPTOR,ThermalDESCRIPTOR> const& converter);
  bool operator() (T output[], const int input[]);
};

} // end namespace olb

#endif
