/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Lukas Baron, Tim Dornieden, Mathias J. Krause,
 *  Albert Mink
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

#ifndef SUPER_LATTICE_LOCAL_F_3D_H
#define SUPER_LATTICE_LOCAL_F_3D_H

#include<vector>

#include "superBaseF3D.h"
#include "superCalcF3D.h"
#include "blockLatticeLocalF3D.h"
#include "indicator/indicatorBaseF3D.h"
#include "core/superLattice3D.h"

/** Note: Throughout the whole source code directory genericFunctions, the
 *  template parameters for i/o dimensions are:
 *           F: S^m -> T^n  (S=source, T=target)
 */

namespace olb {


////////////////////////////////////////////////////////////////////////////////
//////if globIC is not on the local processor, the returned vector is empty/////
////////////////////////////////////////////////////////////////////////////////

/// functor to get pointwise f population on local lattices
template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticeFpop3D final : public SuperLatticeF3D<T,DESCRIPTOR> {
public:
  SuperLatticeFpop3D(SuperLattice3D<T,DESCRIPTOR>& sLattice);
  bool operator() (T output[], const int input[]) override;
};

/// functor to get pointwise dissipation density on local lattices
template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticeDissipation3D final : public SuperLatticeF3D<T,DESCRIPTOR> {
private:
  const UnitConverter<T,DESCRIPTOR>& _converter;
public:
  SuperLatticeDissipation3D(SuperLattice3D<T,DESCRIPTOR>& sLattice,
                            const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};

/// functor to get pointwise dissipation density on local lattices
template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticePhysDissipation3D final : public SuperLatticePhysF3D<T,DESCRIPTOR> {
public:
  SuperLatticePhysDissipation3D(SuperLattice3D<T,DESCRIPTOR>& sLattice,
                                const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};

/// functor to get pointwise turbulent dissipation density on local lattices
template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticeEffevtiveDissipation3D final : public SuperLatticeF3D<T,DESCRIPTOR> {
private:
  const UnitConverter<T,DESCRIPTOR>& _converter;
public:
  SuperLatticeEffevtiveDissipation3D(SuperLattice3D<T,DESCRIPTOR>& sLattice,
                                     const UnitConverter<T,DESCRIPTOR>& converter, T smagoConst,
                                     LESDynamics<T, DESCRIPTOR>& LESdynamics);
  bool operator() (T output[], const int input[]) override;
};

/// functor to get pointwise turbulent dissipation density on local lattices
template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticePhysEffevtiveDissipation3D final : public SuperLatticePhysF3D<T,DESCRIPTOR> {
public:
  SuperLatticePhysEffevtiveDissipation3D(SuperLattice3D<T,DESCRIPTOR>& sLattice,
                                         const UnitConverter<T,DESCRIPTOR>& converter, T smagoConst,
                                         LESDynamics<T, DESCRIPTOR>& LESdynamics);
  bool operator() (T output[], const int input[]) override;
};

/// functor to get pointwise density rho on local lattices
template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticeDensity3D final : public SuperLatticeF3D<T,DESCRIPTOR> {
public:
  SuperLatticeDensity3D(SuperLattice3D<T,DESCRIPTOR>& sLattice);
  bool operator() (T output[], const int input[]) override;
};


/// functor to get pointwise velocity on local lattice
template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticeVelocity3D final : public SuperLatticeF3D<T,DESCRIPTOR> {
public:
  SuperLatticeVelocity3D(SuperLattice3D<T,DESCRIPTOR>& sLattice);
  bool operator() (T output[], const int input[]) override;
};

/// functor to get pointwise velocity on local lattice
template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticeExternalVelocity3D final : public SuperLatticeF3D<T,DESCRIPTOR> {
public:
  SuperLatticeExternalVelocity3D(SuperLattice3D<T,DESCRIPTOR>& sLattice);
  bool operator() (T output[], const int input[]);
};

/// functor to get pointwise strain rate on local lattice
/// s_ij = 1/2*(du_idr_j + du_jdr_i)
template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticeStrainRate3D final : public SuperLatticeF3D<T,DESCRIPTOR> {
private:
  const UnitConverter<T,DESCRIPTOR>& _converter;
public:
  SuperLatticeStrainRate3D(SuperLattice3D<T,DESCRIPTOR>& sLattice,
                           const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};

/// functor to get pointwise phys strain rate on local lattice
/// s_ij = 1/2*(du_idr_j + du_jdr_i)
template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticePhysStrainRate3D final : public SuperLatticePhysF3D<T,DESCRIPTOR> {
public:
  SuperLatticePhysStrainRate3D(SuperLattice3D<T,DESCRIPTOR>& sLattice,
                               const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]);
};

/// functor to get pointwise the material no. presenting the geometry on local lattice
template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticeGeometry3D final : public SuperLatticeF3D<T,DESCRIPTOR> {
private:
  SuperGeometry3D<T>& _superGeometry;
  const int _material;
public:
  SuperLatticeGeometry3D(SuperLattice3D<T,DESCRIPTOR>& sLattice,
                         SuperGeometry3D<T>& superGeometry, const int material = -1);
  bool operator() (T output[], const int input[]) override;
};


/// functor to get pointwise the rank no. + 1 on local lattice
template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticeRank3D final : public SuperLatticeF3D<T,DESCRIPTOR> {
public:
  SuperLatticeRank3D(SuperLattice3D<T,DESCRIPTOR>& sLattice);
  bool operator() (T output[], const int input[]) override;
};


/// functor to get pointwise the cuboid no. + 1 on local lattice
template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticeCuboid3D final : public SuperLatticeF3D<T,DESCRIPTOR> {
public:
  SuperLatticeCuboid3D(SuperLattice3D<T,DESCRIPTOR>& sLattice);
  bool operator() (T output[], const int input[]) override;
};


/// functor to get pointwise phys pressure from rho on local lattices
template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticePhysPressure3D final : public SuperLatticePhysF3D<T,DESCRIPTOR> {
public:
  SuperLatticePhysPressure3D(SuperLattice3D<T,DESCRIPTOR>& sLattice,
                             const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};


/// functor to get pointwise phys velocity on local lattice
template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticePhysVelocity3D final : public SuperLatticePhysF3D<T,DESCRIPTOR> {
public:
  SuperLatticePhysVelocity3D(SuperLattice3D<T,DESCRIPTOR>& sLattice,
                             const UnitConverter<T,DESCRIPTOR>& converter, bool print=false);
  bool operator() (T output[], const int input[]) override;
private:
  bool _print;
};

template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticePhysExternalPorosity3D final : public SuperLatticePhysF3D<T,DESCRIPTOR> {
public:
  SuperLatticePhysExternalPorosity3D(SuperLattice3D<T,DESCRIPTOR>& sLattice,
                                     const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]);
};

template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticePhysExternalVelocity3D final : public SuperLatticePhysF3D<T,DESCRIPTOR> {
public:
  SuperLatticePhysExternalVelocity3D(SuperLattice3D<T,DESCRIPTOR>& sLattice,
                                     const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]);
};

template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticePhysExternalParticleVelocity3D final : public SuperLatticePhysF3D<T,DESCRIPTOR> {
public:
  SuperLatticePhysExternalParticleVelocity3D(SuperLattice3D<T,DESCRIPTOR>& blockLattice,
      const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]);
};

template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticeExternal3D final : public SuperLatticeF3D<T,DESCRIPTOR> {
public:
  SuperLatticeExternal3D(SuperLattice3D<T,DESCRIPTOR>& sLattice, int start, int size);
  bool operator() (T output[], const int input[]);
};

template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticePhysExternal3D final : public SuperLatticeF3D<T,DESCRIPTOR> {
public:
  SuperLatticePhysExternal3D(SuperLattice3D<T,DESCRIPTOR>& sLattice,
                             T convFactorToPhysUnits,
                             int offset, int size);
  bool operator() (T output[], const int input[]) override;
};

/// functor to get pointwise phys force acting on a boundary with a given material on local lattice
template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticePhysBoundaryForce3D final : public SuperLatticePhysF3D<T,DESCRIPTOR> {
private:
  SuperGeometry3D<T>& _superGeometry;
  const int _material;
public:
  SuperLatticePhysBoundaryForce3D(SuperLattice3D<T,DESCRIPTOR>& sLattice,
                                  SuperGeometry3D<T>& superGeometry, const int material,
                                  const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};

/// functor to get pointwise phys wall shear stress with a given material on local lattice
template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticePhysWallShearStress3D final : public SuperLatticePhysF3D<T,DESCRIPTOR> {
private:
  SuperGeometry3D<T>& _superGeometry;
  const int _material;
public:
  SuperLatticePhysWallShearStress3D(SuperLattice3D<T,DESCRIPTOR>& sLattice,
                                    SuperGeometry3D<T>& superGeometry, const int material,
                                    const UnitConverter<T,DESCRIPTOR>& converter,
                                    IndicatorF3D<T>& indicator);
  bool operator() (T output[], const int input[]) override;
};

/// functor to get pointwise phys force acting on a boundary with a given indicator on local lattice
template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticePhysBoundaryForceIndicator3D final : public SuperLatticePhysF3D<T,DESCRIPTOR> {
private:
  SuperGeometry3D<T>& _superGeometry;
  ParticleIndicatorF3D<T,T>& _indicator;
//  const UnitConverter<T,DESCRIPTOR>& _converter;
public:
  SuperLatticePhysBoundaryForceIndicator3D(SuperLattice3D<T,DESCRIPTOR>& sLattice,
      SuperGeometry3D<T>& superGeometry,
      ParticleIndicatorF3D<T,T>& indicator,
      const UnitConverter<T,DESCRIPTOR>& converter);
//  SuperLatticePhysBoundaryForceIndicator3D(const SuperLatticePhysBoundaryForceIndicator3D<T,DESCRIPTOR>& rhs){_superGeometry = rhs._superGeometry; _indicator =rhs._indicator; _converter = rhs._converter;};
  bool operator() (T output[], const int input[]) override;
};

/// functor to get pointwise phys torque acting on a boundary with a given indicator on local lattice
template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticePhysBoundaryTorqueIndicator3D final : public SuperLatticePhysF3D<T,DESCRIPTOR> {
private:
  SuperGeometry3D<T>& _superGeometry;
  ParticleIndicatorF3D<T,T>& _indicator;
public:
  SuperLatticePhysBoundaryTorqueIndicator3D(SuperLattice3D<T,DESCRIPTOR>& sLattice,
      SuperGeometry3D<T>& superGeometry,
      ParticleIndicatorF3D<T,T>& indicator,
      const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};


/// functor to get pointwise phys force acting on a boundary with a given material on local lattice
/// see: Caiazzo, Junk: Boundary Forces in lattice Boltzmann: Analysis of MEA
template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticePhysCorrBoundaryForce3D final : public SuperLatticePhysF3D<T,DESCRIPTOR> {
private:
  SuperGeometry3D<T>& _superGeometry;
  const int _material;
public:
  SuperLatticePhysCorrBoundaryForce3D(SuperLattice3D<T,DESCRIPTOR>& sLattice,
                                      SuperGeometry3D<T>& superGeometry,
                                      const int material, const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};


/// functor to get pointwise, lattice-dependent external field
template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticeExternalField3D final : public SuperLatticeF3D<T,DESCRIPTOR> {
private:
  int _beginsAt;
  int _sizeOf;
public:
  SuperLatticeExternalField3D(SuperLattice3D<T,DESCRIPTOR>& sLattice, int beginsAt, int sizeOf);
  bool operator() (T output[], const int input[]) override;
};

/// functor to get pointwise, lattice-dependent porosity values in [0,1]
/// in combination with (Extended)PorousBGKdynamics: 0->solid, 1->fluid
template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticePorosity3D final : public SuperLatticeF3D<T,DESCRIPTOR> {
private:
public:
  SuperLatticePorosity3D(SuperLattice3D<T,DESCRIPTOR>& sLattice);
  bool operator() (T output[], const int input[]);
};

/// functor to get pointwise mesh-independent permeability values in (0,inf) in combination with (Extended)PorousBGKdynamics
/// note: result is cropped to 999999
template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticePhysPermeability3D final : public SuperLatticePhysF3D<T,DESCRIPTOR> {
public:
  SuperLatticePhysPermeability3D(SuperLattice3D<T,DESCRIPTOR>& sLattice, const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]);
};

/// functor to get pointwise mesh-independent permeability values in (0,inf) in combination with (Extended)PorousBGKdynamics
/// note: result is cropped to 1
template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticePhysCroppedPermeability3D final : public SuperLatticePhysF3D<T,DESCRIPTOR> {
public:
  SuperLatticePhysCroppedPermeability3D(SuperLattice3D<T,DESCRIPTOR>& sLattice, const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]);
};


/// computes pointwise -nu/K*u on the lattice, can be used with SuperSum3D as objective
template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticePhysDarcyForce3D final : public SuperLatticePhysF3D<T,DESCRIPTOR> {
private:
  SuperGeometry3D<T>& _superGeometry;
  const int _material;
public:
  SuperLatticePhysDarcyForce3D(SuperLattice3D<T,DESCRIPTOR>& sLattice,
                               SuperGeometry3D<T>& superGeometry,
                               const int material, const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]);
};


/// functor to get a pointwise local average of a passed functor with a given material and radius on local lattice
/// the output data must be of the same size and dimension like f
template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticeAverage3D final : public SuperLatticeF3D<T,DESCRIPTOR> {
private:
  SuperLatticeF3D<T,DESCRIPTOR>& _f;
  SuperGeometry3D<T>& _superGeometry;
  const int _material;
  T _radius;
public:
  SuperLatticeAverage3D(SuperLatticeF3D<T,DESCRIPTOR>& f,
                        SuperGeometry3D<T>& superGeometry, const int material, T radius);
  bool operator() (T output[], const int input[]) override;
};


/// functor that returns pointwise the l2-norm, e.g. of a velocity
template <typename T, template <typename U> class DESCRIPTOR>
class SuperEuklidNorm3D final : public SuperLatticeF3D<T,DESCRIPTOR> {
private:
  SuperLatticeF3D<T,DESCRIPTOR>& _f;
public:
  SuperEuklidNorm3D(SuperLatticeF3D<T,DESCRIPTOR>& f);
  bool operator() (T output[], const int input[]) override;
};


template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticeInterpPhysVelocity3D final : public SuperLatticePhysF3D<T,DESCRIPTOR> {
public:
  SuperLatticeInterpPhysVelocity3D(SuperLattice3D<T,DESCRIPTOR>& sLattice, UnitConverter<T,DESCRIPTOR> const& converter);
  bool operator()(T output[], const int input[]) override;
  void operator()(T output[], const T input[], const int iC);
};

/// functor to get pointwise phys temperature from rho on local lattices
template <typename T, template <typename U> class DESCRIPTOR, template <typename V> class ThermalDESCRIPTOR>
class SuperLatticePhysTemperature3D final : public SuperLatticeThermalPhysF3D<T,DESCRIPTOR,ThermalDESCRIPTOR> { // templatename before <
public:
  SuperLatticePhysTemperature3D(SuperLattice3D<T,ThermalDESCRIPTOR>& sLattice,
                                ThermalUnitConverter<T,DESCRIPTOR,ThermalDESCRIPTOR> const& converter);
  bool operator() (T output[], const int input[]);
};

//template <typename T, template <typename U> class DESCRIPTOR>
//class SuperLatticeInterpPhysVelocity3Degree3Dfinal : public SuperLatticeF3D<T,DESCRIPTOR> {
//private:
//  std::vector<BlockLatticeInterpPhysVelocity3Degree3D<T,DESCRIPTOR>* > _bLattices;
//public:
//  SuperLatticeInterpPhysVelocity3Degree3D(SuperLattice3D<T,DESCRIPTOR>& sLattice, UnitConverter<T,DESCRIPTOR>& conv);
//  bool operator() (T output[], const int input[]) {
//    return 0;
//  }
//  void operator()(T output[], const T input[], const int iC);
//};
//
//template <typename T, template <typename U> class DESCRIPTOR>
//class SuperLatticeInterpDensity3Degree3Dfinal : public SuperLatticeF3D<T,DESCRIPTOR> {
//private:
//  std::vector<BlockLatticeInterpDensity3Degree3D<T,DESCRIPTOR>* > _bLattices;
//public:
//  SuperLatticeInterpDensity3Degree3D(SuperLattice3D<T,DESCRIPTOR>& sLattice, UnitConverter<T,DESCRIPTOR>& conv);
//  bool operator() (T output[], const int input[]) {
//    return 0;
//  }
//  void operator()(T output[], const T input[], const int iC);
//};

/// functor that returns forces acting on a particle surface, returns for every particle 0...2: translational force, 3...5: torque, 7: number of voxels
template <typename T, template <typename U> class DESCRIPTOR>
class SuperLatticeMomentumExchangeForce3D final : public SuperLatticePhysF3D<T,DESCRIPTOR> {
private:
  SuperGeometry3D<T>& _superGeometry;
  std::vector<ParticleIndicatorF3D<T,T>* >& _vectorOfIndicator;
public:
  SuperLatticeMomentumExchangeForce3D(SuperLattice3D<T,DESCRIPTOR>& sLattice,
                                      SuperGeometry3D<T>& superGeometry,
                                      std::vector<ParticleIndicatorF3D<T,T>* >& indicator,
                                      const UnitConverter<T,DESCRIPTOR>& converter);
  bool operator() (T output[], const int input[]) override;
};

} // end namespace olb

#endif
