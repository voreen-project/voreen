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

#include "latticeTimeStepScale3D.h"
#include "latticeGuoZhaoPhysBodyForce3D.h"
#include "latticeGuoZhaoPhysK3D.h"
#include "latticeGuoZhaoEpsilon3D.h"
#include "latticeIndicatorSmoothIndicatorIntersection3D.h"
#include "latticeTauFromBoundaryDistance3D.h"
#include "latticePhysPoreSizeDistribution3D.h"
#include "latticePhysBoundaryDistance3D.h"
#include "latticePhysHeatFlux3D.h"
#include "latticeThermalComfort3D.h"
#include "latticePhysTemperature3D.h"
#include "latticePorousMomentumLossForce3D.h"
#include "latticeMomentumExchangeForce3D.h"
#include "latticeInterpPhysVelocity3D.h"
#include "euklidNorm3D.h"
#include "latticePhysDarcyForce3D.h"
#include "latticePhysCroppedPermeability3D.h"
#include "latticePhysPermeability3D.h"
#include "latticeVolumeFractionApproximation3D.h"
#include "latticePorosity3D.h"
#include "latticeField3D.h"
#include "latticePhysCorrBoundaryForce3D.h"
#include "latticePhysHeatFluxBoundary3D.h"
#include "latticePhysWallShearStress3D.h"
#include "latticePSMPhysForce3D.h"
#include "latticePhysBoundaryForce3D.h"
#include "latticePhysExternalScalar3D.h"
#include "latticePhysExternal3D.h"
#include "latticePhysExternalParticleVelocity3D.h"
#include "latticePhysExternalVelocity3D.h"
#include "latticePhysExternalPorosity3D.h"
#include "latticePhysVelocity3D.h"
#include "latticePhysViscosity3D.h"
#include "latticePhysPressure3D.h"
#include "latticeCuboid3D.h"
#include "latticeRank3D.h"
#include "latticeGeometry3D.h"
#include "latticePhysStrainRate3D.h"
#include "latticePhysShearRateMag3D.h"
#include "latticeStrainRate3D.h"
#include "latticeFlux3D.h"
#include "latticeExternalVelocity3D.h"
#include "latticeVelocity3D.h"
#include "latticeDensity3D.h"
#include "latticePhysEffevtiveDissipation3D.h"
#include "latticeEffevtiveDissipation3D.h"
#include "latticePhysDissipation3D.h"
#include "latticeDissipation3D.h"
#include "latticeFpop3D.h"
#include "particleIndicatorF3D.h"

#include "latticeTimeStepScale3D.hh"
#include "latticeGuoZhaoPhysBodyForce3D.hh"
#include "latticeGuoZhaoPhysK3D.hh"
#include "latticeGuoZhaoEpsilon3D.hh"
#include "latticeIndicatorSmoothIndicatorIntersection3D.hh"
#include "latticeTauFromBoundaryDistance3D.hh"
#include "latticePhysPoreSizeDistribution3D.hh"
#include "latticePhysBoundaryDistance3D.hh"
#include "latticePhysHeatFlux3D.hh"
#include "latticeThermalComfort3D.hh"
#include "latticePhysTemperature3D.hh"
#include "latticePorousMomentumLossForce3D.hh"
#include "latticeMomentumExchangeForce3D.hh"
#include "latticeInterpPhysVelocity3D.hh"
#include "euklidNorm3D.hh"
#include "latticePhysDarcyForce3D.hh"
#include "latticePhysCroppedPermeability3D.hh"
#include "latticePhysPermeability3D.hh"
#include "latticeVolumeFractionApproximation3D.hh"
#include "latticePorosity3D.hh"
#include "latticeField3D.hh"
#include "latticePhysCorrBoundaryForce3D.hh"
#include "latticePhysHeatFluxBoundary3D.hh"
#include "latticePhysWallShearStress3D.hh"
#include "latticePSMPhysForce3D.hh"
#include "latticePhysBoundaryForce3D.hh"
#include "latticePhysExternalScalar3D.hh"
#include "latticePhysExternal3D.hh"
#include "latticePhysExternalParticleVelocity3D.hh"
#include "latticePhysExternalVelocity3D.hh"
#include "latticePhysExternalPorosity3D.hh"
#include "latticePhysVelocity3D.hh"
#include "latticePhysViscosity3D.hh"
#include "latticePhysPressure3D.hh"
#include "latticeCuboid3D.hh"
#include "latticeRank3D.hh"
#include "latticeGeometry3D.hh"
#include "latticePhysStrainRate3D.hh"
#include "latticePhysShearRateMag3D.hh"
#include "latticeStrainRate3D.hh"
#include "latticeFlux3D.hh"
#include "latticeExternalVelocity3D.hh"
#include "latticeVelocity3D.hh"
#include "latticeDensity3D.hh"
#include "latticePhysEffevtiveDissipation3D.hh"
#include "latticeEffevtiveDissipation3D.hh"
#include "latticePhysDissipation3D.hh"
#include "latticeDissipation3D.hh"
#include "latticeFpop3D.hh"
#include "latticeKineticEnergy3D.hh"

#include "dynamics/latticeDescriptors.h"

namespace olb {

template class SuperLatticeFpop3D<double,descriptors::D3Q19<>>;
template class SuperLatticeDissipation3D<double,descriptors::D3Q19<>>;
template class SuperLatticePhysDissipation3D<double,descriptors::D3Q19<>>;
template class SuperLatticeEffevtiveDissipation3D<double,descriptors::D3Q19<>>;
template class SuperLatticePhysEffevtiveDissipation3D<double,descriptors::D3Q19<>>;
template class SuperLatticeDensity3D<double,descriptors::D3Q19<>>;
template class SuperLatticeVelocity3D<double,descriptors::D3Q19<>>;
template class SuperLatticeFlux3D<double,descriptors::D3Q19<>>;
template class SuperLatticeStrainRate3D<double,descriptors::D3Q19<>>;
template class SuperLatticeGeometry3D<double,descriptors::D3Q19<>>;
template class SuperLatticeRank3D<double,descriptors::D3Q19<>>;
template class SuperLatticeCuboid3D<double,descriptors::D3Q19<>>;
template class SuperLatticePhysPressure3D<double,descriptors::D3Q19<>>;
template class SuperLatticePhysVelocity3D<double,descriptors::D3Q19<>>;
template class SuperLatticePhysBoundaryForce3D<double,descriptors::D3Q19<>>;
template class SuperLatticePhysCorrBoundaryForce3D<double,descriptors::D3Q19<>>;
template class SuperEuklidNorm3D<double,descriptors::D3Q19<>>;
template class SuperLatticeInterpPhysVelocity3D<double,descriptors::D3Q19<>>;
template class SuperLatticeIndicatorSmoothIndicatorIntersection3D<double,descriptors::D3Q19<>,false>;
template class BlockLatticeFpop3D<double,descriptors::D3Q19<>>;
template class BlockLatticeDissipation3D<double,descriptors::D3Q19<>>;
template class BlockLatticePhysDissipation3D<double,descriptors::D3Q19<>>;
template class BlockLatticeEffevtiveDissipation3D<double,descriptors::D3Q19<>>;
template class BlockLatticePhysEffevtiveDissipation3D<double,descriptors::D3Q19<>>;
template class BlockLatticeDensity3D<double,descriptors::D3Q19<>>;
template class BlockLatticeVelocity3D<double,descriptors::D3Q19<>>;
template class BlockLatticeFlux3D<double,descriptors::D3Q19<>>;
template class BlockLatticeGeometry3D<double,descriptors::D3Q19<>>;
template class BlockLatticeRank3D<double,descriptors::D3Q19<>>;
template class BlockLatticeCuboid3D<double,descriptors::D3Q19<>>;
template class BlockLatticePhysPressure3D<double,descriptors::D3Q19<>>;
template class BlockLatticePhysVelocity3D<double,descriptors::D3Q19<>>;
template class BlockLatticeStrainRate3D<double,descriptors::D3Q19<>>;
template class BlockLatticePhysBoundaryForce3D<double,descriptors::D3Q19<>>;
template class BlockLatticePhysCorrBoundaryForce3D<double,descriptors::D3Q19<>>;
template class BlockEuklidNorm3D<double,descriptors::D3Q19<>>;
template class BlockLatticeInterpPhysVelocity3D<double,descriptors::D3Q19<>>;
template class BlockLatticeIndicatorSmoothIndicatorIntersection3D<double,descriptors::D3Q19<>,false>;

} // end namespace olb
