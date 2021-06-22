/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2014 Albert Mink, Mathias J. Krause
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
#include "latticeDissipation2D.h"
#include "latticeDissipation2D.hh"
#include "latticePhysDissipation2D.h"
#include "latticePhysDissipation2D.hh"
#include "latticeDensity2D.h"
#include "latticeDensity2D.hh"
#include "latticeVelocity2D.h"
#include "latticeVelocity2D.hh"
#include "latticePhysStrainRate2D.h"
#include "latticePhysStrainRate2D.hh"
#include "latticePhysWallShearStress2D.h"
#include "latticePhysWallShearStress2D.hh"
#include "latticeGeometry2D.h"
#include "latticeGeometry2D.hh"
#include "latticeRank2D.h"
#include "latticeRank2D.hh"
#include "latticeCuboid2D.h"
#include "latticeCuboid2D.hh"
#include "latticePhysPressure2D.h"
#include "latticePhysPressure2D.hh"
#include "latticePhysVelocity2D.h"
#include "latticePhysVelocity2D.hh"
#include "latticePhysBoundaryForce2D.h"
#include "latticePhysBoundaryForce2D.hh"
#include "latticePhysCorrBoundaryForce2D.h"
#include "latticePhysCorrBoundaryForce2D.hh"
#include "latticePorosity2D.h"
#include "latticePorosity2D.hh"
#include "latticePhysPermeability2D.h"
#include "latticePhysPermeability2D.hh"
#include "latticePhysDarcyForce2D.h"
#include "latticePhysDarcyForce2D.hh"
#include "euklidNorm2D.h"
#include "euklidNorm2D.hh"
#include "latticeIndicatorSmoothIndicatorIntersection2D.h"
#include "latticeIndicatorSmoothIndicatorIntersection2D.hh"


#include "latticeAverage2D.h"
#include "latticeAverage2D.hh"
#include "latticeGuoZhaoEpsilon2D.h"
#include "latticeGuoZhaoEpsilon2D.hh"
#include "latticeGuoZhaoPhysBodyForce2D.h"
#include "latticeGuoZhaoPhysBodyForce2D.hh"
#include "latticeGuoZhaoPhysK2D.h"
#include "latticeGuoZhaoPhysK2D.hh"
#include "latticePhysExternalParticleVelocity2D.h"
#include "latticePhysExternalParticleVelocity2D.hh"
#include "latticePhysExternalPorosity2D.h"
#include "latticePhysExternalPorosity2D.hh"
#include "latticePhysExternalVelocity2D.h"
#include "latticePhysExternalVelocity2D.hh"
#include "latticePhysExternalZeta2D.h"
#include "latticePhysExternalZeta2D.hh"
#include "latticePhysHeatFlux2D.h"
#include "latticePhysHeatFlux2D.hh"
#include "latticePhysTemperature2D.h"
#include "latticePhysTemperature2D.hh"
#include "latticePorousMomentumLossForce2D.h"
#include "latticePorousMomentumLossForce2D.hh"
#include "latticeMomentumExchangeForce2D.h"
#include "latticeMomentumExchangeForce2D.hh"
#include "latticePSMPhysForce2D.h"
#include "latticePSMPhysForce2D.h"
#include "latticeVolumeFractionApproximation2D.h"
#include "latticeVolumeFractionApproximation2D.hh"
#include "latticeVolumeFractionPolygonApproximation2D.h"
#include "latticeVolumeFractionPolygonApproximation2D.hh"
#include "latticeStrainRate2D.h"
#include "latticeStrainRate2D.hh"

#include "dynamics/latticeDescriptors.h"

namespace olb {

template class SuperLatticeDissipation2D<double,descriptors::D2Q9<>>;
template class SuperLatticePhysDissipation2D<double,descriptors::D2Q9<>>;
template class SuperLatticeDensity2D<double,descriptors::D2Q9<>>;
template class SuperLatticeVelocity2D<double,descriptors::D2Q9<>>;
template class SuperLatticePhysStrainRate2D<double,descriptors::D2Q9<>>;
template class SuperLatticePhysWallShearStress2D<double,descriptors::D2Q9<>>;
template class SuperLatticeGeometry2D<double,descriptors::D2Q9<>>;
template class SuperLatticeRank2D<double,descriptors::D2Q9<>>;
template class SuperLatticeCuboid2D<double,descriptors::D2Q9<>>;
template class SuperLatticePhysPressure2D<double,descriptors::D2Q9<>>;
template class SuperLatticePhysVelocity2D<double,descriptors::D2Q9<>>;
template class SuperLatticePhysBoundaryForce2D<double,descriptors::D2Q9<>>;
template class SuperLatticePhysCorrBoundaryForce2D<double,descriptors::D2Q9<>>;
template class SuperLatticePorosity2D<double,descriptors::D2Q9<>>;
template class SuperLatticePhysPermeability2D<double,descriptors::D2Q9<>>;
template class SuperLatticePhysDarcyForce2D<double,descriptors::D2Q9<>>;
template class SuperEuklidNorm2D<double,descriptors::D2Q9<>>;
template class SuperLatticeIndicatorSmoothIndicatorIntersection2D<double,descriptors::D2Q9<>,false>;

template class BlockLatticeDissipation2D<double,descriptors::D2Q9<>>;
template class BlockLatticePhysDissipation2D<double,descriptors::D2Q9<>>;
template class BlockLatticeDensity2D<double,descriptors::D2Q9<>>;
template class BlockLatticeVelocity2D<double,descriptors::D2Q9<>>;
template class BlockLatticePhysStrainRate2D<double,descriptors::D2Q9<>>;
template class BlockLatticePhysWallShearStress2D<double,descriptors::D2Q9<>>;
template class BlockLatticeGeometry2D<double,descriptors::D2Q9<>>;
template class BlockLatticeRank2D<double,descriptors::D2Q9<>>;
template class BlockLatticeCuboid2D<double,descriptors::D2Q9<>>;
template class BlockLatticePhysPressure2D<double,descriptors::D2Q9<>>;
template class BlockLatticePhysVelocity2D<double,descriptors::D2Q9<>>;
//template class BlockLatticePhysExternalPorosity2D<double,descriptors::D2Q9<>>;
//template class BlockLatticePhysExternalVelocity2D<double,descriptors::D2Q9<>>;
template class BlockLatticePhysBoundaryForce2D<double,descriptors::D2Q9<>>;
template class BlockLatticePhysCorrBoundaryForce2D<double,descriptors::D2Q9<>>;
template class BlockLatticePorosity2D<double,descriptors::D2Q9<>>;
template class BlockLatticePhysPermeability2D<double,descriptors::D2Q9<>>;
template class BlockLatticePhysDarcyForce2D<double,descriptors::D2Q9<>>;
template class BlockLatticeAverage2D<double,descriptors::D2Q9<>>;
template class BlockEuklidNorm2D<double,descriptors::D2Q9<>>;
template class BlockLatticeIndicatorSmoothIndicatorIntersection2D<double,descriptors::D2Q9<>,false>;

} // end namespace olb
