/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Lukas Baron, Mathias J. Krause, Albert Mink0
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
 * Groups all include files for the directory genericFunctions.
 */

#include "blockBaseF2D.hh"
#include "blockBaseF3D.hh"
#include "blockCalcF3D.hh"
#include "blockLatticeIntegralF3D.hh"
#include "blockMin3D.hh"
#include "blockMax3D.hh"
#include "blockAverage3D.hh"
#include "blockGeometryFaces3D.hh"
#include "blockLocalAverage3D.hh"
#include "latticeFrameChangeF3D.hh"
#include "reductionF3D.hh"
#include "blockReduction3D2D.hh"
#include "superBaseF3D.hh"
#include "superCalcF3D.hh"
#include "superConst3D.hh"
#include "superLatticeIntegralF3D.hh"
#include "superLatticeRefinementMetricF3D.hh"
#include "blockLatticeRefinementMetricF3D.hh"
#include "superLocalAverage3D.hh"
#include "superMin3D.hh"
#include "superMax3D.hh"
#include "superAverage3D.hh"
#include "superGeometryFaces3D.hh"
#include "turbulentF3D.hh"
#include "superErrorNorm3D.hh"
#include "indicator/indicator3D.hh"
#include "integral/integral3D.hh"
#include "timeAveraged/superLatticeTimeAveraged3D.hh"
#include "blockRoundingF3D.hh"
#include "superRoundingF3D.hh"
#include "blockDiscretizationF3D.hh"
#include "superDiscretizationF3D.hh"
#include "latticeExternal3D.hh"
#include "latticeExternalScalarField3D.hh"
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
#include "latticeMomentumExchangeForce.hh"
#include "latticeStokesDragForce.hh"
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
#include "latticeDiscreteNormal3D.hh"
#include "blockStatisticF3D.hh"
#include "superStatisticF3D.hh"
//#include "blockStdDeviationF3D.hh"
//#include "superStdDeviationF3D.hh"
