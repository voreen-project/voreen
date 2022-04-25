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

/** \file
 * Groups all include files for the directory genericFunctions.
 */

#include "blockBaseF2D.h"
#include "blockBaseF3D.h"
#include "blockCalcF3D.h"
#include "blockLatticeIntegralF3D.h"
#include "blockMin3D.h"
#include "blockMax3D.h"
#include "blockAverage3D.h"
#include "blockGeometryFaces3D.h"
#include "blockLocalAverage3D.h"
#include "latticeFrameChangeF3D.h"
#include "reductionF3D.h"
#include "blockReduction3D2D.h"
#include "blockLatticeRefinementMetricF3D.h"
#include "superLatticeRefinementMetricF3D.h"
#include "superBaseF3D.h"
#include "superCalcF3D.h"
#include "superConst3D.h"
#include "superLatticeIntegralF3D.h"
#include "superMin3D.h"
#include "superMax3D.h"
#include "superAverage3D.h"
#include "superGeometryFaces3D.h"
#include "superLocalAverage3D.h"
#include "turbulentF3D.h"
#include "superErrorNorm3D.h"
#include "indicator/indicator3D.h"
#include "integral/integral3D.h"
#include "timeAveraged/superLatticeTimeAveraged3D.h"
#include "superRoundingF3D.h"
#include "blockRoundingF3D.h"
#include "superDiscretizationF3D.h"
#include "blockDiscretizationF3D.h"
#include "latticeExternal3D.h"
#include "latticeExternalScalarField3D.h"
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
#include "latticeMomentumExchangeForce.h"
#include "latticeStokesDragForce.h"
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
#include "latticePhysExternalParticleVelocity3D.h"
#include "latticePhysExternalVelocity3D.h"
#include "latticePhysExternalPorosity3D.h"
#include "latticePhysField3D.h"
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
#include "latticeKineticEnergy3D.h"
#include "particleIndicatorF3D.h"
#include "latticeDiscreteNormal3D.h"
#include "blockStatisticF3D.h"
#include "superStatisticF3D.h"
//#include "blockStdDeviationF3D.h"
//#include "superStdDeviationF3D.h"
