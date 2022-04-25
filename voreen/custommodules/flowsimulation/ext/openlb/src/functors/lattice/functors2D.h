/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Lukas Baron, Mathias J. Krause
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
#include "blockCalcF2D.h"
#include "blockGeometryFaces2D.h"
#include "blockLatticeIntegralF2D.h"
#include "blockMin2D.h"
#include "blockMax2D.h"
#include "blockAverage2D.h"
#include "blockLatticeRefinementMetricF2D.h"
#include "blockLocalAverage2D.h"
#include "blockReduction2D2D.h"
#include "blockReduction2D1D.h"
#include "reductionF2D.h"
#include "superBaseF2D.h"
#include "superCalcF2D.h"
#include "superConst2D.h"
#include "superGeometryFaces2D.h"
#include "superLatticeIntegralF2D.h"
#include "superMin2D.h"
#include "superMax2D.h"
#include "superAverage2D.h"
#include "superLatticeRefinementMetricF2D.h"
#include "superLocalAverage2D.h"
#include "superErrorNorm2D.h"
#include "indicator/indicator2D.h"
#include "integral/integral2D.h"
#include "timeAveraged/superLatticeTimeAveraged2D.h"
#include "blockRoundingF2D.h"
#include "superRoundingF2D.h"
#include "blockDiscretizationF2D.h"
#include "superDiscretizationF2D.h"
#include "latticePhysDissipation2D.h"
#include "latticePhysField2D.h"
#include "latticeDensity2D.h"
#include "latticeExternal2D.h"
#include "latticeVelocity2D.h"
#include "latticePhysStrainRate2D.h"
#include "latticePhysWallShearStress2D.h"
#include "latticeGeometry2D.h"
#include "latticeRank2D.h"
#include "latticeCuboid2D.h"
#include "latticeExternalScalarField2D.h"
#include "latticePhysPressure2D.h"
#include "latticePhysVelocity2D.h"
#include "latticePhysBoundaryForce2D.h"
#include "latticePhysCorrBoundaryForce2D.h"
#include "latticePorosity2D.h"
#include "latticePhysPermeability2D.h"
#include "latticePhysDarcyForce2D.h"
#include "euklidNorm2D.h"
#include "latticeIndicatorSmoothIndicatorIntersection2D.h"
#include "latticeField2D.h"
#include "latticeAverage2D.h"
#include "latticeGuoZhaoEpsilon2D.h"
#include "latticeGuoZhaoPhysBodyForce2D.h"
#include "latticeGuoZhaoPhysK2D.h"
#include "latticePhysExternalParticleVelocity2D.h"
#include "latticePhysExternalPorosity2D.h"
#include "latticePhysExternalVelocity2D.h"
#include "latticePhysVelocity2D.h"
#include "latticePhysExternalZeta2D.h"
#include "latticePhysHeatFlux2D.h"
#include "latticePhysTemperature2D.h"
#include "latticePorousMomentumLossForce2D.h"
#include "latticeMomentumExchangeForce.h"
#include "latticeStokesDragForce.h"
#include "latticePSMPhysForce2D.h"
#include "latticeVolumeFractionApproximation2D.h"
#include "latticeVolumeFractionPolygonApproximation2D.h"
#include "latticeStrainRate2D.h"
#include "latticeDiscreteNormal2D.h"
