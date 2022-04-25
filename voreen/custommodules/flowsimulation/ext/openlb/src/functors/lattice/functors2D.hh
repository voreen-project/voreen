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

#include "blockBaseF2D.hh"
#include "blockCalcF2D.hh"
#include "blockAverage2D.hh"
#include "blockGeometryFaces2D.hh"
#include "blockLatticeIntegralF2D.hh"
#include "blockMin2D.hh"
#include "blockMax2D.hh"
#include "blockLatticeRefinementMetricF2D.hh"
#include "blockLocalAverage2D.hh"
#include "latticeFrameChangeF3D.hh"
#include "blockReduction2D2D.hh"
#include "blockReduction2D1D.hh"
#include "reductionF2D.hh"
#include "superBaseF2D.hh"
#include "superCalcF2D.hh"
#include "superConst2D.hh"
#include "superGeometryFaces2D.hh"
#include "superLatticeIntegralF2D.hh"
#include "superMin2D.hh"
#include "superMax2D.hh"
#include "superAverage2D.hh"
#include "superLatticeRefinementMetricF2D.hh"
#include "superLocalAverage2D.hh"
#include "superErrorNorm2D.hh"
#include "indicator/indicator2D.hh"
#include "integral/integral2D.hh"
#include "timeAveraged/superLatticeTimeAveraged2D.hh"
#include "blockRoundingF2D.hh"
#include "superRoundingF2D.hh"
#include "blockDiscretizationF2D.hh"
#include "superDiscretizationF2D.hh"
#include "latticeExternal2D.hh"
#include "latticeExternalScalarField2D.hh"
#include "latticePhysDissipation2D.hh"
#include "latticeDensity2D.hh"
#include "latticeVelocity2D.hh"
#include "latticePhysStrainRate2D.hh"
#include "latticePhysWallShearStress2D.hh"
#include "latticeGeometry2D.hh"
#include "latticeRank2D.hh"
#include "latticeCuboid2D.hh"
#include "latticePhysPressure2D.hh"
#include "latticePhysVelocity2D.hh"
#include "latticePhysBoundaryForce2D.hh"
#include "latticePhysCorrBoundaryForce2D.hh"
#include "latticePorosity2D.hh"
#include "latticePhysPermeability2D.hh"
#include "latticePhysDarcyForce2D.hh"
#include "euklidNorm2D.hh"
#include "latticeIndicatorSmoothIndicatorIntersection2D.hh"
#include "latticeField2D.hh"
#include "latticeAverage2D.hh"
#include "latticeGuoZhaoEpsilon2D.hh"
#include "latticeGuoZhaoPhysBodyForce2D.hh"
#include "latticeGuoZhaoPhysK2D.hh"
#include "latticePhysExternalParticleVelocity2D.hh"
#include "latticePhysExternalPorosity2D.hh"
#include "latticePhysExternalVelocity2D.hh"
#include "latticePhysVelocity2D.hh"
#include "latticePhysExternalZeta2D.hh"
#include "latticePhysHeatFlux2D.hh"
#include "latticePhysTemperature2D.hh"
#include "latticePorousMomentumLossForce2D.hh"
#include "latticeMomentumExchangeForce.hh"
#include "latticeStokesDragForce.hh"
#include "latticePSMPhysForce2D.hh"
#include "latticeVolumeFractionApproximation2D.hh"
#include "latticeVolumeFractionPolygonApproximation2D.hh"
#include "latticeStrainRate2D.hh"
#include "latticeDiscreteNormal2D.hh"
