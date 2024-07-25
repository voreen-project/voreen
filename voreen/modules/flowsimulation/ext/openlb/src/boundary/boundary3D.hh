/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007 the OpenLB project
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
 * Groups all the generic 3D template files in the boundaryConditions directory.
 */
#include "advectionDiffusionBoundaryPostProcessor3D.hh"
#include "boundaryPostProcessors3D.hh"
#include "extendedFiniteDifferenceBoundary3D.hh"
#include "inamuroAnalyticalDynamics.hh"
#include "inamuroNewtonRaphsonDynamics.hh"
#include "offBoundaryPostProcessors3D.hh"
#include "rtlbmBoundaryDynamics.hh"
#include "wallFunctionBoundaryPostProcessors3D.hh"
#include "zouHeDynamics.hh"
#include "setLocalVelocityBoundary3D.hh"
#include "setInterpolatedVelocityBoundary3D.hh"
#include "setLocalPressureBoundary3D.hh"
#include "setInterpolatedPressureBoundary3D.hh"
#include "setPartialSlipBoundary3D.hh"
#include "setLocalConvectionBoundary3D.hh"
#include "setInterpolatedConvectionBoundary3D.hh"
#include "setAdvectionDiffusionConvectionBoundary3D.hh"
#include "setAdvectionDiffusionTemperatureBoundary3D.hh"
#include "setExtFieldBoundary3D.hh"
#include "setZeroDistributionBoundary3D.hh"
#include "setWallFunctionBoundary3D.hh"
#include "setFreeEnergyInletBoundary3D.hh"
#include "setFreeEnergyOutletBoundary3D.hh"
#include "setFreeEnergyWallBoundary3D.hh"
#include "setBounceBackVelocityBoundary3D.hh"
#include "defineU3D.hh"
#include "setSlipBoundary3D.hh"
#include "setSlipBoundaryWithDynamics3D.hh"
#include "setZouHeVelocityBoundary3D.hh"
#include "setZouHePressureBoundary3D.hh"
#include "setRtlbmDiffuseTemperatureBoundary3D.hh"
#include "setRtlbmDiffuseConstTemperatureBoundary3D.hh"
#include "setRtlbmDirectedTemperatureBoundary3D.hh"

