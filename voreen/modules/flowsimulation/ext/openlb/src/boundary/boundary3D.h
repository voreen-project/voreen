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
 * Groups all the 3D include files in the boundaryConditions directory.
 */
#include "advectionDiffusionBoundaryPostProcessor3D.h"
#include "boundaryPostProcessors3D.h"
#include "extendedFiniteDifferenceBoundary3D.h"
#include "inamuroAnalyticalDynamics.h"
#include "inamuroNewtonRaphsonDynamics.h"
#include "offBoundaryPostProcessors3D.h"
#include "rtlbmBoundaryDynamics.h"
#include "wallFunctionBoundaryPostProcessors3D.h"
#include "zouHeDynamics.h"
#include "setLocalVelocityBoundary3D.h"
#include "setInterpolatedVelocityBoundary3D.h"
#include "setLocalPressureBoundary3D.h"
#include "setInterpolatedPressureBoundary3D.h"
#include "setPartialSlipBoundary3D.h"
#include "setLocalConvectionBoundary3D.h"
#include "setInterpolatedConvectionBoundary3D.h"
#include "setAdvectionDiffusionConvectionBoundary3D.h"
#include "setAdvectionDiffusionTemperatureBoundary3D.h"
#include "setExtFieldBoundary3D.h"
#include "setZeroDistributionBoundary3D.h"
#include "setWallFunctionBoundary3D.h"
#include "setFreeEnergyInletBoundary3D.h"
#include "setFreeEnergyOutletBoundary3D.h"
#include "setFreeEnergyWallBoundary3D.h"
#include "setBouzidiVelocityBoundary3D.hh"
#include "setBouzidiZeroVelocityBoundary3D.hh"
#include "setBounceBackVelocityBoundary3D.h"
#include "defineU3D.h"
#include "setSlipBoundary3D.h"
#include "setSlipBoundaryWithDynamics3D.h"
#include "setZouHeVelocityBoundary3D.h"
#include "setZouHePressureBoundary3D.h"
#include "setRtlbmDiffuseTemperatureBoundary3D.h"
#include "setRtlbmDiffuseConstTemperatureBoundary3D.h"
#include "setRtlbmDirectedTemperatureBoundary3D.h"
#include "helper.h"

#include "boundary3D.hh"
