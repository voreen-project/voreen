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
 * Groups all the 2D include files in the boundaryConditions directory.
 */
#include "boundaryPostProcessors2D.h"
#include "extendedFiniteDifferenceBoundary2D.h"
#include "inamuroAnalyticalDynamics.h"
#include "inamuroNewtonRaphsonDynamics.h"
#include "offBoundaryPostProcessors2D.h"
#include "defineU2D.h"
#include "setLocalVelocityBoundary2D.h"
#include "setInterpolatedVelocityBoundary2D.h"
#include "setLocalPressureBoundary2D.h"
#include "setInterpolatedPressureBoundary2D.h"
#include "setPartialSlipBoundary2D.h"
#include "setFreeEnergyWallBoundary2D.h"
#include "setFreeEnergyInletBoundary2D.h"
#include "setFreeEnergyOutletBoundary2D.h"
#include "setRegularizedTemperatureBoundary2D.h"
#include "setAdvectionDiffusionTemperatureBoundary2D.h"
#include "setRegularizedHeatFluxBoundary2D.h"
#include "setLocalConvectionBoundary2D.h"
#include "setInterpolatedConvectionBoundary2D.h"
#include "setSlipBoundary2D.h"
#include "setSlipBoundaryWithDynamics2D.h"
#include "setPartialSlipBoundary2D.h"
#include "zouHeDynamics.h"
#include "setZouHePressureBoundary2D.h"
#include "setZouHeVelocityBoundary2D.h"
#include "helper.h"
#include "setBouzidiVelocityBoundary2D.hh"
#include "setBouzidiZeroVelocityBoundary2D.hh"

#include "boundary2D.hh"
