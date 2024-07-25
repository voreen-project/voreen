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
 * Groups all the generic 2D template files in the boundaryConditions directory.
 */

#include "boundaryPostProcessors2D.hh"
#include "extendedFiniteDifferenceBoundary2D.hh"
#include "inamuroAnalyticalDynamics.hh"
#include "inamuroNewtonRaphsonDynamics.hh"
#include "offBoundaryPostProcessors2D.hh"
#include "defineU2D.hh"
#include "setLocalVelocityBoundary2D.hh"
#include "setInterpolatedVelocityBoundary2D.hh"
#include "setLocalPressureBoundary2D.hh"
#include "setInterpolatedPressureBoundary2D.hh"
#include "setPartialSlipBoundary2D.hh"
#include "setFreeEnergyWallBoundary2D.hh"
#include "setFreeEnergyInletBoundary2D.hh"
#include "setFreeEnergyOutletBoundary2D.hh"
#include "setRegularizedTemperatureBoundary2D.hh"
#include "setAdvectionDiffusionTemperatureBoundary2D.hh"
#include "setRegularizedHeatFluxBoundary2D.hh"
#include "setLocalConvectionBoundary2D.hh"
#include "setInterpolatedConvectionBoundary2D.hh"
#include "setSlipBoundary2D.hh"
#include "setSlipBoundaryWithDynamics2D.hh"
#include "setPartialSlipBoundary2D.hh"
#include "zouHeDynamics.hh"
#include "setZouHePressureBoundary2D.hh"
#include "setZouHeVelocityBoundary2D.hh"


