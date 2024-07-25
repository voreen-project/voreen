/*  This file is part of the OpenLB library
 *
*  Copyright (C) 2006-2015 Mathias J. Krause
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
 * Groups all the generic 3D template files in the complexDynamics directory.
 */

#include "advectionDiffusionForces.hh"
#include "entropicDynamics.hh"
#include "freeEnergyPostProcessor3D.hh"
#include "freeSurfacePostProcessor3D.hh"
#include "interactionPotential.hh"
#include "navierStokesAdvectionDiffusionCouplingPostProcessor3D.hh"
#include "legacy/porousBGKdynamics.hh"
#include "rtlbmDynamics.hh"
#include "shanChenDynOmegaForcedPostProcessor3D.hh"
#include "shanChenForcedPostProcessor3D.hh"
#include "shanChenForcedSingleComponentPostProcessor3D.hh"
#include "legacy/smagorinskyBGKdynamics.hh"
#include "stochasticSGSdynamics.hh"

#include "porousAdvectionDiffusionDynamics.hh"
#include "porousForcedBGKDynamics.hh"
#include "guoZhaoDynamics.hh"
