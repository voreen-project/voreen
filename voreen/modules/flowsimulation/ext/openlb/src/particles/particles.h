/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Nicolas Hafen, Thomas Henn, Mathias J. Krause, Marie-Luise Maier
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
 * Groups all the include .h-files for 3D particles in
 * the particles directory.
 */

#include "particles/descriptor/particleDescriptors.h"
#include "particles/descriptor/particleDescriptorAlias.h"
#include "particle.h"
#include "particleSystem.h"
#include "particles/functions/bodyMotionFunctions.h"
#include "particles/functions/particleDynamicsFunctions.h"
#include "particles/dynamics/particleDynamics.h"
#include "functions/particleTasks.h"
#include "particleManager.h"
#include "particles/functions/particleMotionFunctions.h"
#include "particles/functions/particleCreatorFunctions.h"
#include "particles/functions/particleUtilities.h"
#include "resolved/smoothIndicatorInteraction.h"
#include "resolved/blockLatticeInteraction.h"
#include "resolved/superLatticeInteraction.h"
#include "resolved/momentumExchangeForce.h"
