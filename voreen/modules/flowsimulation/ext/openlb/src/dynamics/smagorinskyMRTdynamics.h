/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Patrick Nathen, Mathias J. Krause
 *                2021 Adrian Kummerlaender
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
 * MRT Dynamics with adjusted omega -- header file.
 */
#ifndef SMAGORINSKY_MRT_DYNAMICS_H
#define SMAGORINSKY_MRT_DYNAMICS_H

#include "interface.h"
#include "collisionMRT.h"

namespace olb {

/// Smagorinsky MRT collision step
template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using SmagorinskyMRTdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  collision::SmagorinskyEffectiveOmega<collision::MRT>
>;

/// Smagorinsky MRT collision step with Ladd-Verberg forcing
template<typename T, typename DESCRIPTOR, typename MOMENTA=momenta::BulkTuple>
using SmagorinskyForcedMRTdynamics = dynamics::Tuple<
  T, DESCRIPTOR,
  MOMENTA,
  equilibria::SecondOrder,
  collision::SmagorinskyEffectiveOmega<collision::MRT>,
  forcing::LaddVerberg
>;

}

#endif
